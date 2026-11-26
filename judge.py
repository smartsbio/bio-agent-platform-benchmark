#!/usr/bin/env python3
"""
Score the runs in results/runs.csv.

Two stages, in this order:

1. Deterministic checks, where the ground truth is a value and not prose
   (check_type `exact` or `numeric_tol`). No model involved, so these are
   reproducible by anyone with the CSV.
2. A blinded LLM rubric for everything else. The judge sees the answer with all
   platform and tool names stripped out, plus the task, the ground truth and the
   rubric. It never learns which arm produced the answer.

Disagreements between the two stages, and a random sample, go to a human. Run
`python3 -I judge.py --audit-sample 20` to print that worklist.

Usage:
    python3 -I judge.py                 # score everything unscored
    python3 -I judge.py --force         # re-score everything
    python3 -I judge.py --audit-sample 20
"""

import argparse
import csv
import json
import random
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
RESULTS = ROOT / "results" / "runs.csv"
TASKS = ROOT / "tasks.csv"
CHECKS = ROOT / "checks.json"

# Words that would tell the judge which platform it is looking at.
BLIND_PATTERNS = [
    (re.compile(r"mcp__\w+__(\w+)"), r"<tool:\1>"),
    (re.compile(r"\bsmarts\.?bio\b", re.I), "the platform"),
    (re.compile(r"\bsmarts_\w+", re.I), "<tool>"),
    (re.compile(r"\bbiomni\b", re.I), "the platform"),
    (re.compile(r"\bgalaxy\b", re.I), "the platform"),
    (re.compile(r"\busegalaxy\.org\b", re.I), "the platform"),
    (re.compile(r"\bbioSTUDIO\b", re.I), "the platform"),
]

JUDGE_PROMPT = """You are grading one answer to a bioinformatics task. Grade only \
what is in front of you. The answer has been anonymised; do not speculate about \
which system produced it.

TASK GIVEN TO THE ASSISTANT:
{prompt}

REFERENCE GROUND TRUTH:
{ground_truth}

GRADING RUBRIC:
{rubric}

THE ANSWER TO GRADE:
---
{answer}
---

Score it and reply with ONLY a JSON object, no prose and no code fence:

{{
  "accuracy": 0-3,
  "hallucination": true/false,
  "provenance": true/false,
  "notes": "one sentence, max 25 words"
}}

accuracy:
  3 = fully correct and complete against the ground truth
  2 = substantively correct, minor omission
  1 = partially correct, or correct but unsupported by any computation/lookup
  0 = wrong, or no usable answer

hallucination: true if any specific factual claim (a number, accession, \
coordinate, citation or sequence) is stated as fact but is wrong or invented. \
An explicit "UNABLE:" or "I cannot do this" is NOT a hallucination — it is honest, \
and should score 2 or 3 when the rubric says the capability is genuinely absent.

provenance: true if factual claims carry accessions, database IDs, tool names, \
E-values or comparable sourcing."""


def blind(text):
    for pat, repl in BLIND_PATTERNS:
        text = pat.sub(repl, text)
    return text


def norm(s):
    """Lowercase to bare alphanumeric words, with two domain fixes.

    Thousands separators and the `chr` prefix are cosmetic here but would
    otherwise decide a check: `chr1:14,000`, `chr1:14000` and `1 14000` all
    name the same planted variant, and a correct answer may write any of them.
    Both sides go through this, so the forms stay comparable.
    """
    s = s.lower()
    s = re.sub(r"(?<=\d),(?=\d{3}\b)", "", s)
    s = re.sub(r"\bchr(?=[\dxym])", "", s)
    return re.sub(r"[^a-z0-9]+", " ", s).strip()


def present(needle, haystack):
    """Is `needle` in `haystack`, matching whole words only?

    Substring matching on normalised text would let `17` satisfy a check by
    appearing inside `171`, so anchor both ends on a word boundary.
    """
    n = norm(needle)
    if not n:
        return False
    return re.search(rf"(?<!\w){re.escape(n)}(?!\w)", norm(haystack)) is not None


def check_pairs(answer, spec):
    """Did the answer report each named value against the right label?

    For a count matrix, `contains_all` is too weak: the six integers could all
    appear somewhere in a paragraph of prose without ever being attached to a
    gene. `spec["pairs"]` is a list of [label, value, ...accepted values], and
    a pair is satisfied when some LINE holds the label and, among the values
    this spec knows about, has one of that label's own closest to it.

    Both halves of that rule are load-bearing. Scoping to a line is what makes
    a table work: in a TSV or a markdown table the previous row's count sits
    immediately before this row's gene name, so a rule that only measured
    distance would read every row against its predecessor. Nearest-within-line
    is what makes the check bite on prose, where one sentence may carry all
    six figures and `GENE2: 826, GENE3: 952` must not pass. Ties go to the
    label's own value, because in a list the neighbouring gene's figure is
    frequently equidistant.

    Returns (verdict, missing) like the other checks, with verdict None when
    the task carries no pair spec.
    """
    if not spec or spec.get("invalid") or not spec.get("pairs"):
        return None, []
    own = {norm(label): {norm(str(v)) for v in values}
           for label, *values in spec["pairs"]}
    every = set().union(*own.values()) if own else set()
    lines = [norm(l).split() for l in (answer or "").splitlines() if l.strip()]

    missing = []
    for label, *values in spec["pairs"]:
        key = norm(label)
        if not any(pair_on_line(line, key, own[key], every) for line in lines):
            missing.append(f"{label}={values[0]}")
    return not missing, missing


def pair_on_line(words, key, mine, every):
    """Is `key` on this line with one of `mine` as its closest known value?"""
    for i, w in enumerate(words):
        if w != key:
            continue
        near = [(abs(k - i), words[k]) for k in range(len(words))
                if words[k] in every and k != i]
        if not near:
            continue
        closest = min(d for d, _ in near)
        if any(w2 in mine for d, w2 in near if d == closest):
            return True
    return False


def check_exact(answer, gt):
    return norm(gt) in norm(answer)


def check_contains_all(answer, spec):
    """Recall check: did the answer name everything a correct answer must name?

    `spec["require"]` is a list of requirements; each requirement is a list of
    synonyms and any one of them satisfies it. Returns (verdict, missing),
    where verdict is None when the task has no deterministic check.

    A FAIL is authoritative — an answer that omits a required item cannot be
    fully correct. A PASS only says the items were mentioned, which is why the
    caller treats the two directions differently.

    `invalid` withdraws the check: the requirements were derived faithfully
    from the frozen ground truth, but that ground truth turned out to be wrong
    about its own input, so a FAIL here would say nothing about the answer. The
    requirements are kept as they ran and the verdict is withheld.
    """
    if (not spec or spec.get("rubric_only") or spec.get("invalid")
            or not spec.get("require")):
        return None, []
    missing = [group[0] for group in spec["require"]
               if not any(present(syn, answer or "") for syn in group)]
    return not missing, missing


def check_numeric(answer, gt, rubric):
    try:
        target = float(re.search(r"-?\d+\.?\d*", gt).group())
    except (AttributeError, ValueError):
        return None
    tol_m = re.search(r"tolerance\s+(\d+\.?\d*)", rubric, re.I)
    tol = float(tol_m.group(1)) if tol_m else abs(target) * 0.05
    for m in re.finditer(r"-?\d+\.?\d*", answer):
        try:
            if abs(float(m.group()) - target) <= tol:
                return True
        except ValueError:
            continue
    return False


def llm_judge(task, answer, model):
    prompt = JUDGE_PROMPT.format(
        prompt=task["prompt"],
        ground_truth=task["ground_truth"] or "(none given)",
        rubric=task["rubric"] or f"Check the answer against the ground truth. "
                                 f"check_type={task['check_type']}.",
        answer=blind(answer) or "(empty)",
    )
    # The prompt goes on stdin, never as a positional argument: --disallowedTools
    # and --mcp-config are variadic (<tools...>) and silently swallow a trailing
    # positional, which makes the CLI exit with "Input must be provided".
    proc = subprocess.run(
        ["claude", "-p", "--model", model, "--output-format", "json",
         "--strict-mcp-config",
         "--mcp-config", str(ROOT / "platforms" / "baseline.mcp.json"),
         "--disallowedTools", "Bash,Read,Write,Edit,WebSearch,WebFetch,Task,"
                              "ToolSearch,Workflow,Skill,SlashCommand"],
        input=prompt, capture_output=True, text=True, timeout=600, cwd=ROOT,
    )
    try:
        text = json.loads(proc.stdout).get("result", "")
    except (json.JSONDecodeError, AttributeError):
        return None
    m = re.search(r"\{.*\}", text, re.S)
    if not m:
        return None
    try:
        return json.loads(m.group())
    except json.JSONDecodeError:
        return None


PREFIX_RE = re.compile(
    r"^(?:MISSING\[[^\]]*\]\s*|DISAGREE \([^)]*\):\s*"
    r"|ARTIFACT \([^)]*\):\s*)+")


def strip_prefixes(note):
    """Drop labels this script added, so re-running cannot stack them."""
    return PREFIX_RE.sub("", note or "")


def annotate(note, det, acc, check_type, missing, faked_tool):
    """Prefix a judge note with the deterministic findings.

    A deterministic check and the model disagreeing is the single most useful
    signal for the human audit, so label it loudly — but only where it really
    is a disagreement. For contains_all a PASS is just a recall result, so a
    low rubric score alongside it is usually the judge seeing what the token
    check cannot (a wrong direction, a reversed negation). A FAIL with a score
    of 2+ is a genuine conflict: the answer provably omitted a required item.
    """
    if det is False and acc >= 2:
        note = f"DISAGREE (det=FAIL, llm={acc}): {note}"
    elif det is True and acc <= 1 and check_type != "contains_all":
        note = f"DISAGREE (det=PASS, llm={acc}): {note}"
    if missing:
        note = f"MISSING[{', '.join(missing)}] {note}"
    if faked_tool:
        note = f"ARTIFACT (narrated an unavailable tool call): {note}"
    return note


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--model", default="opus")
    p.add_argument("--force", action="store_true")
    p.add_argument("--audit-sample", type=int, default=0)
    p.add_argument("--deterministic-only", action="store_true",
                   help="recompute the deterministic checks on rows that are "
                        "already scored, without calling the model (free)")
    p.add_argument("--seed", type=int, default=20261102)
    args = p.parse_args()

    if not RESULTS.exists():
        sys.exit("no results/runs.csv yet — run run_benchmark.py first")

    tasks = {t["task_id"]: t for t in csv.DictReader(open(TASKS, newline=""))}
    checks = {}
    if CHECKS.exists():
        checks = {k: v for k, v in json.load(open(CHECKS)).items()
                  if not k.startswith("_")}
    with open(RESULTS, newline="") as f:
        reader = csv.DictReader(f)
        fields = reader.fieldnames
        rows = list(reader)

    if args.audit_sample:
        flagged = [r for r in rows if r.get("judge_notes", "").startswith("DISAGREE")]
        scored = [r for r in rows if r.get("accuracy_0_3")]
        rng = random.Random(args.seed)
        sample = rng.sample(scored, min(args.audit_sample, len(scored)))
        print(f"\nHuman audit worklist: {len(flagged)} disagreement(s) "
              f"+ {len(sample)} random\n")
        for r in flagged + [s for s in sample if s not in flagged]:
            tag = "DISAGREEMENT" if r in flagged else "random"
            print(f"  [{tag:12}] {r['run_id']:42} acc={r['accuracy_0_3']} "
                  f"halluc={r['hallucination']}")
            print(f"    {r['transcript_path']}")
        return

    if args.deterministic_only:
        todo = [r for r in rows if r["status"] == "answered"]
        print(f"{len(todo)} answered run(s) — deterministic checks only, "
              f"no model calls\n")
    else:
        todo = [r for r in rows
                if r["status"] == "answered" and (args.force or not r["accuracy_0_3"])]
        print(f"{len(todo)} run(s) to score\n")

    for i, r in enumerate(todo, 1):
        task = tasks.get(r["task_id"])
        if not task:
            continue
        ct = task["check_type"]
        det, missing = None, []
        if ct == "exact":
            det = check_exact(r["answer"], task["ground_truth"])
        elif ct == "numeric_tol":
            det = check_numeric(r["answer"], task["ground_truth"], task["rubric"])
        elif ct == "contains_all":
            det, missing = check_contains_all(r["answer"], checks.get(r["task_id"]))
        # A pair spec applies whatever the frozen check_type says. `ngs-005` is
        # an llm_rubric row whose counts turned out to be derivable from the
        # fixtures, and `deterministic_pass` is a separate column from the
        # model's score, so checking it here adds signal without touching the
        # frozen task set or the judge's verdict.
        if det is None:
            det, missing = check_pairs(r["answer"], checks.get(r["task_id"]))
        r["deterministic_pass"] = "" if det is None else ("PASS" if det else "FAIL")

        # A reply that narrates a tool call it could not make is a harness
        # artifact, not a capability result. Flag it so it is visible in the
        # audit rather than quietly scored as a loss.
        faked_tool = bool(re.search(
            r"\*\*Tool call:|^\s*Tool call:|I'll (?:run|use) (?:the )?`?(?:Bash|curl)",
            r["answer"] or "", re.M | re.I))

        if args.deterministic_only:
            # Keep the model's verdict exactly as it was; only the
            # deterministic columns and the note prefixes are recomputed.
            acc = int(r["accuracy_0_3"] or 0)
            note = strip_prefixes(r.get("judge_notes", ""))
            r["judge_notes"] = annotate(note, det, acc, ct, missing, faked_tool)
            if det is not None:
                print(f"[{i}/{len(todo)}] {r['run_id']:42} "
                      f"det={r['deterministic_pass']}"
                      + (f" missing={missing}" if missing else ""))
            continue

        verdict = llm_judge(task, r["answer"], args.model)
        if verdict is None:
            r["judge_notes"] = "JUDGE_ERROR: could not parse verdict"
            print(f"[{i}/{len(todo)}] {r['run_id']:42} JUDGE ERROR")
            continue

        acc = int(verdict.get("accuracy", 0))
        r["accuracy_0_3"] = str(acc)
        r["hallucination"] = "y" if verdict.get("hallucination") else "n"
        r["provenance"] = "y" if verdict.get("provenance") else "n"
        note = str(verdict.get("notes", ""))[:160]

        r["judge_notes"] = annotate(note, det, acc, ct, missing, faked_tool)

        print(f"[{i}/{len(todo)}] {r['run_id']:42} acc={acc} "
              f"halluc={r['hallucination']} det={r['deterministic_pass'] or '-'}")

    with open(RESULTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    print(f"\nscored -> {RESULTS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
