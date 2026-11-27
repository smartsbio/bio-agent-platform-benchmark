#!/usr/bin/env python3
"""
Reproduce every number reported in the paper from the primary files.

This script implements the two analysis specifications of §2.5. It is the sole
and authoritative implementation: the figures printed here are the figures in
the manuscript.

Inputs, all in this repository:
    tasks.csv                     the frozen task set
    results/runs.csv              one row per run
    results/reference_audit.csv   the audit layer, never merged into runs.csv

No model and no third-party packages: given the same three files, anyone with a
stock Python 3 gets exactly these numbers.

Usage:
    python3 -I analyze.py            # both specifications
    python3 -I analyze.py --csv      # also write results/spec_*.csv

The two specifications (§2.5), in full:

  PRIMARY (audit-corrected), n = 46 per condition
    - the two defective tasks of §2.2 are dropped from every condition;
      they are identified by an `arm == "ALL"` audit row with
      audit_accuracy == "VOID", not by a hard-coded task id
    - audit corrections are applied where the audit disagreed with the judge
    - timeouts score 0 (intention-to-complete)
    - no condition carries an arm-specific exclusion

  SENSITIVITY (raw), n = 48 per condition
    - all task records retained, no audit material used at all
    - timeouts score 0

Rules common to both, each with the clause of the paper that fixes it:

  Repetition filter (§2.5, final paragraph). Only `rep == 1` records enter
  either specification. Three repetition-2 records exist, for three
  control-condition tasks only; they do not form a balanced repeated-measures
  design, so they are excluded rather than averaged in. This is why a run count
  of 46 per condition, not 195, is the denominator.

  Staleness (§2.5). An audit correction applies only to a run that started on
  or before that audit's date, so a cell re-measured after being audited
  carries its own score rather than a stale correction. Flagged in the paper as
  a post-hoc choice; its most visible effect discards a correction that would
  have raised the authoring condition's score on a CRISPR task from 1 to 3.

  What counts as a score correction. Only an audit row whose `audit_accuracy`
  is a digit. Three kinds of audit row therefore do NOT move a score, and the
  distinction is deliberate in each case:

    - "INFRA" (mol-002, smarts.bio). The pre-registered plan records quota
      failures as infra_error. That run nonetheless returned an answer and was
      graded 2, so §2.5 retains it. Dropping it would have removed the single
      task on which smarts.bio was outscored by both comparator conditions,
      while comparator runs that returned nothing scored 0 — an asymmetry the
      paper declines to adopt. Retaining it LOWERS the authoring condition's
      mean. The audit's own verdict is preserved in the audit file unedited;
      this override lives here, in code, where it is visible.
    - blank (des-006, smarts.bio). Same shape: the audit proposed voiding the
      cell as infra_error. §2.5 allows no arm-specific exclusion, so the
      judge's score of 1 stands, again against the authoring condition.
    - a note row with no scores (struct-003 "ALL", and the harness-defect row).
      These record ground-truth and harness defects in the audit trail. The
      task text was not changed after the freeze, so no score moves.

  Judge scores are never edited in place (§2.4). runs.csv holds the judge's
  raw output; corrections exist only in results/reference_audit.csv. Both
  readings stay inspectable, which is the only reason the primary and
  sensitivity specifications can both be computed from the released files.
"""

import argparse
import csv
import math
import statistics
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TASKS = ROOT / "tasks.csv"
RUNS = ROOT / "results" / "runs.csv"
AUDIT = ROOT / "results" / "reference_audit.csv"

# Report order throughout, best-performing first, control last.
ARMS = ["smartsbio-key", "biomni", "galaxy", "baseline"]
LABEL = {"smartsbio-key": "smarts.bio", "biomni": "Biomni",
         "galaxy": "Galaxy", "baseline": "untooled control"}
TOOL_ARMS = ["smartsbio-key", "biomni", "galaxy"]

# Category ids in tasks.csv, with the names the manuscript uses.
CATEGORY = {
    "lookup": "Database lookups", "blast": "Sequence & BLAST",
    "sequence": "Sequence analysis", "ngs": "NGS pipelines",
    "expression": "Expression analysis", "variant": "Variant interpretation",
    "structure": "Structure & annotation", "molecule": "Small molecules",
    "design": "Protein design (GPU)", "literature": "Literature & evidence",
    "everyday": "Everyday multi-step",
}


def num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def load():
    for p in (TASKS, RUNS, AUDIT):
        if not p.exists():
            sys.exit(f"missing {p.relative_to(ROOT)}")
    tasks = {t["task_id"]: t for t in csv.DictReader(open(TASKS, newline=""))}
    runs = list(csv.DictReader(open(RUNS, newline="")))
    audit = list(csv.DictReader(open(AUDIT, newline="")))
    return tasks, runs, audit


def voided_tasks(audit):
    """The defective tasks of §2.2, read from the audit layer rather than
    hard-coded: an `arm == "ALL"` row carrying a VOID verdict."""
    return {a["task_id"] for a in audit
            if a["arm"] == "ALL" and a["audit_accuracy"].strip() == "VOID"}


def score_corrections(audit):
    """(task_id, arm) -> [(audit_date, score)]. Digit verdicts only; see the
    module docstring for why INFRA, blank and note rows are excluded."""
    out = defaultdict(list)
    for a in audit:
        if a["task_id"] == "ALL" or a["arm"] == "ALL":
            continue
        v = a["audit_accuracy"].strip()
        if v.isdigit():
            out[(a["task_id"], a["arm"])].append((a["audit_date"], int(v)))
    return out


def build(runs, audit, *, apply_audit, drop_voided):
    """One score per (task, condition) cell under the chosen specification."""
    void = voided_tasks(audit) if drop_voided else set()
    corr = score_corrections(audit)
    cells, meta = {}, {}
    for r in runs:
        if r["rep"] != "1":                  # §2.5 repetition filter
            continue
        if r["task_id"] in void:
            continue
        key = (r["task_id"], r["arm"])
        # Intention-to-complete: a run that never returned is a zero, not a
        # missing value. Averaging over answered runs only would reward a
        # condition for the runs it failed to finish.
        acc = 0.0 if r["status"] == "timeout" else (num(r["accuracy_0_3"]) or 0.0)
        if apply_audit:
            started = r["started_at"][:10]
            for date, v in corr.get(key, []):
                if started <= date:          # §2.5 staleness rule
                    acc = float(v)
        cells[key] = acc
        meta[key] = r
    return cells, meta


def arm_table(cells):
    by = defaultdict(list)
    for (_, arm), v in cells.items():
        by[arm].append(v)
    rows = []
    for arm in ARMS:
        v = by.get(arm, [])
        if not v:
            continue
        ok = sum(1 for x in v if x >= 2)
        rows.append({
            "arm": LABEL[arm], "n": len(v),
            "mean_score": round(statistics.fmean(v), 2),
            # Sample SD, matching Table 1. §2.6: this describes dispersion
            # across the fixed task set and is not a standard error.
            "sd": round(statistics.stdev(v), 2) if len(v) > 1 else "",
            "perfect_3": sum(1 for x in v if x == 3),
            "zeros": sum(1 for x in v if x == 0),
            "acceptable": ok,
            "acceptable_pct": round(100 * ok / len(v), 1),
        })
    return rows


def category_table(cells, tasks):
    cell = defaultdict(list)
    cats = set()
    for (tid, arm), v in cells.items():
        c = tasks.get(tid, {}).get("category", "")
        cell[(c, arm)].append(v)
        cats.add(c)
    rows = []
    for c in sorted(cats, key=lambda x: CATEGORY.get(x, x)):
        row = {"category": CATEGORY.get(c, c),
               "n_tasks": len(cell.get((c, ARMS[0]), []))}
        pooled = []
        for arm in ARMS:
            v = cell.get((c, arm), [])
            row[LABEL[arm]] = round(statistics.fmean(v), 2) if v else ""
            if arm in TOOL_ARMS:
                pooled += v
        row["pooled tool-enabled"] = round(statistics.fmean(pooled), 2) if pooled else ""
        base = cell.get((c, "baseline"), [])
        row["gap"] = (round(statistics.fmean(pooled) - statistics.fmean(base), 2)
                      if pooled and base else "")
        rows.append(row)
    return rows


def leave_category_out(cells, tasks, category):
    """Re-run the arm table with one category removed (§3.1).

    Reported for `design` because those six tasks are the category most
    likely to favour the authoring condition's managed GPU catalogue, which
    makes the task-set composition objection of §5.1 sharpest there. Note
    that it is not a category the comparators could not enter: on des-001 all
    three tool-enabled conditions scored 3, and on des-006 two did.
    """
    kept = {k: v for k, v in cells.items()
            if tasks.get(k[0], {}).get("category") != category}
    return arm_table(kept)


def sign_test(k, n):
    """Exact two-sided binomial p for k successes in n discordant pairs."""
    if n == 0:
        return 1.0
    probs = [math.comb(n, i) for i in range(n + 1)]
    total = float(sum(probs))
    obs = probs[k] / total
    return min(1.0, sum(p / total for p in probs if p / total <= obs + 1e-12))


def head_to_head(cells):
    tasks_by_arm = defaultdict(dict)
    for (tid, arm), v in cells.items():
        tasks_by_arm[arm][tid] = v
    rows = []
    for i, a in enumerate(ARMS):
        for b in ARMS[i + 1:]:
            shared = set(tasks_by_arm[a]) & set(tasks_by_arm[b])
            win = sum(1 for t in shared if tasks_by_arm[a][t] > tasks_by_arm[b][t])
            loss = sum(1 for t in shared if tasks_by_arm[a][t] < tasks_by_arm[b][t])
            tie = len(shared) - win - loss
            rows.append({
                "pair": f"{LABEL[a]} vs {LABEL[b]}", "n_shared": len(shared),
                "win": win, "tie": tie, "loss": loss,
                "sign_test_p": round(sign_test(win, win + loss), 4),
            })
    # Holm-Bonferroni across the six pairwise tests (§2.6). All six are
    # reported, not only those involving the authoring condition.
    order = sorted(range(len(rows)), key=lambda i: rows[i]["sign_test_p"])
    running = 0.0
    for rank, i in enumerate(order):
        adj = min(1.0, max(running, (len(rows) - rank) * rows[i]["sign_test_p"]))
        running = adj
        rows[i]["holm_p"] = round(adj, 4)
    return rows


def deterministic_layer(cells, tasks, runs):
    """How much of the primary record rests on checks vs the judge (§2.4).

    These counts are published in the text, so they are computed here rather
    than by hand; a figure in the paper that this does not reproduce is wrong.
    """
    keys = set(cells)
    declaring = {t for t in {tid for tid, _ in keys}
                 if tasks[t]["check_type"] != "llm_rubric"}
    prim = [r for r in runs if (r["task_id"], r["arm"]) in keys and r["rep"] == "1"]
    on_decl = [r for r in prim if r["task_id"] in declaring]
    verdict = [r for r in on_decl if r["deterministic_pass"].strip()]
    timeout = [r for r in on_decl if r["status"] == "timeout"]
    unread = [r for r in on_decl
              if not r["deterministic_pass"].strip() and r["status"] != "timeout"]
    n_tasks = len({tid for tid, _ in keys})
    return [
        {"quantity": "tasks in specification", "n": n_tasks},
        {"quantity": "  declaring a checkable type", "n": len(declaring)},
        {"quantity": "  resting on the rubric only", "n": n_tasks - len(declaring)},
        {"quantity": "primary runs", "n": len(prim)},
        {"quantity": "  on declaring tasks", "n": len(on_decl)},
        {"quantity": "    check returned a verdict", "n": len(verdict)},
        {"quantity": "    timed out", "n": len(timeout)},
        {"quantity": "    answer unreadable by the check", "n": len(unread)},
        {"quantity": "  resting on the judge alone", "n": len(prim) - len(verdict)},
    ]


def rerun_substitution(cells, tasks):
    """§5.9: six smarts.bio cells were run again during the study. The audit
    layer recorded a judge score for each before that happened. Substituting
    those earlier scores bounds how much the re-running moved the result.

    The superseded runs themselves were not retained, so this is a bound
    reconstructed from the audit layer, not a measurement of the original
    record. Computed here rather than by hand because the paper quotes it.
    """
    path = ROOT / "results" / "superseded_scores.csv"
    if not path.exists():
        return []
    # score_before is on the same basis as the published cell it replaces: the
    # audited score where the audit returned a digit verdict, the judge score
    # otherwise. Mixing a raw judge score into an audit-corrected comparison
    # would compare two different specifications (day-004 is the case that
    # matters: the audit overturned the judge on its earlier run).
    pre = {}
    for r in csv.DictReader(open(path, newline="")):
        v = num(r["score_before"])
        if v is not None:
            pre[(r["task_id"], r["arm"])] = v
    alt = {k: pre.get(k, v) for k, v in cells.items()}
    rows = []
    for arm in ARMS:
        now = [v for (t, a), v in cells.items() if a == arm]
        sub = [v for (t, a), v in alt.items() if a == arm]
        rows.append({
            "arm": LABEL[arm],
            "published": round(statistics.fmean(now), 2),
            "earlier_scores": round(statistics.fmean(sub), 2),
            "delta": round(statistics.fmean(sub) - statistics.fmean(now), 2),
        })
    for label, c in (("published", cells), ("earlier_scores", alt)):
        pooled = [v for (t, a), v in c.items() if a in TOOL_ARMS]
        base = [v for (t, a), v in c.items() if a == "baseline"]
        rows.append({"arm": f"  execution gap ({label})",
                     "published": "", "earlier_scores": "",
                     "delta": round(statistics.fmean(pooled)
                                    - statistics.fmean(base), 2)})
    for other in ("biomni", "galaxy"):
        out = []
        for label, c in (("published", cells), ("earlier", alt)):
            shared = ({t for t, a in c if a == "smartsbio-key"}
                      & {t for t, a in c if a == other})
            win = sum(1 for t in shared if c[(t, "smartsbio-key")] > c[(t, other)])
            loss = sum(1 for t in shared if c[(t, "smartsbio-key")] < c[(t, other)])
            out.append(f"{win}-{loss} p={sign_test(win, win + loss):.4f}")
        rows.append({"arm": f"  smarts.bio vs {LABEL[other]}",
                     "published": out[0], "earlier_scores": out[1], "delta": ""})
    return rows


def audit_direction(audit):
    """Net effect of the audit layer, by condition (§5.1).

    Reported because the paper makes a claim about which condition the audit
    cost, and a claim of that kind should not be hand-computed.
    """
    net = defaultdict(float)
    up = defaultdict(int)
    down = defaultdict(int)
    for r in audit:
        j, a = num(r.get("judge_accuracy")), num(r.get("audit_accuracy"))
        if j is None or a is None:
            continue
        d = a - j
        net[r["arm"]] += d
        if d > 0:
            up[r["arm"]] += 1
        elif d < 0:
            down[r["arm"]] += 1
    return [{"arm": LABEL.get(k, k), "net_points": round(net[k], 1),
             "upgrades": up[k], "downgrades": down[k]}
            for k in sorted(net)]


def latency_and_cost(cells, meta, runs):
    """Median and p90 wall time (§3.5), and driver cost (§3.6).

    Latency is over runs that returned an answer, on whichever cells the
    current specification retains: a timed-out run has no wall time to report,
    and including its ceiling would measure the timeout rather than the
    platform.

    Driver cost is deliberately NOT restricted by the specification's task
    exclusions. It is the money actually spent operating each condition
    through this harness, and a run on a task later found defective still
    cost what it cost. The totals are therefore over every rep-1 record for
    the condition, which is why they do not change between specifications;
    only the acceptable-answer denominator does.

    `missing_cost_runs` counts records with no cost recorded at all, which is
    an undercount of the total rather than a zero. All of them are timed-out
    runs, but not every timed-out run is missing: one Galaxy run that was
    resumed twenty times before timing out did record its cost. §3.6 states
    this undercount rather than imputing it.
    """
    lat = defaultdict(list)
    ok = defaultdict(int)
    for key, score in cells.items():
        r = meta[key]
        if r["status"] != "timeout":
            w = num(r["wall_time_s"])
            if w is not None:
                lat[r["arm"]].append(w)
        if score >= 2:
            ok[r["arm"]] += 1

    spend = defaultdict(float)
    missing = defaultdict(int)
    timeouts = defaultdict(int)
    for r in runs:
        if r["rep"] != "1":
            continue
        arm = r["arm"]
        if r["status"] == "timeout":
            timeouts[arm] += 1
        c, p = num(r["claude_cost_usd"]), num(r["platform_cost_usd"])
        if c is None and p is None:
            missing[arm] += 1
            continue
        spend[arm] += (c or 0.0) + (p or 0.0)

    rows = []
    for arm in ARMS:
        w = sorted(lat.get(arm, []))
        if not w:
            continue
        idx = min(len(w) - 1, int(math.ceil(0.9 * len(w))) - 1)
        rows.append({
            "arm": LABEL[arm],
            "median_wall_s": round(statistics.median(w), 1),
            "p90_wall_s": round(w[idx], 1),
            "timeouts": timeouts[arm],
            "driver_cost_usd": round(spend[arm], 2),
            "missing_cost_runs": missing[arm],
            "acceptable": ok[arm],
            "cost_per_acceptable_usd": (round(spend[arm] / ok[arm], 2)
                                        if ok[arm] else ""),
        })
    return rows


def show(title, rows, key_order=None):
    print(f"\n{title}")
    if not rows:
        print("  (no rows)")
        return
    cols = key_order or list(rows[0].keys())
    width = {c: max(len(str(c)), *(len(str(r.get(c, ""))) for r in rows)) for c in cols}
    print("  " + "  ".join(str(c).ljust(width[c]) for c in cols))
    for r in rows:
        print("  " + "  ".join(str(r.get(c, "")).ljust(width[c]) for c in cols))


def write_csv(path, rows):
    if not rows:
        return
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"  wrote {path.relative_to(ROOT)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--csv", action="store_true",
                    help="also write results/spec_*.csv")
    args = ap.parse_args()

    tasks, runs, audit = load()
    void = voided_tasks(audit)
    print(f"tasks {len(tasks)} | run records {len(runs)} | audit rows {len(audit)}")
    print(f"defective tasks dropped from the primary specification: "
          f"{', '.join(sorted(void)) or 'none'}")

    specs = [
        ("PRIMARY (audit-corrected)", dict(apply_audit=True, drop_voided=True)),
        ("SENSITIVITY (raw)", dict(apply_audit=False, drop_voided=False)),
    ]
    for name, opts in specs:
        cells, meta = build(runs, audit, **opts)
        print("\n" + "=" * 72)
        print(f"{name} — n = {len(cells) // len(ARMS)} per condition")
        print("=" * 72)
        arms = arm_table(cells)
        show("Mean score and acceptable-answer rate (§3.1)", arms)

        pooled = [v for (_, a), v in cells.items() if a in TOOL_ARMS]
        base = [v for (_, a), v in cells.items() if a == "baseline"]
        print(f"\nExecution gap (§3.2): pooled tool-enabled "
              f"{statistics.fmean(pooled):.2f} vs control "
              f"{statistics.fmean(base):.2f} = {statistics.fmean(pooled) - statistics.fmean(base):.2f}")

        show("Excluding the protein-design category (§3.1 sensitivity)",
             leave_category_out(cells, tasks, "design"))

        cats = category_table(cells, tasks)
        show("Mean score by category (§3.3, Figure 2)", cats)
        show("Paired task-level outcomes (§3.7)", head_to_head(cells))
        show("Latency and driver cost (§3.5, §3.6)",
             latency_and_cost(cells, meta, runs))
        show("Deterministic layer vs judge (§2.4)",
             deterministic_layer(cells, tasks, runs))
        if opts["apply_audit"]:
            show("Net effect of the audit layer, by condition (§5.1)",
                 audit_direction(audit))
            show("Cells that were run again: earlier-score substitution (§5.9)",
                 rerun_substitution(cells, tasks))

        if args.csv:
            tag = "primary" if opts["apply_audit"] else "raw"
            write_csv(ROOT / "results" / f"spec_{tag}_by_arm.csv", arms)
            write_csv(ROOT / "results" / f"spec_{tag}_by_category.csv", cats)


if __name__ == "__main__":
    main()
