#!/usr/bin/env python3
"""
Run the platform benchmark: every task, against every arm, N times.

Each run is one `claude -p` subprocess restricted to a single platform's MCP
server. Nothing in this file reasons about biology; it sets up the arm, times
the call, parses the stream-json transcript and appends a row to results/runs.csv.

Long-running jobs are waited on: when a run hands an asynchronous job to a
platform, the same session is resumed with a fixed nudge every --poll-interval
seconds until it answers or the task's timeout_min is spent. See NUDGE below.

Usage:
    python3 -I run_benchmark.py --arms smartsbio,baseline --tasks look-001,seq-001
    python3 -I run_benchmark.py --dry-run          # print commands, run nothing
    python3 -I run_benchmark.py --repeats 3        # the full run
    python3 -I run_benchmark.py --no-resume-jobs   # one-shot, pre-2026-11-05

Resumable: rows already in results/runs.csv are skipped unless --force.
"""

import argparse
import csv
import json
import os
import random
import re
import subprocess
import threading
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
TASKS = ROOT / "tasks.csv"
RESULTS = ROOT / "results" / "runs.csv"
TRANSCRIPTS = ROOT / "transcripts"
SYSTEM_PROMPT = ROOT / "system_prompt.md"

ARMS = ["baseline", "smartsbio-key", "biomni", "galaxy"]

# Built-in tools are denied on every arm so that a platform cannot be rescued by
# Claude writing its own code, searching the web, or reading local files. This is
# the fairness guarantee of the whole benchmark.
# ToolSearch is deliberately NOT in this list. MCP tools are deferred and get
# loaded through it, so denying it leaves every platform arm with no tools while
# the run still looks healthy -- measured: with ToolSearch denied the smarts.bio
# arm answered a UniProt lookup from memory and reported "I have no
# database-query tools available in this session".
# Allowing it does not weaken isolation: --strict-mcp-config means only the
# arm's own server exists, and a tool named here cannot be loaded through
# ToolSearch either. Both were verified rather than assumed.
# What a run actually had is recorded per row in tools_available, taken from the
# `tools` array of the stream-json init event.
BUILTIN_TOOLS = [
    "Bash", "Read", "Write", "Edit", "NotebookEdit", "Glob", "Grep",
    "WebFetch", "WebSearch", "Task", "Agent", "TodoWrite", "BashOutput",
    "KillShell", "SlashCommand", "Skill", "ListMcpResourcesTool",
    "ReadMcpResourceTool", "Artifact", "Workflow", "TaskOutput",
    "TaskStop", "Monitor", "SendMessage", "ScheduleWakeup", "PushNotification",
    "RemoteTrigger", "DesignSync", "ReportFindings", "EnterWorktree",
    "ExitWorktree", "CronCreate", "CronDelete", "CronList",
]

# --- Long-running jobs: the resume loop -------------------------------------
# A platform that submits an asynchronous job returns a job id, not a result.
# One `claude -p` call therefore ends as soon as the model has reported that id,
# which is typically minutes before the job finishes -- so without this loop the
# row measures "did the job finish before the model stopped talking", not "did
# the platform solve the task". Measured on the 2026-11-05 re-run: rna-001's
# viennarna job completed 2 min, and expr-002's scanpy job 3 min, after their
# sessions had already exited; both scored 0 with correct results sitting in the
# workspace.
#
# The loop resumes the SAME session with a fixed, content-free nudge until the
# model gives an answer or the task's timeout_min budget is spent. It is applied
# identically to every arm; arms whose tools return results inline (Biomni's
# Python REPL) never trigger it, because there is no job to wait on.
NUDGE = (
    "Check on the job you started. If it has finished, give the final answer "
    "now, in full, ending with the '## ANSWER' section as instructed. If it is "
    "still running, reply with exactly one word: PENDING"
)

# A tool that hands back a job id rather than a result. Matched on the bare tool
# name, so it covers mcp__smartsbio__smarts_run_tool / smarts_run_pipeline and
# mcp__galaxy__run_tool, while mcp__biomni__run_python_repl -- which returns its
# output inline -- is correctly excluded.
SUBMIT_TOOL_RE = re.compile(r"(?:^|_)(?:run_tool|run_pipeline)$")

# The model does not always obey the one-word sentinel, so these phrases are a
# fallback: stopping the loop on a leg that says "still running" would record a
# non-answer as the platform's final answer.
PENDING_PHRASES = (
    "pending", "still running", "still in progress", "still queued",
    "currently running", "not yet complete", "not yet finished",
    "once it completes", "once it finishes", "check back",
)

FIELDS = [
    "run_id", "task_id", "arm", "rep", "started_at", "status", "wall_time_s",
    "wait_time_s", "resume_legs",
    "time_to_first_tool_s", "num_turns", "num_tool_calls", "tool_errors",
    "wrong_tool_name_calls", "claude_input_tokens", "claude_output_tokens",
    "claude_cost_usd", "platform_cost_usd", "deterministic_pass", "accuracy_0_3",
    "hallucination", "provenance", "artifacts_produced", "human_intervention_needed",
    "judge_notes", "model", "tools_available", "mcp_status",
    "platform_tool_calls", "session_id",
    "answer", "transcript_path",
]


def load_dotenv():
    """Load .env into the environment so ${VAR} in the MCP configs resolves.

    Values already in the environment win, so a shell export overrides the file.
    """
    env = ROOT / ".env"
    if not env.exists():
        return
    for line in env.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, _, v = line.partition("=")
        k, v = k.strip(), v.strip().strip('"').strip("'")
        if k and k not in os.environ:
            os.environ[k] = v


def load_tasks():
    with open(TASKS, newline="") as f:
        return list(csv.DictReader(f))


def load_done():
    """(task_id, arm, rep) tuples already recorded, so a run can be resumed.

    An infra_error is not a result -- the platform was never fairly asked -- so
    those rows are not counted as done and get another attempt on the next run.
    A timeout is counted: a platform that cannot answer inside its budget has
    told us something, and re-running it would only spend the budget again.
    """
    if not RESULTS.exists():
        return set()
    with open(RESULTS, newline="") as f:
        return {(r["task_id"], r["arm"], r["rep"]) for r in csv.DictReader(f)
                if r.get("status") != "infra_error"}


def append_row(row):
    """Write a row, replacing any existing row for the same task+arm+rep.

    Replacing rather than appending matters for --force: duplicate rows for one
    run_id would silently double-count that task in every aggregate.
    """
    RESULTS.parent.mkdir(parents=True, exist_ok=True)
    existing = []
    if RESULTS.exists():
        with open(RESULTS, newline="") as f:
            existing = [r for r in csv.DictReader(f)
                        if (r["task_id"], r["arm"], r["rep"]) !=
                           (row["task_id"], row["arm"], row["rep"])]
    with open(RESULTS, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=FIELDS, extrasaction="ignore")
        w.writeheader()
        w.writerows(existing + [row])


def server_name(arm):
    """The MCP server key inside an arm's config, which sets the tool prefix.

    This is not always the arm name: `smartsbio-key` is a different config file
    for the same server, whose tools are still called mcp__smartsbio__*.
    """
    cfg = json.loads((ROOT / "platforms" / f"{arm}.mcp.json").read_text())
    names = list(cfg.get("mcpServers", {}))
    return names[0] if names else ""


def build_cmd(task, arm, model, session_id=None, resume=None, skip_perms=False,
              prompt_override=None):
    cfg = ROOT / "platforms" / f"{arm}.mcp.json"
    if prompt_override is not None:
        prompt = prompt_override
    else:
        prompt = task["prompt"]
        if task["input_files"]:
            files = ", ".join(task["input_files"].split(";"))
            prompt += f"\n\n(The relevant files are already available to you: {files})"

    cmd = [
        "claude", "-p",
        "--model", model,
        "--append-system-prompt", SYSTEM_PROMPT.read_text(),
        "--mcp-config", str(cfg),
        "--strict-mcp-config",
        "--disallowedTools", ",".join(BUILTIN_TOOLS),
        "--output-format", "stream-json",
        "--verbose",
    ]
    # ToolSearch is allowed on every arm, including the baseline, because MCP
    # tools are deferred and loaded through it: denying it leaves an arm with no
    # tools at all while looking healthy. It cannot break arm isolation --
    # --strict-mcp-config means only this arm's server exists, and a denied
    # built-in cannot be loaded through it, both verified.
    srv = server_name(arm)
    allowed = ["ToolSearch"] + ([f"mcp__{srv}__*"] if srv else [])
    cmd += ["--allowedTools", ",".join(allowed)]
    if skip_perms:
        cmd += ["--dangerously-skip-permissions"]
    if resume:
        cmd += ["--resume", resume]
    elif session_id:
        cmd += ["--session-id", session_id]
    # The prompt is returned separately and fed on stdin rather than appended as
    # a positional argument: --disallowedTools, --allowedTools and --mcp-config
    # are all variadic (<tools...>) and will swallow a trailing positional.
    return cmd, prompt


def parse_stream(lines, arm):
    """Pull metrics out of a stream-json transcript."""
    srv = server_name(arm)
    m = {
        "num_tool_calls": 0, "tool_errors": 0, "wrong_tool_name_calls": 0,
        "num_turns": "", "claude_input_tokens": "", "claude_output_tokens": "",
        "claude_cost_usd": "", "session_id": "", "answer": "",
        "time_to_first_tool_s": "", "status": "failed",
        "model": "", "tools_available": "", "mcp_status": "",
        "platform_tool_calls": 0, "mcp_never_ready": False,
        # The driver model itself was blocked by Anthropic's safeguards, so
        # this leg never carried the platform's answer back.
        "driver_refusal": "",
        # True when this leg handed a job off to the platform rather than
        # computing an answer, which is what the resume loop waits on.
        "submitted_job": False,
    }
    first_tool_at = None
    t0 = None
    for raw in lines:
        raw = raw.strip()
        if not raw:
            continue
        try:
            ev = json.loads(raw)
        except json.JSONDecodeError:
            continue
        if t0 is None:
            t0 = time.time()
        et = ev.get("type")
        if et == "system" and ev.get("subtype") == "init":
            m["session_id"] = ev.get("session_id", "")
            # Record what this run could actually reach, so arm isolation is
            # auditable after the fact rather than merely asserted.
            m["model"] = ev.get("model", "")
            m["tools_available"] = "|".join(
                sorted(ev.get("tools", []) or []) +
                sorted(s.get("name", "") for s in ev.get("mcp_servers", []) or [])
            )
            # An unauthorised or unreachable server is an infrastructure problem.
            # Without this, every row for that arm would read as the platform
            # being unable to do the task.
            m["mcp_status"] = ";".join(
                f"{s.get('name')}={s.get('status')}"
                for s in ev.get("mcp_servers", []) or []
            )
        elif et == "assistant":
            for block in ev.get("message", {}).get("content", []) or []:
                if block.get("type") == "tool_use":
                    m["num_tool_calls"] += 1
                    if first_tool_at is None:
                        first_tool_at = time.time()
                    # A call to a tool that does not exist on this arm.
                    # ToolSearch is excluded: it is how deferred MCP tools get
                    # loaded, so it is a legitimate discovery step on every arm,
                    # not a mistake. Counting it inflated this metric by the
                    # number of discovery calls -- on look-006 every one of the
                    # four "wrong" calls was a ToolSearch on a run that scored 3.
                    name = block.get("name", "")
                    if srv and name.startswith(f"mcp__{srv}__"):
                        m["platform_tool_calls"] += 1
                    if (srv and name != "ToolSearch"
                            and not name.startswith(f"mcp__{srv}__")):
                        m["wrong_tool_name_calls"] += 1
                    if SUBMIT_TOOL_RE.search(name):
                        m["submitted_job"] = True
        elif et == "user":
            for block in ev.get("message", {}).get("content", []) or []:
                if block.get("type") != "tool_result":
                    continue
                if block.get("is_error"):
                    m["tool_errors"] += 1
                # The CLI says this in a non-error result when an MCP server
                # has not finished connecting. Harmless once, but a run that
                # ends still seeing it never reached the platform at all.
                if srv and "still connecting" in str(block.get("content", "")):
                    m["mcp_never_ready"] = True
        elif et == "system" and ev.get("subtype") == "model_refusal_no_fallback":
            m["driver_refusal"] = ev.get("api_refusal_category") or "unknown"
        elif et == "result":
            m["num_turns"] = ev.get("num_turns", "")
            m["claude_cost_usd"] = ev.get("total_cost_usd", "")
            usage = ev.get("usage", {}) or {}
            m["claude_input_tokens"] = usage.get("input_tokens", "")
            m["claude_output_tokens"] = usage.get("output_tokens", "")
            m["answer"] = ev.get("result", "") or ""
            m["session_id"] = ev.get("session_id", m["session_id"])
            if ev.get("is_error"):
                m["status"] = "failed"
            elif ev.get("subtype") == "success":
                m["status"] = "answered"
            else:
                m["status"] = ev.get("subtype", "failed")
    if first_tool_at and t0:
        m["time_to_first_tool_s"] = round(first_tool_at - t0, 1)
    return m


def extract_answer(text):
    """Return just the ## ANSWER section, which is what the judge scores."""
    if not text:
        return ""
    parts = re.split(r"^##\s*ANSWER\s*$", text, flags=re.MULTILINE | re.IGNORECASE)
    return (parts[-1] if len(parts) > 1 else text).strip()


def run_leg(cmd, prompt, timeout):
    """Run one `claude -p` call. Returns (lines, returncode, timed_out)."""
    # Stream stdout to the transcript as it arrives rather than capturing it at
    # the end: on a timeout the captured buffer comes back empty, which loses the
    # whole transcript and makes a genuine platform timeout look like a run that
    # never started.
    chunks: list[str] = []
    proc = subprocess.Popen(
        cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL, text=True, bufsize=1,
        cwd=ROOT, env={**os.environ},
    )

    def _pump():
        for line in proc.stdout:
            chunks.append(line)

    pump = threading.Thread(target=_pump, daemon=True)
    pump.start()
    try:
        proc.stdin.write(prompt)
        proc.stdin.close()
    except (BrokenPipeError, OSError):
        pass
    timed_out = False
    try:
        proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        proc.kill()
        proc.wait(timeout=30)
    pump.join(timeout=10)
    return "".join(chunks).splitlines(), proc.returncode, timed_out


# Polling backs off, because each resume leg is a fresh process that re-creates
# part of the prompt cache and so costs real money: measured on msa-001, four
# 30 s polls added $2.69 to a $1.25 run. A flat 30 s on a 20-minute budget would
# be ~40 polls. The tiers keep the timing error at or below the base interval for
# jobs that finish in the first few minutes -- where nearly all of them land --
# and let it grow only on long jobs, where two minutes on a forty-minute job is
# under 5%. Multiples of --poll-interval, so that flag still controls the cadence.
POLL_TIER_POLLS = 8        # first 8 polls at the base interval
POLL_TIER_ELAPSED_S = 720  # past 12 minutes, the slowest tier


def poll_delay(polls_done, elapsed_s, base):
    """Seconds to wait before the next poll."""
    if elapsed_s >= POLL_TIER_ELAPSED_S:
        return base * 4
    if polls_done >= POLL_TIER_POLLS:
        return base * 2
    return base


def answer_is_pending(text):
    """True when a leg reported that the job has not finished yet.

    Checked against the ## ANSWER section. The one-word sentinel is the
    requested form; the phrase fallback is length-limited so that a long real
    answer which happens to mention "pending" is not mistaken for a non-answer.
    """
    ans = extract_answer(text).strip()
    if not ans:
        return True
    if ans.upper().rstrip(".") == "PENDING":
        return True
    low = ans.lower()
    return len(ans) < 300 and any(p in low for p in PENDING_PHRASES)


def merge_legs(legs, answer_leg):
    """One row's metrics from all of its legs.

    Counters and cost are summed, because the whole loop is what the task cost.
    The init fields and time_to_first_tool_s come from the first leg, which is
    where the session and its tools were established.
    """
    out = dict(legs[0])
    for k in ("num_tool_calls", "tool_errors", "wrong_tool_name_calls"):
        out[k] = sum(l[k] for l in legs)
    for k in ("num_turns", "claude_input_tokens", "claude_output_tokens"):
        vals = [l[k] for l in legs if l[k] not in ("", None)]
        out[k] = sum(vals) if vals else ""
    costs = [l["claude_cost_usd"] for l in legs
             if l["claude_cost_usd"] not in ("", None)]
    out["claude_cost_usd"] = round(sum(costs), 6) if costs else ""
    out["answer"] = answer_leg["answer"]
    out["status"] = answer_leg["status"]
    return out


def run_one(task, arm, rep, args):
    run_id = f"{task['task_id']}__{arm}__r{rep}"
    tpath = TRANSCRIPTS / f"{run_id}.jsonl"
    TRANSCRIPTS.mkdir(parents=True, exist_ok=True)

    budget = int(task["timeout_min"] or 10) * 60
    cmd, prompt = build_cmd(task, arm, args.model,
                            skip_perms=args.skip_permissions)

    if args.dry_run:
        print(f"  [dry-run] {run_id}")
        shown = [c for c in cmd if "\n" not in c]
        print("    " + " ".join(shown) + "  <<< " + prompt[:60].replace("\n", " "))
        return None

    started = datetime.now(timezone.utc).isoformat()
    t0 = time.time()
    status_override = None
    waited = 0.0
    resume_legs = 0

    lines, rc, timed_out = run_leg(cmd, prompt, budget)
    all_lines = list(lines)
    legs = [parse_stream(lines, arm)]
    m = legs[0]
    if timed_out:
        status_override = "timeout"
    elif rc != 0 and not all_lines:
        status_override = "infra_error"

    # The resume loop. Entered only when the platform was handed an async job,
    # so fast tasks and arms that compute inline are untouched. The first nudge
    # goes out immediately rather than after a sleep: a job that finished during
    # the first leg should not be charged a poll interval it never waited.
    sid = m["session_id"]
    pending = m["submitted_job"] or answer_is_pending(m["answer"])
    while (args.resume_jobs and status_override is None and sid and pending
           and resume_legs < args.max_legs):
        remaining = budget - (time.time() - t0)
        if resume_legs:
            want = poll_delay(resume_legs, time.time() - t0, args.poll_interval)
            nap = min(want, max(0.0, remaining - 5))
            if nap <= 0:
                status_override = "timeout"
                break
            time.sleep(nap)
            waited += nap
            remaining = budget - (time.time() - t0)
        if remaining <= 5:
            status_override = "timeout"
            break
        rcmd, rprompt = build_cmd(task, arm, args.model, resume=sid,
                                  prompt_override=NUDGE,
                                  skip_perms=args.skip_permissions)
        # Recorded in the transcript so the published record shows exactly where
        # the harness intervened and what it said.
        all_lines.append(json.dumps({
            "type": "harness", "subtype": "resume", "leg": resume_legs + 1,
            "waited_s": round(waited, 1), "prompt": NUDGE,
        }))
        lines, rc, timed_out = run_leg(rcmd, rprompt, int(remaining))
        resume_legs += 1
        all_lines += lines
        if timed_out:
            status_override = "timeout"
            break
        if rc != 0 and not lines:
            # A failed resume leg is a harness/CLI problem, not a platform
            # answer. Keep whatever the last good leg said.
            break
        m = parse_stream(lines, arm)
        legs.append(m)
        sid = m["session_id"] or sid
        pending = answer_is_pending(m["answer"])
    else:
        if (args.resume_jobs and status_override is None and pending
                and resume_legs >= args.max_legs):
            status_override = "timeout"

    wall = round(time.time() - t0, 1)
    tpath.write_text("\n".join(all_lines))

    # The answer of record is the last leg that produced one. A later leg only
    # ever supersedes an earlier one when the earlier one had nothing to report,
    # and the nudge asks for the full answer in the required format so that a
    # confirming leg does not replace a good answer with a terser restatement.
    answered = [l for l in legs if not answer_is_pending(l["answer"])]
    m = merge_legs(legs, answered[-1] if answered else legs[-1])

    row = {
        "run_id": run_id, "task_id": task["task_id"], "arm": arm, "rep": str(rep),
        "started_at": started, "wall_time_s": wall,
        "wait_time_s": round(waited, 1), "resume_legs": resume_legs,
        "platform_cost_usd": "", "deterministic_pass": "", "accuracy_0_3": "",
        "hallucination": "", "provenance": "", "artifacts_produced": "",
        "human_intervention_needed": "", "judge_notes": "",
        "transcript_path": str(tpath.relative_to(ROOT)),
        **m,
    }
    row["answer"] = extract_answer(m["answer"])
    if status_override:
        row["status"] = status_override
    # "pending" is the normal init-event status for an HTTP server that finishes
    # connecting a moment later, so it is not an error. Only an authorisation or
    # connection failure is. Whether the tools were actually reachable is
    # recorded in tools_available and visible in the transcript.
    if any(bad in m["mcp_status"] for bad in ("needs-auth", "failed")):
        row["status"] = "infra_error"
    # Anthropic's safeguards blocked the driver model, not the platform. The
    # CLI reports this as is_error, which read as a platform failure: des-004
    # scored smarts.bio 'failed' after 48 clean tool calls and a submitted job.
    # It fires on the arm that gets FURTHEST, because only that arm's context
    # holds real design output, so left uncorrected it penalises competence.
    if m["driver_refusal"]:
        row["status"] = "infra_error"
        row["judge_notes"] = (
            f"driver model refused (category={m['driver_refusal']}) after "
            f"{m['platform_tool_calls']} platform tool call(s); "
            f"job submitted={m['submitted_job']}. Harness limit, not platform.")
    # A server stuck on "pending" is the same outage as one that reports
    # "failed", but it used to be read as the platform having no such tool:
    # des-006 scored smarts.bio for "no prediction run" when the CLI had told
    # it, 16 times, that the server was still connecting. Only an arm that both
    # saw that message and never reached its own server counts, so an honest
    # "this platform cannot do it" after real tool calls is untouched.
    arm_srv = server_name(arm)
    if arm_srv and m["mcp_never_ready"] and m["platform_tool_calls"] == 0:
        row["status"] = "infra_error"
        row["judge_notes"] = (
            f"{arm_srv} MCP never became callable: "
            f"{m['num_tool_calls']} tool "
            f"call(s), none reaching the platform. Never fairly asked.")
    # No init event at all means the CLI never started the session, so the
    # platform was never asked. Scoring that as a timeout would charge a
    # platform for a flake in the harness. Observed once in 12 runs, with an
    # empty transcript, while `uvx galaxy-mcp` itself starts in 2-4 s.
    if not all_lines and status_override != "timeout":
        row["status"] = "infra_error"
        row["judge_notes"] = "no output from the CLI; session never started"
    if row["status"] == "infra_error" and not getattr(args, "_is_retry", False):
        print(f"  {run_id:42} INFRA     retrying once")
        args._is_retry = True
        try:
            retried = run_one(task, arm, rep, args)
        finally:
            args._is_retry = False
        if retried is not None:
            return retried

    append_row(row)

    flag = {"answered": "ok", "timeout": "TIMEOUT",
            "infra_error": "INFRA"}.get(row["status"], row["status"].upper())
    polls = f"  polls={resume_legs} wait={waited:.0f}s" if resume_legs else ""
    print(f"  {run_id:42} {flag:9} {wall:7.1f}s  "
          f"tools={row['num_tool_calls']:<3} err={row['tool_errors']}{polls}")
    return row


def long_running_history():
    """True once at least one long_running row has completed, so estimates for a
    pipeline batch are grounded rather than extrapolated from lookups."""
    if not RESULTS.exists():
        return False
    long_ids = {t["task_id"] for t in load_tasks() if t["long_running"] == "y"}
    with RESULTS.open() as fh:
        return any(r["task_id"] in long_ids and r.get("status") == "answered"
                   for r in csv.DictReader(fh))


def fmt_dur(secs):
    """A duration a human can read at a glance."""
    secs = int(secs)
    if secs < 90:
        return f"{secs}s"
    if secs < 5400:
        return f"{secs // 60}m"
    return f"{secs / 3600:.1f}h"


def history_by_arm():
    """Mean cost and wall time per arm from completed runs, for estimates.

    Returns {arm: (mean_cost_usd, mean_wall_s)}. Empty until the first runs
    exist, in which case no estimate is printed rather than a made-up one.
    """
    if not RESULTS.exists():
        return {}
    acc = {}
    with RESULTS.open() as fh:
        for r in csv.DictReader(fh):
            if r.get("status") != "answered":
                continue
            try:
                c = float(r.get("claude_cost_usd") or 0)
                w = float(r.get("wall_time_s") or 0)
            except ValueError:
                continue
            cs, ws, n = acc.get(r["arm"], (0.0, 0.0, 0))
            acc[r["arm"]] = (cs + c, ws + w, n + 1)
    return {a: (cs / n, ws / n) for a, (cs, ws, n) in acc.items() if n}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--arms", default=",".join(ARMS))
    p.add_argument("--tasks", default="", help="comma-separated task_ids")
    p.add_argument("--categories", default="", help="comma-separated categories")
    p.add_argument("--repeats", type=int, default=1)
    # The installed CLI (2.1.220) rejects full model ids such as claude-opus-5-5
    # with "version 2.1.280 or newer is required", so the alias is the default.
    # The id the alias actually resolved to is recorded per run in the `model`
    # column, which is what reproducibility needs.
    p.add_argument("--model", default="opus")
    p.add_argument("--dry-run", action="store_true")
    p.add_argument("--force", action="store_true", help="re-run completed rows")
    p.add_argument("--skip-long", action="store_true",
                   help="skip long_running rows (fast smoke pass)")
    p.add_argument("--only-long", action="store_true")
    p.add_argument("--skip-permissions", action="store_true",
                   help="pass --dangerously-skip-permissions to the CLI")
    # 30 s rather than 60: the error in a recorded wall time is roughly the poll
    # interval, since the job finishes somewhere inside the sleep and is only
    # noticed at the next wake. Going much below 30 s stops helping, because each
    # resume leg itself costs 5-15 s of model latency, which then dominates.
    p.add_argument("--poll-interval", type=float, default=30.0,
                   help="base seconds between resume polls (default 30); backs "
                        f"off to 2x after {POLL_TIER_POLLS} polls and 4x past "
                        f"{POLL_TIER_ELAPSED_S // 60} minutes")
    p.add_argument("--max-legs", type=int, default=60,
                   help="safety cap on resume polls per run (default 60)")
    p.add_argument("--no-resume-jobs", dest="resume_jobs", action="store_false",
                   help="disable the resume loop (one-shot runs, pre-2026-11-05 "
                        "behaviour)")
    p.add_argument("--seed", type=int, default=20261102)
    args = p.parse_args()
    load_dotenv()

    arms = [a for a in args.arms.split(",") if a]
    for a in arms:
        if not (ROOT / "platforms" / f"{a}.mcp.json").exists():
            sys.exit(f"no MCP config for arm '{a}'")

    # Fail fast on an unauthorised server rather than after N wasted runs.
    if not args.dry_run:
        for a in arms:
            if not server_name(a):
                continue
            probe = subprocess.run(
                ["claude", "-p", "--model", args.model, "--strict-mcp-config",
                 "--mcp-config", str(ROOT / "platforms" / f"{a}.mcp.json"),
                 "--output-format", "stream-json", "--verbose"],
                input="hi", capture_output=True, text=True, timeout=300, cwd=ROOT,
            )
            first = (probe.stdout or "{}").splitlines()[:1]
            try:
                srv = json.loads(first[0]).get("mcp_servers", []) if first else []
            except (json.JSONDecodeError, IndexError):
                srv = []
            # "pending" is not a failure: an HTTP server finishes connecting
            # after the init event is emitted, so a working arm reports pending
            # here and serves its tools moments later. Only an authorisation or
            # connection failure is fatal. Whether the tools were really there
            # is settled per run by the tools_available and mcp_status columns.
            bad = [s for s in srv
                   if s.get("status") not in ("connected", "pending")]
            if bad:
                detail = ", ".join(f"{s.get('name')}={s.get('status')}" for s in bad)
                sys.exit(
                    f"\narm '{a}' is not usable: {detail}\n\n"
                    f"If the status is needs-auth, authorise it once in an "
                    f"interactive session:\n"
                    f"    claude mcp add --transport http {server_name(a)} "
                    f"{json.loads((ROOT / 'platforms' / f'{a}.mcp.json').read_text())['mcpServers'][server_name(a)]['url']}\n"
                    f"    claude      # then run /mcp and sign in\n\n"
                    f"Or use an API key instead: see .env.example and pass "
                    f"--arms {a}-key (if that config exists).\n"
                )
            seen = ", ".join(f"{s.get('name')}={s.get('status')}" for s in srv)
            print(f"arm '{a}': ok ({seen or 'no server'})")

    tasks = load_tasks()
    if args.tasks:
        want = set(args.tasks.split(","))
        tasks = [t for t in tasks if t["task_id"] in want]
    if args.categories:
        want = set(args.categories.split(","))
        tasks = [t for t in tasks if t["category"] in want]
    if args.skip_long:
        tasks = [t for t in tasks if t["long_running"] != "y"]
    if args.only_long:
        tasks = [t for t in tasks if t["long_running"] == "y"]
    if not tasks:
        sys.exit("no tasks selected")

    # Cost gate: GPU and pipeline rows bill real money.
    costed = [t for t in tasks if t["est_cost"].strip() not in ("", "0")]
    if costed and not args.dry_run:
        print(f"\n{len(costed)} selected task(s) have a non-zero compute cost, "
              f"across {len(arms)} arm(s) x {args.repeats} repeat(s):")
        for t in costed:
            print(f"  {t['task_id']:10} {t['est_cost']:>12}  {t['prompt'][:52]}")
        print("\nA compute-quota failure arrives as a tool result inside a healthy "
              "session, so it is recorded 'answered' and scored as a capability "
              "loss. Check the quota before proceeding (21.4/100 hrs on "
              "2026-11-05).")
        if input("Proceed? [y/N] ").strip().lower() != "y":
            sys.exit("aborted at cost gate")

    done = set() if args.force else load_done()
    rng = random.Random(args.seed)

    planned = []
    for rep in range(1, args.repeats + 1):
        for t in tasks:
            # GPU/pipeline rows run once regardless of --repeats.
            if rep > 1 and t["est_cost"].strip() not in ("", "0"):
                continue
            order = arms[:]
            rng.shuffle(order)          # spreads outages evenly across arms
            for arm in order:
                if (t["task_id"], arm, str(rep)) in done:
                    continue
                planned.append((t, arm, rep))

    # What this batch will cost and how long it will take, from whatever this
    # machine has measured before -- stated up front rather than inferred after.
    hist = history_by_arm()
    if hist:
        # Fall back to the mean across measured arms, so an arm with no history
        # is estimated rather than silently counted as costing nothing.
        fb = (sum(c for c, _ in hist.values()) / len(hist),
              sum(w for _, w in hist.values()) / len(hist))
    else:
        fb = (0.0, 0.0)
    est_cost = sum(hist.get(a, fb)[0] for _, a, _ in planned)
    est_secs = sum(hist.get(a, fb)[1] for _, a, _ in planned)
    print(f"\n{len(planned)} run(s) to do "
          f"({len(tasks)} tasks x {len(arms)} arms x {args.repeats} reps, "
          f"{len(done)} already done)")
    if est_cost:
        print(f"estimated: ${est_cost:.2f} driver cost, {fmt_dur(est_secs)} "
              f"wall (from previous runs on this machine)")
        # A pipeline row waits on a job for minutes to hours. Estimating it from
        # lookup-shaped history is worse than not estimating it, so say so.
        if any(t_["long_running"] == "y" for t_, _, _ in planned) and \
                not long_running_history():
            print("  NOTE: this batch contains long_running rows and there is no "
                  "history for them yet. The figures above are extrapolated from "
                  "fast tasks and will understate both cost and wall time, "
                  "probably by a lot. Track runs completed, not percent.")
    print()

    t_batch = time.time()
    for i, (t, arm, rep) in enumerate(planned, 1):
        pct = 100.0 * (i - 1) / len(planned)
        if i > 1:
            per = (time.time() - t_batch) / (i - 1)
            eta = f"  eta {fmt_dur(per * (len(planned) - i + 1))}"
        else:
            eta = f"  eta ~{fmt_dur(est_secs)}" if est_secs else ""
        print(f"[{i}/{len(planned)} {pct:3.0f}%{eta}] "
              f"{t['task_id']} / {arm} / r{rep}")
        try:
            run_one(t, arm, rep, args)
        except KeyboardInterrupt:
            sys.exit("\ninterrupted — rerun to resume")

    if not args.dry_run:
        print(f"\ndone. results -> {RESULTS.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
