# Bioinformatics agent platform benchmark

The same model, the same prompts, four different platforms. Each run may call
exactly one platform's MCP server — or none at all — and is scored for accuracy,
reliability, latency and cost.

**Disclosure: this benchmark was built by [smarts.bio](https://smarts.bio), which
is one of the platforms under test.** The task set is frozen and committed before
the first run, the prompts and transcripts are published, the judge is blinded to
which platform produced an answer, and the categories where smarts.bio loses are
reported alongside the ones where it wins. Read the method and check the
transcripts rather than taking the summary table on trust.

The full write-up, with figures and the complete method, is
[`smartsbio-benchmark-paper.pdf`](smartsbio-benchmark-paper.pdf)
(source: `paper-print.html`).

## What was found

46 tasks per arm, one run each, scored 0–3. This is the audit-corrected primary
specification; `analyze.py` prints it and the raw 48-task sensitivity analysis
side by side, and reproduces every number below.

| Arm | Mean score | Acceptable (≥2) | Perfect (3) | Zeros |
|---|---|---|---|---|
| smarts.bio | **2.70** | 93.5% | 36 | 1 |
| Biomni | 2.30 | 78.3% | 30 | 6 |
| Galaxy | 2.00 | 67.4% | 23 | 8 |
| Untooled control | 1.00 | 28.3% | 6 | 19 |

Paired by task, smarts.bio beats Biomni on 14 tasks and loses 3 (exact binomial
*p* = 0.013), and beats Galaxy on 21 and loses 4 (*p* = 0.0009). Biomni versus
Galaxy is not significant (*p* = 0.13).

The largest effect is not between platforms but between having tools and not:
pooled tool-enabled 2.33 against 1.00 for the same model with no MCP server at
all. On the NGS and expression categories the control scores 0.0 and 0.33.

**Where smarts.bio loses.** Biomni is better on database lookups (3.00 against
2.62) and, with Galaxy, on small molecules (2.67 against 2.33). Galaxy is better
on NGS pipelines (3.00 against 2.83). smarts.bio is also the most expensive per
acceptable answer, $4.08 against Biomni's $1.46.

Two of the 48 authored tasks are voided for defects found in the audit, which is
why the primary specification is 46; `TASK-CHANGES.md` records each one.

## What is compared

Four arms were run, one run per task per arm:

| Arm | What it is | Access |
|---|---|---|
| `smartsbio-key` | smarts.bio hosted MCP | `sk_live_` bearer token from `$SMARTSBIO_API_KEY` |
| `biomni` | Biomni's tool library, in Docker | E1 image at a pinned commit, served over stdio MCP |
| `galaxy` | Galaxy via the official `galaxy-mcp` | `uvx galaxy-mcp` against usegalaxy.org |
| `baseline` | The model alone | No MCP server, no built-in tools |

A fifth config, `platforms/smartsbio.mcp.json`, reaches the same server over
OAuth instead of an API key. It is kept because it is the path a human user
takes, but it produced **no rows in the published run** — the key arm was used
throughout, so that authentication is scripted rather than interactive. Check
`results/runs.csv`: the only `arm` values are the four above.

### The Biomni arm

Biomni is benchmarked from its open-source package (`snap-stanford/Biomni`,
Apache-2.0) rather than through any hosted service, so that the thing being
measured is Biomni itself and the version under test is an exact, citable one.

It ships both an agent and a tool library. `platforms/biomni_server.py` serves
the **tool library** over stdio MCP, so the same Claude model that drives every
other arm does the reasoning here too and the comparison stays a comparison of
tool layers. Biomni's own LLM is never called.

#### How Biomni was installed

Recorded in full, because how a platform is installed decides what it scores.
Biomni runs in a container: nothing is installed on the host.

```bash
docker build --platform linux/amd64 -f platforms/biomni.Dockerfile \
    -t biomni-bench:e1 platforms/        # ~37 GB image
./platforms/fetch_biomni_data.sh         # the 76-file data lake, ~14 GB
```

Set `BIOMNI_DATA_PATH` and `BIOMNI_IMAGE` in `.env`. To reclaim the space
afterwards, remove the image with `docker rmi biomni-bench:e1` and delete the
data lake directory that `BIOMNI_DATA_PATH` points at.

**Why a container, and why `linux/amd64`.** The first attempt installed E1 on
the macOS host and failed: bioconda has no osx-arm64 build for the perl chain
that `blast` depends on, conda solves all-or-nothing, and so two unsatisfiable
packages took down all 66 in `bio_env.yml` — biopython, scanpy, gget, gseapy,
rdkit, samtools, bedtools and BLAST itself. `setup.sh` treats that step as
optional and carried on regardless. Upstream states E1 is tested on Ubuntu
22.04, and on Ubuntu 22.04 the same solve succeeds, so the container is not a
convenience: it is the difference between measuring Biomni and measuring the
host. Apple Silicon runs it through Rosetta rather than qemu, which is fast
enough that the whole build takes about an hour.

**Five fixes were needed to make Biomni's own code run in Biomni's own
environment.** They are listed here because each one, left alone, would have
cost Biomni marks it did not deserve. None of them adds capability.

| Problem at commit `400c1f3` | What was done |
|---|---|
| `tool/genomics.py` imports `esm` and `pybiomart` at module level; E1 installs neither, so the module cannot import and every genomics tool vanishes | `pip install fair-esm pybiomart` |
| E1 installs `mcp 2.3.0`, which renamed `FastMCP` to `MCPServer`, but `create_mcp_server()` imports `mcp.server.fastmcp.FastMCP` — so the MCP server cannot start at all | pin `mcp<2` (resolves to 1.30.0) |
| `setup.sh` step 6 aborts during an image build and is swallowed as optional: `setup_path.sh` is written but `$BIOMNI_TOOLS_DIR` never appears, leaving FastTree, iqtree2, muscle, plink2 and HOMER absent — and IQ-TREE backs the alignment-and-tree task | run `install_cli_tools.sh` directly |
| `A1.__init__` downloads the data lake one file at a time with no timeout. `data_lake.zip` returns 403, so every install takes the per-file path, and a dropped connection leaves it spinning: observed alive and burning CPU for an hour after 2 of 76 files, socket in `CLOSE_WAIT` | `fetch_biomni_data.sh` — parallel, resumable, timeout-bounded, size-checked against `Content-Length` |
| `A1` loads an embedding-based tool retriever on init, which exceeded Docker Desktop's default 7.75 GB memory and got the container OOM-killed | `use_tool_retriever=False` |

The retriever is worth a note. It exists so Biomni's *own* agent can narrow a
large tool set before choosing; here Claude does the choosing and
`create_mcp_server()` reads the registry directly, so it is never consulted.
Switching it off removes an unused component, not a tool.

That Docker Desktop's out-of-box memory limit is too small for Biomni as
shipped is itself a finding about adoption cost, distinct from whether the
tools work.

**The resulting arm**, verified at build time and on first start:

```
11/11 tool modules import
blastn blastp makeblastdb samtools bedtools Rscript macs2 fastqc   (conda)
iqtree2 FastTree muscle plink2                                     (E1 CLI tools)
144 tools served, ~46K tokens of schema, 79 s from cold container
```

The version under test is the commit in `/opt/biomni-commit.txt` inside the
image; cite it with the results.

Which tool modules are served is listed explicitly at the top of
`biomni_server.py`. Serving all of them would overflow the context window with
tool definitions, which is a limit of the client rather than of Biomni, so the
selection is recorded in the file instead of being silently truncated — and it
covers every category in the task set.

`create_mcp_server()` catches a failed module import and merely warns, which
would serve a quietly smaller catalogue. The server therefore imports every
module itself first and refuses to start if any one fails, so a missing
dependency is always a setup error and never a benchmark result.

`BIOMNI_COMMERCIAL_MODE=1` restricts it to tools whose data sources are licensed
for commercial use. Worth setting when comparing against a platform that gates
those sources; the default serves Biomni's full capability.

#### Biomni's tools are of two kinds, and one needs an API key

Of the 144 tools served, **35 take a natural-language `prompt`** and have
Biomni's own model turn it into an API query — `query_uniprot`, `query_clinvar`,
`query_dbsnp`, `query_ensembl`, `query_kegg`, `query_stringdb`, `query_geo`,
`query_opentarget` and the rest of the database lookups. The other **109 take
typed arguments** and call their APIs directly.

So the arm needs `ANTHROPIC_API_KEY` in `.env`. Without it every LLM-mediated
tool returns `Error querying LLM: Could not resolve authentication method`,
which would zero the whole `lookup` category — measuring the configuration
rather than the platform. The two paths, verified:

```
query_uniprot  (prompt → LLM → REST)   ok, success: true, 15 s cold
query_alphafold(uniprot_id="P04637")   ok, AF-P04637-F1, 7 s cold
```

Two consequences for reading the results. Biomni's lookup latency includes its
own model's round trip, so it is not comparable like-for-like with a typed call
on another platform. And an LLM mis-translation of a prompt surfaces as a
*wrong answer*, not an error — a failure mode the typed tools do not have.

Biomni's default model is `claude-sonnet-4-5`, which is **deprecated and retires
30 November 2026**: after that date these 35 tools stop working on the default
configuration. The default is what a user gets, so it is what was benchmarked;
the run predates the retirement. Cost is a few dollars at $3/$15 per MTok.

### Authorising the smarts.bio arm

The hosted server uses OAuth, so no API key is needed. Authorise it once in an
interactive session and the CLI reuses the stored token, keyed by server URL, for
every later `claude -p` subprocess:

```bash
claude mcp add --transport http smartsbio https://mcp.smarts.bio/mcp
claude        # then run /mcp and sign in
```

A claude.ai connector is *not* enough: connectors live in the claude.ai session
layer and a spawned CLI subprocess does not inherit them, so the server would
report `needs-auth`.

`run_benchmark.py` probes every arm before the first run and aborts with
instructions if a server is not `connected`. Any run whose server reports
`needs-auth`, `pending` or `failed` is recorded as `infra_error` and excluded
from the scores, so an authorisation problem can never be published as a
platform limitation.

For a run intended for publication, prefer the `smartsbio-key` arm. An outside
reader can reproduce a documented key; they cannot reproduce an OAuth token
sitting in your keychain.

The `baseline` arm is the control that matters most. If a platform cannot beat a
plain assistant with no tools on a given category, that category is not a reason
to adopt it.

## The fairness rule

Every run is restricted to one platform's MCP tools. All built-in tools — `Bash`,
`Read`, `Write`, `WebSearch`, `WebFetch`, `Task`, `Workflow` and the rest — are
denied on every arm, including the baseline. Without that, the model could
install BLAST locally, `curl` an API, or search the web, and a weak platform
would be rescued by the agent driving it.

`ToolSearch` is the one exception, and it is allowed on every arm. MCP tools are
deferred in this CLI and are loaded through it, so denying it leaves a platform
arm with no tools at all while the run still looks healthy — measured: with
`ToolSearch` denied, the smarts.bio arm answered a UniProt lookup from memory and
reported that it had no database tools. It cannot break arm isolation, because
`--strict-mcp-config` means only that arm's server exists in the session and a
denied built-in cannot be loaded through it either. Both were verified rather
than assumed.

This is verified rather than asserted: every run records the `tools` array from
the CLI's own `init` event in the `tools_available` column, so you can confirm
after the fact what each run could actually reach.

## Long-running jobs

A platform that submits an asynchronous job returns a job id, not a result. One
`claude -p` call therefore ends as soon as the model has reported that id —
often minutes before the job finishes. Measured on 2026-11-05: an RNA-folding
job completed two minutes, and a single-cell clustering job three minutes, after
their sessions had already exited. Both produced correct, complete results that
no one was listening for, and both scored 0.

So when a run hands off an async job, the harness resumes **the same session**
with one fixed, content-free nudge:

> Check on the job you started. If it has finished, give the final answer now,
> in full, ending with the '## ANSWER' section as instructed. If it is still
> running, reply with exactly one word: PENDING

It repeats every 30 seconds (`--poll-interval`) until the model answers or the
task's `timeout_min` budget is spent. Rules:

- **Identical on every arm.** The trigger is architectural, not per-vendor: a
  tool that hands back a job id (`smarts_run_tool`, `smarts_run_pipeline`,
  Galaxy's `run_tool`). Arms whose tools return results inline — Biomni's Python
  REPL — never enter the loop, because there is nothing to wait on. The baseline
  never enters it either.
- **The first nudge is immediate**, with no sleep, so a job that finished during
  the first leg is not charged a poll interval it never waited.
- **Wall time includes the waiting**, since the user is waiting too. The error in
  a recorded wall time is roughly the poll interval — the job finishes somewhere
  inside the sleep and is noticed at the next wake — which is why the base
  interval is 30 s and not 60 s. Below about 30 s it stops helping: each resume
  leg itself costs 5–15 s of model latency, which then dominates.
- **Polling backs off**: 30 s for the first 8 polls, then 60 s, then 120 s past
  12 minutes. Each resume leg is a fresh process that re-creates part of the
  prompt cache, so polls cost real money — measured on `msa-001`, four 30 s polls
  added $2.69 to a $1.25 run. A flat 30 s on a 20-minute budget would be ~31
  polls; the tiers make it 19. The timing error stays at or below 30 s for jobs
  finishing in the first few minutes, where nearly all of them land, and grows
  only on long jobs, where two minutes on a forty-minute job is under 5%.
  `--poll-interval` sets the base and the tiers scale with it.
- **Cost is summed across legs**, because polling re-sends a growing context and
  that is real money. This raises the measured driver cost of async arms. The
  `resume_legs` and `wait_time_s` columns record how much of a row was polling.
- **The answer of record** is the last leg that produced one; a confirming leg
  only ever supersedes a leg that had nothing to report.
- `--no-resume-jobs` restores the old one-shot behaviour, for reproducing runs
  recorded before 2026-11-05.

Without this loop, the measurement is "did the job finish before the model
stopped talking", which systematically penalises asynchronous architectures and
costs nothing for platforms that block in-process. That is a property of the
harness, not of the platforms.

## Running it

```bash
cp .env.example .env         # then fill in the keys
./fetch_fixtures.sh          # pulls the one large fixture

# smoke test: one task, one arm, no credentials needed
python3 -I -u run_benchmark.py --arms baseline --tasks look-001

# a batch: run, judge, analyse, read it, then start the next one
python3 -I -u run_benchmark.py --categories lookup,literature --repeats 3 \
    | tee -a results/batch-A.log
python3 -I -u judge.py
python3 -I analyze.py

# the long pipeline and GPU rows (prompts for cost confirmation first)
python3 -I -u run_benchmark.py --only-long

python3 -I judge.py --audit-sample 20    # human audit worklist
```

`-u` matters: without it Python buffers stdout and a redirected run shows
nothing at all until it finishes, which is indistinguishable from a hang.

### Watching a batch

Each run prints its position, percent and a rolling ETA, and the batch prints
its estimated cost and duration up front, derived from whatever this machine has
already measured rather than from a guess:

```
86 run(s) to do (8 tasks x 4 arms x 3 reps, 15 already done)
estimated: $35.12 driver cost, 68m wall (from previous runs on this machine)

[14/86  15%  eta 59m] look-003 / smartsbio-key / r2
  look-003__smartsbio-key__r2   ok   155.3s  tools=14 err=0
```

From any other terminal, including after the one that started the batch is
closed, `results/runs.csv` is the authoritative record of what has completed,
because runs are resumable: re-issuing the same command skips the rows already
present and picks up the rest.

For the pipeline and GPU batches, **read runs completed rather than percent**: a
single row can take hours while a lookup takes ten seconds, so a percentage
jumps in steps and means little. The runner says so itself when a batch contains
long-running rows it has no history for.

The harness itself needs no third-party Python packages — a stock Python 3.9+ and
the `claude` CLI. The arms bring their own, each in its own environment so none
of it reaches the harness or the system Python: `uv` for Galaxy
(`brew install uv`), and for Biomni a Docker image, built as described above.

Runs are resumable: a `(task, arm, repeat)` already in `results/runs.csv` is
skipped unless you pass `--force`. Arm order is reshuffled per repeat from a
fixed seed, so a platform outage or a slow afternoon spreads evenly across arms
instead of landing on whichever arm happens to run first.

## The task set

48 tasks in `tasks.csv`, weighted toward what smarts.bio users actually run —
database lookups and BLAST dominate real usage — plus the categories where the
platforms genuinely diverge.

| Category | n | |
|---|---|---|
| `lookup` | 8 | gene, protein, variant, expression, interaction, bioactivity records |
| `sequence` | 6 | GC content, ORFs, translation, alignment, RNA folding |
| `blast` | 3 | protein identification, primer placement, organism identification |
| `ngs` | 6 | QC, alignment, germline and somatic calling, RNA-seq, ChIP-seq |
| `expression` | 3 | differential expression, single-cell clustering, enrichment |
| `variant` | 3 | annotation, clinical interpretation, prioritisation |
| `structure` | 4 | secondary structure, pockets, predicted models, metadata |
| `design` | 6 | affinity prediction, de novo binders, nanobodies, enzyme redesign, stability |
| `molecule` | 3 | drug-likeness, analogue generation, bioactivity |
| `literature` | 2 | papers, clinical trials |
| `everyday` | 4 | the questions a bench biologist actually asks |

`day-004` (CRISPR guide design) is an **expected-loss task**, included
deliberately: smarts.bio has no guide-design capability, and the correct answer
is to say so. Full credit goes to the arm that admits the gap; inventing guide
sequences with fabricated off-target scores scores zero and is flagged as a
hallucination.

### Ground truth

Where a task has a checkable answer, it is computed rather than recalled. The
synthetic fixtures carry known answers — 3 planted somatic variants, 2 germline
variants that must *not* be called somatic, 6 ChIP enrichment regions, and a
count matrix with a known 4x fold change — so those rows are scored
deterministically. Ground truth asserted from memory is the single easiest way to
poison a benchmark: during construction, one task's reference answer was wrong
(an enzyme misidentified as a lyase when BLAST showed it is a ligase) and was
only caught by running the search. Every ground-truth value that could be
verified, was.

### Operationalising the frozen ground truth (`checks.json`)

`tasks.csv` was pre-registered and frozen before the first run, and its
`ground_truth` column is prose written for a human reader and for the rubric
judge — `1EMA is NOT wild type; it is the S65T point mutant of GFP`. Prose
cannot be string-matched, which is why the 19 `contains_all` rows initially
carried no deterministic verdict at all and fell through to the judge.

`checks.json` turns that prose into checkable items **without editing the frozen
task set**, so `tasks.csv` keeps the commit hash cited in the article. Each task
maps to a list of requirements; each requirement is a list of synonyms, and any
one of them satisfies it (`F508del`, `delta F508` and `Phe508del` are the same
finding). Coordinates are normalised before matching, so `chr1:14,000`,
`chr1:14000` and `1 14000` compare equal.

Two properties keep this honest, and both matter more than the code:

- **Every item was derived from its own task's frozen `ground_truth` or
  `rubric` text, never from reading an answer.** The checks were written from
  the task definition, not reverse-engineered from the responses they grade.
- **A `PASS` is a recall result, not a correctness result.** It says the answer
  named the things a correct answer must name. A `FAIL` is authoritative — an
  answer that omits a required item cannot be fully correct — so `judge.py`
  treats the two directions asymmetrically: `FAIL` against a rubric score of 2+
  is flagged as a disagreement for the human audit, while `PASS` against a low
  score is not, because the judge can see what tokens cannot.

Where presence genuinely cannot capture the claim, the file records that rather
than pretending. `partial` marks a row whose ground truth makes a claim the
check cannot see — `expr-001` checks that all six genes are named but direction
(up vs down) stays with the rubric, and `ngs-004` checks the three somatic
positions but not the harder negative requirement, since a correct answer
legitimately mentions the two germline coordinates while labelling them
germline. `rubric_only` marks the two rows with no deterministic check at all:
`ngs-006`, because planted peak boundaries are approximate and exact matching
would fail a correct caller, and `blast-001`, which is void pending a decision
because its frozen prompt contains no sequence to search.

Re-running the checks over already-scored rows costs nothing and calls no model:

```bash
python3 -I judge.py --deterministic-only
```

This brought deterministic coverage from 13 rows to 69 of 147. Read the
resulting pass rates with the recall caveat above firmly in mind: all three
platform arms score 100% on it and the Claude-alone baseline 62%, which says
the baseline omits required identifiers, **not** that the platforms are 100%
correct.

### Pre-run truth derivation for the NGS rows

`ngs-003` was scored against variants that were not in its input, which is
recorded under Known limitations. The remaining NGS rows draw their ground
truth from the same fixtures, so each was re-derived from the fixtures
themselves *before* being run rather than audited afterwards. No aligner is
involved and none is installed: the reads are error-free 100 bp cuts, so
locating one by seed-and-extend is exact rather than heuristic.

| script | row | outcome |
|---|---|---|
| `audit/derive_ngs004_truth.py` | `ngs-004` | **truth confirmed.** All three somatic SNVs at VAF 1.00 in the tumour and 0.00 in the normal; both germline SNVs at VAF 1.00 in *both*, so the "must not call these as somatic" half of the task is real. Run as frozen. |
| `audit/derive_ngs004_truth.py` | `ngs-006` | **truth confirmed.** All six planted peaks lie in real sequence at 6.0–7.8× enrichment over the input control. Run as frozen. |
| `audit/derive_ngs005_truth.py` | `ngs-005` | **three of six counts usable, after a correction.** `GENE2`, `GENE4` and `GENE6` come out identical under every placement setting, in both counting conventions: 952/960/961 reads, or 551/561/564 fragments. `GENE1` (712–770 reads) and `GENE5` (745–773) are not determined by the data — they sit against hg19mini's telomeric repeat tracts, where which position a read is given is a property of the aligner. `GENE3` was withdrawn after the run; see below. |

The `ngs-005` sweep is the point of that script rather than an afterthought. An
earlier version printed six confident integers, and three of them moved when
the tie-break changed. Only counts that survive all five settings are used, so
the benchmark never asks an arm for a number two correct aligners would
disagree on.

`ngs-005` is a frozen `llm_rubric` row, and the exact counts are checked
through a `pairs` spec in `checks.json`. A pair spec applies whatever
`check_type` was frozen, because `deterministic_pass` is a separate column from
the model's score: it adds signal without editing the frozen task set or
changing what the judge is shown. The row's rubric still decides whether a real
matrix was produced at all, which is what its ground truth asks.

**The `ngs-005` check was wrong on its first run, and the arms caught it.** It
was frozen with read-level counts only. featureCounts counts each mate
separately by default but counts a pair once under `-p --countReadPairs`, and
all three tool-running arms reported that second column, so the check failed
three correct count matrices. The tell was in the numbers: the gap was a
constant ~396 per gene, not a ratio — the 4 × (readlen − 1) positions where a
mate straddles an exon edge. The script now computes both conventions and the
check accepts either. Two things kept this honest rather than circular: the
corrected values are derived from the fixtures alone, with no reference to any
answer, and the baseline still fails the check, so it did not become a rubber
stamp. The arms' agreement is what prompted the re-derivation; it is not what
the new numbers rest on.

`GENE3` was withdrawn rather than corrected. The derivation calls it exact
(826 reads / 503 fragments) under every setting, but the three arms reported
453, 553 and 556 — agreeing with neither the derivation nor each other. An
unresolved discrepancy must not be scored as ground truth in either direction,
so the gene is dropped from the check and the disagreement is left on the
record. Both corrections are logged per-row in `results/reference_audit.csv`,
and neither touched the `llm_rubric` scores the row was judged on.

## Scoring

1. **Deterministic checks** for `exact`, `numeric_tol` and `contains_all` tasks.
   No model involved, so anyone with `runs.csv` reproduces them exactly. The
   `contains_all` items live in `checks.json`, described above.
2. **A blinded LLM rubric** for everything else. The judge receives the answer
   with all platform and tool names stripped out (`smarts.bio`, `biomni`,
   `galaxy`, `smarts_*`, `mcp__*__*` → "the platform"), plus the task, the ground
   truth and the rubric. It never learns which arm it is grading.
3. **A human audit** of a random 20% plus every disagreement between the two.

Scores are 0-3; success means 2 or better. An explicit "I cannot do this" is
treated as honest and scores well where the capability is genuinely absent — a
fast, clear refusal beats a slow wrong answer, and both beat a fabricated number.

### Things that are not scored as losses

- `infra_error` — quota exhaustion, rate limits, 5xx. Reported separately.
  Compute quota in particular will look exactly like a capability failure if you
  are not careful, so check it before running the GPU rows.
- `ARTIFACT` — a reply that narrates a tool call it could not make. This showed
  up in early baseline runs, where the model wrote a markdown `Tool call: Bash`
  block instead of answering, because it read denied tools as withheld rather
  than absent. The system prompt now states plainly that the visible tools are
  the only ones that exist, and the judge flags any remaining cases for review.

## Known limitations

- **`ngs-003`'s frozen ground truth names the wrong file, and the row is
  scored against variants its own input cannot contain.** The frozen
  `ground_truth` gives `chr1:14000 A>C` and `chr2:11500 A>G` as the germline
  calls to recover from `aligned.bam`. Those two SNVs were planted into the
  tumour/normal FASTQ pair by `fixtures/make_pair_fixtures.py`; `aligned.bam`
  is a separate, pre-existing fixture the generator never touched, and a pileup
  shows the reference base in every read at both loci. The variants the reads
  do carry are six homozygous SNVs at `1:11001 A>G`, `1:14501 C>T`,
  `2:11001 C>T`, `2:14501 C>T`, `4:11001 G>A` and `4:14501 G>A`. Both arms that
  ran the task were therefore failed for being right: Galaxy reported all six
  and smarts.bio five of six, and the blinded judge scored each 0 with a
  hallucination flag. The task, the prompt and the judge scores are left exactly
  as they ran — the correction lives in `results/reference_audit.csv`, derived
  deterministically by `audit/verify_ngs003.py`, which re-scores Galaxy 3 and
  smarts.bio 2 and clears both hallucination flags. smarts.bio's missing
  `2:11001 C>T` is a genuine recall loss and stays on the record. Anyone reading
  `runs.csv` alone will read this row backwards; the audit file is not optional
  here.

- **Biomni's wall time on compute-heavy rows is not comparable to the other
  arms.** The container is `linux/amd64` and the host is `arm64`, so every
  instruction runs through Rosetta translation. This is invisible on rows whose
  work is API calls — Biomni answered `ngs-002` in 59s — but it dominates rows
  that actually compute. `ngs-003` (germline variant calling) hit the 120-minute
  budget as a `timeout` after 14 tool calls, with the container measured at
  99.88% CPU throughout: a single tool call stayed outstanding for over 40
  minutes of emulated execution. Galaxy answered the same row in 172s and
  smarts.bio in 328s. The other arms do not pay this: smarts.bio runs on cloud
  ECS and Galaxy on usegalaxy.org, both native. The penalty is a property of
  *this machine and this packaging*, not of Biomni, so **wall time and any
  timeout attributed to Biomni on a compute-bound row should be read as an
  upper bound, not a measurement of the product.** Accuracy, tool-call counts
  and error rates are unaffected and remain comparable. An amd64 host would remove the distortion;
  the hosted Biomni product, had it been usable here, would not have had it at
  all.

- One tester, one model, one week. Latency depends on time of day and on each
  platform's queue.
- Free and entry tiers are used where possible, as the arms table above records.
  A paid tier may behave differently.
- smarts.bio's design models (Boltz-2, RFdiffusion3) are reachable only through
  its **pipelines**, not as individual tools — `DISABLED_TOOLS` hides them from
  the agent's tool registry while the pipelines call them directly through a
  separate execution path. Design tasks are therefore phrased at pipeline level,
  which is how a scientist would ask anyway.
- The arms are not identically shaped. smarts.bio's MCP offers both individual
  tools and an agent entry point; the Biomni arm serves tools only, and Galaxy
  serves tools only. Where a task is solved by the smarts.bio agent planning
  several steps, part of what is being measured is that planning layer. Runs
  record every tool call, so this is visible per task rather than buried in the
  totals.
- Public benchmark items from LAB-Bench (CC-BY-SA-4.0) and BixBench (Apache-2.0)
  are *not* yet included. They should be, with attribution, to anchor these
  results against an external standard — but their questions must be imported
  from the published datasets, never paraphrased from memory.

## Corrections and right of reply

This benchmark was built by one of the platforms it measures, so the people best
placed to find its mistakes are the maintainers and users of the other three.

**If you think a result here is wrong, we want to hear it, and we will do our
best to address it in the next version.** That applies to anyone: the platforms
under test, their users, and readers with no stake in any of them.

What that means in practice:

- Open an issue, or write to the corresponding author on the paper. Point at the
  task id and the arm; every run has a transcript under `transcripts/` named
  `<task_id>__<arm>__r1.jsonl`, so a disagreement can be argued over the actual
  record rather than over the summary table.
- The most useful objections are that a task is defective, that an arm was
  misconfigured or under-resourced, that a prompt was phrased in a way one
  platform's conventions favour, or that a score is wrong on the evidence. Three
  task defects have already been found this way by internal audit and are
  documented in `TASK-CHANGES.md`; we expect there are more.
- Corrections are published, not quietly applied. A revision will say what
  changed, who raised it, and what it did to the numbers — including if it moves
  them against smarts.bio.
- Scores already published are not rewritten in place. The frozen record stays
  as it is, and a corrected reading is added alongside it, which is the same
  separation the audit layer already uses (`results/reference_audit.csv`).

We do not claim the right to the last word on a benchmark we authored.

## Layout

```
README.md              this file: method, fairness rules, how to reproduce
LICENSE                MIT for code, CC BY 4.0 for data
smartsbio-benchmark-paper.pdf   the full write-up, with figures
paper-print.html       the paper's source; prints to the PDF above

the inputs
  tasks.csv            the frozen task set, exactly as run
  tasks-corrected.csv  three task defects repaired; prepared, never run
  system_prompt.md     identical instructions for every arm
  checks.json          checkable items for the contains_all rows
  fixtures/            input files, with known answers (see Fixture provenance)
  fetch_fixtures.sh    fetches the one fixture too large for git
  requirements.txt     no third-party packages; stdlib plus the `claude` CLI

the arms
  platforms/*.mcp.json           one MCP config per arm (secrets via ${ENV_VAR})
  platforms/biomni.Dockerfile    Biomni's E1 environment, pinned
  platforms/biomni_docker.sh     stdio transport into that container
  platforms/biomni_server.py     serves Biomni's tool library over MCP
  platforms/fetch_biomni_data.sh pre-fetches the data lake
  stage_fixtures.py              uploads fixtures into each arm's workspace
  .env.example                   the credentials each arm needs, with no values

running and scoring
  run_benchmark.py     tasks x arms, one `claude -p` per run
  judge.py             deterministic checks + blinded rubric
  analyze.py           both specifications; reproduces every published figure

the record
  results/runs.csv              one row per run
  results/reference_audit.csv   the audit layer, kept separate from the judge's scores
  results/identifier_audit.csv  registry resolution for every identifier emitted
  transcripts/                  raw stream-json for every run
  audit/                        ground-truth derivation and identifier verification
  TASK-CHANGES.md               every change made to the frozen task set, and when
```

Derived tables (`results/summary.csv`, `by_category.csv`, `head_to_head.csv`)
are deliberately not committed: `analyze.py` regenerates them from `runs.csv`,
so there is one source of truth rather than two that can disagree.

## Fixture provenance

Most fixtures are synthetic: generated by `fixtures/make_pair_fixtures.py`
against the miniature reference, seeded so that re-running reproduces them byte
for byte. They carry no third-party rights. The rest are derived from public
records and are listed here with their accession and licence, so that a reader
can obtain the same inputs from the original source rather than trusting this
copy.

| File(s) | Origin | Terms |
|---|---|---|
| `hg19mini.fasta`, `.fai`, `.dict` | GATK test resources, Broad Institute (`broadinstitute/gatk`). Four 16 kb contigs; the `.dict` still carries the upstream author's local path in its `UR:` field | Apache-2.0 |
| `1UBQ.pdb` | Protein Data Bank entry `1UBQ`, ubiquitin at 1.8 Å (Vijay-Kumar, Bugg & Cook, 1987) | CC0 1.0 — wwPDB releases archive data free of copyright |
| `phix174.fasta` | NCBI RefSeq `NC_001422.1`, *Escherichia* phage phiX174 | Public domain (US Government work; GenBank records are not copyrighted) |
| `hbb_orthologs.fasta` | UniProtKB `P68871`, `P02088`, `P02112` (HBB human, mouse, chicken) | **CC BY 4.0 — attribution required.** The UniProt Consortium, https://www.uniprot.org |
| `pbmc3k_raw.h5ad` | 10x Genomics PBMC 3k, a public demonstration dataset of de-identified donor PBMCs, fetched from the CZ CELLxGENE Discover mirror | 10x Genomics public dataset terms; **not redistributed here** — gitignored and fetched by `fetch_fixtures.sh` |
| `test_variants.vcf`, `test_variants_grch38.vcf` | Positions and rsIDs from dbSNP/Ensembl; the files themselves were written for this benchmark | Public records; files CC BY 4.0 with the task data |

Synthetic, and ours: `reads_R{1,2}`, `tumor_R{1,2}`, `normal_R{1,2}`,
`chip_R{1,2}`, `input_R{1,2}`, `aligned.bam`(`.bai`), `annotation.gtf`,
`counts.csv`, `peaks_hg19mini.bed`, `test_intervals.bed`, `rna.fasta`,
`info_test.csv`, `archive_test.{zip,tar.gz}`, `pair_fixtures_truth.json`.

No fixture contains human subject data beyond the public PBMC 3k demonstration
dataset, which is de-identified and not redistributed here.

## Citing this

```
smarts.bio (2026). A four-arm benchmark of agentic bioinformatics platforms:
smarts.bio, Biomni, Galaxy and an untooled control.
https://github.com/smartsbio/bio-agent-platform-benchmark
```

Cite the commit hash along with it: the hash identifies exactly which version of
the task set, run record and transcripts you read. Note that it dates the
*publication* of those files, not the freeze — `tasks.csv` was frozen before the
first run, but it was not under version control at the time, so the freeze is an
assertion by the author rather than something the history can prove. §2.2 of the
paper says so plainly, and `TASK-CHANGES.md` records every change made since.

A Zenodo DOI will be minted from a tagged release; cite that once it exists.

## Licence

Code MIT, data CC BY 4.0; the split is set out in `LICENSE`. Third-party
fixtures keep the terms in the table above.

Platforms under test keep their own licences: Biomni is Apache-2.0
(`snap-stanford/Biomni`), `galaxy-mcp` is MIT (galaxyproject). Record the commit
or release of each with the results.
