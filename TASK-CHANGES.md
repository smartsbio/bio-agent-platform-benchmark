# Changes to the frozen task set

`tasks.csv` is pre-registered: it is frozen before a run, and it is the set that
produced `results/runs.csv` and every number in the paper. This file records the
changes made to it after that run, so a reader can see which defects the audit
found and what a corrected set would look like.

**Two task files ship in this repository, and the distinction matters:**

| file | what it is |
|---|---|
| `tasks.csv` | the frozen set, exactly as run. `runs.csv`, `judge.py` and the paper all refer to this one. |
| `tasks-corrected.csv` | the same set with the three defects below repaired. **Prepared but never executed:** no run was made against it, and no published number comes from it. |

The benchmark was run once, across all four arms. The corrections below were
written for a repeat run that was not carried out, so they are published as a
starting point for anyone who reruns rather than as part of this result.

**Rule followed here:** nothing was changed while a run was in flight. These
edits were made after the run had completed in full.

Every change below fixes a defect in the *task*, found by audit. None adjusts a
rubric's difficulty, and none was made because of how an arm scored. Two of the
three tasks (`blast-001`, `var-001`) are voided in the analysis for exactly the
defects described here; the third (`ngs-003`) was rescored in the audit layer
against the corrected ground truth, which is where §2.5 of the paper requires
such corrections to live.

---

## 2026-11-10 — three task defects, all found in the audit

Each heading below gives the defect as it stood in `tasks.csv`, and the repair as
it stands in `tasks-corrected.csv`.

### blast-001 — the prompt had no sequence in it

The prompt read *"I pulled this protein sequence out of a structure file and
don't know what it is. Can you identify it?"* and then stopped. No sequence
followed, and `input_files` was empty.

All four arms correctly reported that nothing had arrived and declined to guess.
The rubric required them to name ubiquitin, so all four scored 0, and the task
was voided in the audit layer.

**Change:** the ubiquitin sequence is now inline in the prompt, which is also how
a user would paste it. Ground truth and rubric unchanged.

### ngs-003 — the ground truth named variants from a different fixture

The frozen ground truth was *"chr1:14000 A>C and chr2:11500 A>G"*. Those SNVs do
not exist in `aligned.bam`: an mpileup shows the reference base in every read at
both loci. They belong to the tumour/normal fixtures.

Against the variants the reads actually carry, Galaxy had reported 6 of 6 and
smarts.bio 5 of 6, both with no false positives — and both scored 0.

**Change:** ground truth now names the six SNVs the fixture actually contains,
taken from `samtools mpileup -f hg19mini.fasta aligned.bam`:

| position | REF>ALT | depth |
|---|---|---|
| 1:11001 | A>G | 43 |
| 1:14501 | C>T | 42 |
| 2:11001 | C>T | 49 |
| 2:14501 | C>T | 46 |
| 4:11001 | G>A | 34 |
| 4:14501 | G>A | 42 |

All six are homozygous-alt. The rubric now asks for all six with correct
REF>ALT for full credit. The fixture itself is unchanged.

### var-001 — the fixture was GRCh37 wearing a GRCh38 header

`test_variants.vcf` declared `##reference=hg38`. Its REF alleles match GRCh37
exactly and mismatch GRCh38 at four of five positions:

| chr1 pos | VCF REF | GRCh37 | GRCh38 |
|---|---|---|---|
| 881918 | G | G ✓ | A |
| 888659 | T | T ✓ | G |
| 909238 | G | G ✓ | G |
| 978762 | G | G ✓ | C |
| 1269554 | T | T ✓ | G |

Checked against Ensembl for both builds. Corroborated independently: tool-tests
row `cpu-057` already passed `assembly=GRCh37` for this same file, so the
mislabelling had been noticed before, in another repository, without the header
being corrected.

Two things followed from it. The frozen ground truth promised
missense/synonymous/intronic classifications, but annotated as the header claims,
the file contains **no coding variants at all** — so no arm could reach 3, and
the three 2s are capped rather than ranked. And the judge docked smarts.bio for
HGVSc alleles "contradicting" the stated REF>ALT, when those HGVSc strings are
correct VEP output from GRCh38 bases, reverse-complemented for the minus-strand
transcripts (PERM1, UBE2J2), while the REF>ALT is correct for the VCF. They
disagree because the build is wrong. The agent read the header, passed
`assembly=GRCh38` accordingly, and VEP did exactly as instructed.

**Change:** the task is voided across all four arms (as blast-001 was), and in
`tasks-corrected.csv` the fixture is replaced with five variants whose REF
alleles are the GRCh38 reference base, verified individually against Ensembl:

| variant | rsID | gene | consequence |
|---|---|---|---|
| 1:11796321 G>A | rs1801133 | MTHFR | missense, p.Ala222Val |
| 6:26092913 G>A | rs1800562 | HFE | missense, p.Cys282Tyr |
| 11:5227002 T>A | rs334 | HBB | missense, p.Glu7Val |
| 7:87509329 A>T | rs1045642 | ABCB1 | synonymous, p.Ile1145= |
| 2:135851076 G>A | rs4988235 | MCM6 | intron, c.1917+326C>T |

Consequences are Ensembl VEP's own output on the canonical transcript, so the
ground truth is reproducible rather than asserted. The set now actually contains
the three consequence classes the task asks about, which the old fixture never
did. Three of the five sit on minus-strand transcripts, so a correct answer still
has to handle the strand convention that caused the original confusion — the task
has not been made easier.

**Renamed, deliberately.** The new fixture is `test_variants_grch38.vcf`, not
`test_variants.vcf`. The benchmark stages into `tests/` of the Mors Lab
workspace, which is the same folder the `smarts-bio-mcp` tool-tests suite uses —
18 fixture names are shared between the two suites. Until now every shared file
was a byte-identical copy, so sharing was harmless; this is the first file where
the two suites need *different* contents. tool-tests keeps a GRCh37
`test_variants.vcf` because its rows depend on it (`cpu-046` queries
`1:880000-910000`, `cpu-057` passes `assembly=GRCh37`), so the benchmark takes a
new name rather than overwriting a file another suite relies on.

The old files are kept as `fixtures/test_variants.vcf*.grch37-mislabelled.bak`.

### Also corrected: the tool-tests copy's header

The same mislabelled VCF lives in `smarts-bio-mcp/tool-tests/fixtures/`. Its data
is GRCh37 and its rows already compensate, but the header still said
`##reference=hg38` — and that header is what led the benchmark's agent to choose
GRCh38 in the first place. It now reads `##reference=GRCh37`. The variant records
are byte-identical, so `cpu-046` and `cpu-057` are unaffected; `.gz` and `.tbi`
were regenerated with `bgzip`/`tabix` and the region query still returns its
three records.

### Effect on the headline

Voiding var-001 moves mean accuracy by at most 0.02 and changes no ordering:

| Arm | before | after |
|---|---|---|
| smarts.bio | 2.60 | 2.61 |
| Biomni | 2.30 | 2.30 |
| Galaxy | 2.00 | 2.00 |
| Claude alone | 0.98 | 1.00 |

### Staging of the corrected fixture — done 2026-11-10

Preparation for the repeat run that was not carried out. `test_variants_grch38.vcf`
and its `.gz`/`.tbi` were staged on every arm:

| arm | how | state |
|---|---|---|
| smarts.bio | `stage_fixtures.py --arm smartsbio` → `tests/` | uploaded, read back and verified |
| Galaxy | `stage_fixtures.py --arm galaxy` → history `bbd44e69cb8906b5` | uploaded |
| Biomni | `biomni_docker.sh` bind-mounts `fixtures/` read-only at `/fixtures` | automatic |
| Claude alone | nothing staged, by design | n/a |

The corrected tool-tests fixtures were re-uploaded to `tests/` as well, so the
workspace no longer serves a VCF whose header contradicts its contents.

`stage_fixtures.py` gained two flags for this, because the script is idempotent
by name and so reported the stale fixture as "present" and left it in place:

- `--replace a,b,c` re-uploads named files even when present, deleting the stale
  copy first on Galaxy (a history holds two same-named datasets happily, and a
  tool would then pick between them at random).
- `--verify NAME=MARKER` reads the file back from smarts.bio and requires the
  marker to appear. A 2xx on upload only means the request was accepted; the
  question for a replace is whether the old bytes are gone.

### If the benchmark is rerun

`blast-001`, `ngs-003` and `var-001` must be rerun for **all four arms**, not
only the ones whose scores looked wrong. Nothing in this repository has been
rerun: every number published comes from the single run against `tasks.csv`.
