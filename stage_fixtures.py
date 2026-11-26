#!/usr/bin/env python3
"""Stage the benchmark fixtures onto each platform. Phase 0, untimed.

Batches C, D and E hand the model the user's own data. Those fixtures are
synthetic and have planted ground truth (6 SNPs, 3 somatic + 2 germline
variants, 6 ChIP peak regions, exact featureCounts), so unlike phiX174 they
cannot be fetched from a public database — if they are not staged, the tasks
are unanswerable and every arm scores a loss for the wrong reason.

    python3 -I stage_fixtures.py --check              # report, change nothing
    python3 -I stage_fixtures.py --arm smartsbio
    python3 -I stage_fixtures.py --arm galaxy
    python3 -I stage_fixtures.py --arm all
    python3 -I stage_fixtures.py --arm all --replace test_variants.vcf,...

Per the plan this is untimed setup, logged separately as "data-in friction"
rather than counted against any arm's wall clock.

  smartsbio  uploaded to the workspace in .env under tests/
  galaxy     uploaded into one dedicated history on GALAXY_URL
  biomni     nothing to do: biomni_docker.sh already bind-mounts fixtures/
             read-only at /fixtures

Idempotent: a file already present on a platform is skipped, so re-running
after a partial failure costs nothing. Credentials come from .env and are
never printed.

That skipping is wrong when a fixture's *contents* change while its name stays
the same, which is what happened to test_variants.vcf: the staged copy was the
old file and a plain re-run reported it present and left it there. --replace
names the files to re-upload regardless, deleting the stale copy on platforms
where a same-named upload would otherwise accumulate, and reads one back
afterwards to prove the new bytes actually landed.
"""

import argparse
import json
import mimetypes
import os
import pathlib
import sys
import urllib.error
import urllib.request
import uuid

ROOT = pathlib.Path(__file__).resolve().parent
FIX = ROOT / "fixtures"
GALAXY_HISTORY = "platform-benchmark fixtures"

# Every file any task declares in input_files, plus the index/sidecar files a
# tool needs to open them. Derived from tasks.csv rather than listed by hand.
SIDECARS = {
    "hg19mini.fasta": ["hg19mini.fasta.fai", "hg19mini.dict"],
    "aligned.bam": ["aligned.bam.bai"],
    # Both variant fixtures are listed so this script serves either task file.
    # tasks.csv (the set that ran) declares test_variants.vcf; tasks-corrected.csv
    # declares test_variants_grch38.vcf. The corrected file takes a new name
    # rather than overwriting the old one because the staging folder is shared
    # with another test suite that depends on a GRCh37 test_variants.vcf.
    "test_variants.vcf": ["test_variants.vcf.gz",
                          "test_variants.vcf.gz.tbi"],
    "test_variants_grch38.vcf": ["test_variants_grch38.vcf.gz",
                                 "test_variants_grch38.vcf.gz.tbi"],
}


def wanted():
    """The fixture set the task table actually asks for."""
    import csv
    need = set()
    with (ROOT / "tasks.csv").open() as fh:
        for t in csv.DictReader(fh):
            for f in (t["input_files"] or "").split(";"):
                f = f.strip()
                if f:
                    need.add(f)
    for f in list(need):
        need.update(SIDECARS.get(f, []))
    missing = sorted(f for f in need if not (FIX / f).exists())
    if missing:
        sys.exit(f"fixtures/ is missing {len(missing)} declared file(s): "
                 f"{', '.join(missing)}\nRun fixtures/make_pair_fixtures.py "
                 f"or fetch them before staging.")
    return sorted(need)


def env(name, required=True):
    v = os.environ.get(name, "").strip()
    if required and not v:
        sys.exit(f"{name} is not set. Add it to .env (which is gitignored).")
    return v


def load_env():
    p = ROOT / ".env"
    if not p.exists():
        return
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        k, v = line.split("=", 1)
        os.environ.setdefault(k.strip(), v.strip().strip('"').strip("'"))


def multipart(fields, filefield, filename, payload):
    """Build a multipart/form-data body without pulling in requests."""
    b = uuid.uuid4().hex
    out = bytearray()
    for k, v in fields.items():
        out += (f"--{b}\r\nContent-Disposition: form-data; name=\"{k}\"\r\n\r\n"
                f"{v}\r\n").encode()
    ctype = mimetypes.guess_type(filename)[0] or "application/octet-stream"
    out += (f"--{b}\r\nContent-Disposition: form-data; name=\"{filefield}\"; "
            f"filename=\"{filename}\"\r\nContent-Type: {ctype}\r\n\r\n").encode()
    out += payload + b"\r\n"
    out += f"--{b}--\r\n".encode()
    return bytes(out), f"multipart/form-data; boundary={b}"


def _curl(url, method, headers, data, ctype, timeout):
    """Fallback transport.

    The system python on macOS is 3.9 against LibreSSL 2.8.3, which
    usegalaxy.org rejects with TLSV1_ALERT_PROTOCOL_VERSION. curl on the same
    machine negotiates fine, so rather than require a newer interpreter just
    for staging, shell out. Headers go via stdin to keep the API key off the
    process list, where `ps` would otherwise show it.
    """
    import subprocess, tempfile
    cmd = ["curl", "-sS", "--max-time", str(timeout), "-X", method,
           "-w", "\n%{http_code}", "-K", "-"]
    config = [f'url = "{url}"']
    for k, v in (headers or {}).items():
        config.append(f'header = "{k}: {v}"')
    if ctype:
        config.append(f'header = "Content-Type: {ctype}"')
    tmp = None
    if data:
        tmp = tempfile.NamedTemporaryFile(delete=False)
        tmp.write(data)
        tmp.close()
        config.append(f'data-binary = "@{tmp.name}"')
    try:
        p = subprocess.run(cmd, input="\n".join(config).encode(),
                           capture_output=True, timeout=timeout + 30)
        out = p.stdout.decode(errors="replace")
        body, _, code = out.rpartition("\n")
        code = int(code.strip() or 0)
        if code >= 400:
            raise urllib.error.HTTPError(url, code, body[:300], None, None)
        try:
            return json.loads(body)
        except json.JSONDecodeError:
            return {"raw": body[:400]}
    finally:
        if tmp:
            os.unlink(tmp.name)


def req(url, method="GET", headers=None, data=None, ctype=None, timeout=300):
    r = urllib.request.Request(url, method=method, data=data)
    for k, v in (headers or {}).items():
        r.add_header(k, v)
    if ctype:
        r.add_header("Content-Type", ctype)
    try:
        with urllib.request.urlopen(r, timeout=timeout) as fh:
            body = fh.read()
    except urllib.error.URLError as e:
        if isinstance(e, urllib.error.HTTPError):
            raise
        return _curl(url, method, headers, data, ctype, timeout)
    try:
        return json.loads(body)
    except json.JSONDecodeError:
        return {"raw": body[:400].decode(errors="replace")}


# ---------------------------------------------------------------- smarts.bio

def smartsbio_existing(base, key, ws):
    # Three API quirks, each of which silently produced a wrong answer here:
    #  - the parameter is workspace_id, snake_case. workspaceId returns HTTP
    #    500 "User does not have permission to view files in this workspace",
    #    which reads like an auth failure and is not one.
    #  - the listing is NOT recursive and is capped at 100 entries, so a
    #    top-level call reports files inside tests/ as absent and would
    #    re-upload all of them on every run.
    #  - the subfolder parameter is `path`. `prefix` and `folder` are accepted
    #    and ignored, returning the top level as though the folder were empty.
    try:
        d = req(f"{base}/v1/files?workspace_id={ws}&path=tests",
                headers={"Authorization": f"Bearer {key}"})
    except urllib.error.HTTPError as e:
        sys.exit(f"smarts.bio: listing files failed with HTTP {e.code}: "
                 f"{e.read()[:300].decode(errors='replace')}")
    # An error body must never be read as "no files staged" — that would make
    # this script report a clean slate and re-upload everything on top of it.
    if isinstance(d, dict) and d.get("status") == "error":
        sys.exit(f"smarts.bio: {d.get('message')}")
    inner = d.get("data", d) if isinstance(d, dict) else d
    items = inner if isinstance(inner, list) else (
        inner.get("files") or inner.get("items") or [])
    names = set()
    for f in items:
        if isinstance(f, dict):
            v = f.get("key") or f.get("path") or f.get("name") or ""
            names.add(v.rstrip("/").rsplit("/", 1)[-1])
    return names


def smartsbio_verify(base, key, ws, name, expect):
    """Read a staged file back and confirm a marker string is in it.

    A 2xx on upload only means the request was accepted. For a replace, the
    question is whether the old bytes are gone, and the only way to know is to
    read what the platform now serves.
    """
    # Two things this must not do, both of which produced a false verdict here.
    # The parameter is `key`, not `path`: /v1/files/content takes a file key
    # (workspace-relative is fine), while `path` belongs to the listing endpoint
    # and is a 400 here. And the body cannot go through req(), which truncates a
    # non-JSON response to 400 characters -- long enough to look like a real
    # answer and short enough to miss a marker further down the file, which is
    # how a correctly staged fixture got reported as stale.
    import subprocess
    url = f"{base}/v1/files/content?workspace_id={ws}&key=tests/{name}"
    try:
        out = subprocess.run(
            ["curl", "-sS", "--max-time", "120", "-w", "\n%{http_code}",
             "-H", f"Authorization: Bearer {key}", url],
            capture_output=True, timeout=180)
    except subprocess.TimeoutExpired:
        return None, "timed out"
    text = out.stdout.decode(errors="replace")
    body, _, code = text.rpartition("\n")
    if code.strip() != "200":
        return None, f"HTTP {code.strip() or '?'}"
    return (expect in body), None


def stage_smartsbio(files, check, replace=()):
    base = os.environ.get("SMARTSBIO_BASE_URL", "https://api.smarts.bio").rstrip("/")
    key, ws = env("SMARTSBIO_API_KEY"), env("SMARTSBIO_WORKSPACE_ID")
    have = smartsbio_existing(base, key, ws)
    todo = [f for f in files if f not in have or f in replace]
    print(f"\nsmarts.bio  {base}  workspace ...{ws[-6:]}")
    print(f"  present {len(have & set(files))}/{len(files)}, to upload {len(todo)}")
    if check or not todo:
        return
    for i, name in enumerate(todo, 1):
        payload = (FIX / name).read_bytes()
        body, ctype = multipart({"workspace_id": ws, "path": "tests",
                                 "description": "platform-benchmark fixture"},
                                "file", name, payload)
        try:
            req(f"{base}/v1/files/upload", "POST",
                {"Authorization": f"Bearer {key}"}, body, ctype)
            tag = " (replaced)" if name in replace else ""
            print(f"  [{i}/{len(todo)}] {name}  {len(payload)/1024:.0f} KB  ok{tag}")
        except urllib.error.HTTPError as e:
            print(f"  [{i}/{len(todo)}] {name}  FAILED HTTP {e.code} "
                  f"{e.read()[:200].decode(errors='replace')}")


# -------------------------------------------------------------------- Galaxy

def stage_galaxy(files, check, replace=()):
    url = env("GALAXY_URL").rstrip("/")
    key = env("GALAXY_API_KEY")
    hdr = {"x-api-key": key}
    hs = req(f"{url}/api/histories", headers=hdr)
    hist = next((h for h in hs if h.get("name") == GALAXY_HISTORY), None)
    print(f"\ngalaxy  {url}")
    if hist is None:
        if check:
            print(f"  history '{GALAXY_HISTORY}' does not exist; "
                  f"{len(files)} file(s) would be uploaded")
            return
        body = json.dumps({"name": GALAXY_HISTORY}).encode()
        hist = req(f"{url}/api/histories", "POST", hdr, body, "application/json")
        print(f"  created history {hist['id']}")
    hid = hist["id"]
    contents = req(f"{url}/api/histories/{hid}/contents", headers=hdr)
    have = {c.get("name") for c in contents if not c.get("deleted")}

    # A history happily holds two datasets of the same name, and a tool picking
    # one of them would be picking at random between the old fixture and the new
    # one. So the stale copies go first.
    if replace and not check:
        for c in contents:
            if c.get("name") in replace and not c.get("deleted"):
                try:
                    req(f"{url}/api/histories/{hid}/contents/{c['id']}",
                        "DELETE", hdr)
                    print(f"  deleted stale {c['name']} ({c['id']})")
                    have.discard(c["name"])
                except urllib.error.HTTPError as e:
                    print(f"  could NOT delete stale {c.get('name')}: "
                          f"HTTP {e.code} — upload skipped to avoid a duplicate")
                    replace = tuple(x for x in replace if x != c.get("name"))

    todo = [f for f in files if f not in have or f in replace]
    print(f"  history {hid}: present {len(have & set(files))}/{len(files)}, "
          f"to upload {len(todo)}")
    if check or not todo:
        return
    for i, name in enumerate(todo, 1):
        payload = (FIX / name).read_bytes()
        body, ctype = multipart(
            {"tool_id": "upload1", "history_id": hid,
             "inputs": json.dumps({"dbkey": "?", "file_count": "1",
                                   "files_0|type": "upload_dataset",
                                   "files_0|NAME": name})},
            "files_0|file_data", name, payload)
        try:
            req(f"{url}/api/tools", "POST", hdr, body, ctype)
            print(f"  [{i}/{len(todo)}] {name}  {len(payload)/1024:.0f} KB  ok")
        except urllib.error.HTTPError as e:
            print(f"  [{i}/{len(todo)}] {name}  FAILED HTTP {e.code} "
                  f"{e.read()[:200].decode(errors='replace')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--arm", default="", choices=["smartsbio", "galaxy", "all", ""])
    ap.add_argument("--check", action="store_true",
                    help="report what is staged and change nothing")
    ap.add_argument("--replace", default="",
                    help="comma-separated fixture names to re-upload even if "
                         "already present (use when contents changed but the "
                         "name did not)")
    ap.add_argument("--verify", default="",
                    help="NAME=MARKER — read NAME back from smarts.bio after "
                         "staging and require MARKER to appear in it")
    args = ap.parse_args()
    load_env()
    files = wanted()
    print(f"{len(files)} fixture(s) declared by tasks.csv, all present locally")

    arms = ["smartsbio", "galaxy"] if args.arm in ("all", "") else [args.arm]
    if not args.arm and not args.check:
        sys.exit("Pass --arm smartsbio|galaxy|all to upload, or --check to report.")
    replace = tuple(x.strip() for x in args.replace.split(",") if x.strip())
    unknown = [r for r in replace if r not in files]
    if unknown:
        sys.exit(f"--replace names file(s) no task declares: {', '.join(unknown)}")
    if replace:
        print(f"replacing on every arm: {', '.join(replace)}")

    for a in arms:
        if a == "smartsbio":
            stage_smartsbio(files, args.check, replace)
        else:
            stage_galaxy(files, args.check, replace)

    if args.verify and not args.check:
        name, _, marker = args.verify.partition("=")
        base = os.environ.get("SMARTSBIO_BASE_URL",
                              "https://api.smarts.bio").rstrip("/")
        ok, err = smartsbio_verify(base, env("SMARTSBIO_API_KEY"),
                                   env("SMARTSBIO_WORKSPACE_ID"),
                                   name.strip(), marker.strip())
        if err:
            print(f"\nverify  {name}: could not read it back ({err})")
        elif ok:
            print(f"\nverify  {name} on smarts.bio contains {marker!r} — ok")
        else:
            sys.exit(f"\nverify  {name} on smarts.bio does NOT contain "
                     f"{marker!r}: the stale copy is still being served.")
    print("\nbiomni  nothing to do: biomni_docker.sh bind-mounts fixtures/ "
          "read-only at /fixtures")
    print("baseline  nothing to stage, by design (no tools)")


if __name__ == "__main__":
    main()
