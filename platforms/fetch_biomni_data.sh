#!/bin/bash
# Pre-fetch Biomni's data lake.
#
# A1's own downloader fetches the 76 files one at a time with no timeout, and a
# dropped connection leaves it spinning indefinitely — observed hanging for an
# hour after two files, with the socket in CLOSE_WAIT. There is a data_lake.zip
# on the bucket but it returns 403, which is why their code falls back to the
# per-file path. So fetch the files here instead: parallel, resumable, and with
# timeouts, before the server ever starts. A benchmark run should not have an
# 11 GB download inside its first timed call.
#
# Usage: ./fetch_biomni_data.sh [parallelism]
set -uo pipefail

BUCKET="https://biomni-release.s3.amazonaws.com/data_lake"
IMAGE="${BIOMNI_IMAGE:-biomni-bench:e1}"
DEST="${BIOMNI_DATA_PATH:-$HOME/biomni-data}/biomni_data/data_lake"
JOBS="${1:-6}"

mkdir -p "$DEST"

# Take the file list from the image, so it cannot drift from the pinned commit.
LIST=$(docker run --rm --platform linux/amd64 --entrypoint python "$IMAGE" -c \
    'from biomni.env_desc import data_lake_dict as d; print("\n".join(sorted(d)))')
total=$(echo "$LIST" | grep -c .)
echo "$total files expected in $DEST"

fetch() {
    local name="$1" dest="$2" bucket="$3"
    local out="$dest/$name"
    local remote
    remote=$(curl -sI -m 30 "$bucket/$name" | awk 'tolower($1)=="content-length:"{print $2+0}' | tr -d '\r')
    if [ -f "$out" ] && [ -n "$remote" ] && [ "$(wc -c < "$out")" -eq "$remote" ]; then
        echo "  have  $name"; return 0
    fi
    if curl -fL --retry 5 --retry-delay 3 --retry-all-errors \
            --connect-timeout 20 --max-time 3600 --speed-limit 1024 --speed-time 120 \
            -C - -o "$out" "$bucket/$name" 2>/dev/null; then
        echo "  ok    $name"
    else
        echo "  FAIL  $name"; return 1
    fi
}
export -f fetch

echo "$LIST" | grep . | xargs -P "$JOBS" -I{} bash -c 'fetch "$@"' _ {} "$DEST" "$BUCKET"

have=$(find "$DEST" -type f | wc -l | tr -d ' ')
echo "done: $have/$total files, $(du -sh "$DEST" | cut -f1)"
[ "$have" -eq "$total" ]
