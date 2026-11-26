#!/bin/bash
# stdio MCP transport for the Biomni arm: one container per run, talking the
# protocol over the container's stdin/stdout. `-i` is what makes that work.
#
# Mounts, all of them necessary:
#   /data       the data lake, so the ~11 GB download survives the container
#   /app        the server script, so editing it needs no rebuild
#   /fixtures   the benchmark's input files, read-only
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
IMAGE="${BIOMNI_IMAGE:-biomni-bench:e1}"
DATA="${BIOMNI_DATA_PATH:-$HOME/biomni-data}"

if ! docker image inspect "$IMAGE" >/dev/null 2>&1; then
    echo "Image $IMAGE is not built. See the Biomni arm section of the README." >&2
    exit 1
fi

mkdir -p "$DATA"

exec docker run --rm -i \
    --platform linux/amd64 \
    -v "$DATA:/data" \
    -v "$ROOT/platforms/biomni_server.py:/app/biomni_server.py:ro" \
    -v "$ROOT/fixtures:/fixtures:ro" \
    -e BIOMNI_DATA_PATH=/data \
    -e PYTHONUNBUFFERED=1 \
    -e "ANTHROPIC_API_KEY=${ANTHROPIC_API_KEY:-}" \
    -e "BIOMNI_COMMERCIAL_MODE=${BIOMNI_COMMERCIAL_MODE:-}" \
    "$IMAGE"
