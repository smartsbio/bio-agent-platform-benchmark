# Biomni's documented E1 environment, on the platform its authors test against.
#
# Built for linux/amd64 deliberately. bioconda has no osx-arm64 build for the
# perl chain that `blast` depends on, so on an Apple Silicon host the bio_env.yml
# solve fails and takes all 66 of its packages with it — biopython, scanpy,
# gget, gseapy, rdkit, samtools, bedtools and BLAST itself. Benchmarking that
# environment would measure the host, not Biomni.
#
# Build (hours, under emulation on an arm64 host):
#   docker build --platform linux/amd64 -f platforms/biomni.Dockerfile \
#     -t biomni-bench:e1 platforms/
#
FROM --platform=linux/amd64 ubuntu:22.04

SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
        git curl wget bzip2 unzip xz-utils ca-certificates \
        build-essential procps \
    && rm -rf /var/lib/apt/lists/*

# Miniforge gives conda + mamba on the conda-forge channel.
ENV CONDA_DIR=/opt/conda
RUN curl -fsSL -o /tmp/miniforge.sh \
        "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" \
    && bash /tmp/miniforge.sh -b -p "$CONDA_DIR" \
    && rm /tmp/miniforge.sh
ENV PATH="$CONDA_DIR/bin:$PATH"

# The version under test. Pinned so the image is reproducible and citable.
ARG BIOMNI_COMMIT=400c1f366b96a35ca253e13c9b06c5076af41d65
RUN git clone https://github.com/snap-stanford/Biomni.git /opt/biomni \
    && git -C /opt/biomni checkout "$BIOMNI_COMMIT" \
    && git -C /opt/biomni rev-parse HEAD > /opt/biomni-commit.txt

# E1 proper: base env, core bioinformatics tools, R packages, CLI tools.
# NON_INTERACTIVE=1 answers the script's prompts and installs every component.
ENV BIOMNI_TOOLS_DIR=/opt/biomni_tools \
    NON_INTERACTIVE=1
RUN cd /opt/biomni/biomni_env && bash setup.sh

# setup.sh's step 6 treats the CLI tool installer as optional and continues when
# it aborts, which it does during an image build: setup_path.sh gets written but
# $BIOMNI_TOOLS_DIR never appears, leaving FastTree, iqtree2, muscle, plink2 and
# HOMER absent. One of them backs the alignment-and-tree task, so run E1's own
# installer directly and let the verification step below enforce the result.
RUN cd /opt/biomni/biomni_env \
    && BIOMNI_AUTO_INSTALL=1 bash install_cli_tools.sh

# environment.yml installs Biomni's dependencies but not Biomni itself.
RUN "$CONDA_DIR/envs/biomni_e1/bin/pip" install --no-cache-dir -e /opt/biomni

# DEVIATIONS FROM UPSTREAM. Three packages, each one needed to make Biomni's own
# code run in Biomni's own environment. Nothing here adds capability; without
# them the platform would score worse than it deserves.
#
#   fair-esm    biomni/tool/genomics.py imports `esm` at module level and calls
#               esm.pretrained.load_model_and_alphabet("esm2_t6_8M_UR50D").
#               E1 does not install it, so the module cannot import and every
#               genomics tool disappears from the catalogue.
#   pybiomart   imported by the same module, likewise absent from E1.
#   mcp<2       E1 installs mcp 2.3.0, in which FastMCP was renamed MCPServer,
#               but a1.create_mcp_server() imports mcp.server.fastmcp.FastMCP.
#               As shipped, Biomni's MCP server cannot start at all.
RUN "$CONDA_DIR/envs/biomni_e1/bin/pip" install --no-cache-dir \
        fair-esm==2.0.0 pybiomart==0.2.0 "mcp<2"

ENV BIOMNI_PY="$CONDA_DIR/envs/biomni_e1/bin/python" \
    PATH="$CONDA_DIR/envs/biomni_e1/bin:/opt/biomni_tools/bin:$PATH"

# Fail the build, not the benchmark, if the environment is incomplete: every
# served tool module must import, and the CLI tools must be on PATH.
RUN failed=""; \
    for m in database genomics genetics molecular_biology biochemistry \
             pharmacology literature systems_biology cancer_biology \
             biophysics support_tools; do \
        if err=$("$BIOMNI_PY" -c "import biomni.tool.$m" 2>&1); then \
            echo "ok      biomni.tool.$m"; \
        else \
            echo "FAILED  biomni.tool.$m: $(echo "$err" | tail -1)"; \
            failed="$failed $m"; \
        fi; \
    done; \
    echo "--- CLI tools ---"; \
    for t in blastn blastp makeblastdb samtools bedtools Rscript macs2 fastqc \
             iqtree2 FastTree muscle plink2; do \
        path=$(command -v "$t" || true); \
        printf '%-14s %s\n' "$t" "${path:-MISSING on PATH}"; \
        [ -n "$path" ] || failed="$failed $t"; \
    done; \
    if [ -n "$failed" ]; then \
        echo "Incomplete environment, refusing to tag:$failed"; exit 1; \
    fi

# The data lake is a mounted volume, not an image layer: ~11 GB, and it is data
# rather than environment.
ENV BIOMNI_DATA_PATH=/data
VOLUME /data

# The server script is mounted so it can be edited without a rebuild.
ENTRYPOINT ["/opt/conda/envs/biomni_e1/bin/python", "/app/biomni_server.py"]
