#!/usr/bin/env python3
"""Expose the Biomni tool library over stdio MCP, for the `biomni` arm.

Biomni ships an agent (A1) and a tool library. This serves the *tools*, so the
same Claude model that drives every other arm does the reasoning here too and
the comparison stays a comparison of tool layers. A1's own LLM is never called.

Setup, once:

    pip install biomni
    export BIOMNI_DATA_PATH=~/biomni_data    # ~11 GB, downloaded on first run

Modules are chosen to cover the benchmark's categories. Serving every module
would put a very large tool catalogue in the context window, which is a
limitation of the client rather than of Biomni, so the selection is explicit
and recorded here instead of being silently truncated.
"""

import os
import sys

MODULES = [
    "biomni.tool.database",
    "biomni.tool.genomics",
    "biomni.tool.genetics",
    "biomni.tool.molecular_biology",
    "biomni.tool.biochemistry",
    "biomni.tool.pharmacology",
    "biomni.tool.literature",
    "biomni.tool.systems_biology",
    "biomni.tool.cancer_biology",
    "biomni.tool.biophysics",
    "biomni.tool.support_tools",
]


def main():
    from biomni.agent.a1 import A1

    data_path = os.environ.get("BIOMNI_DATA_PATH")
    if not data_path:
        sys.exit(
            "BIOMNI_DATA_PATH is not set. Point it at a directory with room for "
            "the ~11 GB data lake; see .env.example."
        )

    # commercial_mode=True drops tools whose data sources are not licensed for
    # commercial use. Set BIOMNI_COMMERCIAL_MODE=1 to compare like with like
    # against a platform that gates those sources; default is full capability.
    commercial = os.environ.get("BIOMNI_COMMERCIAL_MODE", "") == "1"

    # create_mcp_server() catches a failed module import and merely warns, which
    # would serve a quietly smaller tool set and understate Biomni's capability.
    # Import every module here first and refuse to start if any is unavailable,
    # so a missing dependency is a setup error and never a benchmark result.
    import importlib

    broken = {}
    for name in MODULES:
        try:
            importlib.import_module(name)
        except Exception as exc:  # noqa: BLE001 - reported, not handled
            broken[name] = f"{type(exc).__name__}: {exc}"
    if broken:
        lines = "\n".join(f"  {m}\n      {e}" for m, e in broken.items())
        sys.exit(
            f"{len(broken)} of {len(MODULES)} Biomni tool modules will not "
            f"import, so their tools would be missing from the catalogue "
            f"without being missing from Biomni:\n{lines}\n\n"
            "Install the dependencies into this venv, or see the Biomni arm "
            "section of the README."
        )

    # use_tool_retriever=False: the retriever is how Biomni's own agent narrows
    # a large tool set before choosing. Here Claude does the choosing and
    # create_mcp_server() reads the registry directly, so the retriever is
    # never consulted — and loading its embedding model cost enough memory to
    # get the container OOM-killed. Switching it off removes an unused
    # component; it does not remove a tool.
    agent = A1(path=data_path, commercial_mode=commercial,
               use_tool_retriever=False)
    mcp = agent.create_mcp_server(tool_modules=MODULES)
    print(
        f"Biomni MCP server: {len(MODULES)} modules, "
        f"commercial_mode={commercial}",
        file=sys.stderr,
    )
    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
