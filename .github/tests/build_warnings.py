#!/usr/bin/env python3
"""List the compiler warnings in a build log that belong to FLARE's own sources.

FLARE CI fails on warnings in FLARE's code; warnings in ERF and its submodules
(``Submodules/``) are checked by ERF's CI. Dropping every line that names
``Submodules/`` is not enough: when GCC inlines library code it reports the warning
against the system header (for example ``/usr/include/c++/11/bits/stl_algobase.h``)
and names the ERF file only in the context lines above it::

    In file included from /usr/include/c++/11/bits/char_traits.h:39,
                     from .../Submodules/ERF/Source/IO/ERF_WriteScalarProfiles.cpp:4:
    In static member function ...,
        inlined from ... at .../Submodules/ERF/Source/IO/ERF_WriteScalarProfiles.cpp:172:17:
    /usr/include/c++/11/bits/stl_algobase.h:431:30: warning: ... [-Wstringop-overflow=]

Each warning or error is therefore attributed to a file:

* a location inside ``--root`` is its own file;
* a location outside ``--root`` (a system or third-party header) is attributed to
  the innermost frame inside ``--root`` among its context lines: the last
  ``inlined from ... at`` or ``required from`` frame, else the translation unit that
  ends the include chain.

The warning is dropped when that file lies under ``Submodules/``. A header warning
with no frame inside ``--root`` is kept, so nothing is dropped on a guess. Any line
that names ``Submodules/`` and linker ``lto-wrapper`` warnings are dropped, as the
grep filter this replaces did.

Usage (prints the kept lines, sorted and unique, then ``Warnings: N``)::

    python3 .github/tests/build_warnings.py --root "$GITHUB_WORKSPACE" build-output.txt
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path
from typing import List, Optional, Sequence

SUBMODULES = "Submodules/"

# A GCC/Clang diagnostic: "path:line[:col]: [fatal ]warning|error: ...".
DIAGNOSTIC_RE = re.compile(r"^(?P<path>[^\s:][^:]*):\d+(?::\d+)?: (?:fatal )?(?:warning|error):")
# Any line the old grep counted.
COUNTED_RE = re.compile(r"warning:|error:", re.IGNORECASE)

# Context lines that GCC prints above a diagnostic, and the path each one names.
INCLUDE_RE = re.compile(r"^(?:In file included from|\s+from) (?P<path>[^\s:][^:]*):\d+[:,]")
FRAME_RE = re.compile(r"^\s*inlined from .* at (?P<path>[^\s:][^:]*):\d+(?::\d+)?[:,]")
REQUIRED_RE = re.compile(r"^(?P<path>[^\s:][^:]*):\d+:\d+:\s+(?:recursively )?required (?:from|by)")
IN_FUNCTION_RE = re.compile(r"^(?:(?P<path>[^\s:][^:]*): )?In (?:static |member |constructor|destructor"
                            r"|function|instantiation|lambda|substitution)")


def inside(path: str, root: str) -> bool:
    """Return whether ``path`` (absolute, or relative to the build) lies in ``root``."""

    return not path.startswith("/") or path == root or path.startswith(root + "/")


def attribute(path: str, context: Sequence[str], root: str) -> Optional[str]:
    """Return the file a diagnostic at ``path`` belongs to, or ``None`` if unknown."""

    if inside(path, root):
        return path
    frames: List[str] = []
    includes: List[str] = []
    for line in context:
        for regex, bucket in ((FRAME_RE, frames), (REQUIRED_RE, frames), (INCLUDE_RE, includes)):
            match = regex.match(line)
            if match:
                bucket.append(match.group("path"))
                break
    for candidate in reversed(frames):
        if inside(candidate, root):
            return candidate
    # The include chain runs from the header out to the translation unit.
    for candidate in reversed(includes):
        if inside(candidate, root):
            return candidate
    return None


def is_context(line: str) -> bool:
    """Return whether ``line`` is a context line GCC prints above a diagnostic."""

    return any(regex.match(line) for regex in (INCLUDE_RE, FRAME_RE, REQUIRED_RE, IN_FUNCTION_RE))


def flare_warnings(lines: Sequence[str], root: str) -> List[str]:
    """Return the sorted, unique warning lines that belong to FLARE's sources."""

    root = root.rstrip("/")
    kept = set()
    context: List[str] = []
    for raw in lines:
        line = raw.rstrip("\n")
        if not COUNTED_RE.search(line):
            context = context + [line] if is_context(line) else []
            continue
        if SUBMODULES in line or "lto-wrapper: warning:" in line:
            context = []
            continue
        match = DIAGNOSTIC_RE.match(line)
        owner = attribute(match.group("path"), context, root) if match else None
        if owner is None or SUBMODULES not in owner:
            kept.add(line)
        context = []
    return sorted(kept)


def main(argv: Optional[Sequence[str]] = None) -> int:
    """Print FLARE's warnings from a build log; the CI step fails on a nonzero count."""

    parser = argparse.ArgumentParser(description="List the warnings in a build log that belong to FLARE.")
    parser.add_argument("--root", required=True, help="repository root as it appears in the log")
    parser.add_argument("log", type=Path, help="build log (compiler output)")
    args = parser.parse_args(argv)
    try:
        text = args.log.read_text(encoding="utf-8", errors="replace")
    except OSError as exc:
        print(f"build_warnings: could not read {args.log}: {exc}", file=sys.stderr)
        return 2
    kept = flare_warnings(text.splitlines(), args.root)
    for line in kept:
        print(line)
    print(f"Warnings: {len(kept)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
