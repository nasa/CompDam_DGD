#!/usr/bin/env python3
"""Check Markdown links in a repository.

Validates:
- Internal links: relative paths and markdown anchors.
- External links: HTTP/HTTPS responses using requests.

Exit codes:
- 0: all links valid
- 1: one or more broken/dead links found
"""

from __future__ import annotations

import argparse
import bisect
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict
from typing import Iterable
from typing import Iterator
from typing import List
from typing import Optional
from typing import Sequence
from typing import Set
from typing import Tuple
from urllib.parse import unquote
from urllib.parse import urlsplit

import requests


SKIP_SCHEMES = {"mailto", "tel", "javascript", "data"}
FALLBACK_TO_GET_CODES = {403, 405, 501}
ACCEPTABLE_RESPONSE_CODES = {401, 403, 429}
SKIP_DIR_NAMES = {".git", ".venv", "venv", "__pycache__", "node_modules"}
HEADING_RE = re.compile(r"^#{1,6}\s+(.+?)\s*$", re.MULTILINE)


@dataclass(frozen=True)
class LinkRef:
    source_file: Path
    line: int
    target: str
    is_image: bool


def _line_starts(text: str) -> List[int]:
    starts = [0]
    for idx, char in enumerate(text):
        if char == "\n":
            starts.append(idx + 1)
    return starts


def _line_for_offset(starts: Sequence[int], offset: int) -> int:
    return bisect.bisect_right(starts, offset)


def _extract_markdown_links(text: str, source_file: Path) -> Iterator[LinkRef]:
    """Yield inline markdown links and image links with line numbers.

    Supports forms like:
    - [text](target)
    - ![alt](target)
    - [text](<target with spaces>)
    """

    starts = _line_starts(text)
    i = 0
    n = len(text)
    while i < n:
        is_image = text.startswith("![", i)
        link_start = i + 2 if is_image else i + 1
        if not (is_image or text.startswith("[", i)):
            i += 1
            continue

        close_bracket = text.find("]", link_start)
        if close_bracket == -1 or close_bracket + 1 >= n or text[close_bracket + 1] != "(":
            i += 1
            continue

        j = close_bracket + 2
        depth = 1
        in_angle = False
        while j < n and depth > 0:
            ch = text[j]
            if ch == "<" and depth == 1:
                in_angle = True
            elif ch == ">" and depth == 1 and in_angle:
                in_angle = False
            elif not in_angle:
                if ch == "(":
                    depth += 1
                elif ch == ")":
                    depth -= 1
            j += 1

        if depth != 0:
            i += 1
            continue

        raw_target = text[close_bracket + 2 : j - 1].strip()
        if raw_target:
            target = _parse_link_target(raw_target)
            if target:
                yield LinkRef(
                    source_file=source_file,
                    line=_line_for_offset(starts, i),
                    target=target,
                    is_image=is_image,
                )
        i = j


def _parse_link_target(raw_target: str) -> Optional[str]:
    if not raw_target:
        return None

    if raw_target.startswith("<"):
        end = raw_target.find(">")
        if end != -1:
            return raw_target[1:end].strip()

    # Split optional title part: [text](target "title")
    # Keep only the first token for classic markdown inline links.
    if " " in raw_target and not raw_target.startswith("#"):
        first = raw_target.split(" ", 1)[0].strip()
        if first:
            return first
    return raw_target.strip()


def _slugify_heading(heading: str) -> str:
    text = heading.strip()
    text = re.sub(r"\s+#*$", "", text)
    text = re.sub(r"`([^`]*)`", r"\1", text)
    text = re.sub(r"\[([^\]]+)\]\([^\)]+\)", r"\1", text)
    text = text.lower()
    text = re.sub(r"[^a-z0-9\s_-]", "", text)
    text = re.sub(r"[\s_]+", "-", text)
    text = re.sub(r"-+", "-", text).strip("-")
    return text


def _anchors_for_markdown(md_path: Path, cache: Dict[Path, Set[str]]) -> Set[str]:
    cached = cache.get(md_path)
    if cached is not None:
        return cached

    text = md_path.read_text(encoding="utf-8", errors="replace")
    counts: Dict[str, int] = {}
    anchors: Set[str] = set()
    for match in HEADING_RE.finditer(text):
        base = _slugify_heading(match.group(1))
        if not base:
            continue
        suffix = counts.get(base, 0)
        counts[base] = suffix + 1
        if suffix == 0:
            anchors.add(base)
        else:
            anchors.add(f"{base}-{suffix}")

    cache[md_path] = anchors
    return anchors


def _is_external(target: str) -> bool:
    scheme = urlsplit(target).scheme.lower()
    return scheme in {"http", "https"}


def _should_skip(target: str) -> bool:
    scheme = urlsplit(target).scheme.lower()
    return scheme in SKIP_SCHEMES


def _resolve_internal_target(source_file: Path, target: str) -> Tuple[Path, str]:
    parts = urlsplit(target)
    rel = unquote(parts.path)
    fragment = unquote(parts.fragment)

    if not rel:
        return source_file, fragment
    return (source_file.parent / rel).resolve(), fragment


def _check_external_link(
    session: requests.Session,
    url: str,
    timeout: float,
    retries: int,
) -> Tuple[bool, str]:
    last_error = "unknown error"
    attempts = max(1, retries + 1)
    for _ in range(attempts):
        try:
            response = session.head(url, allow_redirects=True, timeout=timeout)
            if response.status_code in ACCEPTABLE_RESPONSE_CODES:
                return True, f"HTTP {response.status_code}"
            if response.status_code in FALLBACK_TO_GET_CODES:
                response = session.get(url, allow_redirects=True, timeout=timeout)
                if response.status_code in ACCEPTABLE_RESPONSE_CODES:
                    return True, f"HTTP {response.status_code}"

            if response.status_code >= 400:
                last_error = f"HTTP {response.status_code}"
                continue

            return True, f"HTTP {response.status_code}"
        except requests.RequestException as exc:
            last_error = str(exc)

    return False, last_error


def _iter_markdown_files(root: Path) -> Iterable[Path]:
    for md_file in root.rglob("*.md"):
        if any(part in SKIP_DIR_NAMES for part in md_file.parts):
            continue
        if md_file.is_file():
            yield md_file


def _display_path(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def run(root: Path, timeout: float, retries: int, user_agent: str) -> int:
    md_files = sorted(_iter_markdown_files(root.resolve()))
    anchor_cache: Dict[Path, Set[str]] = {}
    external_cache: Dict[str, Tuple[bool, str]] = {}

    failures: List[str] = []
    total_links = 0

    session = requests.Session()
    session.headers.update({"User-Agent": user_agent})

    for md_file in md_files:
        text = md_file.read_text(encoding="utf-8", errors="replace")
        for link in _extract_markdown_links(text, md_file):
            total_links += 1
            target = link.target.strip()
            if not target or _should_skip(target):
                continue

            if _is_external(target):
                ok, reason = external_cache.get(target, (True, "cached"))
                if target not in external_cache:
                    ok, reason = _check_external_link(session, target, timeout, retries)
                    external_cache[target] = (ok, reason)
                if not ok:
                    failures.append(
                        f"{_display_path(md_file, root)}:{link.line}: {target} -> external link failed ({reason})"
                    )
                continue

            resolved_path, fragment = _resolve_internal_target(md_file, target)
            if not resolved_path.exists():
                failures.append(
                    f"{_display_path(md_file, root)}:{link.line}: {target} -> target path does not exist"
                )
                continue

            if fragment and resolved_path.suffix.lower() == ".md":
                anchors = _anchors_for_markdown(resolved_path, anchor_cache)
                fragment_slug = _slugify_heading(fragment)
                if fragment_slug not in anchors:
                    failures.append(
                        f"{_display_path(md_file, root)}:{link.line}: {target} -> anchor '#{fragment}' not found in {_display_path(resolved_path, root)}"
                    )

    print(f"Markdown files scanned: {len(md_files)}")
    print(f"Links checked: {total_links}")
    print(f"Failures: {len(failures)}")
    if failures:
        print("\nBroken/dead links:")
        for failure in failures:
            print(f"- {failure}")
        return 1
    return 0


def _parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Check links in markdown files.")
    parser.add_argument(
        "root",
        help="Repository root to scan",
    )
    parser.add_argument(
        "--timeout",
        type=float,
        default=10.0,
        help="HTTP timeout in seconds for external links (default: 10)",
    )
    parser.add_argument(
        "--retries",
        type=int,
        default=1,
        help="External link retry count after first attempt (default: 1)",
    )
    parser.add_argument(
        "--user-agent",
        default="CompDam-Markdown-Link-Checker/1.0",
        help="User-Agent header for HTTP requests",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = _parse_args(argv)
    root = Path(args.root).resolve()
    if not root.exists() or not root.is_dir():
        print(f"Error: root directory does not exist: {root}")
        return 2

    return run(
        root=root,
        timeout=args.timeout,
        retries=args.retries,
        user_agent=args.user_agent,
    )


if __name__ == "__main__":
    sys.exit(main())
