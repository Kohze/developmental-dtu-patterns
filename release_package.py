"""Build deterministic, screened manuscript release archives.

Each paper supplies a ``release_config.json`` with an explicit allow-list.  The
builder refuses public packaging when the paper license or a configured source-
rights confirmation is missing, unless the corresponding audit-only override
is supplied.  ZIP entry order and timestamps are fixed so identical content
produces an identical archive hash.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys
import zipfile
from pathlib import Path


FIXED_ZIP_TIME = (2026, 1, 1, 0, 0, 0)
TEXT_SUFFIXES = {
    ".bib", ".cff", ".csv", ".json", ".md", ".ps1", ".py", ".r",
    ".tex", ".tsv", ".txt", ".yml", ".yaml",
}
FORBIDDEN_PATH = re.compile(
    r"(^|/)(?:\.miktex(?:-tmp)?|__pycache__|work)(?:/|$)|"
    r"\.(?:aux|bbl|blg|log|out|pyc|rdata)$",
    re.IGNORECASE,
)
SECRET_PATTERNS = {
    "OpenAI-style secret": re.compile(r"sk-[A-Za-z0-9_-]{12,}"),
    "GitHub token": re.compile(r"gh[opsu]_[A-Za-z0-9]{20,}"),
    "AWS access key": re.compile(r"AKIA[A-Z0-9]{16}"),
    "private key": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH )?PRIVATE KEY-----"),
    "absolute user path": re.compile(r"[A-Za-z]:\\Users\\[^\\\s]+", re.IGNORECASE),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def safe_relative(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    try:
        return resolved.relative_to(root.resolve())
    except ValueError as exc:
        raise ValueError(f"Release source escapes the paper directory: {path}") from exc


def collect_sources(source: Path, config: dict) -> list[tuple[Path, Path]]:
    selected: dict[str, tuple[Path, Path]] = {}
    for relative_text in config.get("include_files", []):
        relative = Path(relative_text)
        path = source / relative
        if not path.is_file():
            raise FileNotFoundError(f"Required release file is missing: {relative_text}")
        selected[relative.as_posix()] = (path, relative)

    for pattern in config.get("include_globs", []):
        matches = [path for path in source.glob(pattern) if path.is_file()]
        if not matches:
            raise FileNotFoundError(f"Release glob matched no files: {pattern}")
        for path in matches:
            relative = safe_relative(path, source)
            selected[relative.as_posix()] = (path, relative)

    for pattern in config.get("optional_globs", []):
        for path in source.glob(pattern):
            if path.is_file():
                relative = safe_relative(path, source)
                selected[relative.as_posix()] = (path, relative)

    return [selected[key] for key in sorted(selected)]


def collect_components(
    source: Path,
    config: dict,
    allow_missing_license: bool,
    allow_unconfirmed_rights: bool,
) -> tuple[list[tuple[Path, Path]], list[dict]]:
    """Collect namespaced sources from sibling release configurations.

    Composite releases keep each companion package under a distinct directory
    and apply that package's existing allow-list and rights gates. Component
    roots must be relative siblings/descendants within the composite source's
    parent directory; arbitrary absolute filesystem paths are rejected.
    """

    sources: list[tuple[Path, Path]] = []
    reports: list[dict] = []
    names: set[str] = set()
    portfolio_root = source.parent.resolve()

    for component in config.get("components", []):
        name = component.get("name", "")
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", name):
            raise ValueError(f"Invalid component name: {name!r}")
        if name in names:
            raise ValueError(f"Duplicate component name: {name}")
        names.add(name)

        source_text = component.get("source", "")
        component_relative = Path(source_text)
        if not source_text or component_relative.is_absolute():
            raise ValueError(f"Component source must be relative: {source_text!r}")
        component_source = (source / component_relative).resolve()
        try:
            component_source.relative_to(portfolio_root)
        except ValueError as exc:
            raise ValueError(
                f"Component source escapes the portfolio directory: {source_text}"
            ) from exc
        if not component_source.is_dir():
            raise FileNotFoundError(f"Component source is missing: {source_text}")

        config_name = component.get("config", "release_config.json")
        config_relative = Path(config_name)
        if config_relative.is_absolute() or ".." in config_relative.parts:
            raise ValueError(f"Invalid component config path: {config_name!r}")
        component_config_path = component_source / config_relative
        if not component_config_path.is_file():
            raise FileNotFoundError(
                f"Component release config is missing: {component_config_path}"
            )
        component_config = json.loads(
            component_config_path.read_text(encoding="utf-8")
        )
        if component_config.get("components"):
            raise ValueError(f"Nested composite components are not supported: {name}")

        license_path = component_source / "LICENSE"
        if not license_path.is_file() and not allow_missing_license:
            raise FileNotFoundError(
                f"Component {name} requires LICENSE; the override is audit-only."
            )
        rights_gate = component_config.get("rights_gate_file")
        rights_confirmed = not rights_gate or (component_source / rights_gate).is_file()
        if not rights_confirmed and not allow_unconfirmed_rights:
            raise FileNotFoundError(
                f"Component {name} requires {rights_gate} before deposit; "
                "the override is audit-only."
            )

        prefix = Path(name)
        component_sources = collect_sources(component_source, component_config)
        component_sources.append((component_config_path, prefix / "release_config.json"))
        if license_path.is_file():
            component_sources.append((license_path, prefix / "LICENSE"))
        if rights_gate and rights_confirmed:
            component_sources.append(
                (component_source / rights_gate, prefix / rights_gate)
            )

        for path, relative in component_sources:
            if relative.parts and relative.parts[0] == name:
                target = relative
            else:
                target = prefix / relative
            sources.append((path, target))

        reports.append({
            "name": name,
            "source": source_text,
            "license_present": license_path.is_file(),
            "source_rights_confirmed": rights_confirmed,
        })

    return sources, reports


def screen_sources(sources: list[tuple[Path, Path]]) -> None:
    problems: list[str] = []
    for path, relative in sources:
        rel = relative.as_posix()
        if FORBIDDEN_PATH.search(rel):
            problems.append(f"forbidden path: {rel}")
            continue
        if path.suffix.lower() not in TEXT_SUFFIXES and path.name not in {
            ".gitignore", "LICENSE"
        }:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        for label, pattern in SECRET_PATTERNS.items():
            if pattern.search(text):
                problems.append(f"{label} in {rel}")
    if problems:
        raise RuntimeError("Release screening failed:\n- " + "\n- ".join(problems))


def write_manifest(destination: Path) -> Path:
    manifest = destination / "RELEASE_CONTENTS.tsv"
    files = sorted(
        (
            path for path in destination.rglob("*")
            if path.is_file() and path != manifest
        ),
        key=lambda path: path.relative_to(destination).as_posix(),
    )
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["relative_path", "bytes", "sha256"])
        for path in files:
            writer.writerow([
                path.relative_to(destination).as_posix(),
                path.stat().st_size,
                sha256(path),
            ])
    return manifest


def deterministic_zip(directory: Path, archive_path: Path) -> None:
    with zipfile.ZipFile(
        archive_path, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9
    ) as archive:
        files = sorted(
            (item for item in directory.rglob("*") if item.is_file()),
            key=lambda path: path.relative_to(directory).as_posix(),
        )
        for path in files:
            relative = path.relative_to(directory).as_posix()
            info = zipfile.ZipInfo(relative, date_time=FIXED_ZIP_TIME)
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o100644 << 16
            archive.writestr(info, path.read_bytes(), compress_type=zipfile.ZIP_DEFLATED)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--destination", type=Path, required=True)
    parser.add_argument("--allow-missing-license", action="store_true")
    parser.add_argument("--allow-unconfirmed-rights", action="store_true")
    parser.add_argument("--skip-zip", action="store_true")
    args = parser.parse_args()

    config_path = args.config.resolve()
    source = config_path.parent
    destination = args.destination.resolve()
    if destination.exists():
        raise FileExistsError(f"Release destination already exists: {destination}")
    if destination == source or source in destination.parents:
        raise ValueError("Release destination must be outside the paper directory.")

    config = json.loads(config_path.read_text(encoding="utf-8"))
    license_path = source / "LICENSE"
    if not license_path.is_file() and not args.allow_missing_license:
        raise FileNotFoundError(
            "LICENSE is required for a depositable release; the override is audit-only."
        )
    rights_gate = config.get("rights_gate_file")
    rights_confirmed = not rights_gate or (source / rights_gate).is_file()
    if not rights_confirmed and not args.allow_unconfirmed_rights:
        raise FileNotFoundError(
            f"{rights_gate} is required before record-level outputs may be deposited; "
            "the override is audit-only."
        )

    sources = collect_sources(source, config)
    component_sources, component_reports = collect_components(
        source,
        config,
        allow_missing_license=args.allow_missing_license,
        allow_unconfirmed_rights=args.allow_unconfirmed_rights,
    )
    sources.extend(component_sources)
    if license_path.is_file():
        sources.append((license_path, Path("LICENSE")))
    screen_sources(sources)

    destinations: dict[str, Path] = {}
    for path, relative in sources:
        key = relative.as_posix()
        previous = destinations.get(key)
        if previous is not None and previous.resolve() != path.resolve():
            raise RuntimeError(
                f"Multiple release sources target {key}: {previous} and {path}"
            )
        destinations[key] = path

    destination.mkdir(parents=True)
    for path, relative in sources:
        target = destination / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    shutil.copy2(Path(__file__).resolve(), destination / "release_package.py")
    shutil.copy2(config_path, destination / "release_config.json")
    manifest = write_manifest(destination)

    archive_path = destination.parent / f"{destination.name}.zip"
    if archive_path.exists():
        raise FileExistsError(f"Release archive already exists: {archive_path}")
    if not args.skip_zip:
        deterministic_zip(destination, archive_path)

    audit_only = (
        (not license_path.is_file())
        or (not rights_confirmed)
        or any(
            (not item["license_present"]) or (not item["source_rights_confirmed"])
            for item in component_reports
        )
    )
    report = {
        "destination": str(destination),
        "files": sum(1 for path in destination.rglob("*") if path.is_file()),
        "bytes": sum(path.stat().st_size for path in destination.rglob("*") if path.is_file()),
        "manifest_sha256": sha256(manifest),
        "zip": None if args.skip_zip else str(archive_path),
        "zip_sha256": None if args.skip_zip else sha256(archive_path),
        "license_present": license_path.is_file(),
        "source_rights_confirmed": rights_confirmed,
        "audit_only": audit_only,
        "components": component_reports,
        "python": sys.version.split()[0],
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
