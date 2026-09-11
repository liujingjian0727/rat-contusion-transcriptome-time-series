#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from pathlib import Path
import csv
import hashlib
import sys


PACKAGE = Path("submission_package_scientific_data")
MANIFEST_DIR = PACKAGE / "06_Manifests"

MANIFEST = MANIFEST_DIR / "submission_package_manifest.tsv"
SHA_FILE = MANIFEST_DIR / "SHA256SUMS.txt"


EXCLUDE_NAMES = {
    MANIFEST.name,
    SHA_FILE.name,
}

EXCLUDE_SUFFIXES = {
    ".bak",
    ".tmp",
    ".swp",
}

EXCLUDE_BASENAMES = {
    ".DS_Store",
    "Thumbs.db",
}


def sha256_file(path, chunk_size=1024 * 1024):
    h = hashlib.sha256()

    with path.open("rb") as f:
        while True:
            chunk = f.read(chunk_size)

            if not chunk:
                break

            h.update(chunk)

    return h.hexdigest()


def category_from_path(rel):
    parts = rel.parts

    if not parts:
        return "Unknown"

    return parts[0]


def main():

    if not PACKAGE.is_dir():
        print(
            f"ERROR: package directory not found: {PACKAGE}",
            file=sys.stderr
        )
        return 1

    MANIFEST_DIR.mkdir(
        parents=True,
        exist_ok=True
    )

    # Remove old self-referential files before rebuilding.
    if MANIFEST.exists():
        MANIFEST.unlink()

    if SHA_FILE.exists():
        SHA_FILE.unlink()

    files = []

    forbidden = []

    for p in sorted(PACKAGE.rglob("*")):

        if not p.is_file():
            continue

        rel = p.relative_to(PACKAGE)

        if p.name in EXCLUDE_NAMES:
            continue

        if p.name in EXCLUDE_BASENAMES:
            forbidden.append(str(rel))
            continue

        if p.suffix.lower() in EXCLUDE_SUFFIXES:
            forbidden.append(str(rel))
            continue

        if p.name.startswith("~$"):
            forbidden.append(str(rel))
            continue

        files.append(p)

    if forbidden:
        print("=" * 78)
        print("ERROR: temporary/backup files found inside submission package")
        print("=" * 78)

        for x in forbidden:
            print(" ", x)

        print()
        print("Move or remove these files before rebuilding the package.")
        return 1

    records = []

    print("=" * 78)
    print("REFRESHING SCIENTIFIC DATA SUBMISSION MANIFEST")
    print("=" * 78)

    for i, p in enumerate(files, 1):

        rel = p.relative_to(PACKAGE)
        digest = sha256_file(p)

        records.append({
            "relative_path": rel.as_posix(),
            "category": category_from_path(rel),
            "filename": p.name,
            "size_bytes": p.stat().st_size,
            "sha256": digest,
        })

        if i % 25 == 0 or i == len(files):
            print(f"Hashed {i}/{len(files)} files")

    # --------------------------------------------------------
    # Write manifest
    # --------------------------------------------------------

    with MANIFEST.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        writer = csv.DictWriter(
            f,
            delimiter="\t",
            fieldnames=[
                "relative_path",
                "category",
                "filename",
                "size_bytes",
                "sha256",
            ],
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(records)

    # --------------------------------------------------------
    # Write GNU-style SHA256 file
    # --------------------------------------------------------

    with SHA_FILE.open(
        "w",
        encoding="utf-8",
        newline=""
    ) as f:

        for r in records:
            f.write(
                f"{r['sha256']}  {r['relative_path']}\n"
            )

    # --------------------------------------------------------
    # Summary by package section
    # --------------------------------------------------------

    counts = {}

    for r in records:
        counts[r["category"]] = (
            counts.get(r["category"], 0) + 1
        )

    print()
    print("=" * 78)
    print("PACKAGE SUMMARY")
    print("=" * 78)

    for k in sorted(counts):
        print(f"{k:<35} {counts[k]:>6} files")

    print("-" * 78)
    print(f"{'TOTAL HASHED FILES':<35} {len(records):>6}")
    print()

    print("Manifest :", MANIFEST)
    print("SHA256   :", SHA_FILE)

    print()
    print("STATUS: PASS")
    print("=" * 78)

    return 0


if __name__ == "__main__":
    sys.exit(main())
