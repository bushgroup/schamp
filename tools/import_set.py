"""Copy a set of .raw acquisitions into an experiment directory, verified.

    uv run tools/import_set.py <source> <set_dir> [--profile NAME|PATH] [--select PATTERN ...]
    uv run tools/import_set.py --verify <set_dir>

`source` is wherever the acquisition PC left the data -- a share, a mapped drive, a USB
stick -- or a single .raw. Every file is hashed as it is read and again after it is
written, and an acquisition is renamed into place only when all of them agree. The set
directory gains MANIFEST.sha256 (sha256sum format), a seeded experiment.toml and
conditions.csv, and a README.md that records each import. Run it again after each new
acquisition: what is already there and unchanged is skipped.

A thin wrapper; `schamp.importer` is the whole of it. Exit status is 1 if a copy failed
to verify, if --verify found anything, or if an acquisition changed at the source and
--replace-changed was not given.
"""

from __future__ import annotations

import argparse
import sys

from schamp import importer


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("source", nargs="?", help="share, drive or folder holding the .raw directories")
    parser.add_argument("set_dir", help="the experiment directory to import into")
    parser.add_argument(
        "--profile",
        default=importer.DEFAULT_PROFILE,
        help="instrument profile written into a new experiment.toml (default %(default)s)",
    )
    parser.add_argument(
        "--select",
        action="append",
        metavar="PATTERN",
        help="import only acquisitions whose name matches (fnmatch, repeatable)",
    )
    parser.add_argument("--title", help="title for a new experiment.toml (default: the set's name)")
    parser.add_argument(
        "--replace-changed",
        action="store_true",
        help="re-import acquisitions whose source no longer matches the manifest",
    )
    parser.add_argument(
        "--verify", action="store_true", help="only check set_dir against its manifest"
    )
    args = parser.parse_args(argv)

    if args.verify:
        if args.source is not None:
            parser.error("--verify takes only the set directory")
        problems = importer.verify(args.set_dir)
        for problem in problems:
            print(problem)
        print(f"{len(problems)} problem(s)" if problems else "manifest verified")
        return 1 if problems else 0

    if args.source is None:
        parser.error("give a source and a set directory")
    report = importer.import_set(
        args.source,
        args.set_dir,
        profile=args.profile,
        select=args.select,
        replace_changed=args.replace_changed,
        title=args.title,
    )
    print(report.describe())
    return 0 if report.ok and not report.changed else 1


if __name__ == "__main__":
    sys.exit(main())
