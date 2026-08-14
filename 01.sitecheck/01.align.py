"""Generate MACSE alignment commands for the target sequence sets."""

import argparse
import shlex
from pathlib import Path


def write_macse_commands(root_dir, macse_jar, output_script):
    """Write one MACSE command for each subdirectory in ``root_dir``."""
    root = Path(root_dir).expanduser().resolve()
    jar = Path(macse_jar).expanduser().resolve()
    output = Path(output_script).expanduser().resolve()

    with output.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("#!/usr/bin/env bash\n")
        handle.write("set -euo pipefail\n\n")

        for target_dir in sorted(path for path in root.iterdir() if path.is_dir()):
            input_file = target_dir / "target.fa"
            command = (
                f"cd {shlex.quote(str(target_dir))} && "
                f"java -jar {shlex.quote(str(jar))} "
                "-prog alignSequences "
                f"-seq {shlex.quote(str(input_file.name))} "
                "-out_NT target.aligned_1.fa "
                "-out_AA target.aligned_2.fa "
                "-gc_def 1"
            )
            handle.write(command + "\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Generate MACSE v2.06 alignment commands for target sequence sets."
    )
    parser.add_argument(
        "root_dir",
        help="Directory containing one subdirectory per target sequence set.",
    )
    parser.add_argument(
        "--macse-jar",
        required=True,
        help="Path to the MACSE v2.06 JAR file.",
    )
    parser.add_argument(
        "--output",
        default="01.align.sh",
        help="Output shell script (default: 01.align.sh).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    write_macse_commands(args.root_dir, args.macse_jar, args.output)
