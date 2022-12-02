#!/usr/bin/python3

import pathlib
import os
import subprocess
from typing import Dict, List

d = pathlib.Path(__file__).parent.parent.resolve() / "src/matrices"


mod = open(d.joinpath("mod.rs"), "w")
mod.write(
    """
/// Matrix is a single alignment matrix used in scoring an alignment.
///
/// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = [[i32; 128]; 128];

/// All the modules beneath here were auto-generated.
"""
)

print("walking", d)

for dirpath, _, filenames in os.walk(d):
    for filename in sorted(filenames):
        if not filename or any(
            p in filename for p in ["README", ".tar.gz", ".rs", ".py"]
        ):
            continue

        src = pathlib.Path(dirpath) / filename
        dst_filename = filename.replace(".", "_").upper()
        dst = d.joinpath(f"{dst_filename}.rs")
        print(f"parsing {src.name} to {dst_filename}")

        # parse src
        matrix: Dict[int, Dict[int, int]] = {c: {} for c in range(128)}
        chars: List[int] = []
        for l in open(src).readlines():
            if l.startswith("#"):
                continue

            if l.startswith(" "):
                chars = [ord(c) for c in l.split()]
                for c in chars:
                    matrix[c] = {}
                continue

            cols = l.split()
            if not cols:
                continue

            char = ord(cols[0])
            vals = [int(v) for v in cols[1:]]

            for c, v in zip(chars, vals):
                if char not in matrix:
                    matrix[char] = {}
                matrix[char][c] = v

        # write to dst
        f = open(dst, "w")
        f.write(
            """
use super::Matrix;

pub static MATRIX: Matrix = [
    """
        )
        for i in range(128):
            f.write("[")
            for j in range(128):
                f.write(
                    f"{matrix[i][j] if i in matrix and j in matrix[i] else 'i32::MIN'}"
                )
                if j < 127:
                    f.write(",")
            f.write("],\n")

        f.write("];")
        f.close()

        mod.write(f"#[allow(non_snake_case)] pub mod {dst_filename};\n")


mod.close()

# autoformat the new files
print("formatting")
subprocess.Popen(["cargo", "fmt", "--all"])
