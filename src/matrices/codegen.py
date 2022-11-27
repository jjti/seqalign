#!/usr/bin/python3

import pathlib
import os
import subprocess
from typing import Dict, List

d = pathlib.Path(__file__).parent.resolve()


mod = open(d.joinpath("mod.rs"), "w")
mod.write(
    """
use std::collections::HashMap;

/// Matrix is a single alignment matrix used in scoring an alignment.
///
/// Maps a char to another char and the corresponding substitution penalty.
pub type Matrix = HashMap<u8, HashMap<u8, i32>>;

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
        matrix: Dict[int, Dict[int, int]] = {}
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
use std::collections::HashMap;

lazy_static! {
    pub static ref MATRIX: Matrix = HashMap::from([
    """
        )
        for left, r in matrix.items():
            f.write(f"({left}, HashMap::from([")
            for k, v in r.items():
                f.write(f"({k},{v}),")
            f.write("])),\n")

        f.write("]); }")
        f.close()

        mod.write(f"#[allow(non_snake_case)] pub mod {dst_filename};\n")


mod.close()

# autoformat the new files
print("formatting")
subprocess.Popen(["cargo", "fmt", "--all"])
