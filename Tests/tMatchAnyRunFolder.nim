import re
import strutils
import unittest

const names = """
data058023_1_143634434.txt
run_001310_data_001230_181017_20-10-14.txt
data019497.txt
"""

suite "Parsing of run files":
  test "parsing with simple regex":
    # The following regex should match any of these filenames
    let reg = re(r".*data.*\.txt")

    var count = 0
    for l in names.splitLines:
      if match(l, reg):
        inc count

    check count == 3

  test "parsing with current strscans impl":
    echo "TODO!"
