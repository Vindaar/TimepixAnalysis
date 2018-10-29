import re
import strutils

const names = """
data058023_1_143634434.txt
run_001310_data_001230_181017_20-10-14.txt
data019497.txt
"""

# The following regex should match any of these filenames
let reg = re(r".*data.*\.txt")

var count = 0
for l in names.splitLines:
  if match(l, reg):
    inc count

doAssert count == 3
