# Package

version       = "0.1.0"
author        = "Sebastian Schmidt"
description   = "Database containing InGrid chip information / procs to retrieve it"
license       = "MIT"
# Dependencies

requires "nim >= 0.18.1"

task koch, "build the binary":
  exec "nim c -d:release --threads:on --out:databaseTool ingridDatabase.nim"
