## this module contains some TimepixAnalysis project related definitions
## and helpers, which allow compile time reflection on the whole project
from macros import getProjectPath
from parseutils import parseUntil
from os import `/`

const OldTosRunlistFname = "Runlist-CAST-D03-W0063.csv"
const AppDir = getProjectPath()
var TpxDir* {.compileTime.} = ""
var OldTosRunList* {.compileTime.} = ""
static:
  discard parseUntil(AppDir, TpxDir, "TimepixAnalysis")
  TpxDir = TpxDir / "TimepixAnalysis"
  OldTosRunList = TpxDir / "resources" / OldTosRunListFname

const OldTosRunListPath* = OldTosRunList
