# executes all tests
import shell

let dir = getCurrentDir()

var toContinue = true

template shellCheck(actions: untyped): untyped =
  if toContinue:
    let res = shellVerbose:
      actions
    toContinue = res[1] == 0
  if not toContinue:
    quit("One of the tests of TimepixAnalysis failed!")

shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/tMatchAnyRunFolder.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/tOldTosFolder.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/io/tReadFadc.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/io/tReadInGrid.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/reconstruction/tInGridGeometry.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/tZipEm.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/raw_data_manipulation/tRawDataTest.nim
shellCheck:
  nim c "-r --threads:on" ($dir)/Tests/reconstruction/tReconstruction.nim
