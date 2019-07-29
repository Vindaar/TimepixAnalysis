import shell, sequtils, algorithm, strformat, strutils


const regions = ["crAll", "crBronze", "crSilver", "crGold"]
const vetoes = ["--fadcveto", "--scintiveto", "--septemveto"]
proc genVetoStr(vetoes: varargs[string]): string =
  for v in vetoes:
    result = result & " " & v
iterator genCombinations(vetoes: varargs[string]): seq[string] =
  yield @[]
  #yield @[vetoes[0]]
  #yield @[vetoes[1]]
  #yield @[vetoes[2]]
  #yield @[vetoes[0], vetoes[1]]
  #yield @[vetoes[0], vetoes[2]]
  #yield @[vetoes[1], vetoes[2]]
  #yield @[vetoes[0], vetoes[1], vetoes[2]]

proc buildFilename(region: string, vetoes: varargs[string]): string =
  result = &"likelihoodOut/likelihood_cdl2018_{region}"
  for v in vetoes:
    result = result & "_" & v.replace("--", "").replace("veto", "")
  result = result & ".h5"
const
  lhood = "ingrid/likelihood"
  data = "/data/CAST/2018_2/DataRuns.h5"
  cdl2018 = "--cdlYear=2018 --altCdlFile=ingrid/calibrationTest-2018.h5 --altRefFile=ingrid/xrayRefTest2018.h5"

proc runCommand(region: string, vetoes: varargs[string]) =
  let vetoStr = genVetoStr(vetoes)
  let fname = buildFilename(region, vetoes)
  let regionStr = &"--region={region}"
  shell:
    `$lhood` `$data` "--h5out" `$fname` `$regionStr` `$cdl2018` `$vetoStr`

for r in regions:
  for c in genCombinations(vetoes):
    runCommand(r, c)
