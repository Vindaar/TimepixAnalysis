import shell, sequtils, strformat, strutils, os
from ingrid / ingrid_types import ChipRegion

#[
Computes all combinations of the following cases
- for year in years of data files given:
  - for region in chip region:
    - no vetoes
    - +scinti veto
    - +FADC veto
    - +septem veto
    - +line veto
All veto combinations are treated additive, with the priority being the one
defined by the local `FlagKind`, but only including the ones you specify!
For that reason it is not possible to in one go generate all *individual* veto
cases at the same time as computing the combined vetoes. Maybe this will be added
in the future.

It will also write the output written by the program to the same output directory
with the same name, but a `.log` extension.
]#

##: XXX: ADD THE `--tracking` FLAG IF DESIRED

type
  ## A simplified reproduction of the `likelihood` `FlagKind` type only containing the
  ## aspects we care about with an additional `fkNoVeto`
  FlagKind = enum
    fkNoVeto   # = "",
    fkScinti   # = "--scintiveto"
    fkFadc     # = "--fadcveto"
    fkSeptem   # = "--septemveto"
    fkLineVeto # = "--lineveto"

proc toStr(fk: FlagKind): string =
  case fk
  of fkNoVeto:   ""
  of fkScinti:   "--scintiveto"
  of fkFadc:     "--fadcveto"
  of fkSeptem:   "--septemveto"
  of fkLineVeto: "--lineveto"

proc genVetoStr(vetoes: set[FlagKind]): string =
  for v in vetoes:
    result = result & " " & (v.toStr())

iterator genCombinations(f2017, f2018: string,
                         regions: set[ChipRegion],
                         vetoes: set[FlagKind]
                        ): tuple[fname: string, region: ChipRegion, vetoes: set[FlagKind]] =
  for fname in @[f2017, f2018].filterIt(it.len > 0):
    for region in regions:
      var vetoSet: set[FlagKind]
      for veto in FlagKind: # iterate over `FlagKind` checking if this veto contained in input
        if veto in vetoes:  # guarantees we return in order of `FlagKind`. Each *additional*
          vetoSet.incl veto # combination is therefore returned
        yield (fname: fname, region: region, vetoes: vetoSet)

proc buildFilename(region: ChipRegion, vetoes: set[FlagKind], outpath: string): string =
  result = &"{outpath}/likelihood_cdl2018_{region}"
  for v in vetoes:
    let vetoStr = (v.toStr).replace("--", "").replace("veto", "")
    if vetoStr.len > 0: # avoid double `_`
      result = result & "_" & vetoStr
  result = result & ".h5"
const
  lhood = "ingrid/likelihood"
  cdl2018 = "--cdlYear=2018 --altCdlFile=ingrid/calibrationTest-2018.h5 --altRefFile=ingrid/xrayRefTest2018.h5"

proc runCommand(fname: string, region: ChipRegion, vetoes: set[FlagKind],
                cdlFile, outpath: string, cdlYear: int, dryRun: bool) =
  let vetoStr = genVetoStr(vetoes)
  let outfile = buildFilename(region, vetoes, outpath)
  let regionStr = &"--region={region}"
  let cdlYear = &"--cdlYear={cdlYear}"
  let cdlFile = &"--cdlFile={cdlFile}"
  if not dryRun:
    let (res, err) = shellVerbose:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile)

    # first write log file
    let logOutput = outfile.extractFilename.replace(".h5", ".log")
    writeFile(&"{outpath}/{logOutput}", res)
    # then check error code. That way  we have the log at least!
    doAssert err == 0, "The last command returned error code: " & $err
  else:
    shellEcho:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile)


proc main(f2017, f2018: string = "", # paths to the Run-2 and Run-3 data files
          regions: set[ChipRegion], # which chip regions to compute data for
          vetoes: set[FlagKind],
          cdlFile: string,
          outpath = "out",
          cdlYear = 2018,
          dryRun = false) =
  for (fname, region, vetoes) in genCombinations(f2017, f2018, regions, vetoes):
    runCommand(fname, region, vetoes, cdlFile, outpath, cdlYear, dryRun)

when isMainModule:
  import cligen
  dispatch main
