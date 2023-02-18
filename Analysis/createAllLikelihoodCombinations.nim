import shell, sequtils, strformat, strutils, os
from std / times import epochTime
from ingrid / ingrid_types import ChipRegion

import cligen / [procpool, mslice, osUt]

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
    fkExclusiveLineVeto # line veto *without* septem veto & lvRegular, ecc cut 1.0
                        # important this is last!

proc toStr(fk: FlagKind): string =
  case fk
  of fkNoVeto:   ""
  of fkScinti:   "--scintiveto"
  of fkFadc:     "--fadcveto"
  of fkSeptem:   "--septemveto"
  of fkLineVeto: "--lineveto"
  of fkExclusiveLineVeto: "--lineveto"


proc genVetoStr(vetoes: set[FlagKind]): string =
  for v in vetoes:
    result = result & " " & (v.toStr())

iterator genCombinations(f2017, f2018: string,
                         regions: set[ChipRegion],
                         vetoes: set[FlagKind]
                        ): tuple[fname: string, year: int, region: ChipRegion, vetoes: set[FlagKind]] =
  for fname in @[f2017, f2018].filterIt(it.len > 0):
    let year = if fname == f2017: 2017 else: 2018
    for region in regions:
      var vetoSet: set[FlagKind]
      for veto in FlagKind: # iterate over `FlagKind` checking if this veto contained in input
        if veto in vetoes:  # guarantees we return in order of `FlagKind`. Each *additional*
          vetoSet.incl veto # combination is therefore returned
          if veto == fkExclusiveLineVeto:
            # remove septem veto
            vetoSet.excl fkSeptem
            vetoSet.excl fkLineVeto # don't need line veto anymore
        yield (fname: fname, year: year, region: region, vetoes: vetoSet)

proc buildFilename(year: int, region: ChipRegion, vetoes: set[FlagKind], outpath: string): string =
  let runPeriod = if year == 2017: "Run2" else: "Run3"
  result = &"{outpath}/likelihood_cdl2018_{runPeriod}_{region}"
  for v in vetoes:
    let vetoStr = (v.toStr).replace("--", "").replace("veto", "")
    if vetoStr.len > 0: # avoid double `_`
      result = result & "_" & vetoStr
  result = result & ".h5"
const
  lhood = "ingrid/likelihood"
  cdl2018 = "--cdlYear=2018 --altCdlFile=ingrid/calibrationTest-2018.h5 --altRefFile=ingrid/xrayRefTest2018.h5"

proc runCommand(fname: string, year: int, region: ChipRegion, vetoes: set[FlagKind],
                cdlFile, outpath: string, cdlYear: int, dryRun: bool, readOnly: bool) =
  let vetoStr = genVetoStr(vetoes)
  let outfile = buildFilename(year, region, vetoes, outpath)
  let regionStr = &"--region={region}"
  let cdlYear = &"--cdlYear={cdlYear}"
  let cdlFile = &"--cdlFile={cdlFile}"
  let readOnly = if readOnly: "--readOnly" else: ""
  if not dryRun:
    let (res, err) = shellVerbose:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile) ($readOnly)
    # first write log file
    let logOutput = outfile.extractFilename.replace(".h5", ".log")
    writeFile(&"{outpath}/{logOutput}", res)
    # then check error code. That way  we have the log at least!
    doAssert err == 0, "The last command returned error code: " & $err
  else:
    shellEcho:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile) ($readOnly)

type
  InputData = object
    fname: array[512, char] # fixed array for the data filename
    year: int
    region: ChipRegion
    vetoes: set[FlagKind]

proc toArray(s: string): array[512, char] = # could mem copy, but well
  doAssert s.len < 512
  for i in 0 ..< s.len:
    result[i] = s[i]

proc fromArray(ar: array[512, char]): string =
  result = newStringOfCap(512)
  for i in 0 ..< 512:
    if ar[i] == '\0': break
    result.add ar[i]

proc `$`(id: InputData): string = $(fname: id.fname.fromArray(), year: id.year, region: id.region, vetoes: id.vetoes)

proc main(f2017, f2018: string = "", # paths to the Run-2 and Run-3 data files
          regions: set[ChipRegion], # which chip regions to compute data for
          vetoes: set[FlagKind],
          cdlFile: string,
          outpath = "out",
          cdlYear = 2018,
          dryRun = false,
          multiprocessing = false) =
  if not multiprocessing: # run all commands in serial
    for (fname, year, region, vetoes) in genCombinations(f2017, f2018, regions, vetoes):
      runCommand(fname, year, region, vetoes, cdlFile, outpath, cdlYear, dryRun, readOnly = false)
  else:
    var cmds = newSeq[InputData]()
    for (fname, year, region, vetoes) in genCombinations(f2017, f2018, regions, vetoes):
      cmds.add InputData(fname: fname.toArray(),
                         year: year,
                         region: region,
                         vetoes: vetoes)

    for cmd in cmds:
      echo "Command: ", cmd
      echo "As filename: ", buildFilename(cmd.year, cmd.region, cmd.vetoes, outpath)
    if not dryRun:
      # run them using a procpool
      let t0 = epochTime()
      let jobs = 28
      # We use a cligen procpool to handle running all jobs in parallel
      var pp = initProcPool((
        proc(r, w: cint) =
          let i = open(r)
          var o = open(w, fmWrite)
          var cmd: InputData
          while i.uRd(cmd):
            echo "Running value: ", cmd
            runCommand(cmd.fname.fromArray(), cmd.year, cmd.region, cmd.vetoes, cdlFile, outpath, cdlYear, dryRun, readOnly = true)
            discard w.wrLine "INFO: Finished input pair: " & $cmd
      ), framesLines, jobs)

      proc prn(m: MSlice) = echo m
      pp.evalOb cmds, prn
      echo "Running all likelihood combinations took ", epochTime() - t0, " s"

when isMainModule:
  import cligen
  dispatch main
