import std / [sequtils, strformat, strutils, os, sets, algorithm]
from std / times import epochTime

import shell
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

  Veto = object
    kind: FlagKind
    additional: bool # whether this veto is treated as an additional or also standalone

  Combination = object
    fname: string
    calib: string
    year: int
    region: ChipRegion
    signalEff: float # the lnL signal efficiency
    vetoes: set[FlagKind]
    vetoPercentile: float # FADC veto percentile

proc toStr(fk: FlagKind): string =
  case fk
  of fkNoVeto:   ""
  of fkScinti:   "--scintiveto"
  of fkFadc:     "--fadcveto"
  of fkSeptem:   "--septemveto"
  of fkLineVeto: "--lineveto"
  of fkExclusiveLineVeto: "--lineveto"

proc `$`(v: Veto): string =
  if v.additional:
    result = "+"
  result.add v.kind.toStr()

proc contains(s: HashSet[Veto], f: FlagKind): bool =
  for x in s:
    if x.kind == f: return true

iterator items(s: HashSet[Veto]): Veto =
  # yield elements in order of `FlagKind`
  var vSeq: seq[Veto]
  var ms = s
  while ms.len > 0:
    vSeq.add ms.pop
  vSeq = vSeq.sortedByIt(it.kind)
  for x in vSeq:
    yield x

proc genVetoStr(vetoes: set[FlagKind]): string =
  for v in vetoes:
    result = result & " " & v.toStr()

template yieldVetoes(): untyped {.dirty.} =
  for vSet in vetoSets:
    var flags: set[FlagKind]
    for veto in vSet:      # iterate over `HashSet[Veto]`, adding this veto to the current flags
      flags.incl veto.kind # combination is therefore returned
      if veto.kind == fkExclusiveLineVeto:
        # remove septem veto
        flags.excl fkSeptem
        flags.excl fkLineVeto # don't need line veto anymore
      if veto.additional: # in this case we *do not* yield, e.g. `+fkScinti` won't run with scinti only, accumulate
        continue
      var comb = Combination(fname: fname, calib: calib, year: year,
                             region: region,
                             signalEff: eff,
                             vetoes: flags,
                             vetoPercentile: -1.0)
      if fkFadc in flags: # if FADC contained, yield all percentiles after another
        for perc in fadcVetoPercentiles:
          comb.vetoPercentile = perc
          yield comb
      else:
        yield comb


iterator genCombinations(f2017, f2018: string,
                         c2017, c2018: string,
                         regions: set[ChipRegion],
                         signalEff: seq[float],
                         vetoSets: seq[HashSet[Veto]],
                         fadcVetoPercentiles: seq[float]
                        ): Combination =
  for tup in zip(@[f2017, f2018].filterIt(it.len > 0), @[c2017, c2018]):
    let (fname, calib) = (tup[0], tup[1])
    let year = if fname == f2017: 2017 else: 2018
    for region in regions:
      if signalEff.len > 0:
        for eff in signalEff:
          yieldVetoes()
      else:
        let eff = 0.0 # local zero so that it acts as none given
        yieldVetoes()

proc buildFilename(comb: Combination, outpath: string): string =
  let runPeriod = if comb.year == 2017: "Run2" else: "Run3"
  result = &"{outpath}/likelihood_cdl2018_{runPeriod}_{comb.region}"
  if comb.signalEff > 0.0:
    result.add &"_signalEff_{comb.signalEff}"
  for v in comb.vetoes:
    let vetoStr = (v.toStr).replace("--", "").replace("veto", "")
    if vetoStr.len > 0: # avoid double `_`
      result = result & "_" & vetoStr
  if comb.vetoPercentile >= 0.0:
    result = result & "_vetoPercentile_" & $comb.vetoPercentile
  result = result & ".h5"

proc runCommand(comb: Combination, cdlFile, outpath: string,
                cdlYear: int, dryRun: bool, readOnly: bool) =
  let vetoStr = genVetoStr(comb.vetoes)
  let outfile = buildFilename(comb, outpath)
  let regionStr = &"--region={comb.region}"
  let cdlYear = &"--cdlYear={cdlYear}"
  let cdlFile = &"--cdlFile={cdlFile}"
  let calibFile = if fkFadc in comb.vetoes: &"--calibFile={comb.calib}"
                  else: ""
  let vetoPerc = if comb.vetoPercentile > 0.0: &"--vetoPercentile={comb.vetoPercentile}" else: ""
  let readOnly = if readOnly: "--readOnly" else: ""
  let signalEff = if comb.signalEff > 0.0: &"--signalEfficiency={comb.signalEff}" else: ""
  let fname = comb.fname
  if not dryRun:
    let (res, err) = shellVerbose:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile) ($readOnly) ($calibFile) ($vetoPerc) ($signalEff)
    # first write log file
    let logOutput = outfile.extractFilename.replace(".h5", ".log")
    writeFile(&"{outpath}/{logOutput}", res)
    # then check error code. That way  we have the log at least!
    doAssert err == 0, "The last command returned error code: " & $err
  else:
    shellEcho:
      "likelihood -f" ($fname) "--h5out" ($outfile) ($regionStr) ($cdlYear) ($vetoStr) ($cdlFile) ($readOnly) ($calibFile) ($vetoPerc) ($signalEff)

type
  InputData = object
    fname: array[512, char] # fixed array for the data filename
    calib: array[512, char] # fixed array for filename of calibration file
    year: int
    region: ChipRegion
    signalEff: float # the lnL signal efficiency
    vetoes: set[FlagKind]
    vetoPercentile: float

proc toArray(s: string): array[512, char] = # could mem copy, but well
  doAssert s.len < 512
  for i in 0 ..< s.len:
    result[i] = s[i]

proc fromArray(ar: array[512, char]): string =
  result = newStringOfCap(512)
  for i in 0 ..< 512:
    if ar[i] == '\0': break
    result.add ar[i]

proc `$`(id: InputData): string =
  $(fname: id.fname.fromArray(), calib: id.calib.fromArray(),
    year: id.year, region: id.region,
    signalEff: id.signalEff,
    vetoes: id.vetoes,
    vetoPercentile: id.vetoPercentile)

proc toInputData(comb: Combination): InputData =
  result = InputData(fname: comb.fname.toArray(), calib: comb.calib.toArray(),
                     year: comb.year, region: comb.region,
                     signalEff: comb.signalEff,
                     vetoes: comb.vetoes,
                     vetoPercentile: comb.vetoPercentile)

proc toCombination(data: InputData): Combination =
  result = Combination(fname: data.fname.fromArray(), calib: data.calib.fromArray(),
                       year: data.year, region: data.region,
                       signalEff: data.signalEff,
                       vetoes: data.vetoes,
                       vetoPercentile: data.vetoPercentile)

proc main(f2017, f2018: string = "", # paths to the Run-2 and Run-3 data files
          c2017, c2018: string = "", # paths to the Run-2 and Run-3 calibration files (needed for FADC veto)
          regions: set[ChipRegion], # which chip regions to compute data for
          vetoSets: seq[HashSet[Veto]] = @[],
          cdlFile: string,
          outpath = "out",
          cdlYear = 2018,
          dryRun = false,
          multiprocessing = false,
          fadcVetoPercentiles: seq[float] = @[],
          signalEfficiency: seq[float] = @[],
         ) =
  if vetoSets.anyIt(fkFadc in it) and ( # stop if FADC veto used but calibration file missing
     (f2017.len > 0 and c2017.len == 0) or
     (f2018.len > 0 and c2018.len == 0)):
    doAssert false, "When using the FADC veto the corresponding calibration file to the background " &
      "data file is required."
  if not multiprocessing: # run all commands in serial
    for comb in genCombinations(f2017, f2018, c2017, c2018, regions, signalEfficiency, vetoSets, fadcVetoPercentiles):
      runCommand(comb, cdlFile, outpath, cdlYear, dryRun, readOnly = false)
  else:
    var cmds = newSeq[InputData]()
    for comb in genCombinations(f2017, f2018, c2017, c2018, regions, signalEfficiency, vetoSets, fadcVetoPercentiles):
      cmds.add comb.toInputData()

    for cmd in cmds:
      echo "Command: ", cmd
      echo "As filename: ", buildFilename(cmd.toCombination(), outpath)
    if not dryRun:
      # run them using a procpool
      let t0 = epochTime()
      let jobs = 8 # running with 28 jobs _definitely_ runs out of RAM on a machine with 64GB. 10 seems to work fine.
                    # However, most of the jobs are done very quickly anyway. The crAll (esp incl septem/line veto)
                    # are by far the slowest. So while 10 is slower than 28, the difference is small.
      ## See note at the bottom of the file.
      # We use a cligen procpool to handle running all jobs in parallel
      var pp = initProcPool((
        proc(r, w: cint) =
          let i = open(r)
          var o = open(w, fmWrite)
          var cmd: InputData
          while i.uRd(cmd):
            echo "Running value: ", cmd
            runCommand(cmd.toCombination(), cdlFile, outpath, cdlYear, dryRun, readOnly = true)
            discard w.wrLine "INFO: Finished input pair: " & $cmd
      ), framesLines, jobs)

      proc prn(m: MSlice) = echo m
      pp.evalOb cmds, prn
      echo "Running all likelihood combinations took ", epochTime() - t0, " s"

when isMainModule:
  import cligen/argcvt
  from strutils import parseEnum
  proc argParse(dst: var HashSet[Veto], dfl: HashSet[Veto],
                a: var ArgcvtParams): bool =
    var vals = a.val.strip(chars = {'{', '}'}).split(',')
    for val in vals:
      let valC = val.strip
      let additional = if valC[0] == '+': true else: false
      try:
        let kind = parseEnum[FlagKind](valC.strip(chars = {'+'}))
        dst.incl(Veto(kind: kind, additional: additional))
      except ValueError:
        raise newException(Exception, "Invalid flag given: " & $valC)
    result = true

  import cligen
  dispatch main


#[
A test run of:
```
likelihood -f ~/CastData/data/DataRuns2017_Reco.h5 --h5out /tmp/blabla.h5 --region crAll --cdlYear 2018 --cdlFile ~/CastData/data/CDL_2019/calibration-cdl-2018.h5 --lineveto --scintiveto --fadcveto --calibFile ~/CastData/data/CalibrationRuns2017_Reco.h5
```
peaked at 10.1GB of used memory according to `htop`.

Interestingly it seemed to increase only at the _end_ of each run, probably during the
writing portion of the code.

It's a bit annoying that the memory seems to increase significantly each time. As expected run
186 (the 11 day run) caused the largest spike.

Let's test by compiling `likelihood` with `-d:useMalloc` to see if we actually give memory
back to the system.

Update: This seems to have done the trick! When we reached run 186 in the initial run without
using malloc we were already at 7GB and after at over 9. This time we are at 3.3GB with a very
short peak to 4.1GB during the (likely?) writing period, but a direct drop down to about 3GB again!
So instead of increasing we are actually giving memory back.

]#
