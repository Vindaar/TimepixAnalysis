import shell, cligen, parsetoml
import std / [strutils, os, sugar, sequtils]

type
  OptKind = enum
    okTpx3, okRaw, okReco, okEnergy, okPlot

  Config = object
    tpxPrefix: string
    rawPrefix: string
    recoPrefix: string
    rawRecoConfig: string
    outpath: string    # outpath for the currently processed file
    fileSuffix: string # suffix for the currently processed file
    plotDataSuffix: string # suffix controlling which plots to create etc

const SourcePath = currentSourcePath().parentDir

proc readConfig(path: string, outpath = "", rawRecoConfig = "", plotDataSuffix = ""): Config =
  let configPath = if path.len == 0:
                     SourcePath / "config_tpx3.toml"
                   else: path
  echo "Reading config: ", configPath
  let cfg = parseToml.parseFile(configPath)
  let outpath = if outpath.len > 0: outpath else: cfg["General"]["outpath"].getStr
  let recoCfg = if rawRecoConfig.len > 0: rawRecoConfig else: cfg["General"]["rawRecoConfigFile"].getStr
  let plotSuffix = if plotDataSuffix.len > 0: plotDataSuffix
                   else: cfg["General"]["plotDataSuffix"].getStr
  result = Config(tpxPrefix: cfg["General"]["tpx3Prefix"].getStr,
                  rawPrefix: cfg["General"]["rawPrefix"].getStr,
                  recoPrefix: cfg["General"]["recoPrefix"].getStr,
                  outpath: outpath,
                  rawRecoConfig: recoCfg,
                  plotDataSuffix: plotSuffix)

proc setDefaults(cfg: var Config, fname, suffix: string) =
  ## Updates the outpath & suffix fields to defaults based on the input filename,
  ## if not overridden.
  cfg.fileSuffix = if suffix.len > 0: suffix
                   else: dup(fname.extractFilename, removePrefix("DataTake_"), removeSuffix(".h5"))
  if cfg.outpath.len == 0:
    cfg.outpath = parentDir(fname)

proc toName(path, prefix, suffix: string): string = path / prefix & suffix & ".h5"

proc main(
  fnames: seq[string], names: seq[string] = @[],
  tpx3 = false, raw = false, reco = false, energy = false, plot = false,
  outpath = "", config = "", runType = "rtCalibration",
  rawRecoConfig = "",
  plotDataSuffix = ""
     ) =
  ## Performs data parsing, reconstruction and plotting of the input Timepix3 files.
  ##
  ## If more than one file is given, please also hand a name using the `names` argument
  ## for each input file to easier disambiguate between the inputs.
  ##
  ## By default all steps of the analysis are performed one after another. Individual
  ## steps can also be done separately (but then it's up to you to make sure the required
  ## steps are done / the required files exist according to the suffixes given in the
  ## config file!).
  ##
  ## The equivalent of no flags is `--tpx3 --raw --reco --energy --plot`, which
  ## first performs parsing of the raw Tpx3 data (equivalent to the existing Python
  ## code). By default this generates a new `tpx3_<foo>.h5` output file.
  ## From there `--raw` performs raw data parsing, which splits the clusters by
  ## ToA (using the `raw_data_manipulation` command) and outputs a `raw_<foo>.h5` file.
  ## Next is `--reco` via the `reconstruction`  command, which performs (geometric) cluster
  ## finding and computation of geometric cluster properties and generates a `reco_<foo>.h5`
  ## output. The `--energy` step simply adds energy information based on the number of pixels
  ## times 26eV. Finally `--plot` generates all plos by calling `plotData`. The command handed
  ## to the `plotData` tool is configured by the `plotDataSuffix` or in the config file using
  ## the same field.
  ##
  ## The individual steps should be considered in case either something went wrong on some
  ## step and the previous steps are to be avoided on a succeeding run or in the most likely
  ## case the plots should be regenerated with a different command. Check the commands
  ## available for `plotData` by running `plotData -h` to get an idea (or ask me).
  ##
  ## If the input is a (or multiple) ⁵⁵Fe run(s), the default `runType` of `rtCalibration` is
  ## correct. If a background run is to be analyzed, hand `--runType rtBackground`.
  ##
  ## `outpath` (either via CL argument or in config file) controls the path where the H5
  ## files will be stored.
  var cfg = readConfig(config, outpath, rawRecoConfig, plotDataSuffix)
  if fnames.len > 1 and names.len != fnames.len:
    echo "Error: Please give one name for each input file using the `--names` argument!"
    echo "\tGiven input files: ", fnames
    echo "\tGiven names: ", names
    return
  else:
    template walkFiles(inPrefix, outPrefix: string, body: untyped): untyped =
      for i in 0 ..< fnames.len:
        let fname {.inject.} = fnames[i]
        if names.len > 0:
          cfg.setDefaults(fname, names[i])
        let inName {.inject.} = toName(cfg.outpath, inPrefix, cfg.fileSuffix)
        let outName {.inject.} = toName(cfg.outpath, outPrefix, cfg.fileSuffix)
        body

    # If all arguments are `false`, set `all` flag
    let all = not tpx3 and not raw and not reco and not energy and not plot

    if tpx3 or all:
      walkFiles("", cfg.tpxPrefix):
        shell:
          parse_raw_tpx3 -f ($fname) "-o" ($outname)
    let cfgPath = if cfg.rawRecoConfig.len > 0: "--config " & cfg.rawRecoConfig
                  else: ""
    if raw or all:
      walkFiles(cfg.tpxPrefix, cfg.rawPrefix):
        shell:
          raw_data_manipulation -p ($inName) --runType ($runType) --tpx3 "--out" ($outName) ($cfgPath)
    if reco or all:
      walkFiles(cfg.rawPrefix, cfg.recoPrefix):
        shell:
          reconstruction ($inName) "--out" ($outName) ($cfgPath)
    if energy or all:
      walkFiles(cfg.recoPrefix, ""):
        shell:
          reconstruction ($inName) "--only_energy 26.0" ($cfgPath)
    if plot or all:
      var files = newSeq[string]()
      walkFiles(cfg.recoPrefix, ""):
        files.add inName
      let inf = files[0]
      let cmps = files[1 .. ^1].mapIt("--h5Compare " & $it).join(" ")
      let plotParams = cfg.plotDataSuffix
      shell:
        plotData --h5file ($inf) ($cmps) --runType ($runType) ($plotParams)


when isMainModule:
  dispatch(main, help = {
    "fnames" : "The paths to the raw Tpx3 files to be processed",
    "names" : "The names to be associated with each file. Assign the same name to " &
      "multiple files to combine them into a single output.",
    "tpx3" : "If set performs parsing of raw data.",
    "raw" : "If set performs ToA clustering and converting to TPA format.",
    "reco" : "If set performs geomtric clustering and calculation of geometric properties.",
    "energy" : "If set performs calculation of the energy of each cluster based on # pixels.",
    "plot" : "If set generates all plots as described desired via `plotDataSuffix`.",

    "outpath" : "The path in which all generated H5 files will be stored.",
    "config" : "The path to the config file for this tool.",
    "runType" : "The run type of the input files {rtCalibration : ⁵⁵Fe, rtBackground : background data}.",
    "plotDataSuffix" : "The command handed to `plotData` to control which plots are generated."})
