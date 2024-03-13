import shell
import cligen
import strutils, os, sugar

type
  OptKind = enum
    okTpx3, okRaw, okReco, okEnergy, okPlot

const tpRaw = "tpx3_"
const rawPrefix = "raw_"
const recoPrefix = "reco_"

proc toName(path, prefix, suffix: string): string = path / prefix & suffix & ".h5"

proc main(fname: string, suffix = "",
          tpx3 = false, raw = false, reco = false, energy = false, plot = false,
          outpath = "") =
  var suffix = suffix
  if suffix.len == 0:
    suffix = dup(fname.extractFilename, removePrefix("data_take_"), removeSuffix(".h5"))
  let path = if outpath.len == 0: parentDir(fname) else: outpath
  let out1 = if tpx3: toName(path, tpRaw, suffix) else: fname
  let rawf = if raw: toName(path, rawPrefix, suffix) else: out1
  let recof = if reco: toName(path, recoPrefix, suffix) else: rawf
  if tpx3:
    shell:
      "readTpx3RawTest -f" ($fname) "-o" ($out1)
  if raw:
    shell:
      "raw_data_manipulation --tpx3" ($out1) "--runType calib" "--out" ($rawf) "--config raw_reco_config.toml"
  if reco:
    shell:
      reconstruction ($rawf) "--out" ($recof) "--config raw_reco_config.toml"
  if energy:
    shell:
      reconstruction ($recof) "--only_energy 26.0" "--config raw_reco_config.toml"
  if plot:
    shell:
      plotData --h5file ($recof) --runType rtCalibration --ingrid --occupancy --backend bGgPlot --config plotData.toml


when isMainModule:
  dispatch main
