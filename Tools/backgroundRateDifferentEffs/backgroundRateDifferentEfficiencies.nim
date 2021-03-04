import shell, strformat, strutils

const effs = [0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 0.975, 0.9875, 0.99, 0.995]
const Data2017 = "/mnt/1TB/CAST/2017/DataRuns2017_Reco.h5"
const Data2018 = "/mnt/1TB/CAST/2018_2/DataRuns2018_Reco.h5"

proc rewriteToml(path: string, eff: float) =
  ## rewrites the given TOML file in the `path` to use the `interval`
  ## instead of the existing value
  var data = readFile(path).splitLines
  for l in mitems(data):
    if l.startsWith("signalEfficiency"):
      l = "signalEfficiency = " & $eff
  writeFile(path, data.join("\n"))

proc computeLikelihood(f, outName: string, eff: float) =
  let (res, err, code) = shellVerboseErr:
    likelihood ($f) "--h5out" ($outName) "--altCdlFile /mnt/1TB/CAST/CDL_2019/calibration-cdl-2018.h5 --altRefFile /mnt/1TB/CAST/CDL_2019/XrayReferenceFile2018.h5 --cdlYear=2018"
  if code != 0:
    raise newException(Exception, "Error computing likelihood cuts for eff " & $eff)

proc plotBackgroundRate(f1, f2: string, eff: float) =
  let suffix = &"_eff_{eff}"
  let (res, err, code) = shellVerboseErr:
    one:
      cd ~/CastData/ExternCode/TimepixAnalysis/Plotting/plotBackgroundRate
      ./plotBackgroundRate ($f1) ($f2) "--suffix" ($suffix)
      ./plotBackgroundRate ($f1) ($f2) "--separateFiles --suffix" ($suffix)
  if code != 0:
    raise newException(Exception, "Error plotting background rate for eff " & $eff)

for eff in effs:
  ## rewrite toml file
  rewriteToml("/home/basti/CastData/ExternCode/TimepixAnalysis/Analysis/ingrid/config.toml", eff)
  ## compute new likelihood for both data files
  let f1 = &"out/lhood_2017_eff_{eff}.h5"
  let f2 = &"out/lhood_2018_eff_{eff}.h5"
  computeLikelihood(Data2017, f1, eff)
  computeLikelihood(Data2018, f2, eff)
  ## plot background rate
  plotBackgroundRate(f1, f2, eff)
