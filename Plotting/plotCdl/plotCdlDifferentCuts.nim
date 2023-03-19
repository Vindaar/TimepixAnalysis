import std / [os, strformat, strutils, sequtils]
import pkg / [nimhdf5, ggplotnim, seqmath]

import ingrid / [tos_helpers, ingrid_types]

#[
The intent of this tool is to investigate the different cuts we use throughout
the code base for the CDL data.
In the context of generating the CDL related fits etc (`cdl_spectrum_creation`)
we use an internal cutting proc defined in `cdl_cuts`.
For the rest of the code (application of CDL data in `likelihood`) we use the
logic from `likelihood_utils`, namely the templates based on `withCdlData` (which takes
as input the `calibration-cdl` H5 file instead of the regular reconstruction H5
file that is the input to `cdl_spectrum_creation`.
While looking at the code I noticed a likely bug in the cutting templates, hence
this script was born.
]#

var Dsets = newSeq[InGridDsetKind]()
for dset in InGridDsetKind:
  if dset notin {igInvalid, igEnergyFromPixel, igLikelihood, igNumClusters,
              igFractionInHalfRadius, igRadiusDivRmsTrans, igRadius, igBalance, igLengthDivRadius, igEventNumber}:
    Dsets.add dset

proc readRaw*[T](h5f: H5File, runNumber, chip: int, dset: string,
                 _: typedesc[T]): seq[T] =
  ## Reads the raw data of the chip & run (no cuts)
  result = h5f.readAs(recoDataChipBase(runNumber) & $chip / dset, T)

proc readXrayRefCuts*(cdlFile: string, tfKind: TargetFilterKind, dset: InGridDsetKind): seq[float] =
  ## Reads the data using the `withXrayReferenceCuts` template
  withXrayRefCuts(cdlFile, $tfKind, yr2018, igEnergyFromCharge, Dsets):
    result.add data[dset][i]

proc readLogLCuts*(cdlFile: string, tfKind: TargetFilterKind, dset: InGridDsetKind): seq[float] =
  ## Reads the data using the `withLogLFilterCuts` template
  withLogLFilterCuts(cdlFile, $tfKind, yr2018, igEnergyFromCharge, Dsets):
    result.add data[dset][i]

import ingrid / projectDefs
const filename = TpxDir / "resources/cdl_runs_2019.org"
proc readDifferentCuts(h5f: H5File, cdlFile: string, tfKind: TargetFilterKind): DataFrame =
  var dfAll = newDataFrame()
  var dsets = newSeq[string]()
  for dset in Dsets:
    var df = newDataFrame()
    let dsetStr = dset.toDset()
    for (run, grp) in tfRuns(h5f, tfKind, filename):
      let dfRaw = toDf({ dsetStr:  h5f.readRaw(run, 3, dsetStr, float),
                         "Type" : "Raw",
                         "run" : run })
      let dfCdlCut = toDf({ dsetStr : h5f.readCutCDL(run, 3, dsetStr, tfKind, float),
                            "Type" : "CDLCut",
                            "run" : run })
      df.add dfRaw
      df.add dfCdlCut
    let dfXray = toDf({ dsetStr : readXrayRefCuts(cdlFile, tfKind, dset),
                        "Type" : "XrayRef",
                        "run" : 0 })
    let dfLogL = toDf({ dsetStr : readLogLCuts(cdlFile, tfKind, dset),
                        "Type" : "LogLFilter",
                        "run" : 0 })
    df.add dfXray
    df.add dfLogL
    dsets.add dsetStr
    dfAll[dsetStr] = df[dsetStr, float]
    dfAll["run"] = df["run", int]
    dfAll["Type"] = df["Type", string]

  dfAll = dfAll.gather(dsets, "Dset", "Value")
  result = dfAll

proc plotDifferentCuts(df: DataFrame, tfKind: TargetFilterKind, plotPath: string) =
  for (tup, subDf) in groups(df.group_by("Dset")):
    let dset = tup[0][1].toStr
    echo subDf
    let dfP = subDf.filter(f{float -> bool: `Value` >= percentile(col("Value"), 3) and `Value` <= percentile(col("Value"), 97)})

    ggplot(dfP, aes("Value", fill = factor("Type"))) +
      geom_histogram(bins = 100, hdKind = hdOutline, alpha = 0.5, position = "identity", density = false) +
      ggtitle("Dataset " & $dsetStr & " for " & $tfKind & " by different cut approaches") +
      ggsave(&"{plotPath}/{$tfKind}_{dset}_histogram_by_different_cut_approaches.pdf")

proc main(fname: string, # path to CDL Reco H5 file
          cdlFile: string) = # path to `calibration-cdl` H5 file
  ## Compares the different ways to cut the CDL data, as done in `cdl_spectrum_creation` as well
  ## as in the context of the `likelihood` program.
  var h5f = H5open(fname, "r")
  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  for tfkind in TargetFilterKind:
    let df = h5f.readDifferentCuts(cdlFile, tfKind)
    plotDifferentCuts(df, tfKind, plotPath)

when isMainModule:
  import cligen
  dispatch main
