## This tool is used to extract the rate of scintillator triggers
## based on a certain run given in a H5 file.
## What we do is read all scintillator triggers. Then filter to those
## larger than 200 and check the count, calculate rate from
## - \delta t sensitive to trigger
## - total run time of events with FADC trigger?

import ingrid / tos_helpers
import seqmath, nimhdf5, strutils, tables
import os

import ggplotnim

proc getScintiFadcDataFrame(h5f: H5File, runNumber: int): DataFrame =
  ## creates a dataframe for run number `runNumber` with the following columns
  ## - fadc triggered
  ## - scinti 1 counts
  ## - scinti 2 counts
  ## - eventDuration
  const names = ["eventDuration", "fadcReadout", "szint1ClockInt", "szint2ClockInt"]
  result = DataFrame()
  let fileInfo = h5f.getFileInfo()
  for run in fileInfo.runs:
    if runNumber >= 0 and run != runNumber: continue
    var grp = h5f[(recoBase() & $run).grp_str]
    var df = newDataFrame()
    for n in names:
      let data = h5f.readAs(grp.name / n, float)
      df[n.replace("ClockInt", "")] = data
    df["runNumber"] = run
    result.add df

proc calcScinti(df: DataFrame, scinti: string, runNumber: int,
                plotPath: string): float =
  ##  returns the rate of triggers from this scinti for random events
  ggplot(df.filter(f{idx(scinti) < 4095}),
         aes(x = scinti)) +
    geom_histogram(bins = 50) +
    ggsave(plotPath / scinti & "_full_run" & $runNumber & ".pdf")

  let nonTrivial = df.filter(f{float: idx(scinti) > 0 and
                               idx(scinti) < 4095})
  let nonMain = df.filter(f{float: idx(scinti) >= 300 and
                            idx(scinti) < 4095})
  echo nonTrivial
  echo nonMain
  echo "non trivial ", scinti, ":", nonTrivial.len
  echo "outside of main peak ", scinti, ":", nonMain.len
  echo "time of non trivial ", scinti, ":",
    nonTrivial["eventDuration", float].sum
  echo "time of non trivial non main ", scinti, ":",
    nonMain["eventDuration", float].sum
  ggplot(nonMain, aes(x = scinti)) +
    geom_histogram() +
    ggsave(plotPath / scinti & "_nonMain_run" & $runNumber & ".pdf")

  ggplot(df.filter(f{float: idx(scinti) > 0 and idx(scinti) < 300}), aes(scinti, fill = factor("runNumber"))) +
    geom_histogram(binWidth = 1.0) +
    ggsave(plotPath / scinti & "_main_run" & $runNumber & ".pdf")

  # finally a facet plot of both together
  let dfG = df.rename(f{"SiPM" <- "szint1"}, f{"Paddle" <- "szint2"})
    .gather(["SiPM", "Paddle"], "Scintillator", "Clocks")
    .filter(f{float: idx("Clocks") > 0 and idx("Clocks") < 300})
  ggplot(dfG, aes("Clocks", fill = "Scintillator")) +  #, fill = factor("Number"))) +
    facet_wrap("Scintillator", scales = "free") +
    facetMargin(0.5) +
    margin(right = 3.5, bottom = 1) +
    geom_histogram(binWidth = 1.0, position = "identity") +
    xlab("Clock cycles", margin = 0.25) +
    legendPosition(0.83, 0.0) +
    ggtitle("Clock cycles at which scintillators triggered in 2018 data") +
    ggsave(plotPath / "scintillators_facet_main_run" & $runNumber & ".pdf", width = 800, height = 480)


  let nonMainTriggers = nonMain.len
  let mainTriggers = df.filter(f{float: idx(scinti) > 0 and idx(scinti) < 300}).len
  echo "Main triggers: ", mainTriggers
  let fTriggers = df.len
  echo "FADC triggers ", fTriggers
  # calc rate by considering number of FADC triggers - main triggers
  # (since those cannot reliable be counted towards the open shutter time
  # for randoms)
  let triggers = fTriggers - mainTriggers
  echo "Relevant triggers ", triggers

  # time of open shutter
  let deltaT = triggers.float * (4094.0 - 300.0) * 25e-9
  echo "Open shutter time in s: ", deltaT
  result = nonMainTriggers.float / deltaT

proc main(fname: string, runNumber: int = -1) =

  var h5f = H5open(fname, "r")
  defer: discard h5f.close()

  let df = getScintiFadcDataFrame(h5f, runNumber)
  echo df

  # first calculate total time of run
  let runTime = df["eventDuration", float].sum
  echo "Open shutter time in hours: ", runTime / 3600.0

  let dfTriggered = df.filter(f{float: `fadcReadout` > 0})
  let runTimeFadc = dfTriggered["eventDuration", float].sum
  echo "Number of FADC triggers: ", dfTriggered.len
  echo "Open shutter time w/ FADC trigger in hours: ", runTimeFadc / 3600.0
  echo dfTriggered

  let plotPath = h5f.attrs[PlotDirPrefixAttr, string]
  let rate1 = calcScinti(dfTriggered, "szint1", runNumber, plotPath)
  let rate2 = calcScinti(dfTriggered, "szint2", runNumber, plotPath)

  echo "Scinti1 rate: ", rate1, " s^-1"
  echo "Scinti2 rate: ", rate2, " s^-1"


when isMainModule:
  import cligen
  dispatch(main, help = { "fname" : "Input H5 file name of file containing scintillator triggers",
                          "runNumber" : "Optional run number to restrict to this run number" })
