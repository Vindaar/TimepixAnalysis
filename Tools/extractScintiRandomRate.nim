## This tool is used to extract the rate of scintillator triggers
## based on a certain run given in a H5 file.
## What we do is read all scintillator triggers. Then filter to those
## larger than 200 and check the count, calculate rate from
## - \delta t sensitive to trigger
## - total run time of events with FADC trigger?

import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables
import os, sugar

import ggplotnim
import helpers / utils

proc getScintiFadcDataFrame(h5f: var H5FileObj, runNumber: int): DataFrame =
  ## creates a dataframe for run number `runNumber` with the following columns
  ## - fadc triggered
  ## - scinti 1 counts
  ## - scinti 2 counts
  ## - eventDuration
  var grp = h5f[(recoBase() & $runNumber).grp_str]
  const names = ["eventDuration", "fadcReadout", "szint1ClockInt", "szint2ClockInt"]
  result = DataFrame()
  for n in names:
    let data = h5f.readAs(grp.name / n, float)
    result[n.replace("ClockInt", "")] = toVector(%~ data)
    result.len = data.len

proc sumEv(s: seq[float]): float =
  result = s.foldl(a + b, 0.0)
liftVectorFloatProc(sumEv)

proc calcScinti(df: DataFrame, scinti: string, runNumber: int): float =
  ##  returns the rate of triggers from this scinti for random events
  ggplot(df.filter(f{scinti < 4095}),
         aes(x = scinti)) +
    geom_histogram() +
    ggsave(scinti & "_full_run" & $runNumber & ".pdf")

  let nonTrivial = df.filter(f{scinti > 0 and
                               scinti < 4095})
  let nonMain = df.filter(f{scinti >= 300 and
                            scinti < 4095})
  echo nonTrivial
  echo nonMain
  echo "non trivial ", scinti, ":", nonTrivial.len
  echo "outside of main peak ", scinti, ":", nonMain.len
  echo "time of non trivial ", scinti, ":",
    nonTrivial.summarize(f{"time" ~ sumEv("eventDuration")})
  echo "time of non trivial non main ", scinti, ":",
    nonMain.summarize(f{"time" ~ sumEv("eventDuration")})
  ggplot(nonMain, aes(x = scinti)) +
    geom_histogram() +
    ggsave(scinti & "_nonMain_run" & $runNumber & ".pdf")

  ggplot(df.filter(f{scinti > 0 and scinti < 300}), aes(scinti)) +
    geom_histogram() +
    ggsave(scinti & "_main_run" & $runNumber & ".pdf")

  let nonMainTriggers = nonMain.len
  let mainTriggers = df.filter(f{scinti > 0 and scinti < 300}).len
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

proc main(fname: string, runNumber: int) =
  var h5f = H5open(fname, "r")
  defer: discard h5f.close()

  let df = getScintiFadcDataFrame(h5f, runNumber)
  echo df

  # first calculate total time of run
  let runTime = df.summarize(f{"runTime" ~ sumEv("eventDuration")})["runTime"][0].toFloat
  echo "Open shutter time in hours: ", runTime / 3600.0

  let dfTriggered = df.filter(f{"fadcReadout" > 0})
  let runTimeFadcDf = dfTriggered.summarize(f{"runTimeFADC" ~ sumEv("eventDuration")})
  echo "Number of FADC triggers: ", dfTriggered.len
  let runTimeFadc = runTimeFadcDf["runTimeFADC"][0].toFloat
  echo "Open shutter time w/ FADC trigger in hours: ", runTimeFadc / 3600.0
  echo dfTriggered

  let rate1 = calcScinti(dfTriggered, "szint1", runNumber)
  let rate2 = calcScinti(dfTriggered, "szint2", runNumber)

  echo "Scinti1 rate: ", rate1, " s^-1"
  echo "Scinti2 rate: ", rate2, " s^-1"


when isMainModule:
  if paramCount() == 2:
    main(paramStr(1), paramStr(2).parseInt)
  else:
    echo "Please hand a filename and a run number!"
