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
    result[n] = toVector(%~ data)
    result.len = data.len

proc sumEv(s: seq[float]): float =
  result = s.foldl(a + b, 0.0)
liftVectorFloatProc(sumEv)

proc main(fname: string, runNumber: int) =
  var h5f = H5file(fname, "r")
  defer: discard h5f.close()

  let df = getScintiFadcDataFrame(h5f, runNumber)



  # first calculate total time of run
  let runTime = df.summarize(f{"runTime" ~ sumEv("eventDuration")})["runTime"][0].toFloat
  echo "Open shutter time in hours: ", runTime / 3600.0

  let dfTriggered = df.filter(f{"fadcReadout" > 0})
  let runTimeFadcDf = dfTriggered.summarize(f{"runTimeFADC" ~ sumEv("eventDuration")})
  echo "Number of FADC triggers: ", dfTriggered.len
  let runTimeFadc = runTimeFadcDf["runTimeFADC"][0].toFloat
  echo "Open shutter time w/ FADC trigger in hours: ", runTimeFadc / 3600.0
  echo dfTriggered
  echo dfTriggered.summarize(f{"test" ~ max("szint2ClockInt")})
  ggplot(dfTriggered, aes(x = "szint1ClockInt")) +
    geom_histogram() +
    geom_histogram(aes("szint2ClockInt")) +
    ggsave("scinti.pdf")

  let nonTrivialSzint1 = dfTriggered.filter(f{"szint1ClockInt" > 0 and
                                              "szint1ClockInt" < 4095})
  let nonMainSzint1 = dfTriggered.filter(f{"szint1ClockInt" > 200 and
                                           "szint1ClockInt" < 4095})
  let nonTrivialSzint2 = dfTriggered.filter(f{"szint2ClockInt" > 0 and
                                              "szint2ClockInt" < 4095})
  let nonMainSzint2 = dfTriggered.filter(f{"szint2ClockInt" > 200 and
                                           "szint2ClockInt" < 4095})
  echo nonTrivialSzint1
  echo nonMainSzint1
  echo nonTrivialSzint2
  echo nonMainSzint2

  echo "non trivial szint 1: ", nonTrivialSzint1.len
  echo "outside of main peak 1: ", nonMainSzint1.len
  echo "non trivial szint 2: ", nonTrivialSzint2.len
  echo "outside of main peak 2: ", nonMainSzint2.len

  echo "time of non trivial 1 ",
    nonTrivialSzint1.summarize(f{"time" ~ sumEv("eventDuration")})
  echo "time of non trivial non main 1 ",
    nonMainSzint1.summarize(f{"time" ~ sumEv("eventDuration")})
  echo "time of non trivial 2 ",
    nonTrivialSzint2.summarize(f{"time" ~ sumEv("eventDuration")})
  echo "time of non trivial non main 2 ",
    nonMainSzint2.summarize(f{"time" ~ sumEv("eventDuration")})

when isMainModule:
  if paramCount() == 2:
    main(paramStr(1), paramStr(2).parseInt)
  else:
    echo "Please hand a filename and a run number!"
