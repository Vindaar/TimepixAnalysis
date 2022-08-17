import ingrid / private / hdf5_utils
import ggplotnim, nimhdf5, sequtils
import scinim / signals
import times, os

proc calcDf(data: Tensor[float], window, order: int): DataFrame =
  let t0 = cpuTime()
  let svg = savitzkyGolayFilter(data, window, order)
  echo "Computing SVG took ", cpuTime() - t0
  result = toDf({ "time" : toSeq(0 ..< 2560), "fadc" : data,
                  "filtered" : svg })
  #  .gather(["fadc", "filtered"], "Type", "Value")

#defColumn(uint16)

proc readAndCalc(f: string, window, order: int): tuple[event: DataFrame, all: DataFrame] = # Table[colType(uint16)]] =
  let h5f = H5open(f, "r")

  let info = getFileInfo(h5f)
  doAssert info.runs.len == 1
  let run = info.runs[0]
  echo "RUn ", run
  let path = fadcDataBasename(run)
  let fadcDset = h5f[path.dset_str]
  let data = fadcDset[float].toTensor.reshape(fadcDset.shape)

  let dat = data[0, _].squeeze
  result.event = calcDf(dat, window, order)

  ggplot(result.event, aes("time")) +
    geom_point(aes = aes(y = "fadc"), size = 1.5, alpha = 0.2) +
    geom_line(aes = aes(y = "fadc")) +
    geom_line(aes = aes(y = "filtered"), size = 1.0, color = "purple") +
    ggsave("/t/test_fadc.pdf")

  # now read all aux data
  result.all = toDf({ "fallTime" : h5f.readAs(fallTimeBasename(run), int),
                      "riseTime" : h5f.readAs(riseTimeBasename(run), int),
                      "minvals" : h5f[minvalsBasename(run), float] })
  doAssert h5f.close() >= 0

proc main(files: seq[string], window = 15, order = 4) =

  var df = newDataFrame()
  var dfO = newDataFrame() #colType(uint16).newDataTable()
  for f in files:
    var (ev, all) = readAndCalc(f, window, order)
    ev = ev.mutate(f{"File" <- f.extractFileName()})
    all = all.mutate(f{"File" <- f.extractFileName()})
    df.add ev
    dfO.add all
  echo df
  echo dfO

  dfO = dfO.filter(f{float -> bool: abs(`minvals`) < 10.0})
    .gather(["fallTime", "riseTime", "minvals"], "key", "value")
  for (tup, subDf) in groups(dfO.group_by("key")):
    ggplot(subDf, aes("value", fill = "File")) +
      facet_wrap("key", scales = "free") +
      geom_histogram(bins = 100, hdKind = hdOutline, position = "identity", density = true, alpha = 0.5) +
      ggsave("/t/fadc_comparison" & $tup[0][1].toStr & ".pdf")
  df.writeCsv("/t/data_fadc.csv")



when isMainModule:
  import cligen
  dispatch main
