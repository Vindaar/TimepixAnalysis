import ggplotnim, seqmath, sequtils, os, sugar, strscans, strformat
import ingrid / [tos_helpers]
from arraymancer import tensor

import nimhdf5, numericalnim

import cligen

const Data2014 = "../../resources/background-rate-gold.2014+2015.dat"
const Ecol = "Energy"
const Ccol = "Counts"
const Rcol = "Rate"

type
  LogLFile = object
    name: string
    year: int
    totalTime: float
    df: DataFrame

proc scaleDset(data: Column, totalTime, factor: float): Column =
  ## scales the data in `data` according to the area of the gold region,
  ## total time and bin width. The data has to be pre binned of course!
  let area = pow(0.95 - 0.45, 2) # area of gold region!
  const bin_width = 0.2 # 0.392
  const shutter_open = 1.0 # Does not play any role, because totalDuration takes
                           # this into account already!
  let scale = factor / (totalTime * shutter_open * area * bin_width) #* 1e5
  result = toColumn data.toTensor(float).map_inline(x * scale)

proc readTime(h5f: H5FileObj): float =
  let lhGrp = h5f["/likelihood".grp_str]
  result = lhGrp.attrs["totalDuration", float]

proc scaleDset(h5f: H5FileObj, data: Column, factor: float): Column =
  result = scaleDset(data, h5f.readTime(), factor)

proc extractYear(f: string): int =
  ## We assume the filenames
  ## - always contain a year
  ## - always separated by underscores
  let fname = f.extractFilename
  var dummy: string
  if fname.scanf("$*_$i_$*", dummy, result, dummy):
    echo "found a year ", result

proc readFiles(files: seq[string]): seq[LogLFile] =
  ## reads all H5 files given and stacks them to a single
  ## DF. An additional column is added, which corresponds to
  ## the filename. That way one can later differentiate which
  ## events belong to what and decide if data is supposed to
  ## be accumulated or not
  for file in files:
    let h5f = H5open(file, "r")
    var df = h5f.readDsets(likelihoodBase(), some((3, @["energyFromCharge"])))
      .rename(f{Ecol <- "energyFromCharge"})
    let fname = file.extractFilename
    df["File"] = constantColumn(fname, df.len)
    result.add LogLFile(name: fname,
                        totalTime: readTime(h5f),
                        df: df,
                        year: extractYear(fname))
    discard h5f.close()

proc histogram(df: DataFrame): DataFrame =
  ## Calculates the histogam of the energy data in the `df` and returns
  ## a histogram of the binned data
  ## TODO: allow to do this by combining different `File` values
  let (hist, bins) = histogram(df[Ecol].toTensor(float).toRawSeq,
                               range = (0.0, 20.0), bins = 100)
  result = seqsToDf({ Ecol : bins, Ccol : concat(hist, @[0]) })

template sumIt(s: seq[typed], body: untyped): untyped =
  var res: float
  for it {.inject.} in s:
    res += body
  res

proc flatScale(files: seq[LogLFile], factor: float): DataFrame =
  var count {.global.} = 0
  var df: DataFrame
  for f in files:
    df.add f.df
  result = df.histogram()
  result = result.mutate(f{float: "CountErr" ~ sqrt(`Counts`)})
  result[Rcol] = result[Ccol].scaleDset(files.sumIt(it.totalTime), factor)
  result["RateErr"] = result["CountErr"].scaleDset(files.sumIt(it.totalTime), factor)
  result["Dataset"] = constantColumn("2017/18_" & $count, result.len)
  inc count
  result = result.mutate(f{"yMin" ~ `Rate` - `RateErr`}, f{"yMax" ~ `Rate` + `RateErr`})
  result["yMin"] = result["yMin"].toTensor(float).map_inline:
    if x >= 0.0: x
    else: 0.0
  result.drop("Counts")

proc calcIntegratedBackgroundRate(df: DataFrame, factor: float): float =
  ## returns the integrated background rate given by the Rate
  ## stored in the given DataFrame, integrated over the energy
  ## range stored in the DF (typically that should be 0 - 10 keV)
  ##
  ## It is assumed that the energy is given in keV and the rate in
  ## keV⁻¹ cm⁻² s⁻¹.
  let energies = df[Ecol].toTensor(float)
  let rate = df[Rcol].toTensor(float)
  result = simpson(rate.toRawSeq, energies.toRawSeq) / factor

proc main(files: seq[string], log = false, title = "", show2014 = false,
          separateFiles = false,
          suffix = ""
         ) =
  discard existsOrCreateDir("plots")
  let logLFiles = readFiles(files)
  let factor = if log: 1.0 else: 1e5
  var df = newDataFrame()
  if separateFiles:
    for logL in logLFiles:
      df.add flatScale(@[logL], factor)
  else:
    df =  flatScale(logLFiles, factor)
  echo df

  ## NOTE: this has to be calculated before we add 2014 data if we do, of course,
  ## because otherwise we have everything duplicated!
  let intBackRate = calcIntegratedBackgroundRate(df, factor)
  echo &"Integrated background rate: {intBackRate:.4e} cm⁻² s⁻¹"
  echo &"Integrated background rate/keV: {intBackRate / 10.0:.4e} keV⁻² cm⁻² s⁻¹"

  if show2014:
    df.drop(Ccol)
    var df2014 = toDf(readCsv(Data2014, sep = ' ', header = "#"))
      .rename(f{Rcol <- "Rate[/keV/cm²/s]"}, f{"yMax" <- "dRateUp"},
              f{"yMin" <- "dRateDown"}, f{Ecol <- "E[keV]"})
    let lastLine = seqsToDf({ Ecol : @[10.1], Rcol : @[0.0], "yMin" : @[0.0], "yMax" : @[0.0] })
    df2014.add lastLine
    if not log:
      df2014 = df2014.mutate(f{Rcol ~ 1e5 * `Rate`}, f{"yMin" ~ 1e5 * `yMin`},
                             f{"yMax" ~ 1e5 * `yMax`})
    df2014 = df2014.mutate(f{"yMin" ~ `Rate` - `yMin`}, f{"yMax" ~ `Rate` + `yMax`},
                           f{Ecol ~ `Energy` - 0.1})
    df2014["Dataset"] = constantColumn("2014/15", df2014.len)
    echo df2014
    df.add df2014
  echo df
  df = df.filter(f{c"Energy" < 12.0})
  let fnameSuffix = logLFiles.mapIt($it.year).join("_") & "_show2014_" & $show2014 & "_separate_" & $separateFiles & suffix
  let titleSuff = if show2014: " compared to 2014/15" else: ""
  let transparent = color(0, 0, 0, 0)
  let fname = &"plots/background_rate_{fnameSuffix}.pdf"
  echo "INFO: storing plot in ", fname
  ggplot(df, aes(Ecol, Rcol, fill = "Dataset")) +
    geom_histogram(stat = "identity", position = "identity", alpha = some(0.5),
                   color = some(transparent),
                   hdKind = hdOutline) +
    geom_point(binPosition = "center", position = "identity") +
    geom_errorbar(binPosition = "center",
                  aes = aes(yMin = "yMin", yMax = "yMax"),
                  errorBarKind = ebLines) +
    scale_y_continuous() +
    xlab("Energy [keV]") +
    ylab("Rate [10⁻⁵ keV⁻¹ cm⁻² s⁻¹]") +
    xlim(0, 12.0) +
    ggtitle(&"Background rate of Run 2 & 3 (2017/18){titleSuff}, {suffix}") +
    ggsave(fname, width = 800, height = 480)

when isMainModule:
  dispatch main
