import std / [os, strutils, strformat, sequtils, algorithm, math]
import pkg / [mpfit, zero_functional, seqmath, chroma, ggplotnim, unchained]
import ingrid / [ingrid_types, tos_helpers]
import ingrid / calibration / calib_fitting
import ingrid / calibration

import helpers/utils
import ingridDatabase / [databaseDefinitions, databaseRead]

template canImport(x: untyped): bool =
  compiles:
    import x

when canImport(ingridDatabase):
  # only in this case give option to use plot from database
  import ingridDatabase

when defined(linux):
  const commitHash = staticExec("git rev-parse --short HEAD")
else:
  const commitHash = ""

# get date using `CompileDate` magic
const currentDate = CompileDate & " at " & CompileTime
const docTmpl = """
Version: $# built on: $#
"""
const doc = docTmpl % [commitHash, currentDate]

const
  StopToT = 450.0

const
  GoldenMean = (sqrt(5.0) - 1.0) / 2.0  # Aesthetic ratio
  FigWidth = 1200.0                     # width in inches
  FigHeight = FigWidth * GoldenMean     # height in inches

func parseDbChipArg(dbArg: string): string =
  ## given the argument to the `--db` flag, perform the necessary
  ## conversion to get a valid chip name from the InGrid database
  ## NOTE: If the argument is a valid integer, it is assumed that
  ## the chip is on the Septem H board!
  # first try to parse the argument as an integer
  try:
    let chipNum = dbArg.parseInt
    result = getSeptemHChip(chipNum)
  except ValueError:
    # not a valid chip number. Interpret as chip name
    result = dbArg

proc scurveToDf(s: SCurve): DataFrame =
  result = toDf({ "THL" : s.thl.asType(float),
                  "Counts" : s.hits.asType(float),
                  "Voltage" : s.voltage,
                  "Type" : "Data" })

proc plotSCurves*(df: DataFrame, annotation: string, runPeriod, chip = "",
                  legendX = -1.0, legendY = -1.0,
                  useTeX = true, outpath = "out") =
  let dfData = df.filter(f{`Type` == "Data"})
  let dfFit = df.filter(f{`Type` == "Fit"})
  let lX = if legendX > 0.0: legendX else: 250.0
  let lY = if legendY > 0.0: legendY else: 3900.0
  ggplot(df, aes("THL", "Counts", color = factor("Voltage"), shape = "Type")) +
    geom_line(data = dfData) +
    geom_line(data = dfFit, aes = aes("THL", "Counts"), color = "black", lineType = ltDashed) +
    ggtitle(&"SCurves for run period '{runPeriod}', chip {chip}") +
    annotate(annotation, x = lX, y = lY, font = font(10.0, family = "monospace"),
             backgroundColor = color(0.0, 0.0, 0.0, 0.0)) +
    #xlab(r"$U_\text{injected} / \si{mV}$") +
    #ylab("Counts [\#]") +
    #theme_latex() +
    ggsave(&"{outpath}/s_curves_{chip}.pdf", width = 600, height = 450, useTex = useTeX, standalone = true)

import measuremancer
proc createThlAnnotation*(res: FitResult, charge, thl, thlErr: seq[float]): string =
  ## Inverts the fit parameters given for the THL calibration. We perform the fit as
  ## x = Charge, y = THL, yErr = THL_Err
  ## but plot it as x = THL and y = Charge. Therefore we should invert the fit parameters
  ## to show on the plot. Otherwise it's confusing.
  let m = res.pRes[1] ± res.pErr[1]
  let b = res.pRes[0] ± res.pErr[0]
  result.add &"m = {m} THL/e⁻\n"
  result.add &"1/m = {1.0/m} e⁻/THL\n"
  let fy0 = charge[0] - (thl[0] ± thlErr[0]) * (1.0 / m)
  result.add &"f(0) = {b} e⁻\n"
  result.add &"f⁻¹(y=0) = {fy0} THL\n"
  result.add &"χ²/dof = {res.redChiSq:.2f}"

proc plotThlCalib*(thlCalib: FitResult, charge, thl, thlErr: seq[float], chip = "",
                   legendX = -1.0, legendY = -1.0,
                   useTeX = true, outpath = "out") =
  # flip the plot, i.e. show THL on x instead of y as done for the fit to
  # include the errors
  let dfData = toDf({"THL" : thl, "Charge [e⁻]": charge, "thlError" : thlErr, "Type" : "Data"})
  let dfFit = toDf({"THL" : thlCalib.y, "Charge [e⁻]" : thlCalib.x, "thlError" : 0.0, "Type" : "Fit"})
  var df = newDataFrame()
  df.add dfData
  df.add dfFit
  let lX = if legendX > 0.0: legendX else: 440.0
  let lY = if legendY > 0.0: legendY else: 3900.0
  let annot = createThlAnnotation(thlCalib, charge, thl, thlErr)
  ggplot(dfData, aes("THL", "Charge [e⁻]")) +
    geom_line(data = dfFit, color = parseHex"FF00FF") +
    geom_point() +
    geom_errorbar(aes = aes(
                    xMin = f{`THL` - `thlError`},
                    xMax = f{`THL` + `thlError`}
                  )) +
    annotate(annot, x = lX, y = lY, font = font(10.0, family = "monospace"),
             backgroundColor = color(0,0,0,0)) +
    ggtitle(&"THL calibration of chip {chip}") +
    ggsave(&"{outpath}/thl_calibration_chip_{chip}.pdf", width = 600, height = 450,
            useTeX = useTeX, standalone = true)

proc plotToTCalib*(totCalib: FitResult, tot: Tot, chip = 0, chipName = "",
                   useTeX = false) =
  let dfData = toDf({ "U / mV" : tot.pulses.mapIt(it.float),
                      "ToT" : tot.mean,
                      "std" : tot.std.mapIt(if classify(it) == fcNaN: 0.0 else: it) })
  let dfFit = toDf({"U / mV" : totCalib.x, "ToT" : totCalib.y})
  let df = bind_rows([("ToT", dfData), ("Fit", dfFit)], "by")
  var title = ""
  if chipName.len > 0:
    title = &"ToT calibration of Chip {chipName}"
  else:
    title = &"ToT calibration of Chip {chip}"

  ggplot(dfData, aes("U / mV", "ToT")) +
    geom_point() +
    geom_errorbar(aes = aes(yMin = f{`ToT` - `std`}, yMax = f{`ToT` + `std`})) +
    geom_line(data = dfFit, color = some(color(1.0, 0.0, 1.0))) +
    annotate(totCalib.resText, x = 50.0, y = 150.0, font = font(10.0, family = "monospace"),
             backgroundColor = color(0.0, 0.0, 0.0, 0.0)) +
    xlab(r"$U_\text{injected} / \si{mV}$") +
    ylab("ToT / clock cycles") +
    ggtitle(title) +
    #theme_latex() +
    ggsave(&"{outpath}/tot_calib_{chip}.pdf", width = 600, height = 450, useTex = useTeX, standalone = true)

proc plotCharge*(a, b, c, t: float, capacitance: FemtoFarad, chip: int, chipName = "",
                 useTeX = false, outpath = "out") =
  let tots = linspace(0.0, 150.0, 1000)

  let charges = tots.mapIt(calibrateCharge(it, capacitance, a, b, c, t))
  let df = toDf({ "ToT" : tots,
                  "Q" : charges })
  let title = &"Charge calibration of Chip {chipName}"
  ggplot(df, aes("ToT", "Q")) +
    geom_line() +
    xlab(r"ToT [clock cycles]") +
    ylab(r"Charge [e⁻]") +
    ggtitle(title) +
    #theme_latex() +
    ggsave(&"{outpath}/charge_calib_{chip}.pdf", width = 600, height = 450, useTex = true, standalone = true)

iterator sCurves(file, folder, chip, runPeriod: string): SCurve =
  ## yields the traces of the correct argument given to
  ## the program for SCurves
  if chip.len > 0:
    when declared(ingridDatabase):
      let chipName = parseDbChipArg(chip)
      let scurves = getScurveSeq(chipName, runPeriod)
      for curve in scurves.curves:
        yield curve
    else:
      discard
  elif file.len > 0:
    let curve = readScurveVoltageFile(file)
    yield curve
  elif folder.len > 0:
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      let curve = readScurveVoltageFile(f)
      yield curve

proc sCurve(file, folder, chip, runPeriod: string,
            legendX, legendY: float, useTeX: bool,
            outpath: string) =
  ## perform plotting and fitting of SCurves
  var
    voltages: set[int16]
    curve: SCurve

    # calibration seqs
    charge: seq[float] = @[]
    thlMean: seq[float] = @[]
    thlErr: seq[float] = @[]
    # annotations
    annotation = ""


  var dfScurve = newDataFrame()
  for curve in sCurves(file, folder, chip, runPeriod):
    #traces.add trace
    voltages.incl int16(curve.voltage)
    if curve.voltage > 0:
      let df = scurveToDf(curve) #getTrace(curve.thl.asType(float), curve.hits, $curve.voltage)
      # now fit the normal cdf to the SCurve and add its data
      let fitRes = fitSCurve(curve)
      let dfFit = toDf({ "THL" : fitRes.x,
                         "Counts" : fitRes.y,
                         "Voltage" : curve.voltage,
                         "Type" : "Fit" })
      dfScurve.add df
      dfScurve.add dfFit
      annotation.add &"fit {pretty(curve.voltage.mV, 3, short = true):>6} χ²/dof {fitRes.redChiSq:.2f}\n"

      # given `fitRes`, add fit parameters to calibration seqs
      # multiply by 50 to take into account capacitance of test pulse injection
      # to get number of electrons
      ## TODO: replace charge calc by correct calculation using capacitance
      charge.add (charge(getCapacitance(Timepix1), curve.voltage.mV)).float
      thlMean.add (fitRes.pRes[1])
      thlErr.add (fitRes.pErr[1])

  if dfScurve.len > 0:
    plotSCurves(dfScurve, annotation, runPeriod, chip, legendX, legendY, useTeX, outpath)

  # now fit the calibration function
  # however, the voltage data is in a wrong order. need to sort it, before
  # we can fit
  let thlTup = zip(thlMean, thlErr)
  var chThl = zip(charge, thlTup)
  let sortedChThl = chThl.sortedByIt(it[0])
  let
    chSort = sortedChThl.mapIt(it[0])
    thlSort = sortedChThl.mapIt(it[1][0])
    # increase errors artifically by factor 100... Mostly for visualization,
    # but also as a rough guess for typical deviation visible
    thlErrSort = sortedChThl.mapIt(it[1][1] * 100)
  if chSort.len > 0:
    let thlCalib = fitThlCalib(chSort, thlSort, thlErrSort)
    plotThlCalib(thlCalib, chSort, thlSort, thlErrSort, chip, legendX, legendY, useTeX, outpath)

proc parseTotInput(file, folder, chip, runPeriod: string,
                   startTot = 0.0): (int, string, Tot) =
  var tot: Tot
  var chipName = ""
  var chipNum = -1
  if chip.len > 0:
    when declared(ingridDatabase):
      chipName = parseDbChipArg(chip)
      tot = getTotCalib(chipName, runPeriod)
      chipNum = getChipNumber(chipName, runPeriod)
  elif file.len > 0:
    (chipNum, tot) = readToTFile(file, startTot)
  else:
    for f in walkFiles(folder.expandTilde & "/*.txt"):
      # TODO: implement multiple in same plot?
      (chipNum, tot) = readToTFile(f, startTot)

  result = (chipNum, chipName, tot)

proc totCalib(file, folder, chip, runPeriod: string, startFit, startTot: float,
              useTeX: bool, outpath: string) =
  ## perform plotting and analysis of ToT calibration
  var
    bins: seq[float]
    hist: seq[float]

  let (chip, chipName, tot) = parseTotInput(file, folder, chip, runPeriod, startTot)
  let totCalib = fitToTCalib(tot, startFit)
  plotToTCalib(totCalib, tot, chip, chipName, useTeX, outpath)

proc chargeCalib(chip, runPeriod: string, useTeX: bool, outpath: string) =
  ## perform plotting and analysis of charge calibration (inverse of ToT)
  var chipName = ""
  if chip.len > 0:
    when declared(ingridDatabase):
      chipName = parseDbChipArg(chip)
  else:
    quit("Need a chip from the database!")

  let (a, b, c, t) = getTotCalibParameters(chipName, runPeriod)
  let capacitance = getTimepixVersion(chipName, runPeriod).getCapacitance()
  plotCharge(a, b, c, t, capacitance, chip.parseInt, chipName, useTeX, outpath)

import cligen / macUt
proc main(
  scurve = false, tot = false, charge = false,
  chip = "", runPeriod = "", file = "", folder = "", startFit = 0.0, startTot = 0.0,
  legendX = -1.0, legendY = -1.0,
  version = false,
  useTeX = false,
  outpath = "out") =
  docCommentAdd(doc)
  ## A simple tool to plot SCurves or ToT calibrations.

  if chip.len > 0:
    when not declared(ingridDatabase):
      quit("Cannot import InGrid database. --chip option not supported.")

  if scurve:
    sCurve(file, folder, chip, runPeriod, legendX, legendY, useTeX, outpath)
  elif tot:
    totCalib(file, folder, chip, runPeriod, startFit, startTot, useTeX, outpath)
  elif charge:
    chargeCalib(chip, runPeriod, useTeX, outpath)

when isMainModule:
  import cligen
  dispatch(main, help = {
    "scurve" : "If set, perform SCurve analysis",
    "tot" :    "If set, perform ToT calibration analysis",
    "charge" : "If set plot charge calibration (inverse of ToT)",
    "chip" : """If given will read information from InGrid database, if
  available. Either chip number or chip name supported.
  NOTE: if a chip number is used, it is assumed that it
  corresponds to said chip on the Septem H board!""",
    "runPeriod" : "Required if a chip name is given. The run period we want to plot the calibrations for.",
    "file" : "If given will read from a single file",
    "folder" : "If given will read all voltage files from the given folder",
    "chip" : "The number of this chip",
    "startFit" : """Start the TOT fit from this measurement. If differs from
  StartToT constant in source code (or the --startTot value),
  the other datapoints will still be plotted.""",
    "startTot" : "Read the ToT file from this pulse height",
    "useTeX" : "If set will use the TikZ backend of ggplotnim to generate the plot",
    "legendX" : "Used to move the legend left side to X",
    "legendY" : "Used to move the legend left side to X",
    "version" : "Show the version number"})
