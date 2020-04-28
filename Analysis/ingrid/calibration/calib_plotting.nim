import strformat
import hashes
import strutils
import seqmath
import ggplotnim
import sequtils
import ../ingrid_types
import os

#[
Contains all routines that create plots related to InGrid calibration.
]#

proc plotGasGain*[T](charge, counts: seq[T],
                     fitX, fitY: seq[T],
                     G_fit, chiSq: float,
                     chipNumber, runNumber: int,
                     pathPrefix: string) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
  discard existsOrCreateDir(pathPrefix)
  let dfRaw = seqsToDf({ "charge / e-" : charge,
                         "counts" : counts })
  let dfFit = seqsToDf({ "charge / e-" : fitX,
                         "counts" : fitY })
  let df = bind_rows([("Polya", dfRaw),
                      ("Fit", dfFit)],
                     id = "Type")
    # filter to max 2.5e4 electrons
    .filter(fn {c"charge / e-" <= 2.5e4})
  let G = histMean(counts, charge)
  ggplot(df, aes("charge / e-", "counts")) +
    geom_histogram(data = df.filter(f{c"Type" == "Polya"}),
                   stat = "identity") +
    geom_line(data = df.filter(f{c"Type" == "Fit"}),
              color = some(parseHex("FF00FF"))) +
    ggtitle(&"Polya fit of chip {chipNumber}, run {runNumber}: " &
            &"G = {G:.1f}, G_fit = {G_fit:.1f}, χ²/dof = {chiSq:.2f}") +
    ggsave(&"{pathPrefix}/gas_gain_run_{runNumber}_chip_{chipNumber}.pdf")

## TODO: Put "Charge calibration factors vs gas gain. y errors magnified * 100"" into this module!

proc plotFeSpectrum*(feSpec: FeSpecFitData,
                     runNumber: int, chipNumber: int,
                     texts: seq[string],
                     isPixel = true,
                     pathPrefix: string) =
  discard existsOrCreateDir(pathPrefix)
  doAssert feSpec.binning == feSpec.xFit
  let df = seqsToDf({ "hist" : feSpec.hist,
                      "bins" : feSpec.binning,
                      "fit" : feSpec.yFit })
  var
    xLabel: string
    yLabel: string
    suffix: string
  if isPixel:
    xLabel = "# of pixels"
    yLabel = "counts"
  else:
    xLabel = "charge / 10^3 e¯"
    yLabel = "counts"
    suffix = "_charge"

  ggplot(df, aes("bins")) +
    geom_histogram(aes(y = "hist"), stat = "identity") +
    geom_line(aes("bins", y = "fit"),
              color = some(parseHex("FF00FF"))) +
    xlab(xlabel) +
    ylab(ylabel) +
    annotate(texts.join("\n"),
             left = 0.02,
             bottom = 0.15) +
    ggtitle(&"Fe spectrum for run: {runNumber}") +
    ggsave(&"{pathPrefix}/fe_spec_run_{runNumber}_chip_{chipNumber}{suffix}.pdf")

proc plotFeEnergyCalib*(ecData: EnergyCalibFitData,
                        runNumber: int,
                        isPixel = true,
                        pathPrefix: string) =
  discard existsOrCreateDir(pathPrefix)
  let dfEnergy = seqsToDf({ "E" : ecData.energies,
                            "H" : ecData.peaks,
                            "H_err" : ecData.peaksErr })
  let dfEFit = seqsToDf({ "E" : ecData.xFit, "H" : ecData.yFit })
  let dfEC = bind_rows(("Data", dfEnergy), ("Fit", dfEFit), id = "type")
  var
    yLabel: string
    suffix: string
  if isPixel:
    yLabel = "# of primary electrons"
  else:
    yLabel = "Total charge / 10^3 e-"
    suffix = "_charge"

  ggplot(dfEC, aes("E", "H", color = "type")) +
     geom_line(data = dfEC.filter(fn {`type` == "Fit"})) +
     geom_errorbar(data = dfEC.filter(fn {`type` == "Data"}),
                   aes = aes(x = "E", yMin = fn {`H` - `H_err`}, yMax = fn {`H` + `H_err`})) +
     geom_point(data = dfEC.filter(fn {`type` == "Data"})) +
     scale_y_continuous() +
     xlab("E / keV") +
     ylab(yLabel) +
     ggtitle(&"Detector response to X-rays of energies `E` for run: {runNumber}") +
     ggsave(&"{pathPrefix}/energy_calib_run_{runNumber}{suffix}.pdf")

proc plotGasGainVsChargeCalib*(gainVals, calib, calibErr: seq[float],
                               fitResult: FitResult,
                               pathPrefix: string) =
  # now that we have all, plot them first
  discard existsOrCreateDir(pathPrefix)
  let dfData = seqsToDf({ "Gain" : gainVals,
                          "Calib" : calib,
                          "CalibErr" : calibErr })

  # TODO: refactor the following by creating a function which takes care of
  # boilerplate in the whole file here
  let dfFit = seqsToDf({ "Gain" : fitResult.x,
                         "Calib" : fitResult.y })
  let df = bind_rows(("Data", dfData), ("Fit", dfFit), id = "Type")
  # write results to ingrid database

  # use the data to which we fit to create a hash value. That way we only overwrite the file
  # in case we plot the same data
  let fnameHash = concat(gainVals, calib).hash
  ggplot(df, aes("Gain", "Calib")) +
    geom_point(data = df.filter(f{`Type` == "Data"})) +
    geom_errorbar(data = df.filter(f{`Type` == "Data"}),
                  aes = aes(yMin = f{`Calib` - `CalibErr`},
                            yMax = f{`Calib` + `CalibErr`})) +
    geom_line(data = dfFit, color = some(parseHex("FF00FF"))) +
    annotate(&"χ²/dof = {fitResult.redChiSq:.2f}", left = 0.75, bottom = 0.1) +
    xlab("Gas gain `G`") +
    ylab("Calibration factor `a^{-1}` [1e-6 keV / e]") +
    ggtitle("Charge calibration factors vs gas gain. y errors magnified * 100") +
    ggsave(&"{pathPrefix}/gasgain_vs_calibration_charge_{fnameHash}.pdf")
