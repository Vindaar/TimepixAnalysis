import strformat
import strutils
import seqmath
import ggplotnim
import ../ingrid_types

#[
Contains all routines that create plots related to InGrid calibration.
]#

proc plotGasGain*[T](charge, counts: seq[T],
                     fitX, fitY: seq[T],
                     G_fit, chiSq: float,
                     chipNumber, runNumber: int) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
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
    ggsave(&"out/gas_gain_run_{runNumber}_chip_{chipNumber}.pdf")

## TODO: Put "Charge calibration factors vs gas gain. y errors magnified * 100"" into this module!

proc plotFeSpectrum*(feSpec: FeSpecFitData,
                     runNumber: int, chipNumber: int,
                     texts: seq[string],
                     isPixel = true) =
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
    ggsave(&"out/fe_spec_run_{runNumber}_chip_{chipNumber}{suffix}.pdf")

proc plotFeEnergyCalib*(ecData: EnergyCalibFitData,
                        runNumber: int,
                        isPixel = true) =
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
     ggtitle("Detector response to X-rays of energies `E`") +
     ggsave(&"out/energy_calib_run_{runNumber}{suffix}.pdf")
