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
                     xMin, xMax: float,
                     G_fit, chiSq: float,
                     chipNumber, runNumber: int,
                     pathPrefix: string,
                     gasGainInterval = none[GasGainIntervalData]()) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
  let dfRaw = toDf({ "charge / e-" : charge,
                         "counts" : counts })
  let dfFit = toDf({ "charge / e-" : fitX,
                         "counts" : fitY })
  let df = bind_rows([("Polya", dfRaw),
                      ("Fit", dfFit)],
                     id = "Type")
    # filter to max 2.5e4 electrons, more irrelevant for polya
    .filter(f{c"charge / e-" <= 2.5e4})
    .mutate(f{float -> bool: "FitRange" ~ c"charge / e-" >= xMin})
  let G = histMean(counts, charge)
  var
    intTitle: string
    suffix: string
  if gasGainInterval.isSome:
    let g = gasGainInterval.get
    intTitle = " " & $g
    suffix = toPathSuffix g

  createDir(&"/tmp/{path_prefix}/")
  when false:
    df.write_csv(&"/tmp/{path_prefix}/gas_gain_run_{runNumber}_chip_{chipNumber}{suffix}.csv")
    echo xMin
    echo df.filter(f{Value -> bool: c"Type" == %~ "Fit" and c"FitRange" == %~ false})
    echo df.filter(f{string -> bool: c"Type" == "Fit"})
  ggplot(df, aes("charge / e-", "counts")) +
    geom_histogram(data = df.filter(f{c"Type" == "Polya"}),
                   stat = "identity") +
    geom_line(data = df.filter(f{Value -> bool: c"Type" == %~ "Fit" and c"FitRange" == %~ true}),
              color = some(parseHex("FF00FF"))) +
    geom_line(data = df.filter(f{Value -> bool: c"Type" == %~ "Fit" and c"FitRange" == %~ false}),
              color = some(parseHex("FF00FF")),
              lineType = some(ltDashed)) +
              #alpha = some(0.8)) +
    margin(top = 2) +
    ggtitle(&"Polya chip {chipNumber}, run {runNumber}: " &
            &"G = {G:.1f}, G_fit = {G_fit:.1f}, χ²/dof = {chiSq:.2f}" &
            intTitle) +
    ggsave(&"{pathPrefix}/gas_gain_run_{runNumber}_chip_{chipNumber}{suffix}.pdf")

proc escapeLatex(s: string): string =
  result = s.multiReplace([("e^-", r"$e^-$"), ("\n", r"\\")])

proc plotFeSpectrum*(feSpec: FeSpecFitData,
                     runNumber: int, chipNumber: int,
                     texts: seq[string],
                     isPixel = true,
                     isFadc = false,
                     pathPrefix: string,
                     useTeX: bool) =
  discard existsOrCreateDir(pathPrefix)
  doAssert feSpec.binning == feSpec.xFit
  let df = toDf({ "hist" : feSpec.hist,
                  "bins" : feSpec.binning,
                  "fit" : feSpec.yFit })
  var
    xLabel: string
    yLabel: string
    suffix: string
    titleSuffix: string
  if isPixel:
    xLabel = "# of pixels"
    yLabel = "Counts"
  elif isFadc:
    xLabel = if not useTeX: "FADC signal U [V]" else: r"FADC signal U [$\si{V}$]"
    yLabel = "Counts"
    suffix = "_fadc"
    titleSuffix = " for FADC data"
  else:
    xLabel = if not useTeX: "Charge [10³ e¯]" else: r"Charge [$\SI{1e3}{e^-}$]"
    yLabel = "Counts"
    suffix = "_charge"

  let annot = if not useTeX: texts.join("\n")
              else: texts.join("\n").escapeLatex()
  ggplot(df, aes("bins")) +
    geom_histogram(aes(y = "hist"), stat = "identity",
                   hdKind = hdOutline) +
    geom_line(aes("bins", y = "fit"),
              color = some(parseHex("FF00FF"))) +
    xlab(xlabel) +
    ylab(ylabel) +
    annotate(annot,
             left = 0.02,
             bottom = 0.25,
             font = font(12.0, family = "monospace")) +
    ggtitle(&"Fe spectrum for run: {runNumber}{titleSuffix}") +
    ggsave(&"{pathPrefix}/fe_spec_run_{runNumber}_chip_{chipNumber}{suffix}.pdf",
           width = 600, height = 360,
           useTeX = useTeX, standalone = useTeX)

proc plotFeEnergyCalib*(ecData: EnergyCalibFitData,
                        runNumber: int,
                        isPixel = true,
                        isFadc = false,
                        pathPrefix: string,
                        useTeX: bool) =
  discard existsOrCreateDir(pathPrefix)
  let dfEnergy = toDf({ "E" : ecData.energies,
                        "H" : ecData.peaks,
                        "H_err" : ecData.peaksErr })
  let dfEFit = toDf({ "E" : ecData.xFit, "H" : ecData.yFit })
  let dfEC = bind_rows(("Data", dfEnergy), ("Fit", dfEFit), id = "type")
  var
    yLabel: string
    suffix: string
    titlePrefix: string = "Detector"
  if isPixel:
    yLabel = "# of primary electrons"
  elif isFadc:
    yLabel = if not useTeX: "FADC signal U [V]" else: r"FADC signal U [$\si{V}$]"
    suffix = "_fadc"
    titlePrefix = "FADC"
  else:
    yLabel = if not useTeX: "Total charge [10³ e⁻]" else: r"Total charge [$\SI{1e3}{e^-}$]"
    suffix = "_charge"
  let xLabel = if not useTeX: "E [keV]" else: r"E [$\si{keV}$]"

  ggplot(dfEC, aes("E", "H", color = "type")) +
     geom_line(data = dfEC.filter(fn {`type` == "Fit"})) +
     geom_errorbar(data = dfEC.filter(fn {`type` == "Data"}),
                   aes = aes(x = "E", yMin = fn {`H` - `H_err`}, yMax = fn {`H` + `H_err`})) +
     geom_point(data = dfEC.filter(fn {`type` == "Data"})) +
     scale_y_continuous() +
     xlab(xLabel) +
     ylab(yLabel) +
     margin(right = 3) +
     ggtitle(&"{titlePrefix} response to X-rays of energies `E` for run: {runNumber}") +
     ggsave(&"{pathPrefix}/energy_calib_run_{runNumber}{suffix}.pdf",
            width = 600, height = 360,
            useTeX = useTeX, standalone = useTeX)

proc plotGasGainVsChargeCalib*(gainVals, calib, calibErr: seq[float],
                               fitResult: FitResult,
                               pathPrefix: string, useTeX = false) =
  # now that we have all, plot them first
  discard existsOrCreateDir(pathPrefix)
  let dfData = toDf({ "Gain" : gainVals,
                      "Calib" : calib,
                      "CalibErr" : calibErr })

  # TODO: refactor the following by creating a function which takes care of
  # boilerplate in the whole file here
  let dfFit = toDf({ "Gain" : fitResult.x,
                     "Calib" : fitResult.y })
  let df = bind_rows(("Data", dfData), ("Fit", dfFit), id = "Type")
  # write results to ingrid database

  proc buildAnnotate(): seq[string] =
    result.add "f(x) = m·x + b"
    result.add &"m = {fitResult.pRes[1]:.2e} ± {fitResult.pErr[1]:.2e}"
    result.add &"b = {fitResult.pRes[0]:.2f} ± {fitResult.pErr[0]:.2f}"
    result.add &"χ²/dof = {fitResult.redChiSq:.2f}"
  let annotation = buildAnnotate()

  # use the data to which we fit to create a hash value. That way we only overwrite the file
  # in case we plot the same data
  let fnameHash = concat(gainVals, calib).hash
  let ylabel = if useTeX: r"Calibration factor $a⁻¹$ [$\SI{1e-6}{keV.e^-1}$]"
               else: "Calibration factor a⁻¹ [10⁻⁶ keV / e]"
  let width = getEnv("WIDTH", "800").parseFloat
  let height = getEnv("HEIGHT", "480").parseFloat
  let fontScale = getEnv("FONT_SCALE", "1.0").parseFloat
  ggplot(df, aes("Gain", "Calib")) +
    geom_point(data = df.filter(f{`Type` == "Data"})) +
    geom_errorbar(data = df.filter(f{`Type` == "Data"}),
                  aes = aes(yMin = f{`Calib` - `CalibErr`},
                            yMax = f{`Calib` + `CalibErr`})) +
    geom_line(data = dfFit, color = some(parseHex("FF00FF"))) +
    annotate(annotation.join("\n"), left = 0.65, bottom = 0.2, font = font(family = "monospace")) +
    xlab("Gas gain 'G'") +
    ylab(ylabel) +
    theme_font_scale(fontScale) +
    ggtitle("Energy calibration factors vs gas gain") +
    ggsave(&"{pathPrefix}/gasgain_vs_energy_calibration_factors_{fnameHash}.pdf",
           width = width, height = height,
           useTeX = useTeX, standalone = useTeX)

proc plotFeSpectrumInfoFacet*(pos_x, pos_y, ecc, rms_trans: seq[float],
                              hits: seq[int64],
                              runNumber: int,
                              chipNumber: int,
                              pathPrefix: string) =
  ## plots a helper overview of the data going into the Fe spectrum cut. Useful to get a
  ## look at the run
  let hitsf = hits.mapIt(it.float64)
  let df = toDf(pos_x, pos_y, ecc, rms_trans, hitsf)
  #df.write_csv("/tmp/run_305_tpa_data.csv")
  #echo "ELEMENTS ", ecc.len
  #let x_dset = h5f[(group.name / "x").dset_str]
  #let y_dset = h5f[(group.name / "y").dset_str]
  #let xdata = x_dset[special_type(uint8), uint8]
  #let ydata = y_dset[special_type(uint8), uint8]
  #for i in 0 ..< df.len:
  #  echo "I ", i
  #  let dfEv = toDf({"x" : xdata[i], "y" : ydata[i]})
  #  ggplot(dfEv, aes("x", "y")) +
  #    geom_point() + ggsave("/tmp/event_" & $i & ".pdf")
  #  copyFile("/tmp/event_" & $i & ".pdf", "/tmp/event.pdf")

  # gather and filter NaN
  let dfGath = df.gather(getKeys(df), key = "Property", value = "Value")
    .filter(f{float -> bool: classify(`Value`) notin {fcNaN, fcNegInf, fcInf} })
  ggplot(dfGath, aes("Value")) +
    facet_wrap("Property", scales = "free") +
    geom_histogram(bins = 100, position = "identity", binBy = "subset") +
    ggtitle(&"Facet overview of variables used for Fe spectrum, run {runNumber}") +
    ggsave(&"{pathPrefix}/fe_spec_facet_run_{runNumber}_chip_{chipNumber}.pdf")
