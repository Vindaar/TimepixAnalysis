import std / [strformat, strutils, math]
import ggplotnim, seqmath
import numericalnim except linspace
import xrayAttenuation

from std / os import getEnv
let UseTeX = getEnv("USE_TEX", "false").parseBool
let Newline = if UseTeX: r"\\" else: ""

proc readLLNL(fname: string, sep = ','): DataFrame =
  const areaBore = π * (2.15.cm)^2
  var df = readCsv(fname, sep = sep)
  proc renameIf(df: DataFrame, f, t: string): DataFrame =
    result = df
    if f in df:
      result = df.rename(f{t <- f})
  # maybe rename some columns
  echo df
  df = df.renameIf("Energy[keV]", "Energy [keV]")
  df = df.renameIf("E(keV)", "Energy [keV]")
  df = df.renameIf("Area(cm^2)", "Effective Area [cm²]")
  df = df.renameIf("EffectiveArea[cm²]", "Effective Area [cm²]")
  # compute efficiency
  df = df.mutate(f{float -> float: "Efficiency" ~ idx("Effective Area [cm²]").cm² / areaBore})
  result = df

proc interpLLNL(llnl: DataFrame, energies: seq[float]): DataFrame =
  ## interpolates the LLNL telescope efficiencies to the given energies
  let interp = newLinear1D(llnl["Energy [keV]", float].toSeq1D,
                           llnl["Efficiency", float].toSeq1D)
  var effs = newSeq[float](energies.len)
  let eMin = llnl["Energy [keV]", float].min
  let eMax = llnl["Energy [keV]", float].max
  for idx in 0 ..< effs.len:
    effs[idx] = if energies[idx] < eMin:
                  0.0
                elif energies[idx] > eMax:
                  0.0
                else:
                  interp.eval(energies[idx])
  result = toDf({"Energy [keV]" : energies, "LLNL" : effs})

proc plotEfficiency(df: DataFrame, outpath: string) =
  block:
    let df = df
      .gather(["300nm SiN", "Efficiency", "full Eff.", "Eff • SB • ε", "30mm Ar Abs.",
               "200μm Si", "20nm Al", "LLNL", "Eff • ε", "Eff • ε • LLNL"],
              key = "Type", value = "Efficiency")
    ggplot(df, aes("Energy [keV]", "Efficiency", color = "Type")) +
      geom_line() +
      ggtitle("Full detector efficiencies, including window, gas, ε, window SB, telescope") +
      margin(top = 1.75) +
      theme_font_scale(1.0, family = "serif") +
      ggsave(&"{outpath}/window_plus_argon_efficiency.pdf", width = 600, height = 380)
  block:
    echo df
    let df = df.drop(["200μm Si", "SB", "Eff • ε • LLNL", "Eff • ε", "Eff • SB • ε", "full Eff."])
      .rename(f{"Combined" <- "Efficiency"})
      .gather(["300nm SiN", "30mm Ar Abs.", "20nm Al", "LLNL", "Combined"],
              key = "Type", value = "Efficiency")
    echo df
    ggplot(df, aes("Energy [keV]", "Efficiency", color = "Type")) +
      geom_line() +
      ggtitle(&"Detection efficiencies of window, software eff., LLNL efficiency{Newline} and Argon absorption") +
      margin(top = 1.75) +
      themeLatex(fWidth = 0.9, width = 600, baseTheme = singlePlot, useTeX = UseTeX, useWithoutTeX = false) +
      ggsave(&"{outpath}/detection_efficiency.pdf", width = 600, height = 380)

proc castGas(): GasMixture =
  let arC = compound((Ar, 1)) # need Argon gas as a Compound
  let isobutane = compound((C, 4), (H, 10))
  # define the gas mixture
  result = initGasMixture(293.K, 1050.mbar, [(arC, 0.977), (isobutane, 0.023)])

proc calculateEfficiencies(energies: seq[float], gas: GasMixture, Si₃N₄: Compound;
                           Al: Aluminium; Si: Silicon,
                           ε = 0.8): DataFrame =
  let cols = ["300nm SiN", "200μm Si", "30mm Ar Abs.", "20nm Al"]
  result = newDataFrame()
  for c in cols:
    result[c] = newColumn(colFloat, energies.len)
  for i, E in energies:
    result[cols[0], i] = Si₃N₄.transmission(300.nm, E.keV)
    result[cols[1], i] = Si.transmission(200.μm, E.keV)
    result[cols[2], i] = 1.0 - gas.transmission(3.cm, E.keV) # we want absorption for gas!
    result[cols[3], i] = Al.transmission(20.nm, E.keV)
  result["SB"] = 1.0 - 0.222 # coverage of SB in gold region
  result["ε"] = ε
  result["Energy [keV]"] = energies

proc combineEfficiencies(df, llnl: DataFrame): DataFrame =
  result = df
  result["LLNL"] = llnl["LLNL"]
  result = result.mutate(f{"DetectorEfficiency" ~ idx("30mm Ar Abs.") * idx("300nm SiN") * idx("20nm Al")},
                         f{"Efficiency" ~ `DetectorEfficiency` * `LLNL`},
                         f{"Eff • SB • ε" ~ `DetectorEfficiency` * `SB` * `ε`},
                         f{"Eff • ε" ~ `DetectorEfficiency` * `ε`},
                         f{"Eff • ε • LLNL" ~ `Efficiency` * `ε`},
                         f{"full Eff." ~ idx("Eff • SB • ε") * `LLNL`}) # strongback occlusion of 22% and ε = 80%

proc main(llnlEff: string = "~/org/resources/llnl_xray_telescope_cast_effective_area_extended.csv",
          sep = ',',
          plotPath = "~/org/Figs/statusAndProgress/detector",
          outpath = "~/org/resources/detector/",
          ε = 0.8,
          showBrowser = false
         ) = # desired software efficiency to optionall include
  ## Calculates the detection efficiency of the septemboard detector.
  ## Includes the
  ## - 300 nm SiN window with 20 nm Al coating
  ## - (optional) 500 μm Si strongback
  ## - Ar/Isobutane 97.7/2.3 gas mixture X-ray absorption
  ## - Effective area of the LLNL telescope (via `llnlEff`)
  ##   Currently all our attempts to compute the effective area ourselves
  ##   yields numbers that are larger than the existing numbers.
  ## - (optional) software efficiency
  # 1. define CAST gas mixture
  let gas = castGas()
  # 2. SiN window
  let Si₃N₄ = compound((Si, 3), (N, 4))
  let Al = Aluminium.init(ρ = 2.7.g•cm⁻³)
  # 3. Silicon strongback
  let Si = Silicon.init(ρ = 2.336.g•cm⁻³)
  # 4. define common energies to use for all
  let energies = linspace(0.0, 10.0, 1000)
  # 5. read LLNL DF and interpolate
  let llnl = readLLNL(llnlEff, sep).interpLLNL(energies)
  # 6. calculate efficiencies of all detector components
  var df = calculateEfficiencies(energies, gas, Si₃N₄, Al, Si, ε)
  # 7. finish DF
  df = combineEfficiencies(df, llnl)
  # 8. write CSV file
  df.writeCsv(&"{outpath}/combined_detector_efficiencies.csv")
  if showBrowser:
    showBrowser(df, "df_initial.html")
  # 9. plot them!
  df.plotEfficiency(plotPath)

when isMainModule:
  import cligen
  dispatch main
