import ggplotnim, unchained, measuremancer, nimpy
import std / [strformat, json, times]

# for multiprocessing
import cligen / [procpool, mslice, osUt]

## This script uses `PyBoltz`, a (I believe semi-automatic?) port from `Magboltz` to
## `Cython`. It's pretty convenient as it's pretty simple to get up and running
## and iirc its performance is comparable.
## Thanks to `nimpy` we can easily call it from Nim (in contrast to `Magboltz` itself,
## which doesn't provide any kind of library interface I believe).
##
## This script simply simulates the Argon/Isobutane 97.7/2.3% gas mixture used in the
## Septemboard detector at CAST. As the drift field and pressure is fixed the entire
## time, the interesting variable is the temperature (as we knew the temperature
## underwent significant fluctuations). More importantly though, this also just
## serves as a reference to having the expected numbers for the drift velocity,
## longitudinal and transverse diffusion coefficients.

defUnit(V•cm⁻¹)
defUnit(mm•μs⁻¹)
type
  MagRes = object
    E: V•cm⁻¹
    T: K ## XXX: FIX CONVERSION OF C TO K
    v: Measurement[mm•μs⁻¹] ## XXX: FIX UP UNITS. PyBoltz returns mm•μs⁻¹!!!
    σT: Measurement[float] # μm²•cm⁻¹] # we currently do not support √unit :(
    σL: Measurement[float] # μm²•cm⁻¹]
    σT_σL: Measurement[float]

  FlatRes = object # flat object without type info to write to / read from procpool workers
    E, T, v, Δv, σT, ΔσT, σL, ΔσL, σT_σL, ΔσT_σL: float

## Constants matching the conditions of the CAST Septemboard detector.
const
  EDrift = 500.0.V•cm⁻¹
  P_cast = 1050.mbar.to(Torr)
  Gas    = ["ARGON", "ISOBUTANE"] # syntax used in Magbolzt
  GasFraction = [97.7, 2.3] # CAST gas mixture

proc toFlatRes(m: MagRes): FlatRes =
  result = FlatRes(E: m.E.float, T: m.T.float,
                   v: m.v.value.float, Δv: m.v.error.float,
                   σT: m.σT.value, ΔσT: m.σT.error,
                   σL: m.σL.value, ΔσL: m.σL.error,
                   σT_σL: m.σT_σL.value, ΔσT_σL: m.σT_σL.error)

proc toMagRes(f: FlatRes): MagRes =
  result = MagRes(E: f.E.V•cm⁻¹, T: f.T.Kelvin,
                  v: f.v.mm•μs⁻¹ ± f.Δv.mm•μs⁻¹,
                  σT: f.σT ± f.ΔσT,
                  σL: f.σL ± f.ΔσL,
                  σT_σL: f.σT_σL ± f.ΔσT_σL)

proc toMagRes(res: PyObject, E = EDrift, temp: Kelvin = 0.K): MagRes =
  result = MagRes(E: E, T: temp)

  let v = res["Drift_vel"].val[2].to(float)
  let Δv = res["Drift_vel"].err[2].to(float)
  result.v = v.mm•μs⁻¹ ± Δv.mm•μs⁻¹

  # now get diffusion coefficients for a single centimeter (well √cm)
  let σ_T1 = res["DT1"].val.to(float)
  let Δσ_T1 = res["DT1"].err.to(float)
  result.σT = (σ_T1 ± Δσ_T1)

  let σ_L1 = res["DL1"].val.to(float)
  let Δσ_L1 = res["DL1"].err.to(float)
  result.σL = (σ_L1 ± Δσ_L1)

  result.σT_σL = result.σT / result.σL

proc `$`(m: MagRes): string =
  result.add &"T    = {m.T}"
  result.add &"σ_T1 = {m.σT} μm·cm⁻⁰·⁵"
  result.add &"σ_L1 = {m.σL} μm·cm⁻⁰·⁵"

proc toDf(ms: seq[MagRes]): DataFrame =
  let len = ms.len
  result = newDataFrame()
  for m in ms:
    var df = newDataFrame()
    df.len = 1 # a bit of a hack, but works
    for field, data in fieldPairs(m):
      echo "Data ", data
      when typeof(data) is Measurement:
        let uof = unitOf(data.value)
        let unit = &" [{uof}]"
        df[field & unit] = data.value.float
        df["Δ" & field & unit] = data.error.float
      else:
        let uof = unitOf(data)
        let unit = &" [{uof}]"
        df[field & unit] = data.float
    echo df
    result.add df


const UseTemps = true

proc runMc(): DataFrame =
  let pb = pyImport("PyBoltz.PyBoltzRun")

  # Set up helper object
  let PBRun = pb.PyBoltzRun()

  # Configure settings for our simulation
  var Settings = %* { "Gases"                 : Gas,
                      "Fractions"             : GasFraction,
                      "Max_collisions"        : 2e9,
                      "EField_Vcm"            : EDrift.float, # gotta get rid of units :(
                      "Max_electron_energy"   : 0,
                      "Temperature_C"         : 20,
                      "Pressure_Torr"         : P_cast.float,
                      "BField_Tesla"          : 0,
                      "BField_angle"          : 0,
                      "Angular_dist_model"    : 1,
                      "Enable_penning"        : 0,
                      "Enable_thermal_motion" : 1,
                      "ConsoleOutputFlag"     : 1}

  let t0 = epochTime()
  let jobs = 28

  # We use a cligen procpool to handle running all jobs in parallel
  var pp = initProcPool((
    proc(r, w: cint) =
      let i = open(r)
      var o = open(w, fmWrite)
      var val: float
      while i.uRd(val):
        echo "Running value: ", val
        if UseTemps:
          Settings["Temperature_C"] = % val
        else:
          Settings["EField_Vcm"] = % val

        # commence the run!
        var flatRes: FlatRes
        if UseTemps:
          flatRes = PBRun.Run(Settings).toMagRes(temp = (val + 273.15).K).toFlatRes
        else:
          flatRes = PBRun.Run(Settings).toMagRes(E = val.V•cm⁻¹).toFlatRes
        let t1 = epochTime()
        echo "time taken = ", t1-t0
        echo flatRes

        discard o.uWr(flatRes)
    ),
                        framesOb,
                        jobs,
                        aux = sizeof FlatRes)

  # XXX: Make adjustable
  let temps = arange(16.0, 72.0, 2.0)
  # The electric drift fields were used with a different gas mixture (Ar/Iso 90/10) to see if we
  # can reproduce some results as a sanity check (yes, we did).
  let Es = arange(500.0, 3000.0, 100.0)
  var res = newSeq[MagRes]()
  var readRes = proc(s: MSlice) =
    var f: FlatRes
    s.toOb(f)
    res.add toMagRes(f)
  if UseTemps:
    pp.evalOb temps, readRes
  else:
    pp.evalOb Es, readRes

  let df = res.toDf()
  # finally write the DF to a CSV file
  #df.writeCsv("/tmp/test_ar_iso_90_10_diff_Es.csv")
  df.writeCsv("out/ar_iso_97_3_septemboard_cast_different_temps.csv")
  result = df # and return

proc generatePlots(df: DataFrame) =
  if UseTemps:
    ggplot(df, aes("T [K]", "σT [μm•√cm]", color = "σL [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("σT [μm•√cm]") - idx("ΔσT [μm•√cm]")},
                              yMax = f{idx("σT [μm•√cm]") + idx("ΔσT [μm•√cm]")})) +
      ggsave("out/septemboard_cast_gas_temps_σT.pdf")
    ggplot(df, aes("T [K]", "σL [μm•√cm]", color = "σT [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("σL [μm•√cm]") - idx("ΔσL [μm•√cm]")},
                              yMax = f{idx("σL [μm•√cm]") + idx("ΔσL [μm•√cm]")})) +
      ggsave("out/septemboard_cast_gas_temps_σL.pdf")
    ggplot(df, aes("T [K]", "v [mm•μs⁻¹]", color = "σT [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("v [mm•μs⁻¹]") - idx("Δv [mm•μs⁻¹]")},
                              yMax = f{idx("v [mm•μs⁻¹]") + idx("Δv [mm•μs⁻¹]")})) +
      ggsave("out/septemboard_cast_gas_temps_vel.pdf")
    ggplot(df, aes("T [K]", "σT_σL [UnitLess]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("σT_σL [UnitLess]") - idx("ΔσT_σL [UnitLess]")},
                              yMax = f{idx("σT_σL [UnitLess]") + idx("ΔσT_σL [UnitLess]")})) +
      ggsave("out/septemboard_cast_gas_temps_σT_over_σL.pdf")
  else:
    ggplot(df, aes("E [V•cm⁻¹]", "σT [μm•√cm]", color = "σL [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("σT [μm•√cm]") - idx("ΔσT [μm•√cm]")},
                              yMax = f{idx("σT [μm•√cm]") + idx("ΔσT [μm•√cm]")})) +
      ggsave("out/septemboard_cast_gas_temps_Es_σT.pdf")
    ggplot(df, aes("E [V•cm⁻¹]", "σL [μm•√cm]", color = "σT [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("σL [μm•√cm]") - idx("ΔσL [μm•√cm]")},
                              yMax = f{idx("σL [μm•√cm]") + idx("ΔσL [μm•√cm]")})) +
      ggsave("out/septemboard_cast_gas_temps_Es_σL.pdf")

    ggplot(df, aes("E [V•cm⁻¹]", "v [mm•μs⁻¹]", color = "σT [μm•√cm]")) +
      geom_point() +
      geom_errorbar(aes = aes(yMin = f{idx("v [mm•μs⁻¹]") - idx("Δv [mm•μs⁻¹]")},
                              yMax = f{idx("v [mm•μs⁻¹]") + idx("Δv [mm•μs⁻¹]")})) +
      ggsave("out/septemboard_cast_gas_temps_Es_vel.pdf")


proc main(csvInput = "",
          runPyBoltz = false) =
  ## Generates plots of the diffusion and drift properties of the gas as used in the
  ## septemboard detector at CAST conditions. If the CSV file is not available or
  ## the data should be calculted from the ground up, use `runPyBoltz`. This computes
  ## the data using `PyBoltz` (which must be available in the python installation
  ## that is currently available via your PATH).
  var df = newDataFrame()
  if runPyBoltz:
    df = runMc()
  elif csvInput.len > 0:
    df = readCsv(csvInput)

  # rename columns appropriately
  df = df.rename(f{"σT [μm•√cm]" <- "σT [UnitLess]"},
                 f{"σL [μm•√cm]" <- "σL [UnitLess]"},
                 f{"ΔσT [μm•√cm]" <- "ΔσT [UnitLess]"},
                 f{"ΔσL [μm•√cm]" <- "ΔσL [UnitLess]"})

  if UseTemps:
    echo df.arrange("T [K]").toOrgTable()
  else:
    echo df.arrange("E [V•cm⁻¹]").toOrgTable()

  df.generatePlots()

when isMainModule:
  import cligen
  dispatch main
