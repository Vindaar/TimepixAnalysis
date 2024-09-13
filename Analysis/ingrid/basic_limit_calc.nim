import pkg / [unchained, seqmath, ggplotnim, xrayAttenuation]
#from numericalnim import integrate
import numericalnim except linspace
import std / [sequtils, algorithm, math]

defUnit(GeV⁻¹)
defUnit(GeV⁻²)
defUnit(cm⁻²•s⁻¹)
defUnit(keV⁻¹)
defUnit(keV⁻¹•cm⁻²•s⁻¹)
defUnit(mm²)
defUnit(cm²)
defUnit(g•cm⁻³)

##[
Important note:

The coupling constant is "pulled out" of the conversion probability and solar axion
flux in this code. We leave a unit-ful value of `1 GeV⁻¹` in its place. For the axion
electron coupling it would be a unitless value in the flux of course.

Therefore the `g²` variable that appears in the calculation is a simple float without
units.
]##

type
  Config = object # Stores experimental configuration
    B*: Tesla              = 3.5.T                        ## Magnetic field of the magnet
    L*: Meter              = 20.m                         ## Length of the magnet
    totalTime*: Hour       = 100.h                        ## Total time in hours of solar tracking
    boreRadius*: cm        = 70.cm                        ## Radius of the bore (assumes a circular bore, A = πr²)
    chipArea*: mm²         = 5.mm * 5.mm                  ## Area in which all flux is assumed to be collected and in which
                                                          ## all candidates are detected & background is defined
    areaBore: cm²                                         ## Full bore area, computed from `boreRadius`
    # Detector properties
    gas*: seq[string]      = @["Ar,0.977", "C4H10,0.023"] ## Gas mixture of the detector to use.
    pressure*: mbar        = 1050.mbar                    ## Pressure in the detector chamber.
    T*: Kelvin             = 293.15.K                     ## Temperature in the detector chamber
    chamberHeight*: cm     = 3.cm                         ## Height of the detector chamber
    # Window properties
    window*: string        = "Si3N4"                      ## Compound of the detector window
    windowDensity*: g•cm⁻³ = 3.44.g•cm⁻³                  ## Density of the detector window material
    windowThickness*: μm   = 0.3.μm                       ## Thickness of the window (in μm!)

  Context = object
    cfg: Config
    fluxIntegral: cm⁻²•s⁻¹ ## the base integral of the flux over the entire energy
    # Interpolators
    background: InterpolatorType[keV⁻¹•cm⁻²•s⁻¹]
    gm: GasMixture # Gas mixture of the detector, used to evaluate absorption efficiency
    window: Compound # Window material, used to compute transmission efficiency

## Constants defining the channels and background info
const
  Energies =   @[0.5,    1.5,    2.5,    3.5,    4.5,    5.5,     6.5,    7.5,  8.5,    9.5].mapIt(it.keV)
  Background = @[0.5e-5, 2.5e-5, 4.5e-5, 4.0e-5, 1.0e-5, 0.75e-5, 0.8e-5, 3e-5, 3.5e-5, 2.0e-5]
    .mapIt(it.keV⁻¹•cm⁻²•s⁻¹) # convert to a rate
  ## A possible set of candidates from `Background · chipArea · totalTime · 1 keV`
  ## (1e-5 · 5x5mm² · 100h = 0.9 counts•keV⁻¹)
  Candidates = @[0,      2,      7,     3,      1,      0,       1,      4,    3,      2]

proc signalRate(ctx: Context, E: keV): keV⁻¹•cm⁻²•s⁻¹ # forward declare

proc initConfig(): Config =
  ## This stores all 'experimental configuration' parameter. These are all constant and
  ## don't change for one run of the program.
  ##
  ## The majority of fields are already initialized with default values.
  result = Config()
  result.areaBore = π * result.boreRadius^2

proc initContext(cfg: Config): Context =
  # 1. construct background interpolator
  let bkg = newLinear1D(Energies.mapIt(it.float), Background)
  # 2. intiialize the gas mixture & window
  let gm = parseGasMixture(cfg.gas, cfg.pressure, cfg.T)
  let window = initCompound(cfg.windowDensity, parseCompound(cfg.window))
  # 3. integrate the solar flux, combined with efficiencies
  let xs = linspace(0.0, 10.0, 1000) # keV
  result = Context(cfg: cfg,
                   background: bkg,
                   gm: gm,
                   window: window)
  # calc flux & convert back to float for compatibility with `simpson`
  let fl = xs.mapIt(signalRate(result, it.keV).float)
  result.fluxIntegral = simpson(fl, xs).cm⁻²•s⁻¹ # convert back to units (integrated out `keV⁻¹`!)

proc solarAxionFlux(E: keV): keV⁻¹•cm⁻²•s⁻¹ =
  ## Axion flux produced by the Primakoff effect in solar core
  ## in units of keV⁻¹•m⁻²•yr⁻¹
  ##
  ## This is computed *WITHOUT* `g_aγ` present, i.e. we set `g_aγ = 1 GeV⁻¹`, so that
  ## in `signal` we multiply only by `g_aγ⁴` (to simplify between `g_ae` and `g_aγ`).
  let flux = 2.0 * 1e18.keV⁻¹•m⁻²•yr⁻¹ * (1.GeV⁻¹ / 1e-12.GeV⁻¹)^2 * pow(E / 1.keV, 2.450) * exp(-0.829 * E / 1.keV)
  # convert flux to correct units
  result = flux.to(keV⁻¹•cm⁻²•s⁻¹)

func conversionProbability(ctx: Context): UnitLess =
  ## the conversion probability in the CAST magnet (depends on g_aγ)
  ## simplified vacuum conversion prob. for small masses
  ##
  ## This is computed *WITHOUT* `g_aγ` present, i.e. we set `g_aγ = 1 GeV⁻¹`, so that
  ## in `signal` we multiply only by `g_aγ⁴` (to simplify between `g_ae` and `g_aγ`).
  result = pow( (1.GeV⁻¹ * ctx.cfg.B.toNaturalUnit * ctx.cfg.L.toNaturalUnit / 2.0), 2.0 )

proc signalRate(ctx: Context, E: keV): keV⁻¹•cm⁻²•s⁻¹ =
  ## Calculates the expexted signal rate, given the incoming axion flux, conversion probability
  ## (both independent of coupling constant), telescope effective area, window transmission
  ## and detector gas absorption.
  ## As a result this is the expected rate to be detected *without* `g⁴`.
  #
  ## XXX: ADD effective area for telescope
  ## XXX: Make `solarAxionFlux` depending on config
  result = solarAxionFlux(E) * ctx.gm.transmission(ctx.cfg.chamberHeight, E) *
           ctx.window.transmission(ctx.cfg.windowThickness, E) *
           ctx.conversionProbability()

proc totalFlux(ctx: Context, g²: float): float =
  ## Flux integrated to total time, energy and area
  ##
  ## Multiply by `g²^2`. `g²` can be either `g_ae·g_aγ` or `g_aγ²`. We square that
  ## so that we multiply by the 'missing' terms in both P_aγ and the flux.
  result = ctx.fluxIntegral * ctx.cfg.totalTime * ctx.cfg.areaBore * g²^2

## NOTE: only important that signal and background have the same units!
proc signal(ctx: Context, E: keV, g²: float): keV⁻¹ =
  ## Returns the axion flux based on `g` and energy `E`
  result = ctx.signalRate(E) * ctx.cfg.totalTime.to(Second) * ctx.cfg.areaBore * g²^2

proc background(ctx: Context, E: keV): keV⁻¹ =
  ## Compute an interpolation the background at the energy `E`.
  ##
  ## Note: area of interest is the region on the chip, in which the signal is focused!
  ## This also allows us to see that the "closer" we cut to the expected axion signal on the
  ## detector, the less background we have compared to the *fixed* signal flux!
  result = (ctx.background.eval(E.float).keV⁻¹•cm⁻²•s⁻¹ * ctx.cfg.totalTime * ctx.cfg.chipArea).to(keV⁻¹)

proc likelihood(ctx: Context, g²: float, energies: seq[keV], cs: seq[int]): float =
  ## `energies` = energies corresponding to each channel
  ## `cs` = each element is number of counts in that energy channel
  result = exp(-ctx.totalFlux(g²)) # `e^{-s_tot}`
  for i in 0 ..< cs.len:
    let c = cs[i]       # number of candidates in this channel
    let E = energies[i] # energy of this channel
    let s = ctx.signal(E, g²)
    let b = ctx.background(E)
    result *= pow(1 + s / b, c.float)

proc computeLimit(ctx: Context): float =
  let xLin = linspace(0.0, 1e-20, 1000) # g²
  let yLin = xLin.mapIt(ctx.likelihood(it, Energies, Candidates))
  let yCumSum = yLin.cumSum()          # cumulative sum
  let yMax = yCumSum.max               # maximum of the cumulative sum
  let yCdf = yCumSum.mapIt(it / yMax)  # normalize to get (empirical) CDF
  let limitIdx = yCdf.lowerBound(0.95) # limit at 95% of the CDF
  echo "Limit at : ", xLin[limitIdx]
  result = xLin[limitIdx]

proc main(ctx: Context) =
  echo ctx

  ## Let's plot it from 0 to 3 assuming 4 candidates

  ## define coupling constants
  #let xs = logspace(-13, -10, 300).mapIt(it.GeV⁻¹) # logspace 1e-13 GeV⁻¹ to 1e-8 GeV⁻¹
  #let ys = xs.mapIt(likelihood(it, Energies, Candidates))
  #
  #let df = toDf({"xs" : xs.mapIt(it.float), ys})
  #ggplot(df, aes("xs", "ys")) +
  #  geom_line() +
  #  ggsave("/tmp/energy_bins_likelihood.pdf")
  #
  ### Compute limit, CDF@95%
  # limit needs non logspace x & y data! (at least if computed in this simple way)
  # Code outputs:
  # Limit at : 6.44645e-11 GeV⁻¹

when isMainModule:
  import cligen
  const ConfigPath = "limit.config"
  include mergeCfgEnvLocal
  include unchained/cligenParseUnits

  const dfl* = initConfig() # set defaults!=default for type
  var app = initFromCL(dfl, cmdName = "foo", help = {
    "B": "Magnetic field strength of the magnet",
    "L": "Magnet length in meters",
    "boreRadius" : "Radius of the magnet bore",
    "totalTime" : "Total solar tracking time",
    # Detector properties
    "gas" : "Gas mixture of the detector to use.",
    "pressure" : "Pressure in the detector chamber.",
    "T" : "Temperature in the detector chamber",
    "chamberHeight" : "Height of the detector chamber",
    # Window properties
    "window" : "Materiall of the detector window. Silicone nitried 'Si3N4', Mylar: 'C10H8O4'",
    "windowDensity" : "Density of the detector window material",
    "windowThickness" : "Thickness of the detector window, in μm (!)"
  })
  let ctx = initContext(app)
  ctx.main # Only --help/--version/parse errors cause early exit
