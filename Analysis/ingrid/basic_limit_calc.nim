import unchained, math, seqmath, ggplotnim
#from numericalnim import integrate
import numericalnim except linspace
import std / [sequtils, algorithm]
## Assumptions:
const totalTime = 100.0.h # 100 of "tracking time"
const areaBore = π * (2.15 * 2.15).cm²
const chipArea = 5.mm * 5.mm # assume all flux is focused into an area of 5x5 mm²
                             # on the detector. Relevant area for background!

defUnit(GeV⁻¹)
defUnit(GeV⁻²)
defUnit(cm⁻²•s⁻¹)
defUnit(keV⁻¹)
defUnit(keV⁻¹•cm⁻²•s⁻¹)
defUnit(mm²)
defUnit(cm²)

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
    B*: Tesla ## Magnetic field of the magnet
    L*: Meter ## Length of the magnet
    totalTime*: Hour ## Total time in hours of solar tracking
    boreRadius*: cm ## Radius of the bore (assumes a circular bore, A = πr²)
    chipArea*: mm² ## Area in which all flux is assumed to be collected and in which
                   ## all candidates are detected & background is defined
    areaBore: cm²  ## Full bore area, computed from `boreRadius`

  Context = object
    cfg: Config
    fluxIntegral: cm⁻²•s⁻¹ ## the base integral of the flux over the entire energy
    # Interpolators
    background: InterpolatorType[keV⁻¹•cm⁻²•s⁻¹]

## Constants defining the channels and background info
const
  Energies =   @[0.5,    1.5,    2.5,    3.5,    4.5,    5.5,     6.5,    7.5,  8.5,    9.5].mapIt(it.keV)
  Background = @[0.5e-5, 2.5e-5, 4.5e-5, 4.0e-5, 1.0e-5, 0.75e-5, 0.8e-5, 3e-5, 3.5e-5, 2.0e-5]
    .mapIt(it.keV⁻¹•cm⁻²•s⁻¹) # convert to a rate
  ## A possible set of candidates from `Background · chipArea · totalTime · 1 keV`
  ## (1e-5 · 5x5mm² · 100h = 0.9 counts•keV⁻¹)
  Candidates = @[0,      2,      7,     3,      1,      0,       1,      4,    3,      2]

proc solarAxionFlux(ω: keV): keV⁻¹•cm⁻²•s⁻¹ # forward declare

proc initConfig(): Config =
  ## This stores all 'experimental configuration' parameter. These are all constant and
  ## don't change for one run of the program.
  let bR = 70.cm
  result = Config(B: 3.5.Tesla,
                   L: 20.m,
                   totalTime: 100.h,
                   boreRadius: bR,
                   areaBore: π * bR^2)

proc initContext(cfg: Config): Context =
  # 1. integrate the solar flux
  ## NOTE: in practice this integration must not be done in this proc! Only perform once!
  let xs = linspace(0.0, 10.0, 1000)
  let fl = xs.mapIt(solarAxionFlux(it.keV))
  let integral = simpson(fl.mapIt(it.float), # convert units to float for compatibility
                         xs).cm⁻²•s⁻¹ # convert back to units (integrated out `keV⁻¹`!)
  # 2. construct background interpolator
  let bkg = newLinear1D(Energies.mapIt(it.float), Background)
  result = Context(cfg: cfg,
                   fluxIntegral: integral,
                   background: bkg)

proc solarAxionFlux(ω: keV): keV⁻¹•cm⁻²•s⁻¹ =
  ## Axion flux produced by the Primakoff effect in solar core
  ## in units of keV⁻¹•m⁻²•yr⁻¹
  ##
  ## This is computed *WITHOUT* `g_aγ` present, i.e. we set `g_aγ = 1 GeV⁻¹`, so that
  ## in `signal` we multiply only by `g_aγ⁴` (to simplify between `g_ae` and `g_aγ`).
  let flux = 2.0 * 1e18.keV⁻¹•m⁻²•yr⁻¹ * (1.GeV⁻¹ / 1e-12.GeV⁻¹)^2 * pow(ω / 1.keV, 2.450) * exp(-0.829 * ω / 1.keV)
  # convert flux to correct units
  result = flux.to(keV⁻¹•cm⁻²•s⁻¹)

func conversionProbability(ctx: Context): UnitLess =
  ## the conversion probability in the CAST magnet (depends on g_aγ)
  ## simplified vacuum conversion prob. for small masses
  ##
  ## This is computed *WITHOUT* `g_aγ` present, i.e. we set `g_aγ = 1 GeV⁻¹`, so that
  ## in `signal` we multiply only by `g_aγ⁴` (to simplify between `g_ae` and `g_aγ`).
  result = pow( (1.GeV⁻¹ * ctx.cfg.B.toNaturalUnit * ctx.cfg.L.toNaturalUnit / 2.0), 2.0 )

from numericalnim import simpson # simpson numerical integration routine
proc totalFlux(ctx: Context, g²: float): float =
  ## Flux integrated to total time, energy and area
  # 2. compute final flux by "integrating" out the time and area
  ## Multiply by `g²^2`. `g²` can be either `g_ae·g_aγ` or `g_aγ²`. We square that
  ## so that we multiply by the 'missing' terms in both P_aγ and the flux.
  result = ctx.fluxIntegral * totalTime * ctx.cfg.areaBore * ctx.conversionProbability() * g²^2

## NOTE: only important that signal and background have the same units!
proc signal(ctx: Context, E: keV, g²: float): keV⁻¹ =
  ## Returns the axion flux based on `g` and energy `E`
  result = solarAxionFlux(E) * totalTime.to(Second) * ctx.cfg.areaBore * ctx.conversionProbability() * g²^2

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
    "totalTime" : "Total solar tracking time"
  })
  let ctx = initContext(app)
  ctx.main # Only --help/--version/parse errors cause early exit
