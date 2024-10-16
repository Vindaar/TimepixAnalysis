import pkg / [unchained, seqmath, ggplotnim, xrayAttenuation]
import pkg / random / mersenne
import alea / [core, rng, gauss, poisson]

import numericalnim except linspace
import std / [sequtils, algorithm, math, terminal, stats, strformat, strutils, os, logging]

import ./basic_limit_plotting # contains the ggplotnim command for the toy limits histogram


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
  FluxKind = enum
    fkNone, fkAxionElectron, fkAxionPhoton
  ## NOTE: The default values here are not intended to be changed here. Rather you either
  ## give them as command line arguments, or you add fields in the
  ## `limit.config`
  ## file next to this file. Simply add a line like:
  ## `totalTime = 300.h`
  ## etc.
  Config = object # Stores experimental configuration
    # Experiment setup
    B*: Tesla                = 3.5.T                        ## Magnetic field of the magnet
    L*: Meter                = 20.m                         ## Length of the magnet
    totalTime*: Hour         = 100.h                        ## Total time in hours of solar tracking
    boreDiameter*: cm        = 70.cm                        ## Diameter of the bore (assumes a circular bore, A = πr²)
    chipArea*: mm²           = 5.mm * 5.mm                  ## Area in which all flux is assumed to be collected and in which
                                                            ## all candidates are detected & background is defined
    # Limit setup
    g2_max*: float           = 5e-21                        ## Maximum value for `g²` for the limit. Value should be large enough that
                                                            ## the likelihood function is zero there.
    nmc*: int                = 10000                        ## Number of toy candidate sets to sample for the expected limit
    epsilon*: float          = 1e-9                         ## Maximum value that the likelihood function is allowed to have at `g2_max`.
                                                            ## Ideally L(g2_max) = 0!
    steps*: int              = 1000                         ## Number of coupling constant steps to calculate for each limit.
    # Detector properties
    gas*: seq[string]        = @["Ar,0.977", "C4H10,0.023"] ## Gas mixture of the detector to use.
    pressure*: mbar          = 1050.mbar                    ## Pressure in the detector chamber.
    T*: Kelvin               = 293.15.K                     ## Temperature in the detector chamber
    chamberHeight*: cm       = 3.cm                         ## Height of the detector chamber
    # Window properties
    window*: string          = "Si3N4"                      ## Compound of the detector window, Silicon nitried 'Si3N4', Mylar: 'C10H8O4'",
    windowDensity*: g•cm⁻³   = 3.44.g•cm⁻³                  ## Density of the detector window material
    windowThickness*: μm     = 0.3.μm                       ## Thickness of the window (in μm!)
    windowThroughput*: float = 1.0                          ## Throughput of the window geometry for the flux, due to strongbacks etc.
                                                            ## NOTE: In reality this number may be much higher than a naive occlusion
                                                            ## calculation may suggest, due to the axion image not being occluded.
                                                            ## If unsure, leave it at 1.0.
    # Software parameters
    softwareEff*: float      = 0.9                          ## Software efficiency applied to signal collection.
    # Telescope
    effectiveArea*: string   = "../../resources/llnl_2016_dec_final_efficiency.csv" ## Path to a file containing telescope effective area
    tColE*: string           = "Energy [keV]"               ## Column name of the energy column in the CSV file
    tColEff*: string         = "Efficiency"                 ## Column name of the efficiency column in the CSV file
    tSep*: char              = ','                          ## Separator in the C/TSV file
    # Axions
    axionFlux*: string       = ""                           ## Path to a file containing axion flux (if none given, use analytical Primakoff)
    fluxKind*: FluxKind      = fkNone                       ## Type of flux for the given axion flux CSV ({fkAxionElectron, fkAxionPhoton})
    g2_ref*: float           = 0.0                          ## Reference value `g²` for which flux was computed, e.g. `g²_ae = 1e-26`, `g²_aγ = 1e-24`
    aColE*: string           = "Energy [keV]"               ## Column name of the energy column in the CSV file
    aColFlux*: string        = "Flux / keV⁻¹ m⁻² yr⁻¹"      ## Column name of the efficiency column in the CSV file
    aSep*: char              = ','                          ## Separator in the C/TSV file
    # Background
    backgroundFile*: string  = ""                           ## Path to a file containing background rate (if none given, use a 1e-5 range default)
    bColE*: string           = "Energy [keV]"               ## Column name of the energy column in the CSV file
    bColBkg*: string         = "Background"                 ## Column name of the efficiency column in the CSV file
    bSep*: char              = ','                          ## Separator in the C/TSV file
    bHeader*: string         = ""                           ## Optional indicator for the header in the CSV file, e.g. `"#"`.
    # Random number seed
    rngSeed*: int            = 299_792_458                  ## Random number generator seed
    # Output related parameters
    plotPath*: string        = "plots"                      ## Default location where plots are stored

  Context = ref object
    cfg: Config
    # Experiment
    areaBore: cm²                                          ## Full bore area, computed from `boreRadius`
    fluxIntegral: cm⁻²•s⁻¹ ## the base integral of the flux over the entire energy
    rnd: Random # Random sampler
    # Interpolators
    bkg: InterpolatorType[float] # background rate interpolator, keV⁻¹•cm⁻²•s⁻¹
    bMinE, bMaxE: keV # min and max energis
    meanBkgRate: keV⁻¹•cm⁻²•s⁻¹ # mean background rate over entire range. Used to know how many candidates to sample in total (Poisson)
    backgroundCdf: seq[float] # (E)CDF of the background rate. Used for candidate sampling
    energyForBCdf: seq[float] # energies corresponding to (E)CDF values
    telEff: InterpolatorType[float]     # Effective area of the telescope
    axionFlux: InterpolatorType[float]  # axion flux, can be nil if no CSV file given (then use Primakoff), keV⁻¹•cm⁻²•s⁻¹
    # Detector
    gm: GasMixture # Gas mixture of the detector, used to evaluate absorption efficiency
    window: Compound # Window material, used to compute transmission efficiency

  Candidate = object
    E: keV ## Energy of the candidate

  ## A helper object that stores statistical information about the toy limits
  LimitStatistics = object
    mean: float
    median: float ## The expected limit
    std: float ## The standard deviation
    p5, p16, p25, p75, p84, p95: float ## Different percentiles of the distribution

## Constants defining the default background rate
const
  Energies =   @[0.5,    1.5,    2.5,    3.5,    4.5,    5.5,     6.5,    7.5,  8.5,    9.5].mapIt(it.keV)
  Background = @[0.5e-5, 2.5e-5, 4.5e-5, 4.0e-5, 1.0e-5, 0.75e-5, 0.8e-5, 3e-5, 3.5e-5, 2.0e-5]
    .mapIt(it.keV⁻¹•cm⁻²•s⁻¹) # convert to a rate

# set up a logger
if not dirExists("logs"):
  createDir("logs")
var log = newFileLogger("logs" / "basic_limit_calc.log", fmtStr = "[$date - $time] - $levelname: ")
var L = newConsoleLogger()
addHandler(L)
addHandler(log)
proc info(logger: Logger, msgs: varargs[string, `$`]) =
  logger.log(lvlInfo, msgs)

proc `$`(obj: Config | LimitStatistics): string =
  result = $typeof(obj) & "(\n"
  for f, val in fieldPairs(obj):
    result.add "\t" & alignLeft(f, 20) & " = " & $val & "\n"
  result.add ")"

################################################################################
# (Limit) statistics
################################################################################

proc limitStatistics(limits: seq[float]): LimitStatistics =
  ## Compute differen statistics of the given limits
  result = LimitStatistics(
    mean: limits.mean,
    median: limits.median,
    std: limits.standardDeviation,
    p5: limits.percentile(5),
    p16: limits.percentile(16),
    p25: limits.percentile(25),
    p75: limits.percentile(75),
    p84: limits.percentile(84),
    p95: limits.percentile(95)
  )

template toCDF(data: seq[float], isCumSum = false): untyped =
  ## Computes the CDF of binned data
  var dataCdf = data
  if not isCumSum: seqmath.cumsum(dataCdf)
  let integral = dataCdf[^1]
  dataCdf.mapIt(it / integral)

################################################################################
# Procedures only dealing with input parsing
# (telescope effectiv area, axion flux, background rate)
################################################################################

template errorIfNotFound(c, df): untyped =
  if c notin df:
    raise newException(ValueError, "The required column: " & $c & " does not exist in the dataframe. " &
      "Existing columns are: " & $df.getKeys())

proc parseAxionFlux(cfg: Config): InterpolatorType[float] =
  ## Parses the axion flux, taking into account diving out the reference coupling
  if cfg.axionFlux.len > 0:
    if cfg.fluxKind == fkNone:
      raise newException(ValueError, "When giving a custom solar axion flux CSV file, please provide " &
        "a flux kind via `--fluxKind=fkAxionElectron|fkAxionPhoton`")
    if cfg.g2_ref == 0.0:
      raise newException(ValueError, "When giving a custom solar axion flux CSV file, please provide " &
        "a reference coupling used to compute these values via `--g2_ref=<value>`")
    var df = readCsv(cfg.axionFlux, sep = cfg.aSep)
    errorIfNotfound(cfg.aColE, df)
    errorIfNotfound(cfg.aColFlux, df)
    proc convert(x: float): float =
      result = x.keV⁻¹•m⁻²•yr⁻¹.to(keV⁻¹•cm⁻²•s⁻¹).float
    const fCol = "Flux [keV⁻¹•cm⁻²•s⁻¹]"
    if "Flux / keV⁻¹ m⁻² yr⁻¹" in df:
      df = df.filter(f{`type` == "Total flux"}) # filter to only total flux
        .mutate(f{fCol ~ convert(idx("Flux / keV⁻¹ m⁻² yr⁻¹"))})
    else:
      stdout.styledWrite(fgYellow, "[WARNING]: The given column name is not the default name. We assume " &
        "the flux is given in `keV⁻¹•cm⁻²•s⁻¹`.\n")
    # divide out the reference coupling
    df = df.mutate(f{float: fCol ~ idx(fCol) / cfg.g2_ref})
    result = newLinear1D(df[cfg.aColE, float].toSeq1D,
                         df[fCol, float].toSeq1D)

proc parseFile(f, x, y: string, sep: char, header: string = ""): InterpolatorType[float] =
  var df = readCsv(f, sep = sep, header = header)
  errorIfNotfound(x, df)
  errorIfNotfound(y, df)
  result = newLinear1D(df[x, float].toSeq1D, df[y, float].toSeq1D)

proc parseTelescopeEff(cfg: Config): InterpolatorType[float] =
  ## Parses the efficiency associated with the effective area of the telescope
  if cfg.effectiveArea.len == 0:
    raise newException(ValueError, "Please provide a file for the telescope efficiency via " &
      "`--effectiveArea=<file>`")
  result = parseFile(cfg.effectiveArea, cfg.tColE, cfg.tColEff, cfg.tSep)

proc parseBackgroundRate(cfg: Config): InterpolatorType[float] =
  ## Parses the background rate
  if cfg.backgroundFile.len > 0:
    result = parseFile(cfg.backgroundFile, cfg.bColE, cfg.bColBkg, cfg.bSep, cfg.bHeader)
  else:
    result = newLinear1D(Energies.mapIt(it.float), Background.mapIt(it.float)) # strip type info :(


################################################################################
# Configuration & runtime value logic
################################################################################

proc signalRate(ctx: Context, E: keV): keV⁻¹•cm⁻²•s⁻¹ # forward declare

proc initConfig(): Config =
  ## This stores all 'experimental configuration' parameter. These are all constant and
  ## don't change for one run of the program.
  ##
  ## The majority of fields are already initialized with default values.
  result = Config()

proc initContext(cfg: Config): Context =
  # 1. construct background interpolator
  let bkg = parseBackgroundRate(cfg)
  # 2. intiialize the gas mixture & window
  let gm = parseGasMixture(cfg.gas, cfg.pressure, cfg.T)
  let window = initCompound(cfg.windowDensity, parseCompound(cfg.window))
  # 3. integrate the solar flux, combined with efficiencies
  let Es = linspace(0.0, 10.0, 1000) # keV, energies for flux integration & background interp
  let bkgVals = Es.mapIt(bkg.eval(it))
  result = Context(
    cfg: cfg,
    areaBore: π * (cfg.boreDiameter / 2.0)^2,
    rnd: wrap(initMersenneTwister(cfg.rngSeed.uint32)),
    bkg: bkg,
    bMinE: Es.min.keV, bMaxE: Es.max.keV,
    meanBkgRate: mean(bkgVals).keV⁻¹•cm⁻²•s⁻¹,
    backgroundCdf: toCdf( bkgVals ),
    energyForBCdf: Es,
    gm: gm,
    window: window,
    telEff: parseTelescopeEff(cfg),
    axionFlux: parseAxionFlux(cfg)
  )
  if cfg.axionFlux.len == 0:
    result.cfg.fluxKind = fkAxionPhoton
  # calc flux & convert back to float for compatibility with `simpson`
  let fl = Es.mapIt(signalRate(result, it.keV).float)
  result.fluxIntegral = simpson(fl, Es).cm⁻²•s⁻¹ # convert back to units (integrated out `keV⁻¹`!)

################################################################################
# Random sampling logic
################################################################################

proc sampleCandidates(ctx: Context): seq[Candidate] =
  ## This samples candidates from the background distribution, by using
  ## inverse transform sampling:
  ## https://en.wikipedia.org/wiki/Inverse_transform_sampling
  # 1. compute number of expected candidates
  #    `dN/(dE dt dA) · ΔE · t · A`
  let expNum = ctx.meanBkgRate * (ctx.bMaxE - ctx.bMinE) * ctx.cfg.totalTime * ctx.cfg.chipArea
  # 2. construct Poisson sampler with `λ = expNum`
  let pois = poisson(expNum)
  # 3. draw a sample of candidates, and sample their energies according
  #    to the background CDF
  let uni = uniform(0.0, 1.0) # uniform sampler
  let num = ctx.rnd.sample(pois).int
  result = newSeq[Candidate](num)
  for i in 0 ..< num:
    let u = ctx.rnd.sample(uni) # uniform random number
    let eIdx = ctx.backgroundCDF.lowerBound(u) # find closest value in CDF
    let energy = ctx.energyForBCDF[eIdx].keV # lookup energy, i.e. "invert axis" / inverse transform
    result[i] = Candidate(E: energy)

################################################################################
# Limit ingredients
################################################################################

proc primakoffFlux(E: keV): keV⁻¹•cm⁻²•s⁻¹ =
  ## Axion flux produced by the Primakoff effect in solar core
  ## in units of keV⁻¹•m⁻²•yr⁻¹
  ##
  ## This is computed *WITHOUT* `g_aγ` present, i.e. we set `g_aγ = 1 GeV⁻¹`, so that
  ## in `signal` we multiply only by `g_aγ⁴` (to simplify between `g_ae` and `g_aγ`).
  ##
  ## Analytical Primakoff flux expression from: `doi:10.1088/1475-7516/2013/05/010`
  ## (CAST 2013 axion electron paper)
  let flux = 2.0 * 1e18.keV⁻¹•m⁻²•yr⁻¹ * (1.GeV⁻¹ / 1e-12.GeV⁻¹)^2 * pow(E / 1.keV, 2.450) * exp(-0.829 * E / 1.keV)
  # convert flux to correct units
  result = flux.to(keV⁻¹•cm⁻²•s⁻¹)

proc solarAxionFlux(ctx: Context, E: keV): keV⁻¹•cm⁻²•s⁻¹ =
  ## Solar axion flux. Either analytical if no axion flux input given or interpolated.
  if ctx.axionFlux != nil: result = ctx.axionFlux.eval(E.float).keV⁻¹•cm⁻²•s⁻¹
  else: result = primakoffFlux(E)

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
  ## XXX: Make `solarAxionFlux` depending on config
  result = ctx.solarAxionFlux(E) * ctx.conversionProbability() *
           ctx.telEff.eval(E.float) *
           ctx.window.transmission(ctx.cfg.windowThickness, E) *
           ctx.gm.absorption(ctx.cfg.chamberHeight, E) *
           ctx.cfg.softwareEff *
           ctx.cfg.windowThroughput

proc totalFlux(ctx: Context, g²: float): float =
  ## Flux integrated to total time, energy and area
  ##
  ## Multiply by `g²^2`. `g²` can be either `g_ae·g_aγ` or `g_aγ²`. We square that
  ## so that we multiply by the 'missing' terms in both P_aγ and the flux.
  result = ctx.fluxIntegral * ctx.cfg.totalTime * ctx.areaBore * g²^2

## NOTE: only important that signal and background have the same units!
proc signal(ctx: Context, E: keV, g²: float): keV⁻¹ =
  ## Returns the axion flux based on `g` and energy `E`
  result = ctx.signalRate(E) * ctx.cfg.totalTime.to(Second) * ctx.areaBore * g²^2

proc background(ctx: Context, E: keV): keV⁻¹ =
  ## Compute an interpolation the background at the energy `E`.
  ##
  ## Note: area of interest is the region on the chip, in which the signal is focused!
  ## This also allows us to see that the "closer" we cut to the expected axion signal on the
  ## detector, the less background we have compared to the *fixed* signal flux!
  result = (ctx.bkg.eval(E.float).keV⁻¹•cm⁻²•s⁻¹ * ctx.cfg.totalTime * ctx.cfg.chipArea).to(keV⁻¹)

proc likelihood(ctx: Context, g²: float, cs: seq[Candidate]): float =
  ## `cs` = each candidate with an energy `E` inside the `chipArea`.
  ## Note: In this form of the likelihood function, with s/b
  ## we can get away with using non logs for all interesting
  ## parameter ranges. And we need the actual likelihood (not the χ²)
  ## to compute the limit.
  result = exp(-ctx.totalFlux(g²)) # `e^{-s_tot}`
  for i in 0 ..< cs.len:
    let E = cs[i].E # energy of this candidate
    let s = ctx.signal(E, g²)
    let b = ctx.background(E)
    result *= (1 + s / b)

################################################################################
# Limit calculation
################################################################################

proc computeLimit(ctx: Context, cs: seq[Candidate], toPlot = false): float =
  ## Computes the limit for the given candidates, based on the
  ## parameter range `[0, g²_max]` (defined from configuration).
  let g² = linspace(0.0, ctx.cfg.g2_max, ctx.cfg.steps) # g²
  let Ls = g².mapIt(ctx.likelihood(it, cs))
  if Ls[^1] > ctx.cfg.epsilon:
    raise newException(ValueError, "The coupling constant range for the given candidates / experimental " &
      "parameters is not sufficient. The likelihood function: \n" &
      &"L({ctx.cfg.g2_max}) = {Ls[^1]},\n" &
      "which is larger than the target\n" &
      &"ε = {ctx.cfg.epsilon}.")
  let lCdf = toCdf Ls
  let limitIdx = lCdf.lowerBound(0.95) # limit at 95% of the CDF
  result = g²[limitIdx]

  if toPlot:
    ggplot(toDf(g², Ls), aes("g²", "Ls")) +
      geom_line() + ggsave(ctx.cfg.plotPath / "limit_Ls.pdf")

const Parallel {.booldefine.} = false
when Parallel:
  import pkg / forked
proc monteCarloLimits(ctx: Context): seq[float] =
  ## Computes the expected limit for `nmc` toy candidates and returns all limits.
  let nmc = ctx.cfg.nmc
  result = newSeq[float](nmc)

  template loopBody(): untyped {.dirty.} = # actual loop body, to remove redundancy due to `Parallel`
    const InfoNum = 500
    if i mod InfoNum == 0:
      echo &"Iteration: {i} of {nmc} ({(i.float / nmc.float) * 100.0:.2f}%)"
    let cs = ctx.sampleCandidates()
    let limit = ctx.computeLimit(cs, i mod InfoNum == 0)

  when Parallel:
    for (i, x) in forked(0 ..< nmc):
      loopBody()
      send limit
      join: result[i] = limit
  else:
    for i in 0 ..< nmc:
      loopBody()
      result[i] = limit

proc computeExpectedLimit(ctx: Context): LimitStatistics =
  ## Generates toy candidate sets, computes their limits, calculate the
  ## expected limit and produces a histogram of the limit distribution.
  let limits = ctx.monteCarloLimits()
  result = limitStatistics(limits)
  plotToyLimits(limits, ctx.cfg.plotPath / "toy_limits_histo.pdf")

proc main(ctx: Context) =
  ## A simple limit calculation tool that makes it easy to replace experiment and detector parameters.
  ##
  ## It computes the expected limit based on a given number of toy candidates sampled from the
  ## background distribution. Each limit is computed from a linear scan in the coupling
  ## constant parameter, from 0 to `g2_max` (configuraton / CL argument). As such it is
  ## vital that the parameter range is large enough to show the entire range of the likelihood
  ## function.
  log.info "Computing expected limit for: "
  log.info $ctx.cfg
  let limitStats = ctx.computeExpectedLimit()

  case ctx.cfg.fluxKind
  of fkAxionElectron:
    info "Expected limit: g_aγ·g_ae = ", limitStats.median
  of fkAxionPhoton:
    info "Expected limit: g_aγ = ", sqrt(limitStats.median)
  else: discard

  info "Full statistics (as g²):"
  info $limitStats

when isMainModule:
  import cligen
  const ConfigPath = "limit.config"
  include mergeCfgEnvLocal
  include unchained/cligenParseUnits

  const dfl* = initConfig() # set defaults!=default for type
  var app = initFromCL(dfl, cmdName = "foo", help = {
    "B"               : "Magnetic field of the magnet",
    "L"               : "Length of the magnet",
    "totalTime"       : "Total time in hours of solar tracking",
    "boreDiameter"    : "Diameter of the bore (assumes a circular bore, A = πr²)",
    "chipArea"        : """Area in which all flux is assumed to be collected and in which all candidates
are detected & background is defined.
*NOTE*: This is a *vital* parameter, because it directly affects the signal / background ratio. Our assumption
is always that all flux is contained in this region and the background constant inside it.""",
    "g2_max"          : "Maximum value for `g²` for the limit. Value should be large enough that the likelihood function is zero there.",
    "nmc"             : "Number of toy candidate sets to sample for the expected limit",
    "epsilon"         : "Maximum value that the likelihood function is allowed to have at `g2_max`. Ideally L(g2_max) = 0!",
    "steps"           : "Number of coupling constant steps to calculate for each limit.",
    # Detector properties
    "gas"             : "Gas mixture of the detector to use.",
    "pressure"        : "Pressure in the detector chamber.",
    "T"               : "Temperature in the detector chamber",
    "chamberHeight"   : "Height of the detector chamber",
    # Window properties
    "window"          : "Compound of the detector window, Silicon nitried 'Si3N4', Mylar: 'C10H8O4'",
    "windowDensity"   : "Density of the detector window material",
    "windowThickness" : "Thickness of the window (in μm!)",
    "windowThroughput": """Throughput of the window geometry for the flux, due to strongbacks etc.
NOTE: In reality this number may be much higher than a naive occlusion
calculation may suggest, due to the axion image not being occluded.
If unsure, leave it at 1.0.""",
    # Software parameters
    "softwareEff"     : "Software efficiency applied to signal collection",
    # Telescope
    "effectiveArea"   : "Path to a file containing telescope effective area",
    "tColE"           : "Column name of the energy column in the CSV file",
    "tColEff"         : "Column name of the efficiency column in the CSV file",
    "tSep"            : "Separator in the C/TSV file",
    # Axions
    "axionFlux"       : "Path to a file containing axion flux (if none given, use analytical Primakoff)",
    "fluxKind"        : "Type of flux for the given axion flux CSV ({fkAxionElectron, fkAxionPhoton})",
    "g2_ref"          : "Reference value `g²` for which flux was computed, e.g. `g²_ae = 1e-26`, `g²_aγ = 1e-24`",
    "aColE"           : "Column name of the energy column in the CSV file",
    "aColFlux"        : "Column name of the efficiency column in the CSV file",
    "aSep"            : "Separator in the C/TSV file",
    # Background
    "backgroundFile"  : "Path to a file containing background rate (if none given, use a 1e-5 range default)",
    "bColE"           : "Column name of the energy column in the CSV file",
    "bColBkg"         : "Column name of the efficiency column in the CSV file",
    "bSep"            : "Separator in the C/TSV file",
    "bHeader"         : "Optional indicator for the header in the CSV file, e.g. `\"#\"`.",
    # Random number seed
    "rngSeed"         : "Random number generator seed",
    # Output related parameters
    "plotPath"        : "Default location where plots are stored"
  })
  let ctx = initContext(app)
  ctx.main # Only --help/--version/parse errors cause early exit
