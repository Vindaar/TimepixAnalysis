import std / [math, strformat, sequtils, os]
import ggplotnim, unchained
import arraymancer except linspace
import numericalnim except linspace
from seqmath import gauss, linspace
#from ingrid / tos_helpers import geometry

const TrackingBackgroundRatio* = 19.605 ## TODO: replace by actual time!!
const EnergyCutoff* = 12.0

proc toKDE*(df: DataFrame, toPlot = false,
            outname = ""): DataFrame =
  echo "[KDE] Number of clusters in DF: ", df
  echo "[KDE] Number of clusters in DF < ", EnergyCutoff, " keV: ", df.filter(f{`Energy` <= EnergyCutoff})
  let dfFiltered = df.filter(f{`Energy` <= 12.0},
                         f{float -> bool: `centerX` in 4.5 .. 9.5 and `centerY` in 4.5 .. 9.5}
  )
  let energy = dfFiltered["Energy", float]
  #dfFiltered.showBrowser()
  #if true: quit()

  echo "[KDE] Number of clusters ", energy.size, " normalized to tracking ", energy.size.float / TrackingBackgroundRatio

  let xs = linspace(energy.min, energy.max, 1000)
  echo "MAX ENERGY ", energy.max
  var kde = kde(energy, bw = 0.3, normalize = true)
  defUnit(cm²)
  let scale = energy.size.float /
    ( TrackingBackgroundRatio ) / #3318.h.to(Second) / 190.0.h.to(Second) ) /
    ( (5.0.mm * 5.0.mm).to(cm²) ) #.to(cm²) # normalize to cm⁻²
    #/ # * (175.0 * 3600.0) / #/ # the total background time (the detector was live in)
    #( (5.0.mm * 5.0.mm) / (14.0.mm * 14.0.mm) ) # scale counts up from gold region to equivalent in full chip
    #1.0 / (8359.18367347)  # ratio of pixels in gold region
    #pow(0.95 - 0.45, 2.0) / # area of gold region!
    #pow(1.4, 2.0) /
    #12.0 # * # to get to / keV
    #(pow(1.4, 2.0) / 65535.0)
  echo "[KDE] SCALE: ", scale
  kde = kde.map_inline:
    x * scale.float
  #echo kde[0 .. 500]
  result = toDf({"Energy" : xs, "KDE" : kde})

  let integral = simpson(result["KDE", float].toSeq1D,
                         result["Energy", float].toSeq1D)
  echo "XXX: think about this INTEGRAL OF KDE!!!"
  echo "[KDE] ", integral, " and ratio to naive ", (integral / (energy.size.float / TrackingBackgroundRatio))
  if toPlot:
    let outname = if outname.len == 0: "/home/basti/org/Figs/statusAndProgress/backgroundRates/background_2017_2018_kde_rate.pdf"
                  else: outname
    ggplot(result, aes("Energy", "KDE")) +
      geom_line() +
      ggtitle("KDE of the background clusters, normalized to keV⁻¹•cm⁻² in 190 h of tracking") +
      ggsave(outname)
  #if true: quit()

## XXX: make them compile time variables and only allow modification from static vars?
var Radius* = 40.0 # 100.0 #33.0
var Sigma* = 40.0 / 3.0 #33.3 #11.111
var EnergyRange* = 0.6.keV # Inf.keV

#type
#  MyMetric* = object
#    radius*: float
#    sigma*: float
#    energyRange*: keV
#
### A dummy proc that identifies our `MyMetric` as a custom metric
#proc toCustomMetric(m: MyMetric): CustomMetric = CustomMetric()
#let mym = MyMetric()

proc distance*(metric: typedesc[CustomMetric], v, w: Tensor[float]): float =
  #echo "Metric ", metric
  #doAssert v.squeeze.rank == 1
  #doAssert w.squeeze.rank == 1
  doAssert EnergyRange.float != Inf, "WARNING: You did not set the `EnergyRange` (and likely `Sigma` and `Radius`) "&
    "variables! It is required to set them, as we cannot pass these variables to the `distance` procedure. They " &
    "are defined as globals above the procedure body!"
  #result = Euclidean.distance(v, w)
  #let diff = abs(v -. w)
  #let arg1 = diff[0, 0] #abs(v[0] - w[0])
  #let arg2 = diff[0, 1] #abs(v[1] - w[1])
  #let arg3 = diff[0, 2] #abs(v[2] - w[2])
  # NOTE: this is the fastest way to compute the distance
  # - no squeeze
  # - no temp tensor allocation
  # Argument is non squeezed, so we need to access 2nd axis 0 manually
  let arg1 = abs(v[0, 0] - w[0, 0])
  let arg2 = abs(v[0, 1] - w[0, 1])
  let arg3 = abs(v[0, 2] - w[0, 2])
  let xyDist = arg1*arg1 + arg2*arg2
  ##echo "xy dist ", xyDist, " vs ", Euclidean.distance(v, w)
  let zDist = arg3*arg3
  if zDist <= (EnergyRange * EnergyRange).float:
    result = xyDist
  else:
    result = (2 * Radius * Radius).float # just some value larger than Radius² #pow(sqrt(zDist) * 6.0, 2.0)
  #if xyDist > zDist:
  #  result = xyDist
  #elif xyDist < zDist and zDist <= Radius * Radius:
  #  result = xyDist
  #else:
  #  result = zDist

  ## XXX: In order to treat the energy as pure distance without gaussian behavior, we can do:
  ## - compute distance in both xy and z as currently
  ## - use `zDist` *only* (!!!) as an early "return" so to say. I.e. if larger than a cutoff we
  ##   define, return it. Needs to be a global, as we don't know that cutoff from `v`, `w`, or rather
  ##   hardcode into distance proc
  ## - else *always* return `xyDist`. This guarantees to give us the distance information of the points
  ##   *always* along the xy, which is important for the weighing, but *not* along energy


#proc distance(metric: typedesc[CustomMetric], v, w: Tensor[float]): float =
#  MyMetric.distance(v, w)

import helpers/circle_segments
proc correctEdgeCutoff*(val, radius: float, x, y: int): float {.inline.} =
  ## Corrects the effects of area being cut off for the given `val` if it is
  ## positioned at `(x, y)` and the considered radius is `radius`.
  ##
  ## TODO: for our normal weighted values, this edge cutoff is not correct. We need to
  ## renormalize by the *weighted* area and not the unweighted one...
  let refArea = PI * radius * radius
  let areaLeft = areaCircleTwoLinesCut(radius, min(x, 256 - x).float, min(y, 256 - y).float)
  result = val * refArea / areaLeft

proc correctEdgeCutoff(t: var Tensor[float], radius: float) =
  ## Applies the edge correction for every point in the given tensor
  for y in 0 ..< 256:
    for x in 0 ..< 256:
      t[y, x] = correctEdgeCutoff(t[y, x], radius, x, y)

proc correctEdgeCutoff3D(t: var Tensor[float], radius: float) =
  ## Applies the edge correction for every point in the given tensor
  for y in 0 ..< t.shape[0]:
    for x in 0 ..< t.shape[1]:
      for E in 0 ..< t.shape[2]:
        t[y, x, E] = correctEdgeCutoff(t[y, x, E], radius, x, y)

proc plot2d[T](bl: T) =
  let pix = 256
  var xs = newSeq[int](pix * pix)
  var ys = newSeq[int](pix * pix)
  var cs = newSeq[float](pix * pix)
  var idx = 0
  for y in 0 ..< pix:
    for x in 0 ..< pix:
      xs[idx] = x
      ys[idx] = y
      cs[idx] = bl.eval(y.float, x.float)#t[y, x]
      inc idx
  ggplot(toDf(xs, ys, cs), aes("xs", "ys", fill = "cs")) +
    geom_raster() +
    #scale_fill_continuous(scale = (low: 0.0, high: 10.0)) +
    ggsave("/tmp/test.pdf")

proc plot2dTensor*(t: Tensor[float], outname = "/tmp/test_tensor.pdf",
                   title = "",
                   yMax = 0.0) =
  var xs = newSeq[int](t.size)
  var ys = newSeq[int](t.size)
  var cs = newSeq[float](t.size)
  var idx = 0
  for y in 0 ..< t.shape[0]:
    for x in 0 ..< t.shape[1]:
      xs[idx] = x
      ys[idx] = y
      #if t[y, x] > 5.0:
      #  echo "Noisy pixel: ", x, " and ", y, " have count ", t[y, x]
      #  inc sumNoise, t[y, x].int
      cs[idx] = t[y, x]
      inc idx
  #echo "Total noisy things: ", sumNoise
  template low: untyped = 4.5 / 14.0 * 256.0
  template hih: untyped = 9.5 / 14.0 * 256.0

  let df = toDf(xs, ys, cs)
  ggplot(df, aes("xs", "ys", fill = "cs")) +
    geom_raster() +
    geom_linerange(aes = aes(x = low(), yMin = low(), yMax = hih()), color = some(parseHex("FF0000"))) +
    geom_linerange(aes = aes(x = hih(), yMin = low(), yMax = hih()), color = some(parseHex("FF0000"))) +
    geom_linerange(aes = aes(y = low(), xMin = low(), xMax = hih()), color = some(parseHex("FF0000"))) +
    geom_linerange(aes = aes(y = hih(), xMin = low(), xMax = hih()), color = some(parseHex("FF0000"))) +
    scale_fill_continuous(scale = (low: 0.0, high: yMax)) +
    xlim(0, 256) + ylim(0, 256) +
    margin(top = 1.5) +
    ggtitle(title) +
    ggsave(outname)

proc plot3DTensor(t: Tensor[float], outname = "/tmp/test_tensor_3d.pdf",
                  title = "") =
  var xs = newSeq[int](t.size)
  var ys = newSeq[int](t.size)
  var Es = newSeq[int](t.size)
  var cs = newSeq[float](t.size)
  var idx = 0
  var sumNoise = 0
  for y in 0 ..< t.shape[0]:
    for x in 0 ..< t.shape[1]:
      for E in 0 ..< t.shape[2]:
        xs[idx] = x
        ys[idx] = y
        Es[idx] = E
        #if t[y, x] > 5.0:
        #  echo "Noisy pixel: ", x, " and ", y, " have count ", t[y, x]
        #  inc sumNoise, t[y, x].int
        cs[idx] = t[y, x, E]
        inc idx
  echo "Total noisy things: ", sumNoise
  when false:
    ggplot(toDf(xs, ys, Es, cs), aes("xs", "ys", fill = "cs")) +
      facet_wrap("Es", scales = "free") +
      geom_raster() +
      #scale_fill_continuous(scale = (low: 0.0, high: 10.0)) +
      ggtitle(title) +
      ggsave(outname, width = 1900, height = 1500)
  else:
    for tup, subDf in groups(toDf(xs, ys, Es, cs).group_by("Es")):
      ggplot(subDf, aes("xs", "ys", fill = "cs")) +
        geom_raster() +
        #scale_fill_continuous(scale = (low: 0.0, high: 10.0)) +
        ggtitle(title & " Energy: " & $tup[0][1].toFloat) +
        ggsave(&"/tmp/back_plot_energy_{tup[0][1].toFloat}.pdf")

proc plotDf(df: DataFrame, title, outname: string) =
  ggplot(df, aes("centerX", "centerY")) +
    geom_point() +
    ggtitle(title) +
    ggsave(outname)

template compValue*(tup: untyped, byCount = false, energyConst = false): untyped =
  ## Computes the weighted (`byCount`) / unweighted (`not byCount`) value associated
  ## with a position from the given neighbors (`tup` is a return of `query_ball_point`
  ## on a k-d tree)
  if byCount:
    tup.idx.size.float
  else:
    # weigh by distance using gaussian of radius being 3 sigma
    let dists = tup[0]
    var val = 0.0
    for d in items(dists):
      val += seqmath.gauss(d, mean = 0.0, sigma = Sigma)
    val

proc compDistance(t: var Tensor[float], kd: KDTree[float], radius: float,
                  byCount = false) =
  for y in 0 ..< 256:
    for x in 0 ..< 256:
      let tup = kd.query_ball_point([x.float, y.float].toTensor, radius)
      let val = compValue(tup)
      t[y, x] = val

proc compValueTree(kd: KDTree[float], x, y, E: float,
                   radius: float, metric: typedesc[AnyMetric],
                   byCount = false): float {.inline.} =
  ## Queries the tree at the given coordinate and energy and returns the correctly
  ## weighted value at the point.
  let tup = kd.query_ball_point([x, y, E].toTensor, radius, metric = metric)
  if x == 127 and y == 127:
    toDf({"dists" : tup[0]}).writeCsv("/tmp/distances_127_127.csv")
    #let df = seqsDoDf(dists)
  result = compValue(
    tup,
    byCount = byCount
  )

proc compDistance3D(t: var Tensor[float], Es: seq[float], kd: KDTree[float], radius: float,
                    byCount = false,
                    metric = Euclidean) =
  for y in 0 ..< 256:
    echo "Starting y ", y
    for x in 0 ..< 256:
      for E in 0 ..< Es.len:
        t[y, x, E] = kd.compValueTree(x.float, y.float, Es[E], radius, metric, byCount)

defUnit(keV⁻¹•cm⁻²•s⁻¹, toExport = true)
proc normalizeValue*(x, radius: float, energyRange: keV, backgroundTime: Hour): keV⁻¹•cm⁻²•s⁻¹ =
  let pixelSizeRatio = 65536 / (1.4 * 1.4).cm²
  when false:
    # case for regular circle with weights 1
    let area = π * radius * radius # area in pixel
  else:
    let σ = Sigma
    ## This comes for integration with `sagemath` over the gaussian weighting. See the notes.
    let area = -2*π*(σ*σ * exp(-1/2 * radius*radius / (σ*σ)) - (σ*σ))
  let energyRange = energyRange * 2.0 # we look at (factor 2 for radius)
  ## NOTE: for an *expected limit* this time must be the full background time, as it
  ## is the time that describes the number of clusters we have in the input! Thus,
  ## if we change it to `t_back - t_track`, we artificially increase our background!
  #let backgroundTime = 3318.h.to(Second) #(3318.h - 169.h).to(Second)
  let factor = area / pixelSizeRatio * # area in cm²
    energyRange *
    backgroundTime.to(Second)
  result = x / factor

proc normalizeTensor(t: var Tensor[float], energies: int,
                     radius: float) =
  ## Normalizes the tensor to units of /keV /cm^2 /s
  echo "Normalizing tensor by time: \n\n\n"
  for y in 0 ..< 256:
    for x in 0 ..< 256:
      for E in 0 ..< energies:
        t[y, x, E] = normalizeValue(t[y, x, E], radius, EnergyRange, 3318.Hour).float

proc compNormalized(kd: KDTree[float], x, y: int, E: keV,
                    radius: float,
                    energyRange: keV,
                    backgroundTime: Hour,
                    metric: typedesc[AnyMetric]
                   ): float =
  ## Computes a correctly normalized value for the given position and energy,
  ## using the `radius` from the given tree `kd`.
  result = compValueTree(kd, x.float, y.float, E.float, radius, metric)
    .correctEdgeCutoff(radius, x, y)
    .normalizeValue(radius, energyRange, backgroundTime).float

template fillChip(body: untyped): untyped =
  var t {.inject.} = zeros[float]([256, 256])
  for y {.inject.} in 0 ..< 256:
    for x {.inject.} in 0 ..< 256:
      body
  t

proc compInterEnergy(t: var Tensor[float], kd: KDTree[float],
                     energy: keV,
                     radius: float,
                     energyRange: keV,
                     backgroundTime: Hour,
                     metric: typedesc[AnyMetric],
                     byCount = false) =
  t = fillChip:
    t[y, x] = kd.compNormalized(x, y, energy, radius, energyRange, backgroundTime, metric)
    if x == 128 and y == 128:
      echo "Rate at center: ", t[y, x]

func toIdx*(arg: float): int = (arg / 14.0 * 256.0).round.int.clamp(0, 255)
func toInch*(arg: float|int): float = (arg.float / 256.0 * 14.0).clamp(0.0, 14.0)
proc plotGoldRegionBackgroundRate(kd: KDTree[float], outfile: string,
                                  title: string,
                                  backgroundTime = 3318.Hour) =
  var num = 25
  let coords = linspace(4.5, 9.5, num) # the gold region
  var energies = linspace(0.0, 12.0, 75)
  var rates = newSeq[float](energies.len)
  for i, E in energies:
    var val = 0.0
    for y in coords:
      for x in coords:
        val += compNormalized(kd, x.toIdx, y.toIdx, E.keV, Radius.float, EnergyRange,
                              backgroundTime,
                              metric = CustomMetric)
    rates[i] = val / (num * num).float
    echo "At energy ", E, " of index ", i, " rate: ", rates[i]

  let dfL = toDf(energies, rates)
  ggplot(dfL, aes("energies", "rates")) +
    geom_point() +
    ggtitle(title) +
    ggsave(outfile) #"/tmp/background_gold_region_from_interp.pdf")

template plotEnergySlice*(outfile, title: string, yMax: float, body: untyped): untyped =
  let tr = fillChip:
    body
  tr.plot2dTensor(outfile, title, yMax)

proc plotSingleEnergySlice*(kd: KDTree[float], energy: keV,
                            backgroundTime = 3318.Hour,
                            outfile = "", title = "") =
  let title = if title.len > 0: title else: &"Background interpolation at {energy} keV"
  let outfile = if outfile.len > 0: outfile else: &"/tmp/back_interp_energy_{energy}.pdf"
  var tr = zeros[float]([256, 256])
  tr.compInterEnergy(kd, energy, Radius.float, EnergyRange, backgroundTime, byCount = false, metric = CustomMetric)
  tr.plot2dTensor(outfile, title)

proc toNearestNeighborTree*(df: DataFrame): KDTree[float] =
  ## calls the correct interpolation function and returns the interpolated data
  echo "[INFO]: Building tree based on ", df.len, " background clusters in input"
  let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
                     df["centerY", float].map_inline(toIdx(x).float),
                     df["Energy", float].map_inline(x)], axis = 1)
                     #df["Energy", float].map_inline(x * 25.6)], axis = 1)
  result = kdTree(tTree, leafSize = 16, balancedTree = true)

proc studyBackgroundInterpolation*(df: DataFrame, toPlot = false): DataFrame =
  ## generates a kd tree based on the data and generates multiple plots
  ## we use to study the interpolation and determine good parameters
  var t = zeros[float]([256, 256])
  for idx in 0 ..< df.len:
    let x = toIdx df["centerX", float][idx]
    let y = toIdx df["centerY", float][idx]
    t[y, x] += 1
  t.plot2dTensor()
  #if true: quit()
  block Bilinear:
    var bl = newBilinearSpline(t, (0.0, 255.0), (0.0, 255.0)) # bicubic produces negative values!
    bl.plot2d()

  #block kdTree:
  #  let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
  #                     df["centerY", float].map_inline(toIdx(x).float)],
  #                     axis = 1)
  #  let kd = kdTree(tTree, leafSize = 16, balancedTree = true)
  #  var treeDist = zeros[float]([256, 256])
  #
  #  for radius in [30]: #arange(10, 100, 10):
  #    treeDist.compDistance(kd, radius.float, byCount = true)
  #    treeDist.plot2dTensor("/tmp/background_radius_byenergy_" & $radius & "_bycount.pdf",
  #                          "k-d tree interpolation with radius: " & $radius & " pixels")
  #    treeDist.correctEdgeCutoff(radius.float)
  #    treeDist.plot2dTensor("/tmp/background_radius_byenergy_" & $radius & "_bycount_corrected.pdf",
  #                          "k-d tree interpolation with radius: " & $radius & " pixels")
  #    treeDist.compDistance(kd, radius.float)
  #    treeDist.plot2dTensor("/tmp/background_radius_byenergy_" & $radius & ".pdf",
  #                          "k-d tree interpolation with radius: " & $radius & " pixels")
  #    treeDist.correctEdgeCutoff(radius.float)
  #    treeDist.plot2dTensor("/tmp/background_radius_byenergy_" & $radius & "_corrected.pdf",
  #                          "k-d tree interpolation with radius: " & $radius & " pixels")

    # now plot interpolation based on energy
  echo "3d???\n\n"
  block kdTree3D:
    let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
                       df["centerY", float].map_inline(toIdx(x).float),
                       df["Energy", float].map_inline(x)], axis = 1)
                       #df["Energy", float].map_inline(x * 25.6)], axis = 1)
    let kd = kdTree(tTree, leafSize = 16, balancedTree = true)

    when false:
      kd.plotSingleEnergySlice(1.0.keV)
    when false:
      let Es = @[1.0, 2.0, 4.0, 5.0] #linspace(0.0, 12.0, 10)
      #for (radius, sigma, eSigma) in [(100.0, 33.3333, 0.3),
      #                                (100.0, 15.0, 0.3),
      #                                (75.0, 75.0 / 3.0, 0.3),
      #                                (50.0, 50.0 / 3.0, 0.3),
      #                                (33.333, 11.1111, 0.3),
      #                                (25.0, 25.0 / 3.0, 0.3),
      #                                (100.0, 33.3, 0.5),
      #                                (50.0, 50.0 / 3.0, 0.5)]:
      for (radius, sigma, eSigma) in [(33.0, 11.111, 0.3),
                                      (33.0, 11.111, 0.5),
                                      (25.0, 25.0/3.0, 0.3),
                                      (25.0, 25.0/3.0, 0.5),
                                      (20.0, 20.0/3.0, 0.3),
                                      (20.0, 20.0/3.0, 0.5),
                                      (15.0, 15.0/3.0, 0.3),
                                      (15.0, 15.0/3.0, 0.5)]:
        Radius = radius
        Sigma = sigma
        EnergyRange = eSigma.keV
        let path = "/tmp/plots/"
        let suffix = &"radius_{radius:.0f}_sigma_{sigma:.0f}_energyRange_{eSigma:.1f}"
        let suffixTitle = &"Radius: {radius:.0f}, σ: {sigma:.0f}, ΔE: {eSigma:.1f}"
        echo "Generating plots for: ", suffixTitle
        for E in Es:
          kd.plotSingleEnergySlice(E.keV,
                                   outfile = path / &"back_interp_energy_{E}_{suffix}.pdf",
                                   title = &"Background interp, energy = {E} keV, {suffixTitle}")
          kd.plotGoldRegionBackgroundRate(outfile = path / &"background_gold_from_interp_{suffix}.pdf",
                                          title = &"Interp based gold background rate: {suffixTitle}")

    if true: quit()
    var treeDist = zeros[float]([256, 256, 10])
    echo "Start computationssss"
    for radius in [Radius]: #arange(10, 100, 10):
      let Es = linspace(0.0, 12.0, 10)
      echo "comp 3d dist"
      treeDist.compDistance3D(Es, kd, radius.float, byCount = false, metric = CustomMetric)
      echo "correct edges"
      treeDist.correctEdgeCutoff3D(radius.float)
      #  treeDist[_, _, E].plot2dTensor(
      #    &"/tmp/background_radius_byenergy_{E}_{radius}.pdf",
      #    &"k-d tree interpolation with radius: {radius} pixels, energy {E}")
      #treeDist.correctEdgeCutoff3D(radius.float)
      echo "plot 3d"
      #treeDist.plot3DTensor("/tmp/background_3d_radius_byenergy_" & $radius & ".pdf",
      #                      "k-d tree interpolation with radius: " & $radius & " pixels")
      treeDist.normalizeTensor(10, radius)
      treeDist.plot3DTensor("/tmp/background_3d_radius_byenergy_" & $radius & "_normalized.pdf",
                            "k-d tree interpolation with radius: " & $radius & " pixels, normalized")

      #treeDist.correctEdgeCutoff3D(radius.float)
      #treeDist.plot3DTensor("/tmp/background_radius_byenergy_correccted_" & $radius & ".pdf",
      #                      "k-d tree interpolation with radius: " & $radius & " pixels corrected by edges")


    # now plot interpolation based on energy

  #block kdTreeJustMoreStuff:
  #  let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
  #                     df["centerY", float].map_inline(toIdx(x).float)], axis = 1)
  #  let kd = kdTree(tTree, leafSize = 16, balancedTree = true)
  #  var treeDist = zeros[float]([256, 256, 10])
  #  let radius = 30
  #  let Es = linspace(0.0, 12.0, 10)
  #  for E in 0 ..< Es.high:
  #    let df = df.filter(f{`Energy` >= Es[E] and `Energy` < Es[E+1]})
  #    let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
  #                       df["centerY", float].map_inline(toIdx(x).float)], axis = 1)
  #    let kd = kdTree(tTree, leafSize = 16, balancedTree = true)
  #    for y in 0 ..< 256:
  #      for x in 0 ..< 256:
  #        let tup = kd.query_ball_point([x.float, y.float].toTensor, radius.float, metric = CustomMetric)#, metric = CustomMetric)
  #        let val = compValue(tup, byCount = true)
  #        treeDist[y, x, E] = val
  #
  #  treeDist.correctEdgeCutoff(radius.float)
  #  #  treeDist[_, _, E].plot2dTensor(
  #  #    &"/tmp/background_radius_byenergy_{E}_{radius}.pdf",
  #  #    &"k-d tree interpolation with radius: {radius} pixels, energy {E}")
  #  #treeDist.correctEdgeCutoff(radius.float)
  #  treeDist.plot3DTensor("/tmp/background_3d_radius_byenergy_notreally_" & $radius & ".pdf",
  #                        "k-d tree interpolation with radius: " & $radius & " pixels")
  #  #treeDist.correctEdgeCutoff(radius.float)
  #  #treeDist.plot3DTensor("/tmp/background_radius_byenergy_correccted_" & $radius & ".pdf",
  #  #                      "k-d tree interpolation with radius: " & $radius & " pixels corrected by edges")
  #
  #
  #  # now plot interpolation based on energy

  if true: quit()
  block kdTreeEnergy:
    df.showBrowser()
    discard df.toKDE(toPlot = true, outname = "/tmp/all_data.pdf")
    df.plotDf("all clusters", "/tmp/all_clusters.pdf")
    const radius = 50
    for (l, h) in [(0, 2), (2, 10)]: #0 ..< 10:
      let energy = l
      let df = df.filter(f{`Energy` >= energy.float and `Energy` < (energy + h).float})
      discard df.toKDE(toPlot = true,outname = &"/tmp/range_{energy}.pdf")
      df.plotDf(&"clusters in energy {energy}-{energy+h} keV", &"/tmp/clusters_{energy}.pdf")
      echo "Events left : ", df.len
      let tTree = stack([df["centerX", float].map_inline(toIdx(x).float),
                         df["centerY", float].map_inline(toIdx(x).float)], axis = 1)
      let kd = kdTree(tTree, leafSize = 16, balancedTree = true)

      var treeDist = zeros[float]([256, 256])
      treeDist.compDistance(kd, radius.float, byCount = true)
      treeDist.plot2dTensor(
        &"/tmp/background_energy_{energy}_radius_{radius}_bycount.pdf",
        &"k-d tree interp, radius: {radius} pixels, energy: {energy} - {energy+h} keV. # cluster: {df.len}")
      treeDist.correctEdgeCutoff(radius.float)
      treeDist.plot2dTensor(
        &"/tmp/background_energy_{energy}_radius_{radius}_bycount_corrected.pdf",
        &"k-d tree interp, radius: {radius} pixels, energy: {energy} - {energy+h} keV. # cluster: {df.len}, edge corrected")


      treeDist.compDistance(kd, radius.float)
      treeDist.plot2dTensor(
        &"/tmp/background_energy_{energy}_radius_{radius}.pdf",
        &"k-d tree interp, radius: {radius} pixels, energy: {energy} - {energy+h} keV. # cluster: {df.len}")
      treeDist.correctEdgeCutoff(radius.float)
      treeDist.plot2dTensor(
        &"/tmp/background_energy_{energy}_radius_{radius}_corrected.pdf",
        &"k-d tree interp, radius: {radius} pixels, energy: {energy} - {energy+h} keV. # cluster: {df.len}, edge corrected")

    if true: quit()
