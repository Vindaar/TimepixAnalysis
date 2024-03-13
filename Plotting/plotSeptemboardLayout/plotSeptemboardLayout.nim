import std / os except FileInfo
import std / strutils
import ingrid / [tos_helpers, ingrid_types]
import nimhdf5, ggplotnim, ginger

## The Septemboard layout code is a port of the code used in the python based event
## display for TOS.

const
  Width = 14.1
  Height = 14.1
  BondHeight = 2.0
  FullHeight = Height + BondHeight
  NumChips = 7

  # If this is set to `true` the final plot will only contain the actual raster image. No legend or axes
  OnlyRaster = true

  Run = 241

type
  SeptemRow = object
    left: float
    right: float
    wspace: float
    top: float
    bottom: float

proc initSeptemRow(nChips: int, x_size, y_size, x_dist, x_offset, y_t_offset, y_b_offset, dist_to_row_below: float): SeptemRow =
  # this class implements a single row of chips of the septem board
  # nChips: number of chips in row
  # x_dist: distance in x direction between each chip
  # x_offset: offset of left edge of first chip in row from
  #           left side of center row
  # calculate width and height of row, based on chips and dist
  let width = nChips.float * Width + (nChips - 1).float * x_dist
  let height_active = Height
  let height_full = FullHeight + dist_to_row_below
  # using calc gridspec one calculates the coordinates of the row on
  # the figure in relative canvas coordinates
  # include padding by adding or subtracting from left, right, top, bottom
  result.left      = x_offset / x_size
  result.right     = result.left + width / x_size
  result.wspace    = x_dist / x_size
  result.top       = 1.0 - y_t_offset / y_size
  result.bottom    = result.top - height_active / y_size

proc initSeptemBoard(padding, fig_x_size, fig_y_size, scaling_factor: float): seq[SeptemRow] =
  # implements the septem board, being built from 3 septem rows

  proc initRows(y_size, scaled_x_size, scaled_y_size, y_row1_row2, y_row2_row3, row2_x_dist: float): seq[SeptemRow] =
    # this function creates the row objects for the septem class
    # calculation of row 1 top and bottom (in abs. coords.):
    let
      # (top need to add padding to top of row 1)
      row1_y_top    = y_size - BondHeight - Height
      # bottom in abs. coords.
      row1_y_bottom = 2 * FullHeight + y_row1_row2 + y_row2_row3 - Height
      # offset of left side from septem in abs. coords.
      row1_x_offset = 6.95
    # now create the first row with all absolute coordinates
    result.add initSeptemRow(2, scaled_x_size, scaled_y_size, 0.85, row1_x_offset, row1_y_top, row1_y_bottom, y_row1_row2)
    # calculation of row 2 top and bottom (top & bottom of row2 not affected by padding):
    let
      row2_y_top    = y_size - FullHeight - y_row1_row2 - Height
      row2_y_bottom = FullHeight + y_row2_row3 + BondHeight - Height
      # no offset for row2, defines our left most position in abs. coords.
      row2_x_offset = 0.0 #padding * x_size
    result.add initSeptemRow(3, scaled_x_size, scaled_y_size, row2_x_dist, row2_x_offset, row2_y_top, row2_y_bottom, y_row2_row3)
    # calculation of row 3 top and bottom (add padding to bottom):
    let
      row3_y_top    = y_size - 2 * FullHeight - y_row1_row2 - y_row2_row3 - Height
      row3_y_bottom = BondHeight - Height
      row3_x_offset = 7.22
    result.add initSeptemRow(2, scaled_x_size, scaled_y_size, 0.35, row3_x_offset, row3_y_top, row3_y_bottom, 0)

  # include a padding all around the septem event display of 'padding'
  # use size of figure to scale septem accordingly to have it always properly
  # scaled for the given figure
  # take the inverse of the scaling factor (want 1/2 as input to scale to half size)
  let scaling_factor = 1.0 / scaling_factor
  # first calculate the ratio of the figure
  let fig_ratio = float(fig_x_size) / float(fig_y_size)
  # distances between different rows in absolute coordinates
  let
    y_row1_row2 = 0.38
    y_row2_row3 = 3.1
    # size in y direction of whole septem board in absolute coordinates
    y_size      = 3 * FullHeight + y_row1_row2 + y_row2_row3
    # already define row2_x_dist here (in absolute coordinates) to calculate x_size
    row2_x_dist = 0.35
    # 3 chips * width + 2 * distance between chips (in absolute coordinates)
    x_size      = 3 * Width + (3 - 1) * row2_x_dist
  # calculate the ratio of the septem board
  var ratio = float(x_size) / float(y_size)
  # now calculate the needed ratio to get the correct scaling of the septem on any
  # figure scale. fig_ratio / own ratio
  ratio = fig_ratio / ratio
  let
    # scaled x and y sizes
    scaled_x_size = x_size * ratio * scaling_factor
    scaled_y_size = y_size * scaling_factor
  # and now create the row objects
  result = initRows(y_size, scaled_x_size, scaled_y_size, y_row1_row2, y_row2_row3, row2_x_dist)

proc readVlen(h5f: H5File,
              fileInfo: FileInfo,
              runNumber: int,
              dsetName: string,
              chipNumber = 0,
              dtype: typedesc = float): seq[seq[dtype]] =
  ## reads variable length data `dsetName` and returns it
  ## In contrast to `read` this proc does *not* convert the data.
  let vlenDtype = special_type(dtype)
  let dset = h5f[(fileInfo.dataPath(runNumber, chipNumber).string / dsetName).dset_str]
  result = dset[vlenDType, dtype]

proc calcOccupancy[T](x, y: seq[seq[T]], z: seq[seq[uint16]] = @[]): Tensor[float] =
  ## calculates the occupancy of the given x and y datasets
  ## Either for a `seq[seq[T: SomeInteger]]` in which case we're calculating
  ## the occupancy of a raw clusters or `seq[T: SomeFloat]` in which case
  ## we're dealing with center positions of clusters
  result = newTensor[float]([NPix, NPix])
  # iterate over events
  for i in 0 .. x.high:
    let
      xEv = x[i]
      yEv = y[i]
    var zEv: seq[uint16]
    if z.len > 0:
      zEv = z[i]
    ## continue if full event.
    ## TODO: replace by solution that also works for clusters!!
    #if xEv.len >= 4095: continue
    for j in 0 .. xEv.high:
      if zEv.len > 0:
        result[xEv[j].int, yEv[j].int] += zEv[j].float
      else:
        result[xEv[j].int, yEv[j].int] += 1.0

proc occForChip(h5f: H5File, chip: int, fileInfo: FileInfo): (Tensor[int], Tensor[int], Tensor[float]) =
  const NPix = 256
  let
    xD = h5f.readVlen(fileInfo, Run, "x", chip, dtype = uint8)
    yD = h5f.readVlen(fileInfo, Run, "y", chip, dtype = uint8)
    zD = h5f.readVlen(fileInfo, Run, "ToT", chip, dtype = uint16)
  let occ = calcOccupancy(xD, yD) # , zD)
  var
    x = newTensorUninit[int](NPix * NPix)
    y = newTensorUninit[int](NPix * NPix)
    z = newTensorUninit[float](NPix * NPix)
  var i = 0
  for idx, val in occ:
    x[i] = idx[0]
    y[i] = idx[1]
    z[i] = val
    inc i
  result = (x, y, z)

proc handleOccupancy(h5f: H5File,
                     chip: int,
                     fileInfo: FileInfo,
                     quant: float = 0.0): PlotView =
  # get x and y datasets, stack and get occupancies
  let (x, y, z) = h5f.occForChip(chip, fileInfo)
  let df = toDf(x, y, z)
  var quant = quant
  if quant == 0.0:
    quant = percentile(z, 80)
  result = ggcreate(
    block:
      var plt =
        ggplot(df, aes("x", "y", fill = "z"), backend = bkCairo) +
          geom_raster() +
          scale_fill_continuous(scale = (low: 0.0, high: quant)) + #high: 1600.0)) +
          xlim(0, NPix) + ylim(0, NPix)
      if OnlyRaster:
        plt = plt + theme_void() + hideLegend()
      plt
  )
  #ggplot(df, aes("x", "y", fill = "z"), backend = bkCairo) +
  #    geom_raster() +
  #    scale_fill_continuous(scale = (low: 0.0, high: 1600.0)) +
  #    xlim(0, NPix) + ylim(0, NPix) +
  #    ggsave("/t/test_occ_0_.pdf")

proc drawBounds(v: Viewport) =
  v.drawBoundary(writeName = true)
  for ch in mitems(v.children):
    ch.drawBounds()

proc calcQuantileChip3(h5f: H5File, fileInfo: FileInfo): float =
  let (x, y, z) = h5f.occForChip(3, fileInfo)
  result = percentile(z, 80)

proc addRow(view: Viewport, h5f: H5File, septem: seq[SeptemRow], fileInfo: FileInfo, i, num, chipStart: int, showEmpty = false) =
  let width = septem[i].right - septem[i].left
  let height = septem[i].top - septem[i].bottom

  var row = view.addViewport(left = septem[i].left, bottom = septem[i].bottom,
                             width = width, height = height)
  row.layout(num, 1) #, margin = quant(septem[i].wspace, ukRelative))

  #let quant = calcQuantileChip3()

  for j in 0 ..< num:
    if not showEmpty:
      let plt = handleOccupancy(h5f, chipStart + j, fileInfo) #, quant)
      let v = if OnlyRaster: plt.view[4] else: plt.view
      var pltView = v.relativeTo(row[j])
      row.embedAt(j, pltView)
    row[j].drawBoundary()

  view.children.add row

## Read the data from the reconstructed H5 file of run 241
const path = "/mnt/1TB/CAST/2017/development/reco_$#_sparking.h5" % $Run
let h5f = H5open(path, "r")
let fileInfo = getFileInfo(h5f)

let
  fig_x_size = 10.0
  fig_y_size = 12.04186
  ratio = fig_x_size / fig_y_size

let septem = initSeptemBoard(0.0, fig_x_size, fig_y_size, 1.0)
let ImageSize = fig_x_size * DPI
let view = initViewport(c(0.0, 0.0),
                        quant(fig_x_size, ukInch), quant(fig_y_size, ukInch), backend = bkCairo,
                        wImg = ImageSize, hImg = ImageSize / ratio)

view.addRow(h5f, septem, fileInfo, 0, 2, 5, showEmpty = true)
view.addRow(h5f, septem, fileInfo, 1, 3, 2)
view.addRow(h5f, septem, fileInfo, 2, 2, 0)

view.draw("/home/basti/phd/Figs/detector/sparking/sparking_occupancy_80_quantile_run_$#.pdf" % $Run)
