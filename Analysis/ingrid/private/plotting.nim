#[
This is a module for general plotting functionality that is not directly tied
to specific calibration steps (as the code in `../calibration/calib_plotting` is!)
]#

import ggplotnim, arraymancer, os, strformat

proc plotOccupancy*[T](occ: Tensor[T], path: string, run, chip: int,
                       batchNum: int) =
  ## plots an occupancy plot and puts it to `fname`
  # Need to convert the given occupancy tensor to a DF first
  const NPix = 256
  echo "Creating directory: ", path.parentDir
  createDir(path)
  let fname = path / &"occupancy_run_{run}_chip_{chip}_{batchNum}.pdf"
  var
    x = newTensorUninit[int](NPix * NPix)
    y = newTensorUninit[int](NPix * NPix)
    z = newTensorUninit[float](NPix * NPix)
  var i = 0
  for idx, val in occ:
    x[i] = idx[0].int
    y[i] = idx[1].int
    z[i] = val.float
    inc i
  let df = seqsToDf(x, y, z)
  echo df
  if occ.max > T(0):
    ggplot(df, aes("x", "y", fill = "z")) +
      geom_raster() +
      scale_fill_continuous(scale = (low: 0.0, high: occ.max.float)) +
      xlim(0, NPix) + ylim(0, NPix) +
      ggtitle(&"Occupancy of run {run} for chip {chip}") +
      ggsave(fname, height = 800, width = 800)
