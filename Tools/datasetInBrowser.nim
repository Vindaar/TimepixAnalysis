import datamancer, nimhdf5, cligen
import ingrid / tos_helpers

proc printSummary(df: DataFrame, dataset: string) =
  echo "INFO: ", dataset
  echo "Number of elements: ", df.len
  echo "Min: ", df[dataset, float].min
  echo "Max: ", df[dataset, float].max
  echo "Mean: ", df[dataset, float].mean
  echo "Sum: ", df[dataset, float].sum
  echo "Variance: ", df[dataset, float].variance


proc main(fname: string, chip: int, dataset: string) =
  let h5f = H5open(fname, "r")

  let df = readDsets(h5f, chipDsets = some((chip: chip, dsets: @[dataset])))
  echo df
  printSummary(df, dataset)
  df.showBrowser()

when isMainModule:
  dispatch main
