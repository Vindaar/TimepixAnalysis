template plotHandler(body: untyped): untyped {.dirty.} =
  (proc(h5f: H5FileObj, fileInfo: FileInfo, pd: PlotDescriptor, config: Config): (string, PlotV) =
     body)


proc moreCustom(fileInfo: FileInfo, config: Config): seq[PlotDescriptor] =
  let selector = initSelector(config)
  let customPlot = CustomPlot(kind: cpHistogram, x: "ToA")
  var pd = PlotDescriptor(runType: fileInfo.runType,
                            name: "ToA",
                            xLabel: "ToA",
                            yLabel: "Count",
                            title: "ToA histogramming",
                            selector: selector,
                            isCenterChip: true,
                            chip: -1, ## will be set to center chip when reading data!
                            runs: fileInfo.runs,
                            plotKind: pkCustomPlot,
                            customPlot: customPlot)
  pd.processData = plotHandler:
    var allData = newSeq[int32]()
    echo "Reading all TOA"
    for r in pd.runs:
      echo "Reading toa of run ", r
      let toas = h5f.readVlen(r, "ToA", pd.selector, fileInfo.centerChip, dtype = uint16)
      # shift all ToA values to start at 0 in each cluster
      for toa in toas:
        let toaInt = toa.mapIt(it.int32)
        let minToa = toaInt.min()
        allData.add toaInt.mapIt(it - minToa)
    result[0] = buildOutfile(pd, fileDir, fileType)
    let title = buildTitle(pd)
    echo "Plotting hist now"
    result[1] = plotHist(allData, title, pd.name, result[0], binS = 1.0,
                         binR = (-0.5, 100.0))
  result = @[pd]
