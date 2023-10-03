import std / algorithm
import arraymancer / tensor

import ../ [ingrid_types]
import ./geometry


proc getCenterClusterData*(centerData: CenterChipData, idx: int): CenterClusterData =
  ## Returns the data of the center cluster from the `centerData` and the event index.
  result = CenterClusterData(
    logL: centerData.logL[idx],
    energy: centerData.energies[idx],
    cX: centerData.cX[idx] + 14.0, # assign X, Y but convert to global septem coordinates
    cY: centerData.cY[idx] + 14.0,
    hits: centerData.hits[idx],
    rmsT: centerData.rmsT[idx],
    rmsL: centerData.rmsL[idx]
  )

proc containsOriginal*(cl: ClusterObject[PixInt], septemFrame: SeptemFrame): bool =
  ## Checks whether the original cluster (OC) passing lnL is contained in the
  ## given cluster `cl`, making this a "Hypothetical Larger Cluster" (HLC).
  let testPix = septemFrame.centerCluster[septemFrame.centerCluster.high div 2] # we can take any pixel as we don't drop pixels
  result = testPix in cl.data

proc isOriginal*(cl: ClusterObject[PixInt], septemFrame: SeptemFrame): bool =
  ## Checks whether the given cluster is *exactly* the original cluster
  ## (not very efficient, in the case where they contain exactly the same number
  ## of pixels, but otherwise just a check on pixel count)
  if cl.data.len != septemFrame.centerCluster.len:
    result = false
  else:
    let clSort = cl.data.sorted
    let origSort = septemFrame.centerCluster.sorted
    result = true
    for idx in 0 ..< clSort.len:
      if clSort[idx] != origSort[idx]:
        return false

proc getCenterClusterData*(septemFrame: SeptemFrame,
                           centerData: CenterChipData,
                           recoEv: RecoEvent[PixInt],
                           lineVetoKind: LineVetoKind
                         ): CenterClusterData =
  ## Returns the data of the correct cluster depending on the LineVetoKind
  case lineVetoKind
  of lvRegular, lvRegularNoHLC: result = getCenterClusterData(centerData, septemFrame.centerEvIdx)
  of lvCheckHLC:
    ## In this case we have to identify the HLC from _all_ clusters that were
    ## identified in the Septemevent. This means checking if _any_ pixel from the
    ## original cluster is found in a cluster. If any is found, this is the HLC.
    # done by walking all `recoEv` cluster and checking in any if _a_ pixel from
    # orignial center cluster is contained
    for cl in recoEv.cluster:
      if cl.containsOriginal(septemFrame):
        # found it, fill result
        result = CenterClusterData(
          logL: Inf,        # irrelevant and not used!
          energy: Inf,       # irrelevant and not used!
          cX: cl.centerX,
          cY: cl.centerY,
          hits: cl.data.len,
          rmsT: cl.geometry.rmsTransverse,
          rmsL: cl.geometry.rmsLongitudinal
        )
        return result # can return early here.
    doAssert false, "We cannot be here. This implies we *DID NOT FIND* the original cluster in the reconstructed data!"
  of lvNone:
    doAssert false, "This is not intended as a usable veto kind!"


proc getPixels*(allChipData: AllChipData, chip, idx: int, chargeTensor: var Tensor[float]): PixelsInt =
  ## Use `idx` of this cluster & event to look up all x, y, ToT and charge data of the cluster
  let
    chX = allChipData.x[chip][idx].unsafeAddr.toDataView()
    chY = allChipData.y[chip][idx].unsafeAddr.toDataView()
    chToT = allChipData.ToT[chip][idx].unsafeAddr.toDataView()
    chCh = allChipData.charge[chip][idx].unsafeAddr.toDataView()
  let numPix = allChipData.x[chip][idx].len # data view has no length knowlegde!
  result = newSeq[PixInt](numPix)
  for j in 0 ..< numPix:
    let pix = (x: chX[j], y: chY[j], ch: chToT[j]).chpPixToSeptemPix(chip)
    # add current charge into full septem tensor
    chargeTensor[pix.y, pix.x] += chCh[j]
    result[j] = pix
