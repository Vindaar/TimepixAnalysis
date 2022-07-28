import std / [sequtils, strformat]
import nimhdf5
import ggplotnim
import cligen

type
  ## IMPORTANT: At the moment the order ``matters``! It must be the same order
  ## as in the file, because the data is parsed into this structure from binary
  ## using the data sizes of the types as offsets.
  InterTpx3Data = object
    data_header: uint8
    header: uint8
    hit_index: uint64
    x: uint8
    y: uint8
    TOA: uint16
    TOT: uint16
    EventCounter: uint16
    HitCounter: uint8
    FTOA: uint8
    scan_param_id: uint16
    chunk_start_time: cdouble
    iTOT: uint16
    TOA_Extension: uint64
    TOA_Combined: uint64

  Pixel = object
    x: uint8
    y: uint8
    TOT: uint16
  Cluster = seq[Pixel]

proc main(fname: string) =
  let h5f = H5file(fname, "r")
  let path = "interpreted/hit_data_0"
  let data = h5f[path, InterTpx3Data]
  var clusters = newSeq[Cluster]()
  var cluster = newSeq[Pixel]()
  # add new cluster if diff in time larger than 50 clock cycles
  const cutoff = 50
  var clusterTime = 0
  for i, el in data:
    if el.TOA.int > clusterTime + cutoff:
      if cluster.len > 0:
        clusters.add cluster
      cluster.setLen(0)
    cluster.add Pixel(x: el.x, y: el.y, TOT: el.TOT)
    clusterTime = el.TOA.int

  echo "found ", clusters.len, " cluster!"
  for i in 0 ..< clusters.len:
    let df = toDf({ "x" : clusters[i].mapIt(it.x.int),
                        "y" : clusters[i].mapIt(it.y.int),
                        "TOT" : clusters[i].mapIt(it.TOT.float) })
    if df.len > 5:
      ggplot(df, aes("x", "y", color = "TOT")) +
        geom_point() +
        xlim(0, 256) + ylim(0, 256) +
        theme_opaque() +
        ggsave(&"/tmp/events/event_{i}.png")

  discard h5f.close()

dispatch main
