## Just a playground to play around with `ggplotnim` in combination with TimepixAnalysis.
## Might be an important basis for features required in `ggplotnim` for my work.

import ingrid / tos_helpers
import sequtils, seqmath, nimhdf5, strutils, tables, persvector
import os
import ggplotnim
import helpers / utils

import macros

macro echoType(x: typed): untyped =
  echo x.treeRepr

proc getDf(h5f: var H5FileObj, path: string, keys: varargs[string]): DataFrame =
  ## read the datasets form `path` in the `h5f` file and combine them to a single
  ## DataFrame
  var data: OrderedTable[string, seq[Value]]
  var size = 0
  for key in keys:
    static:
      echoType(key)
    var dsetH5 = h5f[(path / key).dset_str]
    if size == 0:
      size = dsetH5.shape[0]
    withDset(dsetH5):
      data[key] = % dset
  result = toDf(data)
  result.len = size


proc main(fname: string) =
  var h5f = H5file(fname, "r")
  defer: discard h5f.close()
  var
    totDuration = 0.0
    numEvents = 0
    numMax = 0
  for num, group in runs(h5f):
    for grp in items(h5f, group):
      if "fadc" notin grp.name:
        let chipNum = grp.attrs["chipNumber", int]
        if chipNum == 3:
          let df = h5f.getDf(grp.name,
                             concat(@["eventNumber"],
                                    @(getFloatGeometryNames()),
                                    @(getIntClusterNames())))
          echo df
          let fname = "figs/" & num
          echo "Saving as ", fname

          # small sizes
          let dfSmall = df.filter(f{"length" < 1.0})
          let dfOne = df.mutate(f{"smaller" ~ "length" < 1.0})
            #.group_by("smaller")
          echo dfOne


          ggplot(df, aes("length", "eccentricity")) +
            geom_point() +
            ggsave(fname & "le.svg")

          ggplot(dfOne, aes("length", "hits", color = "smaller")) +
            geom_point() +
            ggsave(fname & "length.svg")

          #ggplot(df) +
          #  geom_histogram(aes(x = f{) +
          #  ggsave(fname & "hitsSmall.svg")


when isMainModule:
  if paramCount() > 0:
    main(paramStr(1))
  else:
    echo "Hand a filename to be read from!"
