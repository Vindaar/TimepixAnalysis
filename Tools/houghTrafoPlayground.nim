import ggplotnim
import cligen, sets, hashes, tables

#proc hashtest(num: int): DataFrame =
#  var ids = newSeqOfCap[string](num * num)
#  for i in 0 ..< num:
#    for j in 0 ..< num:
#      if i != j:
#        ids.add $i & "/" & $j
#  let c = toColumn ids
#  var s = newSeq[Hash](c.len)
#  s.hashColumn(c, finish = true)
#  var aset = initTable[int, string]()
#  for i, x in s:
#    if x notin aset:
#      aset[x] = ids[i]
#    else:
#      echo "Already contained! ", x
#      echo "with value ", aset[x]
#      echo "trying to insert ", ids[i]
#      echo "Compute hash manually ", hash(ids[i]), "  vs  ", hash(aset[x])
#  echo "length ", s.len
#  echo "set length ", s.toHashSet().card
#  echo "ids to hash set ", ids.toHashSet().card

proc houghTrafo[T: seq | Tensor](x, y: T): DataFrame =
  ## computes a hough transformation of `x` and `y`
  doAssert x.len == y.len
  var xs = newSeqOfCap[int](x.len * x.len)
  var ys = newSeqOfCap[int](x.len * x.len)
  var ids = newSeqOfCap[string](x.len * x.len)
  var slopes = newSeqOfCap[float](x.len * x.len)
  echo x
  for i in 0 ..< x.len:
    #echo "at ", i
    for j in 0 ..< x.len:
      if i != j:
        xs.add x[j]
        ys.add y[j]
        xs.add x[i]
        ys.add y[i]
        ids.add $i & "/" & $j
        ids.add $i & "/" & $j
        #echo "Ys ", ys[^1], "  ", ys[^2]
        #echo "Xs ", xs[^1], "  ", xs[^2]
        slopes.add (ys[^1] - ys[^2]).float / (xs[^1] - xs[^2]).float
  result = toDf(xs, ys, ids)
  let df = toDf(slopes).filter(f{`slopes` <= 3.0 and `slopes` >= -3.0})
  ggplot(df, aes("slopes")) +
    geom_histogram(bins = 100) +
    ggsave("/tmp/histo_slopes.pdf")
  echo "Mean of slopes ", df["slopes", float].mean

  echo result
  ggplot(result, aes("xs", "ys", color = "ids")) +
    geom_line() +
    ggsave("/tmp/lines.png", width = 4000, height = 2000)

proc main(fname: string) =
  let df = readcsv(fname)
    .arrange(["x", "y"])
  #discard hashtest(df["x", int], df["y", int])
  discard houghTrafo(df["x", int], df["y", int])
  ggplot(df, aes(x, y, color = charge)) +
    geom_point(size = some(1.0)) +
    xlim(0, 768) + ylim(0, 768) + scale_x_continuous() + scale_y_continuous() +
    geom_linerange(aes = aes(y = 0, xMin = 128, xMax = 640)) +
    geom_linerange(aes = aes(y = 256, xMin = 0, xMax = 768)) +
    geom_linerange(aes = aes(y = 512, xMin = 0, xMax = 768)) +
    geom_linerange(aes = aes(y = 768, xMin = 128, xMax = 640)) +
    geom_linerange(aes = aes(x = 0, yMin = 256, yMax = 512)) +
    geom_linerange(aes = aes(x = 256, yMin = 256, yMax = 512)) +
    geom_linerange(aes = aes(x = 512, yMin = 256, yMax = 512)) +
    geom_linerange(aes = aes(x = 768, yMin = 256, yMax = 512)) +
    geom_linerange(aes = aes(x = 128, yMin = 0, yMax = 256)) +
    geom_linerange(aes = aes(x = 384, yMin = 0, yMax = 256)) +
    geom_linerange(aes = aes(x = 640, yMin = 0, yMax = 256)) +
    geom_linerange(aes = aes(x = 128, yMin = 512, yMax = 768)) +
    geom_linerange(aes = aes(x = 384, yMin = 512, yMax = 768)) +
    geom_linerange(aes = aes(x = 640, yMin = 512, yMax = 768)) +
    margin(top = 1.5) +
    ggsave("/tmp/plot_s.pdf")

when isMainModule:
  dispatch main
