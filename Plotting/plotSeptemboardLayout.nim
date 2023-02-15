import strutils, seqmath
import ingrid / tos_helpers
import ggplotnim

## Just a test of plotting the septemboard layout

let UseRealLayout = true
const chipOutlineX = @[
  (x: 0,   yMin: 0, yMax: 255),
  (x: 255, yMin: 0, yMax: 255)
]
const chipOutlineY = @[
  (y: 0  , xMin: 0, xMax: 255),
  (y: 255, xMin: 0, xMax: 255)
]

proc addOutline(plt: var GgPlot) =
  proc toLayout(x: int, isX: bool, chip: int): int =
    result = if UseRealLayout: x.chpPixToRealPix(isX, chip) else: x.chpPixToSeptemPix(isX, chip)
  for chip in 0 .. 6:
    for line in chipOutlineX:
      let val  = line.x.toLayout(true, chip)
      let minV = line.yMin.toLayout(false, chip)
      let maxV = line.yMax.toLayout(false, chip)
      plt = plt + geom_linerange(aes = aes(x = val, yMin = minV, yMax = maxV))
    for line in chipOutlineY:
      let val  = line.y.toLayout(false, chip)
      let minV = line.xMin.toLayout(true, chip)
      let maxV = line.xMax.toLayout(true, chip)
      plt = plt + geom_linerange(aes = aes(y = val, xMin = minV, xMax = maxV))

# some dummy test data
let df = toDf({"x" : @[256, 124], "y": @[341, 211]})

var plt = ggplot(df, aes(x, y)) +
  xlim(0, 800) + ylim(0, 950) +
  scale_x_continuous() + scale_y_continuous()
plt.addOutline()
plt + ggsave("/t/test.pdf", height = 900, width = 800)


template dump(s: untyped): untyped =
  echo astToStr(s), " = ", s
dump Width
dump Height
dump BondHeight
dump FullHeight
# helper distances
dump YRow1Row2
dump YRow2Row3
# distances in X between chips on this row
dump Row1XDist
dump Row2XDist
dump Row3XDist
dump XSize
dump YSize
# offsets of the x positions of each row
dump Row1XOffset
dump Row2XOffset
dump Row3XOffset
# sizes in pixel of the full layout
#           v--- size of tight layout in pixel
#                       v--- size of chips in mm ('real size' of tight layout)
#                                   v--- size in mm of the real layout
dump XSizePix
dump YSizePix
