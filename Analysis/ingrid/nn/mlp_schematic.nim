import ggplotnim, seqmath, sequtils
from ginger import ggColorHue

const PuOr_data = @[
    (0.49803921568627452,  0.23137254901960785,  0.03137254901960784),
    (0.70196078431372544,  0.34509803921568627,  0.02352941176470588),
    (0.8784313725490196 ,  0.50980392156862742,  0.07843137254901961),
    (0.99215686274509807,  0.72156862745098038,  0.38823529411764707),
    (0.99607843137254903,  0.8784313725490196 ,  0.71372549019607845),
    (0.96862745098039216,  0.96862745098039216,  0.96862745098039216),
    (0.84705882352941175,  0.85490196078431369,  0.92156862745098034),
    (0.69803921568627447,  0.6705882352941176 ,  0.82352941176470584),
    (0.50196078431372548,  0.45098039215686275,  0.67450980392156867),
    (0.32941176470588235,  0.15294117647058825,  0.53333333333333333),
    (0.17647058823529413,  0.0                ,  0.29411764705882354)
]

const PiYG_data = @[
    (0.55686274509803924,  0.00392156862745098,  0.32156862745098042),
    (0.77254901960784317,  0.10588235294117647,  0.49019607843137253),
    (0.87058823529411766,  0.46666666666666667,  0.68235294117647061),
    (0.94509803921568625,  0.71372549019607845,  0.85490196078431369),
    (0.99215686274509807,  0.8784313725490196 ,  0.93725490196078431),
    (0.96862745098039216,  0.96862745098039216,  0.96862745098039216),
    (0.90196078431372551,  0.96078431372549022,  0.81568627450980391),
    (0.72156862745098038,  0.88235294117647056,  0.52549019607843139),
    (0.49803921568627452,  0.73725490196078436,  0.25490196078431371),
    (0.30196078431372547,  0.5725490196078431 ,  0.12941176470588237),
    (0.15294117647058825,  0.39215686274509803,  0.09803921568627451)
    ]
const PRGn_data = @[
    (0.25098039215686274,  0.0                ,  0.29411764705882354),
    (0.46274509803921571,  0.16470588235294117,  0.51372549019607838),
    (0.6                ,  0.4392156862745098 ,  0.6705882352941176 ),
    (0.76078431372549016,  0.6470588235294118 ,  0.81176470588235294),
    (0.90588235294117647,  0.83137254901960789,  0.90980392156862744),
    (0.96862745098039216,  0.96862745098039216,  0.96862745098039216),
    (0.85098039215686272,  0.94117647058823528,  0.82745098039215681),
    (0.65098039215686276,  0.85882352941176465,  0.62745098039215685),
    (0.35294117647058826,  0.68235294117647061,  0.38039215686274508),
    (0.10588235294117647,  0.47058823529411764,  0.21568627450980393),
    (0.0                ,  0.26666666666666666,  0.10588235294117647)
    ]

import numericalnim except linspace

template getInterps(data: seq[tuple], num: int): untyped =
  let xs = linspace(0.0, 255.0, num)
  let rs = newLinear1D(xs, data.mapIt(it[0]))
  let gs = newLinear1D(xs, data.mapIt(it[1]))
  let bs = newLinear1D(xs, data.mapIt(it[2]))
  (rs, gs, bs)

#let (rs, gs, bs) = getInterps(PuOr_data, 11)
let (rs, gs, bs) = getInterps(PiYG_data, 11)
#let (rs, gs, bs) = getInterps(PRGn_data, 11)

proc getColor(w: float, scale: (float, float)): Color =
  let at = (w - scale[0]) / (scale[1] - scale[0]) * 255.0
  echo "AT: ", at
  result = color(rs.eval(at), gs.eval(at), bs.eval(at))

proc neurons(xpos, ypos, Δy: float, num: int, typ: string, names: seq[string] = @[]): DataFrame =
  var xs = newSeq[float]()
  var ys = newSeq[float]()
  var ns = newSeq[string]()
  let height = (num - 1).float * Δy
  let yStart = ypos - height / 2.0
  for i in 0 ..< num:
    xs.add xpos
    ys.add (yStart + i.float * Δy)
    if names.len > 0:
      doAssert i < names.len, "Not enough names for the given set of neurons! " & $names.len & " vs " & $num
      ns.add names[i]
    else:
      ns.add "Neuron " & $i
  result = toDf({xs, ys, "typ" : typ, "Names" : ns})

proc weights(plt: GgPlot, ip, op: DataFrame, weights: Tensor[float], typ: string): GgPlot =
  result = plt
  echo weights.shape
  echo weights
  echo ip
  echo op
  let sc = (weights.reshape(weights.size.int).percentile(5), weights.reshape(weights.size.int).percentile(95))
  for y in 0 ..< weights.shape[1]:
    for x in 0 ..< weights.shape[0]:
      let p1 = ip.row(y)
      let xpos1 = p1["xs"].toFloat() + 10.0
      let ypos1 = p1["ys"].toFloat()
      let p2 = op.row(x)
      let xpos2 = p2["xs"].toFloat() - 10.0
      let ypos2 = p2["ys"].toFloat()
      var xs = newSeq[float]()
      var ys = newSeq[float]()
      var ws = newSeq[float]()
      xs.add xpos1
      ys.add ypos1
      xs.add xpos2
      ys.add ypos2
      ws.add weights[x, y]
      ws.add weights[x, y]
      let df = toDf({xs, ys, "typ" : typ, ws})
      let c  = getColor(weights[x, y], sc)
      result = result + geom_line(data = df.clone(), color = c)

proc plotSchematic*(nI, h1, h2, nO: int, w1, w2, w3: Tensor[float],
                    dsets: seq[string],
                    file: string) =
  var df = newDataFrame()
  const
    xInput = 130
    xDelta = 110
    Width = 550
    Height = 400

  let input  = neurons(xInput,              200, 25.0, nI, "input", names = dsets)
  let h1     = neurons(xInput + xDelta,     200, 25.0, h1, "h1")
  let h2     = neurons(xInput + 2 * xDelta, 200, 25.0, h2, "h2")
  let output = neurons(xInput + 3 * xDelta, 200, 25.0, nO, "output")
  df = assignStack(@[input, h1, h2, output])

  let colors = ggColorHue(2, 285)

  # y position of first neuron in each layer
  const YOffset = 25.0
  let
    y0I  = input.row(0)["ys"].toFloat()  - YOffset
    y0H1 = h1.row(0)["ys"].toFloat()     - YOffset
    y0H2 = h2.row(0)["ys"].toFloat()     - YOffset
    y0O  = output.row(0)["ys"].toFloat() - YOffset

  let f = font(size = 8.0, bold = true)
  var plt = ggplot(df, aes("xs", "ys")) +
    geom_point(data = input,  size = 10.0, color = colors[0]) + #"lightsteelblue") +
    geom_point(data = h1,     size = 10.0, color = color(0.60, 0.60, 0.60)) +
    geom_point(data = h2,     size = 10.0, color = color(0.60, 0.60, 0.60)) +
    geom_point(data = output, size = 10.0, color = colors[1]) + # "royalblue") +
    geom_text(data = input, aes = aes(x = f{`xs` - 25}, text = "Names"), size = 3.0, alignKind = taRight) +
    geom_text(data = output, aes = aes(x = f{`xs` + 25}, text = "Names"), size = 3.0, alignKind = taLeft) +
    geom_text(aes = aes(x = xInput + 0 * xDelta, y = y0I , text = "Input layer"), font = f, alignKind = taCenter) +
    geom_text(aes = aes(x = xInput + 1 * xDelta, y = y0H1, text = "Hidden layer 1"), font = f, alignKind = taCenter) +
    geom_text(aes = aes(x = xInput + 2 * xDelta, y = y0H2, text = "Hidden layer 2"), font = f, alignKind = taCenter) +
    geom_text(aes = aes(x = xInput + 3 * xDelta, y = y0O , text = "Output layer"), font = f, alignKind = taCenter)
  # Add the weights lines one by one
  plt = plt.weights(input, h1, w1, "ih1")
  plt = plt.weights(h1,    h2, w2, "h1h2")
  plt = plt.weights(h2, output, w3, "h2o")
  plt +
    scale_x_continuous(breaks = linspace(0.0, Width, 9)) +
    scale_y_continuous(breaks = linspace(0.0, Height, 9)) +
    scale_y_reverse() +
    #scale_color_continuous(scale = (low: 0.0, high: 1.0)) +
    hideLegend() +
    theme_void() +
    xlim(0, Width) + ylim(0, Height) +
    margin(0.2, 0.2, 0.2, 0.2) +
    coord_fixed(1.0) +
    Theme(xtickLabelMargin: some(-2.0), ytickLabelMargin: some(2.0)) +
    ggsave(file, width = Width, height = Height)
