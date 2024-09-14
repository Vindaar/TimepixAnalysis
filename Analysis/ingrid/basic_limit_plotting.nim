import pkg / [ggplotnim, seqmath]
proc plotToyLimits*(limits: seq[float], output: string) =
  const Bins = 100
  let expLim = limits.median
  let (hist, bins) = histogram(limits, bins = Bins) # to know maximum
  let lineTo = hist.max.float
  let df = toDf(limits)
  ggplot(df, aes("limits")) +
    geom_histogram(bins = Bins) +
    geom_linerange(aes = aes(x = expLim, y = 0.0, yMin = 0.0, yMax = lineTo),
                   color = some(parseHex("0000FF"))) +
    scale_y_continuous() +
    ylab("# toy limits") + xlab("g²") +
    annotate(text = "Expected limit",
             x = expLim,
             y = lineTo,
             rotate = -90.0,
             font = font(color = parseHex("0000FF")),
             backgroundColor = color(0.0, 0.0, 0.0, 0.0)) +
    ggsave(output)
