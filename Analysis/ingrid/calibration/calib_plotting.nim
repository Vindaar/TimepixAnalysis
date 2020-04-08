import strformat
import seqmath
import ggplotnim

#[
Contains all routines that create plots related to InGrid calibration.
]#

proc plotGasGain*[T](charge, counts: seq[T],
                     fitX, fitY: seq[T],
                     G_fit, chiSq: float,
                     chipNumber, runNumber: int) =
  ## given a seq of traces (polya distributions for gas gain) plot
  ## the data and the fit, save plots as svg.
  let dfRaw = seqsToDf({ "charge / e-" : charge,
                         "counts" : counts })
  let dfFit = seqsToDf({ "charge / e-" : fitX,
                         "counts" : fitY })
  let df = bind_rows([("Polya", dfRaw),
                      ("Fit", dfFit)],
                     id = "Type")
    # filter to max 2.5e4 electrons
    .filter(fn {c"charge / e-" <= 2.5e4})
  let G = histMean(counts, charge)
  ggplot(df, aes("charge / e-", "counts", color = "Type", fill = "Type")) +
    geom_histogram(data = df.filter(f{"Type" == "Polya"}),
                   stat = "identity") +
    geom_line(data = df.filter(f{"Type" == "Fit"})) +
    ggtitle(&"Polya fit of chip {chipNumber}, run {runNumber}: " &
            &"G = {G:.1f}, G_fit = {G_fit:.1f}, χ²/dof = {chiSq:.2f}") +
    ggsave(&"out/gas_gain_run_{runNumber}_chip_{chipNumber}.pdf")

## TODO: Put "Charge calibration factors vs gas gain. y errors magnified * 100"" into this module!
