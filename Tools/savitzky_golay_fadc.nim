import ingrid / tos_helpers
import ggplotnim, nimhdf5, sequtils
import scinim / signals


proc main(f: string, run: int, window = 15, order = 4) =
  let h5f = H5open(f, "r")

  let fadcDset = h5f[fadcDataBasename(run).dset_str]
  let data = fadcDset[float].toTensor.reshape(fadcDset.shape)

  let dat = data[0, _].squeeze
  let df = toDf({ "time" : toSeq(0 ..< 2560), "fadc" : dat,
                  "filtered" : savitzkyGolayFilter(dat, window, order) })
  #  .gather(["fadc", "filtered"], "Type", "Value")


  ggplot(df, aes("time")) +
    geom_line(aes = aes(y = "fadc")) +
    geom_line(aes = aes(y = "filtered"), color = "purple") +
    #geom_point(size = 1.5) +
    ggsave("/t/test_fadc.pdf")


when isMainModule:
  import cligen
  dispatch main
