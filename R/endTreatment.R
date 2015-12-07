endTreatment <- function(name, begin_info, Signal_data) {
  end_time = proc.time() # record it as soon as possible
  start_time = begin_info[["start"]]
  delta_time = end_time - start_time
  delta = delta_time[]
  cat("End", name, "\n")
  cat("It lasted",
      round(delta["user.self"], 3), "s user time,",
      round(delta["sys.self"] , 3), "s system time and",
      round(delta["elapsed"]  , 3), "s elapsed time.\n")
  if (begin_info[["force.real"]]) {
    # The imaginary part is left untouched
    i <- complex(real=0, imaginary=1)
    Signal_data = Signal_data + i * Im(begin_info[["Original_data"]])
  }
  if (begin_info[["vec"]]) {
    Signal_data = Signal_data[1,]
  }
  return(Signal_data)
}