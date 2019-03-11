#' @export ZeroFilling

ZeroFilling <- function(Fid_data, fn = ncol(Fid_data)) {
  
  # Data initialisation and checks ----------------------------------------------
  begin_info <- beginTreatment("ZeroFilling", Fid_data)
  Spectrum_data <- begin_info[["Signal_data"]]
  checkArg(fn, c("num", "pos"))
  
  
  fn <- 2^(round(log(fn)/log(2))) # round fn to the nearest 2^x
  
  # Construct the matrix filled with 0 and adapt their colnames -------------------
  zero_matrix <- matrix(complex(real = 0, imaginary = 0), nrow = nrow(Fid_data), 
               ncol = fn) 
  sp_descr <- diff(as.numeric(colnames(Fid_data)))[1]
  colnames(zero_matrix) <- seq(from = max(as.numeric(colnames(Fid_data))) + sp_descr,
                      to = max(as.numeric(colnames(Fid_data))) + sp_descr * fn,
                      by = sp_descr)
  
  # Zero filling --------------------------------------
  mat_ZeroFilled <- cbind(Fid_data, zero_matrix)
 
  Fid_data <- mat_ZeroFilled
  
  # Data finalisation ----------------------------------------------
  Fid_data <- endTreatment("ZeroFilling", begin_info, Fid_data)
  
  return(Fid_data)
  
  }

