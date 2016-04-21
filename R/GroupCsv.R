#' @export GroupCsv

GroupCsv = function(Spectrum_data, gr.file, csv2 = TRUE){ #comma as decimal point and a semicolon as field separator
  begin_info <- beginTreatment("GroupCsv", Spectrum_data, force.real=T)
  Spectrum_data <- begin_info[["Signal_data"]]
  checkArg(gr.file, "str",  can.be.null=TRUE)
  checkArg(csv2, "bool")
  
  if (csv2==TRUE){
    gr = read.csv2(gr.file)
  } else {gr = read.csv(gr.file)}
  
  gr[,2] = as.numeric(gr[,2])
  
  if(sum(is.na(gr[,2]))>0) {
    warning("NAs introduced by coercion due to non-numeric arguments in the second column containing the groups")
  }
  rn = data.frame(ID = rownames(Spectrum_data))
  colnames(gr)[1] = "ID"
  
  merged = merge(rn,gr,by="ID", sort=FALSE, by.x="ID", by.y = "ID")
  return(merged[,2])
}
