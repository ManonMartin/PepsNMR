#' @export RepValid

RepValid = function(Spectrum_data, group=NULL, read = c("csv2", "csv", NULL), gr.file=NULL){
  begin_info <-  beginTreatment("RepValid", Spectrum_data, force.real=T)
  Spectrum_data <- Re(begin_info[["Signal_data"]])
  checkArg(group, "num", can.be.null=TRUE)
  read <- match.arg(read)
  checkArg(gr.file, "str",  can.be.null=TRUE)
  
  
  
  if(is.null(group)){
    if(is.null(gr.file) |  is.null(read)){warning("gr.file or read arguments miss- or unspecified")}
    
    switch(read,
           "csv2" = { # csv2
             gr = read.csv2(gr.file)
           },
           "csv" = { # csv
             gr = read.csv(gr.file)
           })
    
    gr[,2] = as.numeric(gr[,2])
    
    if(sum(is.na(gr[,2]))>0) {
      warning("NAs introduced by coercion due to non-numeric arguments in the second column containing the groups")
    }
    rn = data.frame(ID = rownames(Spectrum_data))
    colnames(gr)[1] = "ID"
    
    merged = merge(rn,gr,by="ID", sort=FALSE)
    group = merged[,2]
  }
  
  ### PCA
  res.pca = prcomp(Spectrum_data, retx = TRUE)
  

  scores <- data.frame(group, res.pca$x[,1:4])
  pc1.2 <- ggplot2::qplot(x=PC1, y=PC2, data=scores, colour=factor(group))
  pc3.4 <- ggplot2::qplot(x=PC3, y=PC4, data=scores, colour=factor(group))
  
  plot.pca = list(pc1.2, pc3.4)
  
  ### Inertia
  
  # le calcul du barycentre de chaque groupe :
  #############################################
  y = group
  x = Spectrum_data
  
  nGroup=length(unique(y))
  
  #Creation of group lists for data, y and barycenters
  m=dim(x)[2]
  n=dim(x)[1]
  
  yG=c() # number of the class
  nG=c() # number of elements for each group
  xG=vector("list", nGroup) # data for elements of i's group
  xmG=vector("list", nGroup) # barycenter for elements of i's group
  j=1
  for (i in unique(y)){
    yG[j]=i
    nG[j]=length(which(y==i))
    xG[[j]]=x[which(y==i),]
    xmG[[j]]=apply(xG[[j]],2,mean)
    j=j+1
  }
  
  # Raw mean
  G=apply(x,2,mean)
  
  #### Calculation of Between Inertia  (sum ng*d2(ug,u))
  IB=c()
  for (i in 1:length(unique(y))){
    ib = nG[i]*(sum((xmG[[i]]-G)^2))
    IB=c(IB,ib)
    InertiaB=sum(IB)
  }
  
  #### Calculation of Total Inertia
  xmG=matrix(G,nrow=n,ncol=m,byrow=TRUE)
  it=(x-xmG)^2
  IT=apply(it,1,sum)
  InertiaT = sum(IT)
  iner.inter <- (InertiaB/InertiaT) # Between-group Inertia (% of tot inertia) :
  iner.inter100=iner.inter*100
  
  
  #### Calculation of Within Inertia
  iner.intra = 1-iner.inter
  InertiaW=(iner.intra * InertiaT)
  iner.intra100=iner.intra*100 # Within-group Inertia (% of tot inertia) :
  
  
  PInertiaT = InertiaT*100/InertiaT
  PInertiaW = InertiaW*100/InertiaT
  PInertiaB = InertiaB*100/InertiaT
  
  
  res.inertia=matrix(data = round(c(InertiaB, InertiaW, InertiaT, PInertiaB, PInertiaW, PInertiaT),2), ncol = 3,
                     dimnames = list(c("Value", "Percentage"), c("BI", "WI", "TI")), byrow = TRUE)
  
  return(list(plot.pca = plot.pca, res.pca = res.pca, res.inertia = res.inertia))
  
  
}

