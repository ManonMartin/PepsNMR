######################################
# Analyse en composantes principales #
######################################

PrincipalComponentAnalysis <- function (xnorm,z,titlePlot=c("spectres pretraites amix ACP cpmg"),f=c(''),lenom=c("Spectra"),repertoiregraph=c("./"),afficher.couleur="Volunteer",afficher.pch="Volunteer",couleur=c(grDevices::rgb(0,0,0),grDevices::rgb(0,1,0),grDevices::rgb(1,0,0),grDevices::rgb(0,0,1),grDevices::rgb(0.5,0.5,0.5),grDevices::rgb(0.5,1,1),grDevices::rgb(1,0.5,1),grDevices::rgb(1,1,0.5),grDevices::rgb(0.5,0.5,1),grDevices::rgb(0.5,1,0.5),grDevices::rgb(0.3,0.7,0.3),grDevices::rgb(0.7,0.3,0.3),grDevices::rgb(0.3,0.3,0.7)))
{
  
  # Cette fonction réalise une analyse en composantes principales
  # Elle admet comme entrée:
  #                          - xnorm, les spectres à analyser
  #                          - z, la matrice de design
  #                          - f, le nom de chacun des spectres
  #                          - titlePlot, le nom de l'analyse
  
  # Ses sorties sont:
  #                   - une figure comprenant la projection des spectres sur le plan formé par la première et la deuxième dimension et un autre graphique représentant la projection des spectres sur un plan formé par la quatrième et la troisième composante principale.
  #                   - une figure comprenant quatre graphiques représentant la projection des loading.
  
  # Exemple :
  
  #xnorm=spectres.fusionnes[design[,6]==0,]
  #z=design.fusionnes[design[,6]==0,]
  #f=c('')
  #titlePlot=c("spectres pretraites amix ACP cpmg")
  #lenom=c("AmixCpmgACP")
  #afficher.couleur=c("Jour")
  #afficher.pch="Echantillon"
  #couleur=c(rgb(0,0,0),rgb(0,1,0),rgb(1,0,0),rgb(0,0,1),rgb(0.5,0.5,0.5))
  
  
  if(afficher.pch=="Volunteer"){ndesign.pch=3;valeurfacteur.pch=c(1:12)}
  if(afficher.pch=="Day"){ndesign.pch=4;valeurfacteur.pch=c("1","2","3")}
  if(afficher.pch=="Tube"){ndesign.pch=5;valeurfacteur.pch=c("1","2")}
  if(afficher.pch=="Time"){ndesign.pch=6;valeurfacteur.pch=c("1","2")}
  if(afficher.pch=="Outlier"){ndesign.pch=7;valeurfacteur.pch=c("0","1")}
  
  
  if(afficher.couleur=="Volunteer"){ndesign.couleur=3;valeurfacteur.couleur=c(1:12)}
  if(afficher.couleur=="Day"){ndesign.couleur=4;valeurfacteur.couleur=c("1","2","3")}
  if(afficher.couleur=="Tube"){ndesign.couleur=5;valeurfacteur.couleur=c("1","2")}
  if(afficher.couleur=="Time"){ndesign.couleur=6;valeurfacteur.couleur=c("1","2")}
  if(afficher.couleur=="Outlier"){ndesign.couleur=7;valeurfacteur.couleur=c("0","1")}
  
  
  #library(pcurve)
  cat(afficher.couleur)
  cat(afficher.pch)
  #1. Réalisation de l'ACP
  #_______________________
  
  respca <- stats::prcomp(xnorm) #fonction donnée par R (pas acp que l'on veut...)
  
  #Récupération des pcs :
  
  pcs <- respca$pcs
  
  #Récupération des valeurs propres :
  
  pcd <- diag(respca$d)
  nc <- length(pcd)
  
  #Récupération des vecteurs propres :
  
  pcv <- respca$v
  n <- dim(xnorm)[1]
  
  #Calcul des CP
  
  pcs <- pcs/matrix(pcd,nrow=n,ncol=nc,byrow=T)
  pcs
  
  #Calcul des pourcentages de variance des 4 premières valeurs propres
  
  d2 <- sum(pcd) #somme des valeurs propres
  pourc1 <- (pcd[1]/d2)*100 #valeur propre divisée par la somme (*100 pour le %)
  pourc2 <- (pcd[2]/d2)*100
  pourc3 <- (pcd[3]/d2)*100
  pourc4 <- (pcd[4]/d2)*100
  
  
  
  
  #2. Réalisatoin des projections
  #______________________________
  
  #Projection sur plan des 1-2 et 2-4 CP's
  
  # Projection sur le plan 1-2
  grDevices::dev.new(noRStudioGD = TRUE) 
  graphics::par(mfrow=c(1,3),oma=c(1,1,2,1))
  
  graphics::plot(pcs[,1],pcs[,2],xlab=paste("CP1 -",round(pourc1,2),"%"),
       ylab=paste("CP2 -",round(pourc2,2),"%"),pch='')
  graphics::abline(h=0,v=0)
  
  s <- dim(pcs)[1]
  op=c(1:s)
  # Affichage des couleurs
  for (i in op)
  {  # on boucle sur tous les spectres
    
    # On boucle sur chaque état ( "bubble", "bruker") du facteur qui détermine la couleur (ex: "logiciel")
    for(id.facteur.couleur in 1:length(valeurfacteur.couleur))
    {
      # Si le spectres correspond cet état
      if (T | z[i,ndesign.couleur]==valeurfacteur.couleur[id.facteur.couleur])
      {
        # Si l'option est activée, on affiche le n° du spectre sur la projection
        if(length(f)!=0)
          
        {graphics::text(pcs[i,1],pcs[i,2],f[i],cex=0.7,col=couleur[id.facteur.couleur])}
        
        # Si l'option est activée, l'on affiche des marqueurs sur les spectres
        if(afficher.pch!=0)
        {
          # On boucle sur chaque état( "cpmg","noesy","ste") du facteur qui détermine le marqueur (ex: "METHODE")
          for(id.facteur.pch in 1:length(valeurfacteur.pch))
          {
            if(z[i,ndesign.pch]==valeurfacteur.pch[id.facteur.pch])
            {graphics::points(pcs[i,1],pcs[i,2],pch=id.facteur.pch,col=couleur[id.facteur.couleur])}
          } # On ferme la boucle marqueur
        }
      } # On Passe à la couleur suivante
    } # On ferme la boucle couleur
    
  } # On ferme la boucle sur tous les spectres
  
  
  # La légende des deux graphes
  graphics::plot(1,type="n",axes=FALSE,xlab="",ylab="",ylim=c(0,2),xlim=c(0,2))
  graphics::legend("top",legend=valeurfacteur.pch,pch=1:length(valeurfacteur.pch), cex = 1,title=paste(afficher.pch,": "),bty = "n" )
  graphics::legend("center",legend=valeurfacteur.couleur,col=couleur, cex = 1,lty=1,title=paste(afficher.couleur,": "),bty="n")
  
  
  graphics::plot(pcs[,3],pcs[,4],xlab=paste("CP3 -",round(pourc3,2),"%"),
       ylab=paste("CP4 -",round(pourc4,2),"%"))
  graphics::abline(h=0,v=0)
  s <- dim(pcs)[1]
  op=c(1:s)
  # Affichage des couleurs
  for (i in op)
  {  # on boucle sur tous les spectres
    
    # On boucle sur chaque état ( "bubble", "bruker") du facteur qui détermine la couleur (ex: "logiciel")
    for(id.facteur.couleur in 1:length(valeurfacteur.couleur))
    {
      # Si le spectres correspond cet état
      if (T | z[i,ndesign.couleur]==valeurfacteur.couleur[id.facteur.couleur])
      {
        # Si l'option est activée, on affiche le n° du spectre sur la projection
        if(length(f)!=0)
          
        {graphics::text(pcs[i,3],pcs[i,4],f[i],cex=0.7,col=couleur[id.facteur.couleur])}
        
        # Si l'option est activée, l'on affiche des marqueurs sur les spectres
        if(afficher.pch!=0)
        {
          # On boucle sur chaque état( "cpmg","noesy","ste") du facteur qui détermine le marqueur (ex: "METHODE")
          for(id.facteur.pch in 1:length(valeurfacteur.pch))
          {
            if(z[i,ndesign.pch]==valeurfacteur.pch[id.facteur.pch])
            {graphics::points(pcs[i,3],pcs[i,4],pch=id.facteur.pch,col=couleur[id.facteur.couleur])}
          } # On ferme la boucle marqueur
        }
      } # On Passe à la couleur suivante
    } # On ferme la boucle couleur
    
  } # On ferme la boucle sur tous les spectres
  
  
  
  graphics::title(titlePlot,outer=T)
  grDevices::savePlot(paste(repertoiregraph,lenom,"Projection",".png",sep=""),type=c("png"))
  
  #Projection des loadings
  
  grDevices::dev.new(noRStudioGD = TRUE) 
  graphics::par(mfrow=c(4,1),mar=c(2,4,1,1))
  graphics::plot(pcv[,1],type="l")
  graphics::plot(pcv[,2],type="l")
  graphics::plot(pcv[,3],type="l")
  graphics::plot(pcv[,4],type="l")
  
  list(CP1=pourc1,CP2=pourc2,CP3=pourc3,CP4=pourc4)
  grDevices::savePlot(paste(repertoiregraph,lenom,"Loading",".png",sep=""),type=c("png"))
  
}