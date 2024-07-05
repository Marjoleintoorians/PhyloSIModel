# Multi-host PhyloSI model
# "Host community structure can shape pathogen outbreak dynamics through a phylogenetic dilution effect"
# Scenario 2 - Figures 2 and 3
# July 2024

library(phytools)
library(ape)
library(picante)
library(fields)
#library(apTreeshape)
library(treebalance) # new package for imbalance: collessI(tree, method = "original")
library(stats)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggfortify)
library(ggstatsplot)
library(dplyr)
library(RColorBrewer) # colours!
library(vcd) # mosaic function
library(car)
library(phytools)
library(MASS)
library(MuMIn)
palette <- c("#d73027","#fc8d59","#fef590","#91bfdb","#4575b4")
##########################################################################

#############################################
# 1: Initialize and set parameters
#############################################

nspecvec = c(3,4,5,7,10)
nspecvec<-nspecvec[order(-nspecvec)] # Descending vector, for plotting purposes.
ngen = 100   # number of trees
delta <- 0.1 # disease-induced mortality
mu <- 0.2 # natural mortality rate
psi <- 0.01
bvec <- seq(0.1,1,0.1)

# Initialize dataframe
length_df <- ngen*length(nspecvec)*length(bvec)
out.total_raw <- matrix(0,nrow = length_df, ncol=8)
colnames(out.total_raw) <- c("PD","MPD","Ic","ED", "R0","b","S","varPPd")

set.seed(1)

#############################################
# 2: Simulations host phylogenies and SI model
#############################################

for (k in 1:length(nspecvec)){
  # nr of species
  nspec <- nspecvec[k]
  
  #############################################
  #   Choose DD/FD and Substitition/Addition  
  #############################################
  
  # DD: Substitution
  # Ntot = 100
  # N <- rep(Ntot/nspec,nspec)
  # S <- N #*0.9 # 10% prevalence for each of the species
  
  # DD: Additive
  N <- rep(100,nspec)
  S <-N #rep(90,nspec)
  
  #############################################
  
  # FD: Replacing (Appendix)
  # Ntot = 100
  # N <- rep(Ntot/nspec,nspec)
  # S <- N
  
  #FD: Additive (Appendix)
  # N <- rep(100,nspec)
  # S <- N #rep(90,nspec)
  # Ntot <- sum(N)
  
  #############################################

  # Initialize matrices
  cmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  Rmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  ppd <- array(rep(0, nspec*nspec*ngen), dim=c(nspec, nspec, ngen)) # initialize dataframe for all ngen
  
  # kappa (K) for DD: K(Nj) = xi*Nj. Nj is crossed out in the pmatrix, so K = xi
  # kappa (K) is assumed the same for all species
  kappa <- matrix(1,nrow = nspec, ncol=nspec)
  
  # Create p-matrix
  pmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  for (i in 1:nspec){
    for (j in 1:nspec){
      pmatrix[[i,j]] <-  S[i] # DD (NOT 1/Nj)
      # pmatrix[[i,j]] <-  S[i]/(Ntot)  # FD
    }
  }
  
  # Loop over birth rates
  for (l in 1:length(bvec)){
    
    # set b
    bvalue <- bvec[l]
    
    # Loop to create ngen trees.
    for (x in 1:ngen){
      tree_temp <- pbtree(b=bvalue, d=0, n=nspec) 
      
      # pairwise distances
      ppd<-cophenetic(tree_temp)
      speciesnames <- colnames(ppd)
      
      # Assign a reservoir, random (1-nspec)
      resnumber <- sample.int(nspec,1)
      Res = speciesnames[resnumber]
      
      if(nspec==3){
        # Insert the prob suc trans from WHICH IS T3 --> RESERVOIR
        cmatrix[[1,1]] <- psi/(1+ppd["t1",Res]) #A
        cmatrix[[1,2]] <- psi/(1+ppd["t1",Res]) #A
        cmatrix[[1,3]] <- psi/(1+ppd["t1",Res]) #A
        cmatrix[[2,1]] <- psi/(1+ppd["t2",Res]) #B
        cmatrix[[2,2]] <- psi/(1+ppd["t2",Res]) #B
        cmatrix[[2,3]] <- psi/(1+ppd["t2",Res]) #B
        cmatrix[[3,1]] <- psi/(1+ppd["t3",Res]) #C
        cmatrix[[3,2]] <- psi/(1+ppd["t3",Res]) #C
        cmatrix[[3,3]] <- psi/(1+ppd["t3",Res]) #C
      }
      
      if(nspec==4){
        cmatrix[[1,1]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,2]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,3]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,4]] <- psi/(1+ppd["t1",Res])
        
        cmatrix[[2,1]] <- psi/(1+ppd["t2",Res]) 
        cmatrix[[2,2]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,3]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,4]] <- psi/(1+ppd["t2",Res])
        
        cmatrix[[3,1]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,2]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,3]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,4]] <- psi/(1+ppd["t3",Res])
        
        cmatrix[[4,1]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,2]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,3]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,4]] <- psi/(1+ppd["t4",Res])
      }
      
      if(nspec==5){
        cmatrix[[1,1]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,2]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,3]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,4]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,5]] <- psi/(1+ppd["t1",Res])
        
        cmatrix[[2,1]] <- psi/(1+ppd["t2",Res]) 
        cmatrix[[2,2]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,3]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,4]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,5]] <- psi/(1+ppd["t2",Res])
        
        cmatrix[[3,1]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,2]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,3]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,4]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,5]] <- psi/(1+ppd["t3",Res])
        
        cmatrix[[4,1]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,2]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,3]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,4]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,5]] <- psi/(1+ppd["t4",Res])
        
        cmatrix[[5,1]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,2]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,3]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,4]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,5]] <- psi/(1+ppd["t5",Res])
      }
      
      if(nspec==7){
        cmatrix[[1,1]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,2]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,3]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,4]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,5]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,6]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,7]] <- psi/(1+ppd["t1",Res])
        
        cmatrix[[2,1]] <- psi/(1+ppd["t2",Res]) 
        cmatrix[[2,2]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,3]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,4]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,5]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,6]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,7]] <- psi/(1+ppd["t2",Res])
        
        cmatrix[[3,1]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,2]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,3]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,4]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,5]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,6]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,7]] <- psi/(1+ppd["t3",Res])
        
        cmatrix[[4,1]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,2]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,3]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,4]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,5]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,6]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,7]] <- psi/(1+ppd["t4",Res])
        
        cmatrix[[5,1]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,2]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,3]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,4]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,5]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,6]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,7]] <- psi/(1+ppd["t5",Res])
        
        cmatrix[[6,1]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,2]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,3]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,4]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,5]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,6]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,7]] <- psi/(1+ppd["t6",Res])
        
        cmatrix[[7,1]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,2]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,3]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,4]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,5]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,6]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,7]] <- psi/(1+ppd["t7",Res])
      }
      
      if(nspec==10){
        cmatrix[[1,1]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,2]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,3]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,4]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,5]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,6]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,7]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,8]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,9]] <- psi/(1+ppd["t1",Res])
        cmatrix[[1,10]] <- psi/(1+ppd["t1",Res])
        
        cmatrix[[2,1]] <- psi/(1+ppd["t2",Res]) 
        cmatrix[[2,2]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,3]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,4]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,5]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,6]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,7]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,8]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,9]] <- psi/(1+ppd["t2",Res])
        cmatrix[[2,10]] <- psi/(1+ppd["t2",Res])
        
        cmatrix[[3,1]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,2]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,3]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,4]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,5]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,6]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,7]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,8]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,9]] <- psi/(1+ppd["t3",Res])
        cmatrix[[3,10]] <- psi/(1+ppd["t3",Res])
        
        cmatrix[[4,1]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,2]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,3]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,4]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,5]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,6]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,7]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,8]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,9]] <- psi/(1+ppd["t4",Res])
        cmatrix[[4,10]] <- psi/(1+ppd["t4",Res])
        
        cmatrix[[5,1]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,2]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,3]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,4]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,5]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,6]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,7]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,8]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,9]] <- psi/(1+ppd["t5",Res])
        cmatrix[[5,10]] <- psi/(1+ppd["t5",Res])
        
        cmatrix[[6,1]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,2]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,3]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,4]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,5]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,6]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,7]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,8]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,9]] <- psi/(1+ppd["t6",Res])
        cmatrix[[6,10]] <- psi/(1+ppd["t6",Res])
        
        cmatrix[[7,1]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,2]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,3]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,4]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,5]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,6]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,7]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,8]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,9]] <- psi/(1+ppd["t7",Res])
        cmatrix[[7,10]] <- psi/(1+ppd["t7",Res])
        
        cmatrix[[8,1]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,2]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,3]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,4]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,5]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,6]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,7]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,8]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,9]] <- psi/(1+ppd["t8",Res])
        cmatrix[[8,10]] <- psi/(1+ppd["t8",Res])
        
        cmatrix[[9,1]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,2]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,3]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,4]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,5]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,6]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,7]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,8]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,9]] <- psi/(1+ppd["t9",Res])
        cmatrix[[9,10]] <- psi/(1+ppd["t9",Res])
        
        cmatrix[[10,1]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,2]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,3]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,4]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,5]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,6]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,7]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,8]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,9]] <- psi/(1+ppd["t10",Res])
        cmatrix[[10,10]] <- psi/(1+ppd["t10",Res])
      }
      # Create R-matrix fo each tree
      for (i in 1:nspec){
        for (j in 1:nspec){
          
          Rmatrix[[i,j]] <- (kappa[i,j]*cmatrix[i,j]*pmatrix[i,j])/(delta+mu)
          
        }
      }
      
      # MPD
      ppd_mpd <-ppd[lower.tri(ppd,diag=FALSE)]
      
      # Ic: Colless's I-statistic
      Ic <- collessI(tree_temp, method="original")
      
      # tree depth, check
      tree_depth <- max(node.depth.edgelength(tree_temp))
      
      # Calculate Evolutionary distinctiveness (ED) of the Res
      ED <- evol.distinct(tree_temp, type = "fair.proportion",
                          scale = FALSE, use.branch.lengths = TRUE)   # ED: This returns a measure of species evolutionary distinctiveness
      rownames(ED) <-ED[,1]
     
       # What is the ED of Res?
      EDres_val <- ED[Res,2]
      
      # Total output dataset
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),1] = sum(tree_temp$edge.length) #PD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),2] = mean(ppd_mpd) #MPD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),3] = Ic
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),4] = EDres_val # ED of reservoir
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),5] = max(Re(eigen(Rmatrix)$values)) #R0
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),6] = bvalue # b-value
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),7] = nspec # number of species
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),8] = mean(var(ppd))  #varppd
    } # end of x-loop, for generation ngen trees
    
  } # end b loop (l-loop)
  
} # end nspec (k-loop)

#############################################
# 3: Tweak dataset
#############################################

# Add characters for number
out.total.Sc2 <- data.frame(out.total_raw)
out.total.Sc2$Species <- as.character(out.total.Sc2$S)
# New R0 char : ROUND DOWN R0
out.total.Sc2$R0rounded<-floor(out.total.Sc2$R0) #copy
out.total.Sc2$R0value <- as.character(out.total.Sc2$R0rounded)
# Add log vars
out.total.Sc2$logMPD <- log(out.total.Sc2$MPD)
out.total.Sc2$logED <- log(out.total.Sc2$ED)

#############################################
# 4: PCA
#############################################

pca2 <- prcomp(out.total.Sc2[,(1:5)],scale=TRUE)
# Plot (in manuscript)
pca.sc2<-autoplot(pca2, data = out.total.Sc2, label = F,loadings=T, loadings.label=T,
         loadings.colour = 'black',size = 1,loadings.label.colour = 'black', loadings.label.size = 6,
         colour = 'Species', main = "Scenario 2    -   Option 2") + theme_classic() +
         scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
         theme(legend.position = "none")
#plot both scenario 1 and scenario 2:
require(gridExtra)
grid.arrange(pca.sc1, pca.sc2, ncol=2)


# Sort per ED
ED_list<-sort(out.total.Sc2$ED, decreasing=F)
ED_list<-as.character(round(ED_list, digits = 0)) # sort this char list!
df <- data.frame(PD = round(out.total.Sc2$PD,digits=0), EDchar = ED_list, ED = out.total.Sc2$ED, 
                 R0char = out.total.Sc2$R0rounded, R0 = out.total.Sc2$R0)


#############################################
# 5: Plots
#############################################

R0_ED <- ggplot(out.total.Sc2, aes(log(ED),R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("E) Scenario 2: ED") +
  ylab(bquote(R[0])) + xlab("Log(ED)") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")
R0_MPD_2 <- ggplot(out.total.Sc2, aes(log(MPD),R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("B) Scenario 2: MPD") +
  ylab(bquote(R[0])) + xlab("Log(MPD)") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")
R0_Ic_2 <- ggplot(out.total.Sc2, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("D) Scenario 2: Imbalance") +
  ylab(bquote(R[0])) + xlab("Ic") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position=c(1.2,0.5),legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))

# Plot with Scenario 1 plots
require(gridExtra)
grid.arrange(R0_MPD_1, R0_IC_1, R0_ED, R0_MPD_2, R0_Ic_2, nrow = 2, ncol=3,
              top=textGrob("DD model - Addition"))


#############################################
# 6: Regressions
#############################################

# Univariate models
z<-lm(R0 ~ log(ED),data = out.total.Sc2)
summary(z)
AIC(z)

# Multivariate regressions

# Model A from Sc1
ModelA_Sc2 <- lm(R0 ~ scale(logMPD)+scale(Ic)+scale(S), data = out.total.Sc2,na.action = "na.fail")
summary(ModelA_Sc2)
AIC(ModelA_Sc2)

# Model B - Remove MPD due to collinearity!
ModelB <- lm(R0 ~scale(S)+scale(logED)+scale(Ic),
             data = out.total.Sc2,na.action="na.fail")
summary(ModelB) # Back to R2=0.95 when including log(ED)
AIC(ModelB)

