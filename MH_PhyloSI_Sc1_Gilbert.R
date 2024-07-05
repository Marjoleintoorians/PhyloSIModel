# Multi-host PhyloSI model
# "Host community structure can shape pathogen outbreak dynamics through a phylogenetic dilution effect"
# Application toGilbert et al (2012) crop and pests data
# July 2024

library(outliers)
library(phytools)
library(ape)
library(picante)
library(fields)
library(stats)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggfortify)
library(ggstatsplot)
library(dplyr)
library(car)
library(RColorBrewer) # colours!
library(MASS)
library(MuMIn)
library(binom)
library(bbmle)
library(emmeans)
library(visreg)
library(pwr)
require(gridExtra)
library(treebalance)

palette <- c("#d73027","#fc8d59","#fef590","#91bfdb","#4575b4")
##########################################################################

#############################################
# 1: Initialize and set parameters
#############################################

nspecvec = c(3,4,5,7,10)
nspecvec<-nspecvec[order(-nspecvec)] # Descending vector, for pca plotting purposes.
ngen = 100   # number of trees. FOR REGRESSIONS USE 100
delta <- 0.1 # disease-induced mortality
mu <- 0.2 # natural mortality rate
psi <- 0.01 # c_j: the interspecific probability of succesfull transmission j
bvec <- seq(0.1,1,0.1)

#############################################
# Pick organism
#############################################

# Gilbert et al (2012) parameters
# for bacteria:
A <- 3.26
B <- -2.97
# for viruses:
# A <- 8.44
# B <- -5.2
#### Appendix ####
# Funghi
#A <- 4.3961
#B <- -3.3249
# Oomycetes
# A <- 2.0763
# B <- -2.6679
# Insects
# A<- 3.2441
# B<- -2.7004
# Mites 
# A<- 1.9584
# B<- -2.1681
# Mollusks
# A <- -0.4667
# B <- -1.1391
# Nematodes
# A<-2.7157
# B<- -2.6249
#Plants
# A<-1.9979
# B<- -2.2775

set.seed(1)

# dataframe 
length_df <- ngen*length(nspecvec)*length(bvec)
out.total_raw <- matrix(0,nrow = length_df, ncol=8)
colnames(out.total_raw) <- c("PD","MPD","Ic","R0","b","S","PPd12","varPPd")

#############################################
# 2: Simulations host phylogenies and SI model
#############################################

for (k in 1:length(nspecvec)){
  
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
  
  # Initialize vectors to save all values from ngen trees
  out.eigen <- rep(0,ngen)
  out.pd <- rep(0,ngen)
  out.mpd <- rep(0,ngen)
  out.ppdvar <- rep(0,ngen)
  out.Ic <- rep(0,ngen) #matrix(0,nrow = length(bvec), ncol=ngen)
  # Initialize matrices
  cmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  Rmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  ppd <- array(rep(0, nspec*nspec*ngen), dim=c(nspec, nspec, ngen)) # initialize dataframe for all ngen
  kappa <- matrix(1, nrow=nspec,ncol=nspec)
  
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
    
    # LOOP to create ngen trees
    for (x in 1:ngen){
      tree_temp <- pbtree(b=bvalue, d=0, n=nspec) 
      
      # save data
      ppd_temp<-cophenetic(tree_temp)
      
      ### ADD GILBERT FUNCTION
      
      # Create c-matrix for each tree (ppd_temp[i,j])
      for (i in 1:nspec){
        for (j in 1:nspec){
          cmatrix[[i,j]] <- psi * (exp(A+B*log10(ppd_temp[i,j]+1))) / (1+exp(A+B*log10(ppd_temp[i,j]+1)))
        }
      }
      
      # Create R-matrix fo each tree
      for (i in 1:nspec){
        for (j in 1:nspec){
          Rmatrix[[i,j]] <- (kappa[i,j]*cmatrix[i,j]*pmatrix[i,j])/(delta+mu)
        }
      }
      
      ### Output for all ngen trees (x-loop)
      # PPd
      PPd12 <- ppd_temp["t1","t2"]
      # MPD & varPPd step
      ppd_temp_mpd <-ppd_temp[lower.tri(ppd_temp,diag=FALSE)]
      # Ic: Colless's I-statistic
      Ic <- collessI(tree_temp, method="original")
      # Total output dataset for PCA
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),1] = sum(tree_temp$edge.length) #PD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),2] = mean(ppd_temp_mpd) #MPD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),3] = Ic #Ic
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),4] = max(eigen(Rmatrix)$values) #R0
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),5] = bvalue # b-value
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),6] = nspec # number of species
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),7] = PPd12
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),8] = mean(var(ppd_temp_mpd))  #varppd
      
      # print(x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)))
    } # end of x-loop, for generation ngen trees
    
  } # end b loop (l-loop)
  
} # end nspec (k-loop)


#############################################
# 3: Tweak dataset
#############################################

# extract the values of the potential outliers based on the IQR criterium
out.total.Sc1<-as.data.frame(out.total_raw)
# Make some numbers characters for plotting
out.total.Sc1$bNumber <- as.character(out.total.Sc1$b)
out.total.Sc1$Species <- as.character(out.total.Sc1$S)
# Add log vars
out.total.Sc1$logMPD <- log(out.total.Sc1$MPD)


#############################################
# 4: Plots
#############################################

## Bacteria
# R0_MPD_b <- ggplot(out.total.Sc1, aes(logMPD,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
#   ggtitle("A) Bacteria: MPD") +
#   xlab("log(MPD)") + ylab(bquote(R[0])) +
#   scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
#   theme(legend.position = "none")
# R0_IC_b <- ggplot(out.total.Sc1, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
#   ggtitle("C) Bacteria: Imbalance") + ylab(bquote(R[0])) +
#   xlab("Ic") +
#   scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
#   theme(legend.position = "none")

## Viruses
R0_MPD_v <- ggplot(out.total.Sc1, aes(logMPD,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("B) Viruses: MPD") +
  xlab("log(MPD)") + ylab(bquote(R[0])) +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")
R0_IC_v <- ggplot(out.total.Sc1, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("D) Viruses: Imbalance") + ylab(bquote(R[0])) +
  xlab("Ic") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")

## Add PPd figure
ppd_bac <- exp(3.26-2.97*log10(out.total.Sc1$PPd12+1)) / (1+exp(3.26-2.97*log10(out.total.Sc1$PPd12+1)))
ppd_bac_merged <- data.frame(out.total.Sc1$PPd12,ppd_bac)
colnames(ppd_bac_merged) <- c("PPd","c")
ppd_vir <- exp(8.44-5.2*log10(out.total.Sc1$PPd12+1)) / (1+exp(8.44-5.2*log10(out.total.Sc1$PPd12+1)))
ppd_vir_merged <- data.frame(out.total.Sc1$PPd12,ppd_vir)
colnames(ppd_vir_merged) <- c("PPd","c")

gilfun <- ggplot() +
  geom_line(data = ppd_bac_merged, aes(x = PPd, y = c), color = "#d73027") +
  geom_line(data = ppd_vir_merged, aes(x = PPd, y = c), color = "#91bfdb") +
  ggtitle("E) Effect PPd on c with Gilbert et al (2012) data") + 
  ylab(bquote(c[ij])) + xlab(bquote(PPd[ij])) + theme_classic() +
  annotate(geom="text", x=50, y=0.5, label="viruses") + 
  annotate(geom="text", x=30, y=0.35, label="bacteria")

# Plot all panels
plot<-grid.arrange(R0_MPD_b, R0_IC_b,gilfun,R0_MPD_v, R0_IC_v, ncol=3, nrow=2)

## Appendix figure

# R0_MPD_f <- ggplot(out.total.Sc1, aes(logMPD,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
#   ggtitle("A) Funghi MPD") +
#   xlab("log(MPD)") + ylab(bquote(R[0])) +
#   scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
#   theme(legend.position = "none")
# R0_IC_f <- ggplot(out.total.Sc1, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
#   ggtitle("B) Funghi Ic") + ylab(bquote(R[0])) +
#   xlab("Ic") +
#   scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
#   theme(legend.position = "none")

#fun <- grid.arrange(R0_MPD_f, R0_IC_f, ncol=2)
#Oomycetes<- grid.arrange(R0_MPD_o, R0_IC_o, ncol=2)
#insects<- grid.arrange(R0_MPD_i, R0_IC_i, ncol=2)
#mites<- grid.arrange(R0_MPD_m, R0_IC_m, ncol=2)
#mollusks<- grid.arrange(R0_MPD_mo, R0_IC_mo, ncol=2)
#nematodes<- grid.arrange(R0_MPD_n, R0_IC_n, ncol=2)
#plants<- grid.arrange(R0_MPD_p, R0_IC_p, ncol=2)


# plot All
# grid.arrange(R0_MPD_f, R0_IC_f,
#              R0_MPD_o, R0_IC_o,
#              R0_MPD_i, R0_IC_i,
#              R0_MPD_m, R0_IC_m,
#              nrow=4,ncol=2)
# grid.arrange(R0_MPD_mo, R0_IC_mo,
#              R0_MPD_n, R0_IC_n,
#              R0_MPD_p, R0_IC_p,
#              nrow=3,ncol=2)


#############################################
# 6: Regressions
#############################################

ggplot(out.total.Sc1, aes(x=Ic, y=varPPd)) + geom_violin(trim=FALSE)

# Univariate models
z<-lm(R0 ~ S,data = out.total.Sc1)
AIC(z)
summary(z)

# Multivariate model
ModelA_Sc1 <- lm(R0 ~ scale(logMPD)+scale(Ic)+scale(S), data = out.total.Sc1, na.action = "na.fail")
summary(ModelA_Sc1)
AIC(ModelA_Sc1)



