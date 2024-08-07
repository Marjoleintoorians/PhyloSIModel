# Multi-host PhyloSI model
# "Host community structure can shape pathogen outbreak dynamics through a phylogenetic dilution effect"
# Scenario 1 - Figures 1, 2 and 3
# July 2024

library(outliers)
library(phytools)
library(ape)
library(picante)
library(fields)
library(treebalance) # new package for imbalance: collessI(tree, method = "original")
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
palette <- c("#d73027","#fc8d59","#fef590","#91bfdb","#4575b4")
##########################################################################

#############################################
# 1: Initialize and set parameters
#############################################

nspecvec = c(3,4,5,7,10)
nspecvec<-nspecvec[order(-nspecvec)] # Descending vector, for pca plotting purposes.
ngen = 100 # number of trees. FOR REGRESSIONS USE 100
delta <- 0.1 # disease-induced mortality
mu <- 0.2 # natural mortality rate
psi <- 0.01
bvec <- seq(0.1,1,0.1)

set.seed(1)

# dataframe 
length_df <- ngen*length(nspecvec)*length(bvec)
out.total_raw <- matrix(0,nrow = length_df, ncol=8)
colnames(out.total_raw) <- c("PD","MPD","Ic","R0","b","S","PPd12","varPPd")

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
  
  # Initialize vectors to save all values from ngen trees
  out.eigen <- rep(0,ngen)
  out.pd <- rep(0,ngen)
  out.mpd <- rep(0,ngen)
  out.ppdvar <- rep(0,ngen)
  out.Ic <- rep(0,ngen) 
  # Initialize matrices
  cmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  Rmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  ppd <- array(rep(0, nspec*nspec*ngen), dim=c(nspec, nspec, ngen)) # initialize dataframe for all ngen
  
  # kappa (K) for DD: K(Nj) = xi*Nj. Nj is crossed out in the pmatrix, so K = xi
  # kappa (K) is assumed the same for all species
  kappa <- matrix(1,nrow = nspec, ncol=nspec)
  
  # Create p-matrix for FD (and DD)
  pmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  for (i in 1:nspec){
    for (j in 1:nspec){
      pmatrix[[i,j]] <-  S[i] # DD
      #pmatrix[[i,j]] <-  S[i]/(Ntot)  # FD
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
      
      # Create c-matrix for each tree
      for (i in 1:nspec){
        for (j in 1:nspec){
          cmatrix[[i,j]] <- psi/(1+ppd_temp[i,j])
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
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),3] = Ic
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),4] = max(eigen(Rmatrix)$values) #R0
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),5] = bvalue # b-value
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),6] = nspec # number of species
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),7] = PPd12
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),8] = mean(var(ppd_temp_mpd))  #varppd
      
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
# 4: PCA
#############################################

# PCA 1 (only columns 1:5)
pca1 <- prcomp(out.total.Sc1[,(1:4)],scale=TRUE)
# Plot (in manuscript)
pca.sc1<-autoplot(pca1, data = out.total.Sc1, label = F,loadings=T, loadings.label=T,
        loadings.colour = 'black',size = 1, loadings.label.colour = 'black', loadings.label.size = 6,
        colour = 'Species', show.legend = FALSE, main = "Scenario 1") + theme_classic() +
        scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
        theme(legend.position=c(0.08,0.8),legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
scale_color_brewer(palette = "RdYlBu") +theme(legend.position = "none") #+scale_y_reverse()
# New dataset including the PC1, PC2 and PC3
out.total.Sc1 <- data.frame(cbind(out.total.Sc1, pca1$x[,1:3]))

#############################################
# 5: Plots
#############################################

# Figure 1 manuscript
plot1_12 <- ggplot(out.total.Sc1, aes(x = PPd12, y = R0, color = MPD)) +
  geom_point(size=0.7) + ggtitle("A) Scenario 1") +
  xlab(bquote(PPd["1,2"])) + ylab(bquote(R[0]))  +
  theme_classic() +  theme(legend.position="none") + scale_colour_gradient2(
    low = "#c63927",
    mid = "#fee090",
    high = "#4575b4",
    midpoint =9,
    space = "Lab",
    na.value = "grey50",
    guide = "colourbar",
    aesthetics = "colour")


# Plot Figures 2 and 3 of manuscript
R0_MPD_1 <- ggplot(out.total.Sc1, aes(log(MPD),R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("A) Scenario 1: MPD") +
  xlab("Log(MPD)") + ylab(bquote(R[0])) +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")
R0_IC_1 <- ggplot(out.total.Sc1, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("C) Scenario 1: Imbalance") + ylab(bquote(R[0])) +
  xlab("Ic") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette) +
  theme(legend.position = "none")

# For figures 2 and 3, run MH_PhylogeSI_Scenario2.R to generate the Scenario 2 panels, 
# together with R0_IC_1 and R0_MPD_1 generated in lines 219-228.

# Appendix figure 2: Boxplot PD, MPD, b and R0. Substitution and Addition yield same results
par(mfrow=c(1,3))
boxplot(PD~b, data=out.total.Sc1,xlab="Species evolutionary birth rate, r", ylab="PD (MY)",
        outline=FALSE, main = "A) PD distribution at various speciation rates")
boxplot(MPD~b, data=out.total.Sc1,xlab="Species evolutionary birth rate, r", ylab="MPD (MY)",
        outline=FALSE, main = "B) MPD distribution at various speciation rates")
boxplot(R0~b, data=out.total.Sc1,xlab="Species evolutionary birth rate, r", ylab="R0",
        outline=FALSE, main = "C) R0 distribution at various speciation rates")


#############################################
# 6: Regressions
#############################################

# VarPPd, MPD VERY right-skewed. b and S are uniform, Ic double peek, R0 slightly right-skewed (don't log!)

# Univariate models
z<-lm(R0 ~ S,data = out.total.Sc1)
summary(z)
AIC(z)

# Multivariate model
ModelA_Sc1 <- lm(R0 ~ scale(logMPD)+scale(Ic)+scale(S), data = out.total.Sc1, na.action = "na.fail")
vif(ModelA_Sc1)
summary(ModelA_Sc1)
AIC(ModelA_Sc1)
