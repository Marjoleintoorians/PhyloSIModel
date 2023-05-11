# Phylogenetic community structure shapes outbreak potential
# Marjolein Toorians
# May 2023

library(outliers)
library(phytools)
library(ape)
library(picante)
library(fields)
library(apTreeshape)
library(stats)
library(ggplot2)
library(tidyverse)
library(cowplot)
library(ggfortify)
library(ggstatsplot)
library(dplyr)
library(car)
library(MASS)
library(MuMIn)
library(binom)
library(bbmle)
library(emmeans)
library(visreg)
library(pwr)
palette.Ch2 <- c("#d73027","#fc8d59","#fef590","#91bfdb","#4575b4")
##########################################################################

### 1: Parameters
nspecvec = c(3,4,5,7,10)
nspecvec<-nspecvec[order(-nspecvec)] # Descending vector, for pca plotting purposes.
ngen = 100   # number of trees.
delta <- 0.1 # disease-induced mortality
mu <- 0.2 # natural mortality rate
lambda <- 1 # c_j: the interspecific probability of succesfull transmission j
bvec <- seq(0.1,1,0.1)

set.seed(1)

# dataframe initialization
length_df <- ngen*length(nspecvec)*length(bvec)
out.total_raw <- matrix(0,nrow = length_df, ncol=8)
colnames(out.total_raw) <- c("PD","MPD","Ic","R0","b","S","PPd12","varPPd")

### 2: Simulate communities
for (k in 1:length(nspecvec)){
  
  nspec <- nspecvec[k]
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
  kappa <- matrix(1, nrow=nspec,ncol=nspec)
  
  # Create p-matrix
  pmatrix <- matrix(0,nrow = nspec, ncol=nspec)
  for (i in 1:nspec){
    for (j in 1:nspec){
      
      kappatot <- sum(kappa[,j])
      pmatrix[[i,j]] <- kappa[i,j]/kappatot
      
    }
  }
  
  # 3: Loop over birth rates
  for (l in 1:length(bvec)){
    
    # set b
    bvalue <- bvec[l]
    
    # Create 'ngen' number of trees
    for (x in 1:ngen){
      tree_temp <- pbtree(b=bvalue, d=0, n=nspec) 
      
      # save data
      ppd_temp<-cophenetic(tree_temp)
      
      # Create c-matrix for each tree
      for (i in 1:nspec){
        for (j in 1:nspec){
          
          cmatrix[[i,j]] <- lambda/(1+ppd_temp[i,j])
        }
      }
      
      # Create R-matrix fo each tree
      for (i in 1:nspec){
        for (j in 1:nspec){
          
          Rmatrix[[i,j]] <- (cmatrix[i,j]*pmatrix[i,j])/(delta+mu)
          
        }
      }
      
      ### Output for all 'ngen' trees (x-loop)
      # PPd species 1 and 2
      PPd12 <- ppd_temp["t1","t2"]
      # MPD & varPPd step
      ppd_temp_mpd <-ppd_temp[lower.tri(ppd_temp,diag=FALSE)]
      # Ic: Colless's I-statistic
      i.tree<-as.treeshape(tree_temp, model="yule") # this converts from phylo to treeshape
      # Total output dataset
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),1] = sum(tree_temp$edge.length) #PD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),2] = mean(ppd_temp_mpd) #MPD
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),3] = colless(i.tree, norm = "yule") #Ic
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),4] = max(eigen(Rmatrix)$values) #R0
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),5] = bvalue # b-value
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),6] = nspec # number of species
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),7] = PPd12
      out.total_raw[x+(ngen*(l-1))+(length(bvec)*ngen*(k-1)),8] = mean(var(ppd_temp_mpd))  #varppd

    } # end of x-loop, for generation ngen trees
    
  } # end b loop (l-loop)
  
} # end nspec (k-loop)

### 3: Markup dataset
out.total.Sc1<-as.data.frame(out.total_raw)
# Make some numbers characters for plotting
out.total.Sc1$bNumber <- as.character(out.total.Sc1$b)
out.total.Sc1$Species <- as.character(out.total.Sc1$S)
# Add log vars
out.total.Sc1$logMPD <- log(out.total.Sc1$MPD)

# Test for normality
shapiro.test(log(out.total.Sc1$PD))
# HIGHEST: Log(MPD), log(PD), log(varPPd), Ic(outliers removed), 
grubbs.test(out.total.Sc1$Ic, type = 10, opposite = FALSE, two.sided = FALSE)
#boxplot(out.total.Sc1$Ic)


### 4: PCA

# PCA 1 (only columns 1:5)
pca1 <- prcomp(out.total.Sc1[,(1:4)],scale=TRUE)
# Plot (in manuscript)
pca.sc1<-autoplot(pca1, data = out.total.Sc1, label = F,loadings=T, loadings.label=T, 
         loadings.colour = 'black',size = 1, loadings.label.colour = 'black', loadings.label.size = 6,
         colour = 'Species', show.legend = FALSE, main = "Scenario 1") + theme_classic() + 
         scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette.Ch2) +
         theme(legend.position=c(0.08,0.8),legend.background = element_rect(linetype = 2, size = 0.5, colour = 1))
         #scale_color_brewer(palette = "RdYlBu") +theme(legend.position = "none") #+scale_y_reverse()

# New dataset including the PC1, PC2 and PC3
out.total.Sc1 <- data.frame(cbind(out.total.Sc1, pca1$x[,1:3]))

# Correlations between vars and PCs
loading_scores <- pca1$rotation[,1] # first only look at PCA 1
loading_scores_abs <- abs(loading_scores) # absolute values
loading_scores_ranked <- sort(loading_scores_abs,decreasing = TRUE)
top_pars <- names(loading_scores_ranked)
top_pars

### 5: Plots
# Boxplot b and R0
#PD
par(mfrow=c(2,2))
boxplot(PD~b, data=out.total.Sc1,xlab="Species birth rate, b", ylab="PD (MY)", 
        outline=FALSE, main = "PD of 3-species communities at various speciation rates")
boxplot(R0~b, data=out.total.Sc1,xlab="Species birth rate, b", ylab="R0", outline=FALSE) 
        #,main = "R0 distribution of 100 simulated trees with 3 species")
#MPD
#par(mfrow=c(1,2))
boxplot(MPD~b, data=out.total.Sc1,xlab="Species birth rate, b", ylab="MPD (MY)", 
        outline=FALSE, main = "MPD of 3-species communities at various speciation rates")
boxplot(R0~b, data=out.total.Sc1,xlab="Species birth rate, b", ylab="R0", outline=FALSE) 
        #,main = "R0 distribution of 100 simulated trees with 3 species")




# R0 vs metrics plot
R0_MPD_1 <- ggplot(out.total.Sc1, aes(log(MPD),R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("A: Mean Pairwise Distance") + 
  xlab("Log(MPD)") + ylab(bquote(R[0])) +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette.Ch2) +
  theme(legend.position = "none")
R0_IC_1 <- ggplot(out.total.Sc1, aes(Ic,R0)) + geom_point(aes(colour=Species),size=0.9) + theme_classic() +
  ggtitle("B: Imbalance") + ylab(bquote(R[0])) +
  xlab("Ic") +
  scale_colour_manual(limits =c("3", "4", "5","7","10"), values = palette.Ch2)
# Together
require(gridExtra)
grid.arrange(R0_MPD_1, R0_IC_1, ncol=2)


# 6 Regressions: 
# Violin plot data
ggplot(out.total.Sc1, aes(x=Ic, y=PD)) + 
  geom_violin(trim=FALSE)

# VarPPd, MPD VERY right-skewed. b and S are uniform, Ic double peek, R0 slightly right-skewed (don't log!)
# Models show in Q-Q plot that data is skewed as well! It seems to show "Thin tails", as there are far outliers

# 6.1 Univariate models
z<-lm(R0 ~ Ic,data = out.total.Sc1)
summary(z)


# 6.2 Multivariate model
ModelA_Sc1 <- lm(R0 ~ scale(logMPD)+scale(Ic)+scale(S), data = out.total.Sc1, na.action = "na.fail")
summary(ModelA_Sc1)
AIC(ModelA_Sc1)

# VIF: Variance inflation factor, to explore issues of colinearity
vif(ModelA_Sc1)
