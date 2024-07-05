# Multi-host Phylogenetic SI (Susceptible-Infected) compartmental model
This folder contains the code for the manuscript entitled "Host community structure can shape pathogen outbreak dynamics through a phylogenetic dilution effect".
Codes generate a multi-host SI model, generating a host phylogenetic tree from which the distances of host species are calculated, 
to determine the probability of transmission between species.
These code can be used to generate Figures 1, 2, 3 and 4 of the manuscript. Below I briefly describe each different model (corresponding to the 3 files in this folder).

## Scenario 1
Scenario 1 models assume that the probability of transmission is determined by the distance between the donating and receiving host, 
under the assumption that more closely related species have a higher probability of sharing a disease due to conserved immunological traits.

## Scenario 2
Scenario 2 models assume a reservoir host, and therefore the probability of transmission is determined by the phylogenetic distance between 
the receiving host and the reseroir to which the pathogen has been adapted. Code generates phylogenetic trees for the host community and randomly assigns a reservoir host.

## Gilbert et al. (2012) application for scenario 1 model
Here we apply the plant and pests data from the Gilbert et al. (2012) paper entitled "Evolutionary tools for phytosanitary risk analysis: phylogenetic signal as a predictor of host range of plant pests and pathogens".
Here we assume a scenario 1 model, as there is no information on host reservoirs for these pests.
This code generates Figure 4 of the manuscript.
