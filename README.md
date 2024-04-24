# GW_Bayesian
Galton Watson Process following a Bayesian framework.

The purpose of this repository is to recreate the estimations of the transition
probabilities in the article by Cloez, Daufresne, Kerioui and Fontez (2023).

The main difference is the method of MCMC chosen; the authors decided to use 
Gibbs sampler, in this work, I'll be using HMC through Stan. The model will be
compiled in R thanks to CmdStanR.
