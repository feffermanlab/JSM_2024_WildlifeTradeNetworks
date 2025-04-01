# JSM_2024_WildlifeTradeNetworks
Collection of code to model the bidirectional coupling of pathogenic dynamics and economic factors in wildlife trade networks.
This bidirectional coupling and the associated model results are written about in the manuscript associated with this repository.
The goal of the project is to capture the economic behavior of a trade network, the contagion dynamics on a network with acyclic flow, and the risk of spillover all in a single model. 
We are especially interested in the way that each of these three factors interacts with the others.
The hope is that this model will make it so that the network structure can be given as an independent variable and the spillover risk, total trade volume etc. can be measured as a
dependent varaible. 

FixedPointArgument.m is a script which was helpful in proving the result which ultimately became Theorem 1 in the manuscript. The code is left here for instructive purposes

The next grouping of scripts is used to run the simulation when the stakeholders are grouped into four classes, breeders, distributors, retailers, and consumers. 
WTNsimulateInit.m is a script which sets up all the initial data and parameters required for the following simulations. WTNsimulate.m is the script which takes information
from WTNsimulateInit.m and runs the stochastic simulation to give a payoff and infection probabilities for each stakeholder after a long time. 

The final grouping of scripts was used extensively in the manuscript to generate images and get results from the toy model
AdjMatSelect.m is the funciton which encapselates all of the different adjacency matrices used by the toy model
WTNNumericalEx.m is the code which runs each of the simulated scenarios discussed in section 4 of the manuscript. This script does not rely on the simulation scripts discussed previously
because it assumes simple functional forms for all of the interactions and only separeates stakeholders into three catagories, producers, distributors, and consumers. Much of this script 
is dedicated to the generation of figures and the production of numerical results which are further discussed in the manuscript.  
