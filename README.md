# MiMICS

**M**ulti-scale model of **M**etabolism **I**n **C**ellular **S**ystems

Multi-scale Model of Metabolism In Cellular Systems, abbreviated MiMICS, is a computational framework executed in Python and Java to mechanistically simulate heterogeneous metabolism in microbial communities. MiMICS couples a genome-scale metabolic network reconstruction (GENRE) with the established platform Hybrid Automata Library (HAL) [1], which contains an agent-based model, and a continuum-scale reaction-diffusion model. GENREs are optimized using the COBRApy Python package. The metabolic models (Python) and HAL models (Java) are interfaced using a Py4J Gateway Server. MiMICS simulations can be performed using multiple computing nodes and central processing units to reduce computational time.  

MiMICS can be extended to multiple applications to a) uncover multi-scale mechanisms controlling heterogeneous metabolic states in a microbial biofilm, b) predict the fate of heterogeneous biofilms in various metabolic conditions, or c) test potential biofilm treatment strategies. 

## MiMICS overview

MiMICS simulations predict intracellular and extracellular metabolic processes in both space and time in multicellular communities. At each simulation time step, each MiMICS sub-model is performed to update cellular agent properties and metabolite concentrations. Each cellular agent simulates their own metabolic model, which is constrained on the agent's local metabolite concentrations. Each agent's metabolic model predicts a biomass growth rate, and metabolite secretion and uptake fluxes. The biomass growth rate is used to update an agent’s biomass in the agent-based model. Metabolite secretion and uptake fluxes are passed to the reaction-diffusion model to update the agent's local metabolite concentrations. In the agent-based model, microbial agents can perform behaviors like biomass division, motility, and cell mechanics. Metabolite diffusion is performed using HAL's partial differential equation solver methods. Simulation outputs such as agent locations, agent intracellular metabolic reaction fluxes, and metabolite concentrations are saved at each time step. 


## MiMICS extensibility

Inputs into MiMICS allow the user to easily extend MiMICS to simulate multi-scale metabolism in microbial communities and metabolite conditions of interest. MiMICS can simulate single-species or multi-species microbial communities in 2D or 3D. Individual agents can be represented as a single cell or a population of cells. Metabolite inputs, such as a growth media recipe, define the extracellular metabolite concentrations accessible to an agent to use to grow biomass or respond to. 

A key feature of MiMICS is the ability to incorporate multiple -omics data-integrated GENREs, which can represent unique intracellular metabolic states that yield different parameter values passed to the extracellular models. Incorporation of -omics data-integrated GENREs into MiMICS can result in a heterogeneous distribution of agents simulating different intracellular and extracellular metabolic processes, which may impact the fate of the biofilm. With the generation of spatially-resolved -omics data sets of heterogeneous multicellular communities, MiMICS is advantageous to incorporate this data to a) improve prediction accuracy of multicellular metabolic processes and b) gain new insight of mechanisms controlling metabolic state heterogeneity. 

## References
1. Bravo, R. R. et al. Hybrid Automata Library: A flexible platform for hybrid modeling with real-time visualization. PLoS Comput Biol 16, e1007635 (2020).
