# SOd18O
calculates shifts in the Southern Ocean SST front over deglaciation (DLat_SST) using planktic foraminiferal d18O 
it accompanies Gray et al (submitted) and is an updated version of the method used by Gray et al 2020 (https://doi.org/10.1029/2019GL086328)

SO_d18O_data.csv is a compilation of planktic foram d18O data from the southern ocean spanning the last deglaciation (10-20 ka)

lambeck2014.csv is the sealevel curve of Lambeck et al 2014 (https://doi.org/10.1073/pnas.1411762111)
shakun2012.csv  is the global temperature curve of Shakun et al 2012 (doi:10.1038/nature10915)
pmip_dt.csv is the global mean change in SST betwen LGM and PI forcings in the PMIP3/PMIP4 ensemble (https://esgf-node.llnl.gov/search/cmip5/) calculated as part of this study

SO_d18O_DLat.R calculates DLat_SST
- interpolates each d18O record to uniform timesteps using a generalised additive model, correcting the d18O data for whole ocean effects (sea level, global mean SST change)
- models the d18O data as function of latitude using a generalised additive model at each time step
- calculates the latitudinal shift that minimises the Euclidean distance between the modelled d18O data at each time step, relative to a reference timestep
