#Bethany      28th March 2023
#Simulate continuous trait evolution across a phylogeny using different models

library(ape)
library(phytools)
library(mvMORPH)

tree <- read.tree("data/testtree.txt")
plot(tree)

#Simulate weak BM trait evolution
wBM_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                  param = list(trend = FALSE,
                                 theta = 0, #ancestral state
                                 sigma = 0.1 #strength of drift
                  ))
contMap(tree, wBM_trait[,1])

#Simulate strong BM trait evolution
sBM_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                  param = list(trend = FALSE,
                               theta = 0, #ancestral state
                               sigma = 0.5 #strength of drift
                  ))
contMap(tree, sBM_trait[,1])

#Simulate weak trended BM trait evolution
wtrend_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                    param = list(trend = 0.1, #strength of trend
                               theta = 0, #ancestral state
                               sigma = 0.1 #strength of drift
                    ))
contMap(tree, wtrend_trait[,1])

#Simulate strong trended BM trait evolution
strend_trait <- mvSIM(tree = tree, nsim = 1, model = "BM1",
                     param = list(trend = 0.3, #strength of trend
                                  theta = 0, #ancestral state
                                  sigma = 0.1 #strength of drift
                     ))
contMap(tree, strend_trait[,1])

#Simulate weak OU ("SSP") trait evolution, centred on ancestral state
wOUc_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                   param = list(alpha = 0.01, #strength of selection
                                 theta = 0, #ancestral state
                                 sigma = 0.1 #strength of drift
                  ))
contMap(tree, wOUc_trait[,1])

#Simulate strong OU ("SSP") trait evolution, centred on ancestral state
sOUc_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                   param = list(alpha = 10, #strength of selection
                                theta = 0, #ancestral state
                                sigma = 0.1 #strength of drift
                   ))
contMap(tree, sOUc_trait[,1])

#Simulate weak OU ("SSP") trait evolution, with shifted optimum
wOUs_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                    param = list(root = TRUE,
                                 alpha = 0.01, #strength of selection
                                 theta = c(0, 1), #ancestral state,
                                                  #optimum
                                 sigma = 0.1 #strength of drift
                    ))
contMap(tree, wOUs_trait[,1])

#Simulate strong OU ("SSP") trait evolution, with shifted optimum
sOUs_trait <- mvSIM(tree = tree, nsim = 1, model = "OU1",
                    param = list(root = TRUE,
                                 alpha = 10, #strength of selection
                                 theta = c(0, 1), #ancestral state,
                                                  #optimum
                                 sigma = 0.1 #strength of drift
                    ))
contMap(tree, sOUs_trait[,1])

#Simulate weak AC trait evolution
wAC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                    param = list(theta = 0, #ancestral state
                                 beta = 0.1, #exponential rate
                                 sigma = 0.001 #strength of drift
                    ))
contMap(tree, wAC_trait[,1])

#Simulate strong AC trait evolution
sAC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                   param = list(theta = 0, #ancestral state
                                beta = 0.3, #exponential rate
                                sigma = 0.001 #strength of drift
                  ))
contMap(tree, sAC_trait[,1])

#Simulate weak DC ("Early Burst") trait evolution
wDC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                   param = list(theta = 0, #ancestral state
                                beta = -0.1, #exponential rate
                                sigma = 0.001 #strength of drift
                  ))
contMap(tree, wDC_trait[,1])

#Simulate strong DC ("Early Burst") trait evolution
sDC_trait <- mvSIM(tree = tree, nsim = 1, model = "EB",
                   param = list(theta = 1, #ancestral state
                                beta = -0.3, #exponential
                                sigma = 0.001 #strength of drift
                  ))
contMap(tree, sDC_trait[,1])
