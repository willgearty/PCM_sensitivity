# PCM_sensitivity

Plan:
- Simulate trees with extant and extinct tips
- Simulate trait evolution on those trees with different models
- Subsample trees using different regimes: extant vs. extinct
- Use PCMs to test whether trait evolution models are correctly identified
- Use post-hoc tests to identify whether using tips with higher phylogenetic or functional uniqueness affect model fit

To do:
- Simulate trees (W) -> two sets of births/deaths to test for tree shape (births and deaths both 1, and births 1 and deaths 0.25) and vary proportions of extant/extinct tips
- Model traits using mvMORPH (B) -> use parameter values from Slater et al. (2012), change the lowest value for Brownian motion from 0 to 0.1, try different values for theta for OU (the optimum value)
- Put together code for model fitting (P)


## Papers
Slater et al. 2012 (Integrating fossils with molecular phylogenies improves inference of trait evolution)
https://doi.org/10.1111/j.1558-5646.2012.01723.x

Ho and Ane 2014 (Intrinsic inference difficulties for trait evolution with Ornstein-Uhlenbeck models)
https://besjournals.onlinelibrary.wiley.com/doi/10.1111/2041-210X.12285

Cooper et al. 2016 (A cautionary note on the use of Ornstein Uhlenbeck models in macroevolutionary studies)
https://academic.oup.com/biolinnean/article/118/1/64/2440254

Mongiardino Koch & Parry 2020 (Death is on Our Side: Paleontological Data Drastically Modify Phylogenetic Hypotheses)
https://doi.org/10.1093/sysbio/syaa023

Duchen et al. 2021 (On the Effect of Asymmetrical Trait Inheritance on Models of Trait Evolution)
https://doi.org/10.1093/sysbio/syaa055

Bartoszek et al. 2022 (Model Selection Performance in Phylogenetic Comparative Methods Under Multivariate Ornstein–Uhlenbeck Models of Trait Evolution)
https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syac079/6962281

Beaulieu and O'Meara 2022 (Fossils Do Not Substantially Improve, and May Even Harm, Estimates of Diversification Rate Heterogeneity)
https://doi.org/10.1093/sysbio/syac049

Wilson et al. 2022 (Chronogram or phylogram for ancestral state estimation? Model-fit statistics indicate the branch lengths underlying a binary character's evolution)
https://doi.org/10.1111/2041-210X.13872

Liow et al. 2023 (Cross-disciplinary information for understanding macroevolution)
https://doi.org/10.1016/j.tree.2022.10.013

Mynard et al. 2023 (Impact of Phylogenetic Tree Completeness and Misspecification of Sampling Fractions on Trait Dependent Diversification Models)
https://academic.oup.com/sysbio/advance-article/doi/10.1093/sysbio/syad001/6988090
