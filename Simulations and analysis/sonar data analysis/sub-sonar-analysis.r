#rm(list = ls()) # clear environment

# HEADER ====

#library(groupdata2)
#source("../gicf/gicf.R")

#set.seed(789)

#use.banded.structure <- F


if(use.banded.structure){
  bands.rocks <- 17
  bands.metals <- 31
}else{
  bands.rocks <- 60
  bands.metals <- 60
}


# MODEL SELECTION ====


## Model selection (rocks) ====
print("Computing CV scores for ROCKS...")

cv.scores.rocks.lambda <- nested.model.selection.cv(data.rocks,
                                             N.in.folds,
                                             N.lambdas,
                                             1,
                                             k.max = 0,
                                             adj = banded.adj(60, bands.rocks),
                                             s.num = s.num#Alberto
                                             )

cv.scores.rocks.kappa <- nested.model.selection.cv(data.rocks,
                                             N.in.folds,
                                             1,
                                             N.kappas,
                                             l.max = 0,
                                             adj = banded.adj(60, bands.rocks),
                                             s.num = s.num#Alberto
                                             )
cv.scores.rocks <- nested.model.selection.cv(data.rocks,
                                      N.in.folds,
                                      N.lambdas,
                                      N.kappas,
                                      l.max = 15,
                                      k.max = 5,
                                      adj = banded.adj(60, bands.rocks),
                                      s.num = s.num#Alberto
                                      ) 

cv.scores.rocks <- abind::abind(cv.scores.rocks, 
                                cv.scores.rocks.kappa,
                                cv.scores.rocks.lambda,
                                along = 1)

## Model selection (metals) ====
print("Computing CV scores for METALS...")
cv.scores.metals.lambda <- nested.model.selection.cv(data.metals,
                                             N.in.folds,
                                             N.lambdas,
                                             1,
                                             k.max = 0,
                                             adj = banded.adj(60, bands.metals),
                                             s.num = s.num#Alberto
                                             )


cv.scores.metals.kappa <- nested.model.selection.cv(data.metals,
                                            N.in.folds,
                                            1,
                                            N.kappas,
                                            l.max = 0,
                                            adj = banded.adj(60, bands.metals),
                                            s.num = s.num#Alberto
                                            )
cv.scores.metals <- nested.model.selection.cv(data.metals,
                                       N.in.folds,
                                       N.lambdas,
                                       N.kappas,
                                       l.max = 15,
                                       k.max = 5,
                                       adj = banded.adj(60, bands.metals),
                                       s.num = s.num#Alberto
                                       ) 


cv.scores.metals <- abind::abind(cv.scores.metals, 
                                cv.scores.metals.kappa,
                                cv.scores.metals.lambda,
                                along = 1
                                )

# OUTPUT ====
#source("../gicf/gicf.R")

data.folded <- rbind(data.rocks, data.metals)
error.rate.mle <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         F, F)
error.rate.lambda <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         F, T)
error.rate.kappa <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         T, F)
error.rate.gicf <- QDA.cv(data.folded, 
                         cv.scores.rocks, cv.scores.metals,
                         T, T)

print(error.rate.mle)
print(error.rate.lambda)
print(error.rate.kappa)
print(error.rate.gicf)

errors <- c(error.rate.mle, error.rate.lambda, error.rate.kappa, error.rate.gicf)
names(errors) <- c("MLE", "LAMBDA", "KAPPA", "GICF")

#write.table(errors, paste0("sonar_analysis_band_structure_", use.banded.structure, ".csv"))
