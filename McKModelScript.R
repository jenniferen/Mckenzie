######################################
######All combination of models#######
######################################

#g0 models
#g0 covariates: null, session, effort, Julian day, precipitation, Observers

null.model<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~1, sigma~1),
                     mask = McKMask.comb, binomN = 1)

g0.session<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~Session, sigma~1), mask = McKMask.comb, binomN = 1)


g0.Julian<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~Julian, sigma~1),
                    mask = McKMask.comb, binomN = 1)

g0.effort<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort, sigma~1),
                    mask = McKMask.comb, binomN = 1)

g0.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~precip, sigma~1),
                    mask = McKMask.comb, binomN = 1)

g0.observers<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~observers, sigma~1),
                       mask = McKMask.comb, binomN = 1)

g0.univariate.models<-secrlist(null.model, g0.effort, g0.Julian, g0.precip, g0.session, g0.observers)

AIC(g0.univariate.models)

#Model combinations with two covariates
#g0 covariates: null, session, effort, Julian day, precipitation, Observers

#session combos
g0.session.effort<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session+effort, sigma~1),
                            mask = McKMask.comb, binomN = 1)

g0.session.Julian<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session+effort, sigma~1),
                            mask = McKMask.comb, binomN = 1)

g0.session.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session+precip, sigma~1),
                            mask = McKMask.comb, binomN = 1)

g0.session.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session+observers, sigma~1),
                         mask = McKMask.comb, binomN = 1)

#effort combos

g0.effort.Julian<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort+Julian, sigma~1),
                           mask = McKMask.comb, binomN = 1)

g0.effort.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort+precip, sigma~1),
                           mask = McKMask.comb, binomN = 1)

g0.effort.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort+observers, sigma~1),
                        mask = McKMask.comb, binomN = 1)

#Julian day combos
g0.Julian.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~Julian+precip, sigma~1),
                           mask = McKMask.comb, binomN = 1)

g0.Julian.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~Julian+observers, sigma~1),
                        mask = McKMask.comb, binomN = 1)

#AIC
g0.bivariate.models<-secrlist(null.model, g0.effort, g0.Julian, g0.precip, g0.session, g0.observers, 
                              g0.session.effort, g0.session.Julian, g0.session.precip, g0.session.obs, 
                              g0.effort.Julian, g0.effort.precip, g0.effort.obs, g0.Julian.precip,
                              g0.Julian.obs)

AIC(g0.bivariate.models)

#Model combinations with three covariates
#g0 covariates: null, session, effort, Julian day, precipitation, Observers

#Session

g0.session.effort.Julian<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + Julian, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

g0.session.effort.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + precip, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

g0.session.effort.Julian<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + Julian, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

g0.session.effort.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + observer, sigma~1),
                                mask = McKMask.comb, binomN = 1)

g0.session.Julian.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + Julian + precip, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

g0.session.Julian.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + Julian + observer, sigma~1),
                                mask = McKMask.comb, binomN = 1)

g0.session.precip.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + precip + observer, sigma~1),
                                mask = McKMask.comb, binomN = 1)


#effort

g0.effort.Julian.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort + precip + Julian, sigma~1),
                                  mask = McKMask.comb, binomN = 1)

g0.effort.Julian.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort + precip + Julian, sigma~1),
                               mask = McKMask.comb, binomN = 1)

g0.effort.precip.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort + precip + obs, sigma~1),
                               mask = McKMask.comb, binomN = 1)

#Julian Day

g0.Julian.precip.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~Julian + precip + obs, sigma~1),
                               mask = McKMask.comb, binomN = 1)

#AIC

g0.trivariate.models<-secrlist(null.model, g0.effort, g0.Julian, g0.precip, g0.session, g0.observers, 
                               g0.session.effort, g0.session.Julian, g0.session.precip, g0.session.obs, 
                               g0.effort.Julian, g0.effort.precip, g0.effort.obs, g0.Julian.precip,
                               g0.Julian.obs, g0.session.effort.Julian, g0.session.effort.precip, 
                               g0.session.effort.obs, g0.session.Julian.precip, g0.session.Julian.obs,
                               g0.effort.Julian.precip, g0.effort.Julian.obs, g0.effort.precip.obs,
                               g0.Julian., g0.effort.Julian.obs, g0.Julian.precip.obs)

#Model combinations with four covariates
#g0 covariates: null, session, effort, Julian day, precipitation, Observers

g0.session.effort.Julian.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + Julian + precip, sigma~1),
                                          mask = McKMask.comb, binomN = 1)

g0.session.effort.Julian.obs<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + Julian + observer, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

g0.session.Julian.obs.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + Julian + precip + observer, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

g0.effort.Julian.obs.precip<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~effort + Julian + precip + observer, sigma~1),
                                      mask = McKMask.comb, binomN = 1)

g0.quadvariate.models<-secrlist(null.model, g0.effort, g0.Julian, g0.precip, g0.session, g0.observers, 
                                g0.session.effort, g0.session.Julian, g0.session.precip, g0.session.obs, 
                                g0.effort.Julian, g0.effort.precip, g0.effort.obs, g0.Julian.precip,
                                g0.Julian.obs, g0.session.effort.Julian, g0.session.effort.precip, 
                                g0.session.effort.obs, g0.session.Julian.precip, g0.session.Julian.obs,
                                g0.effort.Julian.precip, g0.effort.Julian.obs, g0.effort.precip.obs,
                                g0.Julian., g0.effort.Julian.obs, g0.Julian.precip.obs, g0.session.effort.Julian.precip,
                                g0.session.effort.Julian.obs, g0.session.Julian.obs.precip, g0.effort.Julian.obs.precip)

AIC(g0.quadvariate.models)

#Global model
g0.global<-secr.fit(McKcapt.comb.telem, model = list(D~1, g0~session + effort + Julian + precip + observer, sigma~1),
                    mask = McKMask.comb, binomN = 1)

g0.models<-secrlist(null.model, g0.effort, g0.Julian, g0.precip, g0.session, g0.observers, 
                    g0.session.effort, g0.session.Julian, g0.session.precip, g0.session.obs, 
                    g0.effort.Julian, g0.effort.precip, g0.effort.obs, g0.Julian.precip,
                    g0.Julian.obs, g0.session.effort.Julian, g0.session.effort.precip, 
                    g0.session.effort.obs, g0.session.Julian.precip, g0.session.Julian.obs,
                    g0.effort.Julian.precip, g0.effort.Julian.obs, g0.effort.precip.obs,
                    g0.Julian., g0.effort.Julian.obs, g0.Julian.precip.obs, g0.session.effort.Julian.precip,
                    g0.session.effort.Julian.obs, g0.session.Julian.obs.precip, g0.effort.Julian.obs.precip,
                    g0.global)


#Density models; Run g0 with top g0 model.
#Density covariates: DDE, DFE, Roads, slope, precip, Crops

#univariate combinations
D.DDE.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE, g0~1, sigma ~1),
                    mask = McKMask.comb, binomN = 1)

D.DFE.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE, g0~1, sigma~1),
                    mask = McKMask.comb, binomN = 1)

D.Roads.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads, g0~1, sigma~1),
                      mask = McKMask.comb, binomN = 1)

D.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~PercentSlope, g0~1, sigma~1),
                      mask = McKMask.comb, binomN = 1)

D.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Precip, g0~1, sigma~1),
                       mask = McKMask.comb, binomN = 1)

D.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Crops, g0~1, sigma~1),
                    mask = McKMask.comb, binomN = 1)

D.univariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.Crops.g0.)

AIC(D.univariate.models)

#bivariate combinations
#Density covariates: DDE, DFE, Roads, slope, precip, Crops

D.DDE.DFE.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE, g0~1, sigma~1),
                        mask = McKMask.comb, binomN = 1)

D.DDE.Roads.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + Roads, g0~1, sigma~1),
                          mask = McKMask.comb, binomN = 1)

D.DDE.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + PercentSlope, g0~1, sigma~1),
                          mask = McKMask.comb, binomN = 1)

D.DDE.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + Precip, g0~1, sigma~1),
                           mask = McKMask.comb, binomN = 1)

D.DDE.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + Crops, g0~1, sigma~1),
                           mask = McKMask.comb, binomN = 1)

D.DFE.Roads.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads, g0~1, sigma~1),
                          mask = McKMask.comb, binomN = 1)

D.DFE.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + PercentSlope, g0~1, sigma~1),
                          mask = McKMask.comb, binomN = 1)

D.DFE.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Precip, g0~1, sigma~1),
                           mask = McKMask.comb, binomN = 1)

D.DFE.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Crops, g0~1, sigma~1),
                           mask = McKMask.comb, binomN = 1)

D.Roads.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + PercentSlope, g0~1, sigma~1),
                            mask = McKMask.comb, binomN = 1)

D.Roads.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + Precip, g0~1, sigma~1),
                             mask = McKMask.comb, binomN = 1)

D.Roads.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + Crops, g0~1, sigma~1),
                             mask = McKMask.comb, binomN = 1)

D.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~PercentSlope + Precip, g0~1, sigma~1),
                             mask = McKMask.comb, binomN = 1)

D.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~PercentSlope + Crops, g0~1, sigma~1),
                             mask = McKMask.comb, binomN = 1)

D.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Precip + Crops, g0~1, sigma~1),
                             mask = McKMask.comb, binomN = 1)

D.bivariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.Crops.g0., D.DDE.DFE.g0., D.DDE.Roads.g0.,
                             D.DDE.slope.g0., D.DDE.precip.g0., D.DDE.Crops.g0., D.DFE.Roads.g0., D.DFE.slope.g0., D.DFE.precip.g0.,
                             D.DFE.Crops.g0., D.Roads.slope.g0., D.Roads.precip.g0., D.slope.precip.g0., D.slope.Crops.g0., D.precip.Crops.g0.)

AIC(D.bivariate.models)

#trivariate models
#Density covariates: DDE, DFE, Roads, slope, precip, Crops

D.DDE.DFE.Roads.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads, g0~1, sigma~1),
                              mask = McKMask.comb, binomN = 1)

D.DDE.DFE.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + PercentSlope, g0~1, sigma~1),
                              mask = McKMask.comb, binomN = 1)

D.DDE.DFE.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Precip, g0~1, sigma~1),
                               mask = McKMask.comb, binomN = 1)

D.DDE.DFE.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Crops, g0~1, sigma~1),
                               mask = McKMask.comb, binomN = 1)

D.DFE.Roads.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + PercentSlope, g0~1, sigma~1),
                                mask = McKMask.comb, binomN = 1)

D.DFE.Roads.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + Precip, g0~1, sigma~1),
                                 mask = McKMask.comb, binomN = 1)

D.DFE.Roads.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + Crops, g0~1, sigma~1),
                                mask = McKMask.comb, binomN = 1)

D.DFE.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + PercentSlope + Precip, g0~1, sigma~1),
                                 mask = McKMask.comb, binomN = 1)

D.DFE.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + PercentSlope + Crops, g0~1, sigma~1),
                                mask = McKMask.comb, binomN = 1)

D.Roads.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + PercentSlope + Precip, g0~1, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

D.Roads.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~Roads + PercentSlope + Crops, g0~1, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

D.slope.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~PercentSlope + Precip + Crops, g0~1, sigma~1),
                                   mask = McKMask.comb, binomN = 1)

D.trivariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.DDE.DFE.g0., D.DDE.Roads.g0.,
                              D.DDE.slope.g0., D.DDE.precip.g0., D.DFE.Roads.g0., D.DFE.slope.g0., D.DFE.precip.g0., D.Roads.slope.g0.,
                              D.Roads.precip.g0., D.slope.precip.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                              D.DDE.DFE.Crops.g0., D.DFE.Roads.Crops.g0., D.DFE.Roads.Crops.g0., 
                              D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.DFE.slope.Crops.g0.,
                              D.Roads.slope.precip.g0., D.Roads.slope.Crops.g0., D.slope.precip.Crops.g0.)

#quadvariate models
#Density covariates: DDE, DFE, Roads, slope, precip, Crops

D.DDE.DFE.Roads.slope.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + PercentSlope, g0~1, sigma~1),
                                    mask = McKMask.comb, binomN = 1)

D.DDE.DFE.Roads.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + Precip, g0~1, sigma~1),
                                     mask = McKMask.comb, binomN = 1)

D.DDE.DFE.Roads.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + Crops, g0~1, sigma~1),
                                     mask = McKMask.comb, binomN = 1)

D.DDE.DFE.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + PercentSlope + Precip, g0~1, sigma~1),
                                     mask = McKMask.comb, binomN = 1)

D.DDE.DFE.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + PercentSlope + Crops, g0~1, sigma~1),
                                     mask = McKMask.comb, binomN = 1)

D.DDE.DFE.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Precip + Crops, g0~1, sigma~1),
                                     mask = McKMask.comb, binomN = 1)

D.DFE.Roads.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + PercentSlope + Precip, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

D.DFE.Roads.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + PercentSlope + Crops, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

D.DFE.Roads.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + Precip + Crops, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

D.DFE.slope.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + PercentSlope + Precip + Crops, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

D.quadvariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.DDE.DFE.g0., D.DDE.Roads.g0.,
                               D.DDE.slope.g0., D.DDE.precip.g0., D.DFE.Roads.g0., D.DFE.slope.g0., D.DFE.precip.g0., D.Roads.slope.g0.,
                               D.Roads.precip.g0., D.slope.precip.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                               D.DDE.DFE.Crops.g0., D.DFE.Roads.Crops.g0., D.DFE.Roads.Crops.g0., 
                               D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.DFE.slope.Crops.g0.,
                               D.Roads.slope.precip.g0., D.Roads.slope.Crops.g0., D.slope.precip.Crops.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                               D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.Roads.slope.precip.g0.,
                               D.DDE.DFE.Roads.slope.g0., D.DDE.DFE.Roads.precip.g0., D.DDE.DFE.Roads.Crops.g0., D.DDE.DFE.slope.precip.g0.,
                               D.DDE.DFE.slope.Crops.g0., D.DDE.DFE.precip.Crops.g0., D.DFE.Roads.slope.precip.g0., D.DFE.Roads.slope.Crops.g0.,
                               D.DFE.Roads.precip.Crops.g0., D.DFE.slope.precip.Crops.g0.)

#quintvariate modles
#Density covariates: DDE, DFE, Roads, slope, precip, Crops


D.DDE.DFE.Roads.slope.precip.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + PercentSlope + Precip, g0~1, sigma~1),
                                    mask = McKMask.comb, binomN = 1)

D.DDE.DFE.Roads.slope.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + PercentSlope + Crops, g0~1, sigma~1),
                                    mask = McKMask.comb, binomN = 1)

D.DDE.Roads.slope.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + Roads + PercentSlope + Precip + Crops, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)

D.DDE.DFE.slope.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + PercentSlope + Precip + Crops, g0~1, sigma~1),
                                             mask = McKMask.comb, binomN = 1)

D.DDE.DFE.Roads.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + Roads + Precip + Crops, g0~1, sigma~1),
                                             mask = McKMask.comb, binomN = 1)

D.DFE.Roads.slope.precip.Crops.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DFE + Roads + PercentSlope + Precip + Crops, g0~1, sigma~1),
                                       mask = McKMask.comb, binomN = 1)


D.quadvariate.models<-secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.DDE.DFE.g0., D.DDE.Roads.g0.,
                               D.DDE.slope.g0., D.DDE.precip.g0., D.DFE.Roads.g0., D.DFE.slope.g0., D.DFE.precip.g0., D.Roads.slope.g0.,
                               D.Roads.precip.g0., D.slope.precip.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                               D.DDE.DFE.Crops.g0., D.DFE.Roads.Crops.g0., D.DFE.Roads.Crops.g0., 
                               D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.DFE.slope.Crops.g0.,
                               D.Roads.slope.precip.g0., D.Roads.slope.Crops.g0., D.slope.precip.Crops.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                               D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.Roads.slope.precip.g0.,
                               D.DDE.DFE.Roads.slope.g0., D.DDE.DFE.Roads.precip.g0., D.DDE.DFE.Roads.Crops.g0., D.DDE.DFE.slope.precip.g0.,
                               D.DDE.DFE.slope.Crops.g0., D.DDE.DFE.precip.Crops.g0., D.DFE.Roads.slope.precip.g0., D.DFE.Roads.slope.Crops.g0.,
                               D.DFE.Roads.precip.Crops.g0., D.DFE.slope.precip.Crops.g0.)

#global model
D.globalmodel.g0.<-secr.fit(McKcapt.comb.telem, model = list(D~DDE + DFE + Roads + PercentSlope + Precip + Crops, g0~1, sigma~1),
                            mask = McKMask.comb, binomN = 1)

D.models<--secrlist(D.DDE.g0., D.DFE.g0., D.Roads.g0., D.slope.g0., D.precip.g0., D.DDE.DFE.g0., D.DDE.Roads.g0.,
                    D.DDE.slope.g0., D.DDE.precip.g0., D.DFE.Roads.g0., D.DFE.slope.g0., D.DFE.precip.g0., D.Roads.slope.g0.,
                    D.Roads.precip.g0., D.slope.precip.g0., D.DDE.DFE.Roads.g0., D.DDE.DFE.slope.g0., D.DDE.DFE.precip.g0.,
                    D.DFE.Roads.slope.g0., D.DFE.Roads.precip.g0., D.DFE.slope.precip.g0., D.Roads.slope.precip.g0.,
                    D.DDE.DFE.Roads.slope.g0., D.DDE.DFE.Roads.precip.g0., D.DFE.Roads.slope.precip.g0., D.globalmodel.g0.)

AIC(D.models)

top.model<-
  
save.image(file = "20210104McKRuns.RData")
load("20210104McKRuns.RData")


model$fit
coef(model)