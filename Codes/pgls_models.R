################## CONSTRUINDO OS MODELOS FILOGENETICOS ##########
#### 19/12/2018 - TALES M. A. PAIVA ######

#PACOTES
library(readxl)
library(dplyr)
library(ape)
library(geiger)
library(sensiPhy)
library(caper)

# CARREGANDO DADOS
data <- read_excel("data/dataset_complete_imputed.xlsx")
dat <- dplyr::select(data, Binomial, life_span, growth_form,
                     sla, max_height, herbivores, tm)

## IMPORTANDO A LISTA DE DATAFRAMES
# HERBVIVOROS PUROS
#load("lists/dataset_lists/lista_dataframes.RData")
# HERBIVOROS = RESIDUALS(HERBIVOROS ~ OCORRENCIAS)
#load("lists/dataset_lists/lista_dataframes_residuals_herb_occ.RData")
# HERBIVOROS = RESIDUALS(HERBIVOROS ~ OCORRENCIAS) - SLA E MAX HEIGHT IMPUTED
load("lists/dataset_lists/lista_dataframes_residuals_herb_occ_imputed.RData")

# NOMEANDO AS LINHAS COM OS NOMES DAS ESPECIES
#obj <- lapply(lp, function(x){ row.names(x)<-as.character(x$Binomial); x})

# LENDO A FILOGENIA
fil <- read.tree("phylogeny/ALLOTB.tre")

## OBJETOS COMPARATIVOS
oc <- match_dataphy(formula = tm ~ herbivores, data = as.data.frame(obj[[1]]), phy = fil)
#oc.sla <- match_dataphy(formula = tm ~ log10(sla), data = obj[[1]], phy = fil)
#oc.mh <- match_dataphy(formula = tm ~ log10(max_height), data = obj[[1]], phy = fil)
#oc.bi <- match_dataphy(formula = tm ~ log10(sla) + log10(max_height), 
#                      data = obj[[1]], phy = fil)

############# BASIC MODEL ######
pgls.bas <- lapply(obj, function(x)phylolm(formula = tm ~ herbivores, 
                                       data = as.data.frame(x),
                                       phy = oc$phy, model = "OUfixedRoot"))
# EXTRACT P-VALUES
for (i in 1:length(pgls.bas)) {
  s <- coef(summary(pgls.bas[[i]]))["herbivores", 4]
}
for (i in 1:length(pgls.bas)) {
  s[i] <- coef(summary(pgls.bas[[i]]))["herbivores", 4]
  
  p_bas<-as.data.frame(s)
}

boxplot(p_bas$s)
abline(a = 0.05, b = 0, col = "blue")
mean(p_bas$s)

############# FULL MODEL ######
pgls.full <- lapply(obj, function(x) phylolm(formula = tm ~ herbivores + 
                                               life_span + growth_form +
                                             log10(sla) + log10(max_height),
                                             data = as.data.frame(x), 
                                             phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["herbivores", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["herbivores", 4]
  
  p_full_bas<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["life_spanBiennial", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["life_spanBiennial", 4]
  
  p_full_bi<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["life_spanPerennial", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["life_spanPerennial", 4]
  
  p_full_pe<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["life_spanVaries", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["life_spanVaries", 4]
  
  p_full_va<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formHemiparasite", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formHemiparasite", 4]
  
  p_full_hem<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formHerb", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formHerb", 4]
  
  p_full_her<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formLiana", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formLiana", 4]
  
  p_full_lia<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formShrub", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formShrub", 4]
  
  p_full_shr<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formTree", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formTree", 4]
  
  p_full_tre<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["growth_formVaries", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["growth_formVaries", 4]
  
  p_full_var<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["log10(sla)", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["log10(sla)", 4]
  
  p_full_sla<-as.data.frame(s)
}
for (i in 1:length(pgls.full)) {
  s <- coef(summary(pgls.full[[i]]))["log10(max_height)", 4]
}
for (i in 1:length(pgls.full)) {
  s[i] <- coef(summary(pgls.full[[i]]))["log10(max_height)", 4]
  
  p_full_max<-as.data.frame(s)
}

############# WHOLE PLANT TRAITS ######

#### LIFE SPAN ----------------------------------------------------------------
pgls.wp.ls <- lapply(obj, function(x) phylolm(formula = tm ~ life_span,
                               data = as.data.frame(x), 
                               phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.wp.ls)) {
  s <- coef(summary(pgls.wp.ls[[i]]))["life_spanBiennial", 4]
}
for (i in 1:length(pgls.wp.ls)) {
  s[i] <- coef(summary(pgls.wp.ls[[i]]))["life_spanBiennial", 4]
  
  p_ls_bi<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.ls)) {
  s <- coef(summary(pgls.wp.ls[[i]]))["life_spanPerennial", 4]
}
for (i in 1:length(pgls.wp.ls)) {
  s[i] <- coef(summary(pgls.wp.ls[[i]]))["life_spanPerennial", 4]
  
  p_ls_pe<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.ls)) {
  s <- coef(summary(pgls.wp.ls[[i]]))["life_spanVaries", 4]
}
for (i in 1:length(pgls.wp.ls)) {
  s[i] <- coef(summary(pgls.wp.ls[[i]]))["life_spanVaries", 4]
  
  p_ls_va<-as.data.frame(s)
}

boxplot(p_ls_bi$s)
boxplot(p_ls_pe$s)
boxplot(p_ls_va$s)
###### GROWTH FORM -------------------------------------------------------------
pgls.wp.gf <- lapply(obj, function(x) phylolm(formula = tm ~ growth_form,
                                           data = as.data.frame(x), 
                                           phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formHemiparasite", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formHemiparasite", 4]
  
  p_gf_hem<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formHerb", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formHerb", 4]
  
  p_gf_her<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formLiana", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formLiana", 4]
  
  p_gf_lia<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formShrub", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formShrub", 4]
  
  p_gf_shr<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formTree", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formTree", 4]
  
  p_gf_tre<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.gf)) {
  s <- coef(summary(pgls.wp.gf[[i]]))["growth_formVaries", 4]
}
for (i in 1:length(pgls.wp.gf)) {
  s[i] <- coef(summary(pgls.wp.gf[[i]]))["growth_formVaries", 4]
  
  p_gf_var<-as.data.frame(s)
}

## BIVARIATE
pgls.wp.bi <- lapply(obj, function(x) phylolm(formula = tm ~ life_span + growth_form,
                                           data = as.data.frame(x), 
                                           phy = oc$phy, model = "OUfixedRoot" ))


############################ VEGETATIVE TRAITS ONLY ##############################


######### SPECIFIC LEAF AREA -------------------------------------------------------
pgls.vt.la <- lapply(obj, function(x) phylolm(formula = tm ~ log10(sla),
                                           data = as.data.frame(x), 
                                           phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.vt.la)) {
  s <- coef(summary(pgls.vt.la[[i]]))["log10(sla)", 4]
}
for (i in 1:length(pgls.vt.la)) {
  s[i] <- coef(summary(pgls.vt.la[[i]]))["log10(sla)", 4]
  
  p_sla<-as.data.frame(s)
}

########## MAX HEIGHT --------------------------------------------------
pgls.vt.mh <- lapply(obj, function(x) phylolm(formula = tm ~ log10(max_height),
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.vt.mh)) {
  s <- coef(summary(pgls.vt.mh[[i]]))["log10(max_height)", 4]
}
for (i in 1:length(pgls.vt.mh)) {
  s[i] <- coef(summary(pgls.vt.mh[[i]]))["log10(max_height)", 4]
  
  p_max<-as.data.frame(s)
}

######### BIVARIATE -------------------------------------
pgls.vt.bi <- lapply(obj, function(x) phylolm(formula = tm ~ log10(sla) + log10(max_height),
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.vt.bi)) {
  s <- coef(summary(pgls.vt.bi[[i]]))["log10(sla)", 4]
}
for (i in 1:length(pgls.vt.bi)) {
  s[i] <- coef(summary(pgls.vt.bi[[i]]))["log10(sla)", 4]
  
  p_bi_sla<-as.data.frame(s)
}
for (i in 1:length(pgls.vt.bi)) {
  s <- coef(summary(pgls.vt.bi[[i]]))["log10(max_height)", 4]
}
for (i in 1:length(pgls.vt.bi)) {
  s[i] <- coef(summary(pgls.vt.bi[[i]]))["log10(max_height)", 4]
  
  p_bi_max<-as.data.frame(s)
}

################################ SUMMARY ##############################################
summary(pgls.bas[[1]])
summary(pgls.full[[1]])

summary(pgls.wp.ls[[1]])
summary(pgls.wp.gf[[1]])
summary(pgls.wp.bi[[1]])

summary(pgls.vt.la[[1]])
summary(pgls.vt.mh[[1]])
summary(pgls.vt.bi[[1]])

############################# MODELOS PUROS + HERBIVORIA #################################

################## WHOLE PLANT TRAITS ####

## LIFE SPAN + HERBIVORES
pgls.wp.ls.h <- lapply(obj, function(x) phylolm(formula = tm ~ life_span + 
                                                  herbivores,
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot"))
## GROWTH FORM + HERBIVORES
pgls.wp.gf.h <- lapply(obj, function(x) phylolm(formula = tm ~ growth_form + 
                                                  herbivores,
                                                data = as.data.frame(x), 
                                                phy = oc$phy, model = "OUfixedRoot"))
## WHOLE PLANT TRAITS + HERBIVORES
pgls.wp.bi.h <- lapply(obj, function(x) phylolm(formula = tm ~ life_span + growth_form + 
                                                  herbivores,
                                                data = as.data.frame(x), 
                                                phy = oc$phy, model = "OUfixedRoot"))
# EXTRACT P-VALUES
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanBiennial", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanBiennial", 4]
  
  p_ls.h_bi<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanPerennial", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanPerennial", 4]
  
  p_ls.h_pe<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanVaries", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["life_spanVaries", 4]
  
  p_ls.h_va<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formHemiparasite", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formHemiparasite", 4]
  
  p_gf.h_hem<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formHerb", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formHerb", 4]
  
  p_gf.h_her<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formLiana", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formLiana", 4]
  
  p_gf.h_lia<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formShrub", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formShrub", 4]
  
  p_gf.h_shr<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formTree", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formTree", 4]
  
  p_gf.h_tre<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formVaries", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["growth_formVaries", 4]
  
  p_gf.h_var<-as.data.frame(s)
}
for (i in 1:length(pgls.wp.bi.h)) {
  s <- coef(summary(pgls.wp.bi.h[[i]]))["herbivores", 4]
}
for (i in 1:length(pgls.wp.bi.h)) {
  s[i] <- coef(summary(pgls.wp.bi.h[[i]]))["herbivores", 4]
  
  p_gf.h_hbv<-as.data.frame(s)
}


################## VEGETATIVE TRAITS ####################################

## SPECIFIC LEAF AREA + HERBIVORES
pgls.vt.la.h <- lapply(obj, function(x) phylolm(formula = tm ~ log10(sla) + herbivores,
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.vt.la.h)) {
  s <- coef(summary(pgls.vt.la.h[[i]]))["log10(sla)", 4]
}
for (i in 1:length(pgls.vt.la.h)) {
  s[i] <- coef(summary(pgls.vt.la.h[[i]]))["log10(sla)", 4]
  
  p_la.h_sla<-as.data.frame(s)
}
for (i in 1:length(pgls.vt.la.h)) {
  s <- coef(summary(pgls.vt.la.h[[i]]))["log10(herbivores)", 4]
}
for (i in 1:length(pgls.vt.la.h)) {
  s[i] <- coef(summary(pgls.vt.la.h[[i]]))["log10(herbivores)", 4]
  
  p_la.h_herb<-as.data.frame(s)
}

## MAX HEIGHT + HERBIVORES
pgls.vt.mh.h <- lapply(obj, function(x) phylolm(formula = tm ~ log10(max_height) +
                                                herbivores,
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot" ))
## VEGETATIVE TRAITS + HERBIVORES
pgls.vt.bi.h <- lapply(obj, function(x) phylolm(formula = tm ~ log10(sla) + 
                                                log10(max_height) +
                                                herbivores,
                                              data = as.data.frame(x), 
                                              phy = oc$phy, model = "OUfixedRoot" ))
# EXTRACT P-VALUES
for (i in 1:length(pgls.vt.bi.h)) {
  s <- coef(summary(pgls.vt.bi.h[[i]]))["log10(sla)", 4]
}
for (i in 1:length(pgls.vt.bi.h)) {
  s[i] <- coef(summary(pgls.vt.bi.h[[i]]))["log10(sla)", 4]
  
  p_la.h_sla<-as.data.frame(s)
}
for (i in 1:length(pgls.vt.bi.h)) {
  s <- coef(summary(pgls.vt.bi.h[[i]]))["herbivores", 4]
}
for (i in 1:length(pgls.vt.bi.h)) {
  s[i] <- coef(summary(pgls.vt.bi.h[[i]]))["herbivores", 4]
  
  p_la.h_herb<-as.data.frame(s)
}
for (i in 1:length(pgls.vt.bi.h)) {
  s <- coef(summary(pgls.vt.bi.h[[i]]))["log10(max_height)", 4]
}
for (i in 1:length(pgls.vt.bi.h)) {
  s[i] <- coef(summary(pgls.vt.bi.h[[i]]))["log10(max_height)", 4]
  
  p_la.h_max<-as.data.frame(s)
}

###############################
pgls.null <- lapply(obj, function(x) phylolm(formula = tm ~ 1,
                                                data = as.data.frame(x), 
                                                phy = oc$phy, model = "OUfixedRoot" ))


####################### SUMMARY ###########################
summary(pgls.wp.ls.h[[1]])
summary(pgls.wp.gf.h[[1]])
summary(pgls.wp.bi.h[[1]])
summary(pgls.vt.la.h[[1]])
summary(pgls.vt.mh.h[[1]])
summary(pgls.vt.bi.h[[1]])




 



