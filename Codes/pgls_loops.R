####################### SCRIPT COMPLETO DA ANALISE POR POPULACOES ########
####### 14/12/18 - TALES M. A. PAIVA --

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

## TIRANDO OS RESÍDUOS DOS HERBIVOROS PELO NUMERO DE OCORRENCIAS
m <- read_excel("data/dataset_complete_residuals.xlsx")
dat$herbivores <- m$residuals
head(dat)

## LOOP PARA CRIAR 1000 DATA FRAMES COM VALORES ALEATORIOS POR ESPECIE
lp <- lapply(1:10000, function(x) dat %>%
               group_by(Binomial, life_span, growth_form,
                        sla, max_height, herbivores) %>%
               sample_n(1))

# NOMEANDO AS LINHAS COM OS NOMES DAS ESPECIES
obj <- lapply(lp, function(x){ row.names(x)<-as.character(x$Binomial); x})

# LENDO A FILOGENIA
fil <- read.tree("phylogeny/ALLOTB.tre")

# VENDO SE TODAS AS ESPECIES ESTAO NA ARVORE
#k <- name.check(phy = fil, data = obj[[1]])
#k$data_not_tree

# CRIANDO O OBJETO COMPARATIVO
oc <- match_dataphy(formula = tm ~ herbivores,
                    data = as.data.frame(obj[[1]]), phy = fil)


######################## IR PARA O SCRIPT "pgls_models" PARA CONSTRUIR OS MODELOS #####
###################################################################################




## LOOP PARA CRIAR OS MODELOS
#mods <- lapply(obj, function(x)phylolm(formula = tm ~ log10(herbivores), 
#                                       data = as.data.frame(x),
#                                       phy = oc$phy, model = "lambda"))
## EXTRAINDO P-VALUES
#for (i in 1:length(mods)) {
#  s <- coef(summary(mods[[i]]))["log10(herbivores)", 4]
#}
#for (i in 1:length(mods)) {
#  s[i] <- coef(summary(mods[[i]]))["log10(herbivores)", 4]
#
#df_p<-as.data.frame(s)
#}

#write.csv(df_p, "C:/Users/Tales/Dropbox/Mestrado/Dados/By_populations/p_values.csv")

# PLOTANDO OS P-VALUES
#plot(df_p$s, ylab = "p-values")
#abline(a = 0.05, b = 0, lwd = 2, col = "blue")
#abline(a = mean(df_p$s), b = 0,lty = 2, col = "red")
#mean(df_p$s)

#save(obj, file = "lists/dataset_lists/lista_dataframes_residuals_herb_occ_imputed.Rdata")
