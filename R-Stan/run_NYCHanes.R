##-----------------------------------------------------------------------------------------------
## TITLE:        run_NYCHanes.R: data analysis of NYC Hanes data using Bayesian shrinkage priors for Variable Selection in an Integrated 
##               Dirichlet-Multinomial model 
##
## VERSION:      1st version (5/19/2020).
##
## AUTHORS:      
##
## DESCRIPTION:  Method runs HMC chains based on Stan to sample from the posterior
##               distribution for the model and uses HS, HS+ and Laplace for the regression parameters.
##
## DEPENDS ON:   'rstan' + scriots for plotting and other utility: "simple_heatmap.R" and "functions_utility_stat.R")
##-----------------------------------------------------------------------------------------------



################################ 
### NYC Hanes data preprocessing
################################ 

rm(list = ls())

setwd("C:/Users/jyotishka/Box Sync/Jyotishka/Microbiome Data")

setwd("R-Stan")
source("simple_heatmap.R")
source("functions_utility_stat.R")

## You can skip the raw-data-reading part and 
## go straight to analysis with reading the csv files
## read in the data

raw_data_read = FALSE

if(raw_data_read){
  source("NYCHanes/read_write_nyc_hanes_raw_data.r")
}

########################
### Run data with R stan
########################

rm(list = ls())
setwd("C:/Users/jyotishka/Box Sync/Jyotishka/Microbiome Data/R-Stan/")
# setwd("/Users/pixushi/Box/R-Stan")
source("simple_heatmap.R")
source("functions_utility_stat.R")

## read in the processed data

phyloY <- read.csv("NYCHanes/phylogeny_genus_50pc.csv", row.names=1, colClasses="character")
Ymat <- read.csv("NYCHanes/otu_genus_50pc.csv", row.names=1)
Xmat <- read.csv("NYCHanes/metadata_matrix.csv", row.names=1)
table(rownames(Xmat)==rownames(Ymat)) # check row name match

# remove samples with NA and pregnancy
Ymat <- Ymat[Xmat$REMOVE==0,]
Xmat <- Xmat[Xmat$REMOVE==0,]
Xmat$REMOVE <- NULL
dim(Xmat)
dim(Ymat)

# scale Xmat

Xmat0 <- Xmat
Xmat <- scale(Xmat)

# ?gsub

cleannames = gsub("New.Reference", "", colnames(Ymat))
cleannames = gsub("New.CleanUp.Reference","", cleannames)
cleannames

colnames(Ymat) = cleannames

colnames(Ymat) <- phyloY[,6]

dim(cor(Xmat, scale(Ymat)))


# simple_heatmap(cor(Xmat, scale(Ymat)))+ coord_flip() + plotTheme()

corr.mat <- cor(Xmat, scale(Ymat),method = "spearman")
colnames(corr.mat) <- phyloY[,6]

(corr.map = simple_heatmap(corr.mat) + coord_flip()+
  labs(title="Spearman's rank correlation") + 
    xlab("Features")+ylab("OTU")+
    plotTheme()
  )
print(corr.map)

ggsave("NYCHanes/rankcorrMap_NYCHanes_2.pdf", corr.map, width = 9, height = 5)

normalize.composition <- function(allData){
  allProp <- sweep(allData, 1, rowSums(allData), "/")
  allProp
}

Yprop <-  normalize.composition(Ymat)

sel.sp = (colSums(Ymat>0)>=quantile(colSums(Ymat>0),0.95))

# sel.obs = (rowSums(chabanCounts)>=100
# sel.sp = (colSums(YData>0)>=10)&(colSums(YData>0)>=10)
topOTU <- Yprop[,sel.sp] %>% dplyr::select(-"NA")


library(ggplot2)
library(viridisLite)
library(viridis)
library(hrbrthemes)
library(reshape2)

melted_otu <- melt(t(topOTU), varnames=c('species', 'samples'))

### Stacked + percent
(stackedbarplot <- ggplot(melted_otu, aes(fill=reorder(species,-value), y=value, x=samples)) + 
  # coord_flip()+
  geom_bar(position="fill", stat="identity") + 
  # scale_x_discrete(labels = row.names(Ymat))+
  ylab("Relative abundance")+xlab("Samples")+
  # theme_light()+ 
  scale_fill_discrete(name= "Genera") +
  labs(title = "NYC Hanes II", subtitle = "top 5% abundance") +
  theme_bw() +
  theme(axis.text.y = element_text(size = 9),
        axis.text.x = element_blank(),
        legend.position = "bottom"))

# print(stackedbarplot)

ggsave("NYCHanesStackedBar.pdf", stackedbarplot, width = 8, height = 5, device = cairo_pdf)


## Compile Stan models 

setwd("C:/Users/jyotishka/Box Sync/Jyotishka/Microbiome Data/R-Stan/stan-codes")

set.seed(12345)
library(pacman)
p_load(MCMCpack,Compositional,dirmult,reshape,philentropy,
       tidyverse, dplyr, readr, magrittr,zoo,reshape2, 
       ggplot2,rstan, bayesm)

rstan_options(auto_write = TRUE)
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

# stan.dir.fit <- stan_model(file = 'multinomial-horseshoe.stan', model_name = "Symmetric Dirichlet")

stan.hs.fit <- stan_model(file = 'multinomial-horseshoe-marg.stan', 
                          model_name = "Dirichlet Horseshoe")

stan.hsplus.fit <- stan_model(file = 'multinomial-hsplus-marg.stan', 
                              model_name = "Dirichlet HS+")

stan.laplace.fit <- stan_model(file = 'multinomial-laplace-marg.stan', 
                               model_name = "Dirichlet Laplace")

# stan.zidm.fit <- stan_model(file = 'ZI-DM-horseshoe.stan', 
#                           model_name = "ZI-DM Horseshoe") ## not using

###################################
## Sampling from DM Horseshoe Model 
###################################

n.iters = 1000
n.chains = 1
rng.seed = 12345

set.seed(rng.seed)
#dirfit <- dirmult(Ymat)
#save(dirfit, file="NYCHanes/dirfit.Rdata")

load("NYCHanes/dirfit.Rdata")
NYC.data = list(N = nrow(Ymat), ncolY = ncol(Ymat), ncolX = ncol(Xmat),
                 X = Xmat, Y = Ymat, psi = dirfit$theta, scale_icept = 2, d=1) 

ptm = proc.time()
smpls.hs.res = sampling(stan.hs.fit, 
                        data = NYC.data, 
                        iter = n.iters,
                        init = 0,
                        seed = rng.seed,
                        cores = 2,
                        warmup = floor(n.iters/2),
                        chains = n.chains,
                        control = list(adapt_delta = 0.85),
                        refresh = 100)
proc.time()-ptm
save(smpls.hs.res, file="NYCHanes/DMHorseshoe.Rdata")

# summarize results

load("NYCHanes/DMHorseshoe.Rdata")
beta.smpls.hs <- rstan::extract(smpls.hs.res, pars=c("beta"), permuted=TRUE)[[1]]
beta.mean.hs <- apply(beta.smpls.hs, c(2,3), mean)
beta.median.hs <- apply(beta.smpls.hs, c(2,3), median)
beta.mode.hs <- apply(beta.smpls.hs, c(2,3), Mode)
beta.sd.hs <- apply(beta.smpls.hs, c(2,3),sd)
beta.hs.ci <- apply(beta.smpls.hs, c(2,3), quantile, probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals

p <- dim(beta.hs.ci)[2]
q <- dim(beta.hs.ci)[3]
beta.hs.select <- matrix(0, p, q)
for(i in 1:p){
  for(j in 1:q){
    beta.hs.select[i,j] = ifelse(beta.hs.ci[1,i,j]*beta.hs.ci[3,i,j] < 0, 0, 1)
  }
}

# heatmaps
library(patchwork)
p1 <-  simple_heatmap(beta.mean.hs) + labs(title = "Horseshoe - mean")
p2 <- simple_heatmap(beta.sd.hs) + labs(title="Horseshoe - sd")
p3 <- simple_heatmap(beta.hs.select) + labs(title="Horseshoe - selected")
pdf('NYCHanes/DMHorseshoe_plot.pdf')
p1
p2
p3
dev.off()

## more maps
dimnames(beta.hs.select) = dimnames(cor(Xmat, scale(Ymat)))

(beta.hs.select.map = simple_heatmap(beta.hs.select) + coord_flip()+
                     labs(title="Variable Selection: HS") + 
                    xlab("Features")+ylab("OTU"))
# print(beta.hs.select.map)

ggsave("NYCHanes/beta-hs-select-dm-2.pdf", beta.hs.select.map, width = 9, height = 5)


library(patchwork)
pcomb <- Y.map / beta.hs.select.map
ggsave("NYCHanes/compare_dmhsvarsel_corr.png", pcomb, width = 8, height = 9)

dimnames(beta.mean.hs) = dimnames(cor(Xmat, scale(Ymat)))
(beta.hs.mean.map = simple_heatmap(beta.mean.hs) + 
                  # scale_fill_distiller(palette = "blues")+
    coord_flip() + #scale_fill_gradient(trans = 'boxcox') #+ plotTheme()
                   labs(title="Posterior mean: HS") +
                    xlab("Features")+ylab("OTU"))
print(beta.hs.mean.map)

library(patchwork)
pcomb <- Y.map / beta.hs.select.map

ggsave("compare_dmhsvarsel_corr.png", pcomb, width = 8, height = 9)

# table of selected variables
tab_sel_dm_hs <- which(beta.hs.select!=0, arr.ind=T)
tab_sel_dm_hs <- as.data.frame(tab_sel_dm_hs)
tab_sel_dm_hs[,3] <- beta.mean.hs[cbind(tab_sel_dm_hs[,1], tab_sel_dm_hs[,2])]
tab_sel_dm_hs[,4] <- beta.sd.hs[cbind(tab_sel_dm_hs[,1], tab_sel_dm_hs[,2])]
tab_sel_dm_hs[,5] <- colnames(Xmat)[tab_sel_dm_hs[,1]]
tab_sel_dm_hs[,6] <- colnames(Ymat)[tab_sel_dm_hs[,2]]
tab_sel_dm_hs[,7] <- phyloY[tab_sel_dm_hs[,2],6]
colnames(tab_sel_dm_hs) <- c("indX", "indY", "mean", "sd", "nameX", "nameY", "nameGenus")
tab_sel_dm_hs[,c(3,5,7)]
write.csv(tab_sel_dm_hs, file="NYCHanes/DMHorseshoe_summary_backup.csv")


########################################
## Sampling from DM-Horseshoe Plus Model 
########################################

ptm = proc.time()
smpls.hsplus.res <- sampling(stan.hsplus.fit, 
                            data = NYC.data, 
                            iter = n.iters,
                            init = 0,
                            seed = rng.seed,
                            cores = 2,
                            warmup = floor(n.iters/2),
                            chains = n.chains,
                            control = list(adapt_delta = 0.85),
                            refresh = 100)
proc.time()-ptm
save(smpls.hsplus.res, file="NYCHanes/DMHorseshoePlus.Rdata")

load("NYCHanes/DMHorseshoePlus.Rdata")
beta.smpls.hsplus <- rstan::extract(smpls.hsplus.res, pars=c("beta"), permuted=TRUE)[[1]]
beta.mean.hsplus <- apply(beta.smpls.hsplus, c(2,3), mean)
beta.median.hsplus <- apply(beta.smpls.hsplus, c(2,3), median)
beta.mode.hsplus <- apply(beta.smpls.hsplus, c(2,3), Mode)
beta.sd.hsplus <- apply(beta.smpls.hsplus, c(2,3), sd)
beta.hsplus.ci <- apply(beta.smpls.hsplus, c(2,3), quantile, probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals


p <- dim(beta.hsplus.ci)[2]
q <- dim(beta.hsplus.ci)[3]
beta.hsplus.select <- matrix(0, p, q)
for(i in 1:p){
  for(j in 1:q){
    beta.hsplus.select[i,j] = ifelse(beta.hsplus.ci[1,i,j]*beta.hsplus.ci[3,i,j] < 0, 0, 1)
  }
}

## more maps
dimnames(beta.hsplus.select) = dimnames(cor(Xmat, scale(Ymat)))

(beta.hsplus.select.map = simple_heatmap(beta.hsplus.select) + coord_flip()+
    labs(title="Variable Selection: HS+") + 
    xlab("Features")+ylab("OTU"))
# print(beta.hsplus.select.map)

ggsave("NYCHanes/beta-hsplus-select-dm-2.pdf", beta.hsplus.select.map, width = 9, height = 5)


# heatmaps
library(patchwork)
p1 <-  simple_heatmap(beta.mean.hsplus) + labs(title = "DMHorseshoePlus - mean")
p2 <- simple_heatmap(beta.sd.hsplus) + labs(title="DMHorseshoePlus - sd")
p3 <- simple_heatmap(beta.hsplus.select) + labs(title="DMHorseshoePlus - selected")
pdf('NYCHanes/HorseshoePlus_plot.pdf')
p1
p2
p3
dev.off()

# table of selected variables
tab_sel_dm_hsplus <- which(beta.hsplus.select!=0, arr.ind=T)
tab_sel_dm_hsplus <- as.data.frame(tab_sel_dm_hsplus)
tab_sel_dm_hsplus[,3] <- beta.mean.hsplus[cbind(tab_sel_dm_hsplus[,1], tab_sel_dm_hsplus[,2])]
tab_sel_dm_hsplus[,4] <- beta.sd.hsplus[cbind(tab_sel_dm_hsplus[,1], tab_sel_dm_hsplus[,2])]
tab_sel_dm_hsplus[,5] <- colnames(Xmat)[tab_sel_dm_hsplus[,1]]
tab_sel_dm_hsplus[,6] <- colnames(Ymat)[tab_sel_dm_hsplus[,2]]
tab_sel_dm_hsplus[,7] <- phyloY[tab_sel_dm_hsplus[,2],6]
colnames(tab_sel_dm_hsplus) <- c("indX", "indY", "mean", "sd", "nameX", "nameY", "nameGenus")
tab_sel_dm_hsplus[,c(3,5,7)]
write.csv(tab_sel_dm_hsplus, file="NYCHanes/DMHorseshoePlus_summary_backup.csv")



#################################
## Sampling from DM-Laplace Model
#################################

ptm = proc.time()
smpls.lap.res <- sampling(stan.laplace.fit, 
                          data = NYC.data, 
                          iter = n.iters,
                          init = 0,
                          seed = rng.seed,
                          cores = 2,
                          warmup = floor(n.iters/2),
                          chains = n.chains,
                          control = list(adapt_delta = 0.85),
                          refresh = 100)
proc.time()-ptm
save(smpls.lap.res, file="NYCHanes/DMLaplace.Rdata")

load("NYCHanes/DMLaplace.Rdata")
beta.smpls.lap <- rstan::extract(smpls.lap.res, pars=c("beta"), permuted=TRUE)[[1]]
beta.mean.lap <- apply(beta.smpls.lap, c(2,3), mean)
beta.median.lap <- apply(beta.smpls.lap, c(2,3), median)
beta.mode.lap <- apply(beta.smpls.lap, c(2,3), Mode)
beta.sd.lap <- apply(beta.smpls.lap, c(2,3),sd)
beta.lap.ci <- apply(beta.smpls.lap, c(2,3), quantile, probs=c(0.025,0.5,0.975)) #the median line with 95% credible intervals

p <- dim(beta.lap.ci)[2]
q <- dim(beta.lap.ci)[3]
beta.lap.select <- matrix(0, p, q)
for(i in 1:p){
  for(j in 1:q){
    beta.lap.select[i,j] = ifelse(beta.lap.ci[1,i,j]*beta.lap.ci[3,i,j] < 0, 0, 1)
  }
}

sum(beta.lap.select)
sum(beta.hs.select)
sum(beta.hsplus.select)

## more maps
dimnames(beta.lap.select) = dimnames(cor(Xmat, scale(Ymat)))

(beta.lap.select.map = simple_heatmap(beta.lap.select) + coord_flip()+
    labs(title="Variable Selection: lap") + 
    xlab("Features")+ylab("OTU"))
print(beta.lap.select.map)

ggsave("NYCHanes/beta-lap-select-dm-2.pdf", beta.lap.select.map, width = 9, height = 5)


# heatmaps
library(patchwork)
p1 <-  simple_heatmap(beta.mean.lap) + labs(title = "DMLaplace - mean")
p2 <- simple_heatmap(beta.sd.lap) + labs(title="DMLaplace - sd")
p3 <- simple_heatmap(beta.lap.select) + labs(title="DMLaplace - selected")
pdf("NYCHanes/DMLaplace_plot.pdf")
p1
p2
p3
dev.off()

# table of selected variables
tab_sel_dm_lap <- which(beta.lap.select!=0, arr.ind=T)
tab_sel_dm_lap <- as.data.frame(tab_sel_dm_lap)
tab_sel_dm_lap[,3] <- beta.mean.lap[cbind(tab_sel_dm_lap[,1], tab_sel_dm_lap[,2])]
tab_sel_dm_lap[,4] <- beta.sd.lap[cbind(tab_sel_dm_lap[,1], tab_sel_dm_lap[,2])]
tab_sel_dm_lap[,5] <- colnames(Xmat)[tab_sel_dm_lap[,1]]
tab_sel_dm_lap[,6] <- colnames(Ymat)[tab_sel_dm_lap[,2]]
tab_sel_dm_lap[,7] <- phyloY[tab_sel_dm_lap[,2],6]
colnames(tab_sel_dm_lap) <- c("indX", "indY", "mean", "sd", "nameX", "nameY", "nameGenus")
tab_sel_dm_lap[,c(3,5,7)]
write.csv(tab_sel_dm_lap, file="NYCHanes/DMLaplace_summary_backup.csv")

save.image(file = "run_NYCHanes_outputs.RData")

load("run_NYCHanes_outputs.RData")
