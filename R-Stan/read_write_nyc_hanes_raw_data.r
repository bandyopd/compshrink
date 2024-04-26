# library(nychanesmicrobiome)
library(dplyr)
library(phyloseq)
NYC_HANES <- loadQiimeData() %>% annotateFactors(.)

biosis.tsv <- system.file("extdata","biosis.tsv", package="nychanesmicrobiome", mustWork = TRUE)
biosis <- read.csv(biosis.tsv,header = TRUE, sep = '\t')

## meta data containing clinical features

metadata <- data.frame(sample_data(NYC_HANES))
# check missingness
dim(metadata)
dim(metadata[,colSums(is.na(metadata))==0])
dim(metadata[,colSums(is.na(metadata))<=1])
dim(metadata[,colSums(is.na(metadata))<=5])
#[1] 282 714
#[1] 282 200
#[1] 282 223
#[1] 282 251
write.csv(metadata, "NYCHanes/metadata.csv")

## features with potentials

meta <- read.csv("NYCHanes/metadata.csv", row.names=1)
selected <- c("SPAGE", "BMI", # continuous
              "GENDER", "RACE", # demographic categories
              "EDU4CAT", "INCOME", "smokingstatus", # multiple categories
              "DBTS_NEW", "OHQ_3", "BPQ_2", "US_NSAID") # binary
meta_filt <- meta[,selected]
summary(meta_filt)
# income in categories
meta_filt$INCOMECAT <- rep(NA, nrow(meta_filt))
meta_filt$INCOMECAT[meta$INCOME<20] <- "<20k"
meta_filt$INCOMECAT[meta$INCOME<45 & meta_filt$INCOME>=20] <- "20k-45k"
meta_filt$INCOMECAT[meta$INCOME>=45] <- ">=45k"
meta_filt$INCOMECAT[is.na(meta$INCOME)] <- "NA"
meta_filt$INCOMECAT <- factor(meta_filt$INCOMECAT, levels=c("NA","<20k", "20k-45k", ">=45k"))
meta_filt$INCOME <- NULL
# race in categories
meta_filt$RACE <- factor(meta_filt$RACE, levels=c("Other", "Asian", "Non-Hispanic Black", "Hispanic", "Non-Hispanic White"))
levels(meta_filt$RACE) <- c("Other", "Asian", "Black", "Hispanic", "White")
# education in categories
meta_filt$EDU4CAT <- factor(meta_filt$EDU4CAT, levels=c("Less than High school diploma", "High school graduate/GED", "Some College or associate's degree", "College graduate or more"))
levels(meta_filt$EDU4CAT) <- c("<HighSchool", "HighSchool", "SomeCollege", "College")
# smoking status in categories
meta_filt$smokingstatus <- factor(meta_filt$smokingstatus, levels=c("Never smoker", "Former smoker", "Secondhand", "Alternative smoker", "Cigarette"))
levels(meta_filt$smokingstatus) <- c("Never", "Former", "SecondHand", "Alternative", "Cigarette")
# blood pressure in categories
meta_filt$BPQ_2 <- factor(meta_filt$BPQ_2, levels=c("2","1"))
levels(meta_filt$BPQ_2) <- c("No", "Yes")
# diabetes in categories
meta_filt$DBTS_NEW <- factor(meta_filt$DBTS_NEW, levels=c("No","Yes"))
# gum disease in categories
meta_filt$OHQ_3 <- factor(meta_filt$OHQ_3, levels=c("No","Yes"))
# use of NSAID in categories
meta_filt$US_NSAID <- factor(meta_filt$US_NSAID, levels=c("2","1"))
levels(meta_filt$US_NSAID) <- c("No", "Yes")

summary(meta_filt)

## create numerical matrix for features

# function for dummy variable generation
make_dummy <- function(x, varname){
  namecat <- levels(x)
  ncat <- length(namecat)
  dummymat <- matrix(0, length(x), ncat)
  colnames(dummymat) <- paste0(varname, ":", namecat)
  for (i in 1:ncat){
    dummymat[,i] <- as.numeric(x==namecat[i])
  }
  return(dummymat)
}

meta_mat <- meta_filt[,c("SPAGE","BMI")]
for (i in 3:ncol(meta_filt)){
  newmat <- make_dummy(meta_filt[,i], colnames(meta_filt)[i])
  meta_mat <- cbind(meta_mat, newmat[,-1])
  colnames(meta_mat)[-(1:(ncol(meta_mat)-ncol(newmat)+1))] <- colnames(newmat)[-1]
}
dim(meta_mat)

# label for samples that should be removed
rm <- is.na(rowSums(meta_mat))  # remove subjects with NAs
rm[meta$PREGNANT==1] <- TRUE  # remove pregnant subjects
table(rm)
meta_mat$REMOVE <- as.numeric(rm)

write.csv(meta_mat, file="NYCHanes/metadata_matrix.csv")


## read count table at genus level

NYC_HANES.genus <- tax_glom(NYC_HANES, taxrank = "Genus")
NYC_HANES.genus.otu <- t(otu_table(NYC_HANES.genus))
NYC_HANES.genus.name <- tax_table(NYC_HANES.genus)[,-7]
dim(NYC_HANES.genus.otu) # [1] 282 148
# read count table n by p, colnames are OTU representatives, all genera
write.csv(NYC_HANES.genus.otu, file="NYCHanes/otu_genus.csv") 
# phylogeny info of OTU representatives, all genera
write.csv(NYC_HANES.genus.name, file="NYCHanes/phylogeny_genus.csv")

## read count table at genus level, limited to genera present in >pct samples

pct <- 50
pct <- 20

NYC_HANES.genus.freq <- merge_taxa(NYC_HANES.genus,
                                   taxa_names(filter_taxa(NYC_HANES.genus, function(x)
                                     mean(x!=0) < pct/100, TRUE))) 
NYC_HANES.genus.freq.otu <- t(otu_table(NYC_HANES.genus.freq))
NYC_HANES.genus.freq.name <- tax_table(NYC_HANES.genus.freq)[,-7]
dim(NYC_HANES.genus.freq.otu) # [1] 282 116
# read count table n by p, colnames are OTU representatives, genera present in >90% samples
write.csv(NYC_HANES.genus.freq.otu, file=paste0("NYCHanes/otu_genus_", pct, "pc.csv"))
# phylogeny info of OTU representatives, genera present in >50% samples
write.csv(NYC_HANES.genus.freq.name, file=paste0("NYCHanes/phylogeny_genus_", pct, "pc.csv"))

