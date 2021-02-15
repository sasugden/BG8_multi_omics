############## NOTES AND CONVENTIONS ################
# All nitrogen source fold changes are calculated as NMS/AMS. 
# Therefore, -/+ log2FC indicates AMS/NMS, respectively.
# All carbon source fold changes are calculated as MeOH/CH4.
# Therefore, -/+ log2FC indicates CH4/MeOH.
# To reuse and simplify code, some annotations refer to "left" or "right" treatments.
# "Left" and "right" refer to "negative" or "positive" log2FC, respectively.

# Define colors for plotting.
treatment_colors <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
names(treatment_colors) <- c("CH4-AMS", "MeOH-AMS", "CH4-NMS", "MeOH-NMS")

######################################## I. INITIAL DATA PROCESSING ##################################
#### Load packages and define functions. ####
setwd("E:/research projects/x_R_RESTARTS/bg8_omics/")

library(vegan)
library(edgeR)
library(magrittr)
library(clusterProfiler)
library(ggplot2)

veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#### Import data from external files. ####
#### Import sample metadata
rnaseq_sample_data <- read.csv("raw_data/rnaseq_metadata.csv") # RNA sample metadata
rownames(rnaseq_sample_data) <- rnaseq_sample_data$SampleID

metab_sample_data <- read.csv("raw_data/metab_metadata.csv") # Metabolite sample metadata
rownames(metab_sample_data) <- metab_sample_data$SampleID

# Ensure that the order of factors in both metadata files is the same.
treatment_order <- c("CH4-AMS", "MeOH-AMS", "CH4-NMS", "MeOH-NMS")
rnaseq_sample_data$Treatment <- factor(rnaseq_sample_data$Treatment, levels=treatment_order)
metab_sample_data$Treatment <- factor(metab_sample_data$Treatment, levels=treatment_order)

# Create lists of sample names associated with each single nutrient, for pairwise comparisons later.
rnaseq_sample_keys <- list(
  rna_CH4_samples = c("S01", "S02", "S03", "S04", "S05", "S06"),
  rna_MeOH_samples = c("S07", "S08", "S09", "S10", "S11", "S12"),
  rna_AMS_samples = c("S04", "S05", "S06", "S10", "S11", "S12"),
  rna_NMS_samples = c("S01", "S02", "S03", "S07", "S08", "S09")
)
metab_sample_names <- list( # identify the metabolome samples in each treatment (methane, methanol, AMS, NMS)
  metab_CH4_samples = c("UABT00017", "UABT00018", "UABT00019", "UABT00020", "UABT00025", "UABT00026", "UABT00027", "UABT00028"),
  metab_MeOH_samples = c("UABT00021", "UABT00022", "UABT00023", "UABT00024", "UABT00029", "UABT00030", "UABT00031", "UABT00032"),
  metab_AMS_samples = c("UABT00017", "UABT00018", "UABT00019", "UABT00020", "UABT00021", "UABT00022", "UABT00023", "UABT00024"),
  metab_NMS_samples = c("UABT00025", "UABT00026", "UABT00027", "UABT00028", "UABT00029", "UABT00030", "UABT00031", "UABT00032")
)

tests <- c("AMSvNMS_Methane", "AMSvNMS_Methanol", "CH4vCH3OH_AMS", "CH4vCH3OH_NMS", "Interactions")

#### Import raw transcript and metabolite counts.
rnaseq_counts <- list()
rnaseq_counts$raw <- read.csv("raw_data/rnaseq_counts.csv")
rownames(rnaseq_counts$raw) <- rnaseq_counts[["raw"]]$Accession
rnaseq_counts[["raw"]]$Accession <- NULL

metab_counts <- list()
metab_counts$raw <- read.csv("raw_data/metab_counts.csv")
rownames(metab_counts$raw) <- metab_counts[["raw"]]$MetaboliteID
metab_counts[["raw"]]$MetaboliteID <- NULL

#### Import annotation information for transcripts and metabolites.
rnaseq_annotations <- read.csv("raw_data/rnaseq_annotations.csv")
rownames(rnaseq_annotations) <- rnaseq_annotations$Accession

metab_annotations <- read.csv("raw_data/metab_annotations.csv")
rownames(metab_annotations) <- metab_annotations$MetaboliteID

metab_annotations$Superfamily <- factor(metab_annotations$Superfamily)
metab_annotations$Subfamily <- forcats::fct_inorder(metab_annotations$Subfamily)

#### Create a mapping file that gives 1:1 relationships for KEGG ontology to gene IDs ####
temp_sub <- rnaseq_annotations[,c("Accession", "KEGG_Enzymes")]
temp_sub <- tidyr::separate(temp_sub, col=KEGG_Enzymes, into=c("KEGG1", "KEGG2", "KEGG3", "KEGG4", "KEGG5", "KEGG6",
                                                               "KEGG7", "KEGG8", "KEGG9", "KEGG10", "KEGG11", "KEGG12",
                                                               "KEGG13", "KEGG14", "KEGG15", "KEGG16", "KEGG17", "KEGG18",
                                                               "KEGG19", "KEGG20", "KEGG21", "KEGG22", "KEGG23", "KEGG24",
                                                               "KEGG25", "KEGG26", "KEGG27", "KEGG28", "KEGG29", "KEGG30"),
                            sep=",", remove=TRUE, convert=FALSE)

vectors <- list()
for(i in 1:nrow(temp_sub)){vectors[[i]] <- as.vector(as.character(temp_sub[i,][c(2:31)]))}
vectors2 <- list()
for(i in 1:nrow(temp_sub)){vectors2[[i]] <- rep(as.character(temp_sub[i,][1,1]), 30)}

temp1 <- c(vectors[[1]], vectors[[2]])
for(i in 3:length(vectors)){temp1 <- c(temp1, vectors[[i]])}
temp2 <- c(vectors2[[1]], vectors2[[2]])
for(i in 3:length(vectors2)){temp2 <- c(temp2, vectors2[[i]])}

temp <- data.frame(
  Accession=temp2,
  KEGG_Enzymes=temp1)
temp <- subset(temp, is.na(KEGG_Enzymes)=="FALSE")
temp[] <- data.frame(lapply(temp, function(x) {gsub("ko:", "", x)}))

temp <- merge(temp, rnaseq_annotations[,c("Accession", "UniProtKB")], by="Accession", all=FALSE)

rnaseq_mappings <- temp

rm(temp_sub, vectors, temp, temp1, temp2, vectors2, i)

#### Scale/normalize transcript and metabolite abundances ####
# Calculate TPM values for each transcript.
rnaseq_counts$tpm <- rnaseq_counts$raw
for(i in 1:ncol(rnaseq_counts$tpm)){
  vector <- rnaseq_counts$tpm[,i]/(rnaseq_annotations$DNA_Length/1000)
  scale <- sum(vector)/1000000
  rnaseq_counts$tpm[,i] <- vector/scale
  rm(vector, scale)
}

# Calculate the mean TPM value per treatment.
temp <- as.data.frame(t(rnaseq_counts$tpm))
temp$Treatment <- rnaseq_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(mean=mean))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(rnaseq_counts$tpm)
temp$Accession <- rownames(temp)
rnaseq_counts$avg <- temp
rm(temp)

# Calculate the standard deviation of TPM values per treatment.
temp <- as.data.frame(t(rnaseq_counts$tpm))
temp$Treatment <- rnaseq_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(sd=sd))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(rnaseq_counts$tpm)
temp$Accession <- rownames(temp)
rnaseq_counts$sds <- temp
rm(temp)

# Calculate median-scaled metabolite abundances.
metab_counts$scaled <- as.data.frame(t(metab_counts$raw))
for(i in 1:ncol(metab_counts$scaled )){
  metab_counts$scaled[,i] <- metab_counts$scaled [,i]/median(metab_counts$scaled[,i], na.rm=TRUE)
}
metab_counts$scaled <- as.data.frame(t(metab_counts$scaled))

# Impute missing values as one-half the minimum abundance of each metabolite.
metab_counts$imputed <- as.data.frame(metab_counts$scaled)
for(i in 1:nrow(metab_counts$imputed)){
  rowmin <- 0.5*min(metab_counts$imputed[i,], na.rm=TRUE)
  for(j in 1:ncol(metab_counts$imputed)){
    if(is.na(metab_counts$imputed[i,j])=="TRUE"){
      metab_counts$imputed[i,j] <- rowmin
    }
  }
  rm(rowmin)
}

# Calculate the mean metabolite abundance value per treatment.
temp <- as.data.frame(t(metab_counts$imputed))
temp$Treatment <-  metab_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(mean=mean))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(metab_counts$imputed)
temp$MetaboliteID <- rownames(temp)
metab_counts$avg <- temp
rm(temp)

# Calculate standard deviations of metabolite abundances per treatment.
temp <- as.data.frame(t(metab_counts$imputed))
temp$Treatment <-  metab_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(sd=sd))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(metab_counts$imputed)
temp$MetaboliteID <- rownames(temp)
metab_counts$sds <- temp
rm(temp, i, j)

#### Identify undetected transcripts and metabolites in each sample. ####
## Transcripts ####
# Define transcripts as either present (1) or absent (0).
temp <- rnaseq_counts$tpm
temp[temp > 0] <- 1

# Count the number of times each transcript was detected in a treatment.
temp <- as.data.frame(t(temp))
temp$Treatment <- as.character(rnaseq_sample_data$Treatment)
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(list(sum=sum)) %>%
  as.data.frame()
rownames(temp) <- temp$Treatment
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
rownames(temp) <- rnaseq_annotations$Accession

# Test for significant differences in transcript detection rates.
test <- rnaseq_counts[["raw"]]
test[test != 0] <- 1
test[is.na(test)] <- 0

rnaseq_sample_data$Missing <- nrow(test)-colSums(test)
ggplot(rnaseq_sample_data, aes(x=Treatment, y=Missing)) + geom_boxplot()
summary(aov(Missing ~ Treatment, rnaseq_sample_data)) # F=1.721, df=3, p=0.24
TukeyHSD(aov(Missing ~ Treatment, rnaseq_sample_data))

rm(test)

# Create a data frame where undetected transcripts are rows and the treatments are columns.
temp2 <- list(CH4.AMS = data.frame(Accession=rownames(subset(temp, `CH4-AMS`==0)),
                                   CH4.AMS=rep("Not detected in CH4:AMS")),
              CH4.NMS = data.frame(Accession=rownames(subset(temp, `CH4-NMS`==0)),
                                   CH4.NMS=rep("Not detected in CH4:NMS")),
              MeOH.AMS = data.frame(Accession=rownames(subset(temp, `MeOH-AMS`==0)),
                                    MeOH.AMS=rep("Not detected in MeOH:AMS")),
              MeOH.NMS = data.frame(Accession=rownames(subset(temp, `MeOH-NMS`==0)),
                                    MeOH.NMS=rep("Not detected in MeOH:NMS")),
              All = data.frame(Accession=rownames(subset(temp, rowSums(temp)==0)),
                               All=rep("Not detected in experiment")))

rnaseq_undetected <- merge(temp2[[1]], temp2[[2]], by="Accession", all=TRUE)
rnaseq_undetected <- merge(rnaseq_undetected, temp2[[3]], by="Accession", all=TRUE)
rnaseq_undetected <- merge(rnaseq_undetected, temp2[[4]], by="Accession", all=TRUE)
rnaseq_undetected <- merge(rnaseq_undetected, temp2[[5]], by="Accession", all=TRUE)
rnaseq_undetected <- merge(rnaseq_undetected, rnaseq_annotations[,c("Accession", "COG", "GenBank_Name")],
                           by="Accession", all=FALSE)
rm(temp, temp2)

## Metabolites ####
# Test for significant differences in unobserved metabolites among treatments
test <- metab_counts[["raw"]]
test[test != 0] <- 1
test[is.na(test)] <- 0

metab_sample_data$Missing <- nrow(test)-colSums(test)
ggplot(metab_sample_data, aes(x=Treatment, y=Missing)) + geom_boxplot()
summary(aov(Missing ~ Treatment, metab_sample_data)) # F=227.1, df=3, p<0.001
TukeyHSD(aov(Missing ~ Treatment, metab_sample_data))

#                     diff        lwr       upr     p adj
# MeOH-AMS-CH4-AMS   71.75  62.689986  80.81001 0.0000000 ***
# CH4-NMS-CH4-AMS     6.00  -3.060014  15.06001 0.2533186
# MeOH-NMS-CH4-AMS   23.75  14.689986  32.81001 0.0000257 ***
# CH4-NMS-MeOH-AMS  -65.75 -74.810014 -56.68999 0.0000000 ***
# MeOH-NMS-MeOH-AMS -48.00 -57.060014 -38.93999 0.0000000 ***
# MeOH-NMS-CH4-NMS   17.75   8.689986  26.81001 0.0004149

rm(test)

# Obtain a list of undetected metabolites (same method as for transcripts)
temp <- metab_counts$raw
temp[is.na(temp)] <- 0 # Convert undetected metabolites to 0
temp[temp > 0] <- 1 # Convert detected metabolites to 1
temp <- as.data.frame(t(temp))
temp$Treatment <- as.character(metab_sample_data$Treatment) # Add metadata
temp <- temp %>% dplyr::group_by(Treatment) %>% 
  dplyr::summarise_all(list(sum=sum)) %>% # Count number of detected metabolites per treatment.
  as.data.frame()
rownames(temp) <- temp$Treatment
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
rownames(temp) <- metab_annotations$MetaboliteID

temp2 <- list(CH4.AMS = data.frame(MetaboliteID=rownames(subset(temp, `CH4-AMS`==0)), # Summarise undetected metabolites.
                                   CH4.AMS=rep("Not detected in CH4:AMS")),
              CH4.NMS = data.frame(MetaboliteID=rownames(subset(temp, `CH4-NMS`==0)),
                                   CH4.NMS=rep("Not detected in CH4:NMS")),
              MeOH.AMS = data.frame(MetaboliteID=rownames(subset(temp, `MeOH-AMS`==0)),
                                    MeOH.AMS=rep("Not detected in MeOH:AMS")),
              MeOH.NMS = data.frame(MetaboliteID=rownames(subset(temp, `MeOH-NMS`==0)),
                                    MeOH.NMS=rep("Not detected in MeOH:NMS")))

metab_undetected <- merge(temp2[[1]], temp2[[2]], by="MetaboliteID", all=TRUE)
metab_undetected <- merge(metab_undetected, temp2[[3]], by="MetaboliteID", all=TRUE)
metab_undetected <- merge(metab_undetected, temp2[[4]], by="MetaboliteID", all=TRUE)
metab_undetected <- merge(metab_annotations[,c("MetaboliteID", "Superfamily", "Biochemical")], metab_undetected, 
                          by="MetaboliteID", all=FALSE)
rm(temp, temp2)

############## CLUSTERING ANALYSIS ###############
########## Principal components analysis ####
#### Transcriptome ####
rnaseq_pca_data <- as.data.frame(t(rnaseq_counts$tpm)) 
rnaseq_pca_data <- log(rnaseq_pca_data+0.01) # Natural log transformation, including a pseudocount

# Calculate the PCA using the 'rda' function.
rnaseq_pca_results <- list()
rnaseq_pca_results$pca <- rda(rnaseq_pca_data)

# Obtain the percentage of variation on each axis.
rnaseq_pca_results[["pca"]][["CA"]][["eig"]][[1]] / rnaseq_pca_results[["pca"]][["CA"]]$tot.chi # 36.67%
rnaseq_pca_results[["pca"]][["CA"]][["eig"]][[2]] / rnaseq_pca_results[["pca"]][["CA"]]$tot.chi # 21.37%

# Store the axis loadings for each sample as a separate object for plotting later.
rnaseq_pca_results$scores <- as.data.frame(
  merge(scores(rnaseq_pca_results$pca, choices = c(1:3), display="sites"),
        rnaseq_sample_data, by="row.names", all=TRUE))

# Calculate standard ellipses for plotting later.
plot.new()
pca.ellipse <- ordiellipse(rnaseq_pca_results$pca, 
                            rnaseq_pca_results[["scores"]]$Treatment, 
                            display="sites", 
                            kind="sd", 
                            conf=0.95, 
                            label=T)

rnaseq_pca_results$ellipse <- data.frame()

    
for(g in levels(rnaseq_pca_results[["scores"]]$Treatment)){
  rnaseq_pca_results$ellipse <- rbind(rnaseq_pca_results$ellipse,
                                      cbind(as.data.frame(with(rnaseq_pca_results$scores[rnaseq_pca_results[["scores"]]$Treatment==g,],
                                                               veganCovEllipse(pca.ellipse[[g]]$cov,
                                                                               pca.ellipse[[g]]$center,
                                                                               pca.ellipse[[g]]$scale)))
                                            ,Treatment=g))
}
rm(pca.ellipse, rnaseq_pca_data)

#### Metabolome ####
metab_pca_data <- as.data.frame(t(metab_counts$imputed)) # Extract data.
metab_pca_data <- log(metab_pca_data) # Natural log transformation (no need for pseudocount)

# Calculate the PCA using the 'rda' function.
metab_pca_results <- list()
metab_pca_results$pca <- rda(metab_pca_data)

# Obtain the percentage of variation on each axis.
metab_pca_results[["pca"]][["CA"]][["eig"]][[1]] / metab_pca_results[["pca"]][["CA"]]$tot.chi # 68.75%
metab_pca_results[["pca"]][["CA"]][["eig"]][[2]] / metab_pca_results[["pca"]][["CA"]]$tot.chi # 14.14%

# Store the axis loadings for each sample as a separate object for plotting later.
metab_pca_results$scores <- as.data.frame(
  merge(scores(metab_pca_results$pca, choices = c(1:3), display="sites"),
        metab_sample_data, by="row.names", all=TRUE))

# Calculate standard ellipses for plotting later.
plot.new()
pca.ellipse <- ordiellipse(metab_pca_results$pca, 
                           metab_pca_results[["scores"]]$Treatment, 
                           display="sites", 
                           kind="sd", 
                           conf=0.95, 
                           label=T)

metab_pca_results$ellipse <- data.frame()

for(g in levels(metab_pca_results[["scores"]]$Treatment)){
  metab_pca_results$ellipse <- rbind(metab_pca_results$ellipse,
                                      cbind(as.data.frame(with(metab_pca_results$scores[metab_pca_results[["scores"]]$Treatment==g,],
                                                               veganCovEllipse(pca.ellipse[[g]]$cov,
                                                                               pca.ellipse[[g]]$center,
                                                                               pca.ellipse[[g]]$scale)))
                                            ,Treatment=g))
}
rm(pca.ellipse, metab_pca_data)

########## Sparse partial least squares discriminant analysis ####
#### Transcriptome ####
class <- as.character(rnaseq_sample_data$Treatment)
summary(class) # Ensure the treatments are appropriately specified.

data <- as.data.frame(t(log(rnaseq_counts$tpm+0.01)))
dim(data) # Ensure samples are rows and observations are columns.

list.keepX <- c(1:9,  seq(10, 300, 10))
tune.splsda.srbct <- tune.splsda(data, class, ncomp = 3, validation = 'Mfold', folds = 2, 
                                 progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                 test.keepX = list.keepX, nrepeat = 100, cpus = 3)
error <- tune.splsda.srbct$error.rate  # error rate per component for the keepX grid
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
plot(tune.splsda.srbct, col = color.jet(3))

rnaseq_splsda <- list()
rnaseq_splsda[["model"]] <- splsda(data, class, ncomp = ncomp, keepX = c(5,20,30)) 

plotIndiv(rnaseq_splsda[["model"]], comp = c(1,2),
          group = class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'Initial sPLS-DA', theme=theme(panel.backgound=element_blank())) 

rnaseq_splsda[["loadings"]] <- list(Component1=data.frame(Accession=selectVar(rnaseq_splsda[["model"]], comp = 1)$name,
                                                          Value=selectVar(rnaseq_splsda[["model"]], comp = 1)$value),
                                    Component2=data.frame(Accession=selectVar(rnaseq_splsda[["model"]], comp = 2)$name,
                                                          Value=selectVar(rnaseq_splsda[["model"]], comp = 2)$value),
                                    Component3=data.frame(Accession=selectVar(rnaseq_splsda[["model"]], comp = 3)$name,
                                                          Value=selectVar(rnaseq_splsda[["model"]], comp = 3)$value))
for(i in 1:3){
  colnames(rnaseq_splsda[["loadings"]][[i]]) <- c("Accession", "Axis.Score")
  rnaseq_splsda[["loadings"]][[i]] <- merge(rnaseq_splsda[["loadings"]][[i]], rnaseq_annotations[,c("Accession","COG","Chromosome","GenBank_Name")], by="Accession", all=FALSE)
  rnaseq_splsda[["loadings"]][[i]] <- merge(rnaseq_splsda[["loadings"]][[i]], rnaseq_counts$avg, by="Accession", all=FALSE)
  rnaseq_splsda[["loadings"]][[i]] <- rnaseq_splsda[["loadings"]][[i]][order(rnaseq_splsda[["loadings"]][[i]]$Axis.Score), ]
}

rm(select.keepX, ncomp, error, tune.splsda.srbct, list.keepX, class, data)

#### Metabolome ####
class <- as.character(metab_sample_data$Treatment)
summary(class) # Ensure the treatments are appropriately specified.

data <- as.data.frame(t(log(metab_counts$imputed)))
dim(data) # Ensure samples are rows and observations are columns.

list.keepX <- c(1:9,  seq(10, 300, 10))
tune.splsda.srbct <- tune.splsda(data, class, ncomp = 3, validation = 'Mfold', folds = 2, 
                                 progressBar = TRUE, dist = 'max.dist', measure = "BER",
                                 test.keepX = list.keepX, nrepeat = 100, cpus = 3)
error <- tune.splsda.srbct$error.rate  # error rate per component for the keepX grid
ncomp <- tune.splsda.srbct$choice.ncomp$ncomp # optimal number of components based on t-tests
select.keepX <- tune.splsda.srbct$choice.keepX[1:ncomp]  # optimal number of variables to select
plot(tune.splsda.srbct, col = color.jet(3))

metab_splsda <- list()
metab_splsda[["model"]] <- splsda(data, class, ncomp = ncomp, keepX = c(5,30,5)) 

plotIndiv(metab_splsda[["model"]], comp = c(1,3),
          group = class, ind.names = FALSE, 
          ellipse = TRUE, legend = TRUE, title = 'Initial sPLS-DA', theme=theme(panel.backgound=element_blank())) 

metab_splsda[["loadings"]] <- list(Component1=data.frame(Accession=selectVar(metab_splsda[["model"]], comp = 1)$name,
                                                         Value=selectVar(metab_splsda[["model"]], comp = 1)$value),
                                   Component2=data.frame(Accession=selectVar(metab_splsda[["model"]], comp = 2)$name,
                                                         Value=selectVar(metab_splsda[["model"]], comp = 2)$value),
                                   Component3=data.frame(Accession=selectVar(metab_splsda[["model"]], comp = 3)$name,
                                                         Value=selectVar(metab_splsda[["model"]], comp = 3)$value))
for(i in 1:3){
  colnames(metab_splsda[["loadings"]][[i]]) <- c("Accession", "Axis.Score")
  metab_splsda[["loadings"]][[i]]$MetaboliteID <- rownames(metab_splsda[["loadings"]][[i]])
  metab_splsda[["loadings"]][[i]] <- merge(metab_splsda[["loadings"]][[i]],
                                           metab_annotations[,c("MetaboliteID","Superfamily","Subfamily","Biochemical")],
                                           by="MetaboliteID", all=FALSE)
  metab_splsda[["loadings"]][[i]] <- merge(metab_splsda[["loadings"]][[i]],
                                           metab_counts$avg, by="MetaboliteID", all=FALSE)
  metab_splsda[["loadings"]][[i]] <- metab_splsda[["loadings"]][[i]][order(metab_splsda[["loadings"]][[i]]$Axis.Score), ]
}

rm(select.keepX, ncomp, error, tune.splsda.srbct, list.keepX, class, data)

######################################## II. DOWNSTREAM RNASEQ ANALYSIS ################
#### Differential abundance ####
rnaseq_diff_abund <- list()

for(p in 1:4){ # For each pairwise comparison among treatments
  dgList <- DGEList(counts=rnaseq_counts$raw[,c(rnaseq_sample_keys[[p]])], # Subset to the appropriate samples
                    genes=rownames(rnaseq_counts$raw))
  if(p %in% c(1:2)){designMat <- model.matrix(~Nitrogen, data=rnaseq_sample_data[c(rnaseq_sample_keys[[p]]),])} # Compare by carbon source
  if(p %in% c(3:4)){designMat <- model.matrix(~Carbon, data=rnaseq_sample_data[c(rnaseq_sample_keys[[p]]),])} # Compare by nitrogen source
  
  # Differential abundance testing using edgeR. See the edgeR tutorial for details.
  dgList <- calcNormFactors(dgList, method="TMM")
  plotMDS(dgList)
  dgList <- estimateGLMCommonDisp(dgList, design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
  plotBCV(dgList)
  
  fit <- glmFit(dgList, designMat)
  lrt <- glmLRT(fit, coef=2)
  edgeR_result <- topTags(lrt)
  
  # Extract the edgeR results and calculate BH-adjusted p-values.
  diff.abund.results <- lrt[["table"]]
  diff.abund.results$p.adj <- p.adjust(diff.abund.results$PValue, method="BH")
  diff.abund.results$log10q <- -1*log(diff.abund.results$p.adj, base=10)
  
  # Add information on the mean abundance of each gene.
  if(p==1){diff.abund.results <- cbind(diff.abund.results, rnaseq_counts$avg[,c("CH4-AMS", "CH4-NMS")])}
  if(p==2){diff.abund.results <- cbind(diff.abund.results, rnaseq_counts$avg[,c("MeOH-AMS", "MeOH-NMS")])}
  if(p==3){diff.abund.results <- cbind(diff.abund.results, rnaseq_counts$avg[,c("CH4-AMS", "MeOH-AMS")])}
  if(p==4){diff.abund.results <- cbind(diff.abund.results, rnaseq_counts$avg[,c("CH4-NMS", "MeOH-NMS")])}
  
  # Process and store the data frame.
  diff.abund.results$Accession <- rownames(diff.abund.results)
  diff.abund.results <- merge(diff.abund.results, 
                              rnaseq_annotations[,c("Accession","COG","GenBank_Name")], 
                              by="Accession", all=TRUE)
  
  rownames(diff.abund.results) <- diff.abund.results$Accession
  diff.abund.results$Accession <- NULL
  
  diff.abund.results <- diff.abund.results[,c(1,4,5,6,9,10,7,8)]
  
  rnaseq_diff_abund[[p]] <- diff.abund.results
  
  rm(diff.abund.results, fit, lrt, edgeR_result, dgList, designMat)
}
names(rnaseq_diff_abund) <- tests[c(1:4)]

#### Interaction effects ####
# Create a new model design testing for interaction effects in all four samples.
dgList <- DGEList(counts=rnaseq_counts$raw,
                  genes=rownames(rnaseq_counts$raw))
designMat <- model.matrix(~Nitrogen*Carbon, data=rnaseq_sample_data)

dgList <- calcNormFactors(dgList, method="TMM")
plotMDS(dgList)
dgList <- estimateGLMCommonDisp(dgList, design=designMat)
dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
plotBCV(dgList)

fit <- glmFit(dgList, designMat)
lrt <- glmLRT(fit, coef=4)
edgeR_result <- topTags(lrt)

diff.abund.results <- lrt[["table"]]
diff.abund.results$p.adj <- p.adjust(diff.abund.results$PValue, method="BH")
diff.abund.results$log10q <- -1*log(diff.abund.results$p.adj, base=10)

diff.abund.results <- cbind(diff.abund.results, rnaseq_counts$avg[,c("CH4-AMS", "CH4-NMS", "MeOH-AMS", "MeOH-NMS")])

diff.abund.results$Accession <- rownames(diff.abund.results)
diff.abund.results <- merge(diff.abund.results, 
                            rnaseq_annotations[,c("Accession","COG","GenBank_Name")], 
                            by="Accession", all=TRUE)
rownames(diff.abund.results) <- diff.abund.results$Accession

rnaseq_diff_abund[[5]] <- diff.abund.results
names(rnaseq_diff_abund)[5] <- "Interactions"

rm(fit, lrt, edgeR_result, dgList, designMat, diff.abund.results)

#### COG over-representation analysis ####
## Count the total number of genes assigned to each COG. ####
# Define a function to do this.
temp <- subset(rnaseq_counts$tpm, rowSums(rnaseq_counts$tpm) > 0) # Only use detected transcripts.
rnaseq_cog_counts <- subset(rnaseq_annotations, Accession %in% rownames(temp)) # Subset annotations data.

cog_counter <- function(table){
  temp <- tidyr::separate(table, col=COG, into=c("COG1", "COG2", "COG3"),
                          sep=", ", remove=TRUE, convert=FALSE)
  temp$Blank <- NULL
  temp1 <- temp %>% dplyr::group_by(COG1) %>% dplyr::count()
  temp2 <- temp %>% dplyr::group_by(COG2) %>% dplyr::count()
  temp3 <- temp %>% dplyr::group_by(COG3) %>% dplyr::count()
  colnames(temp1) <- c("COG", "Count1")
  colnames(temp2) <- c("COG", "Count2")
  colnames(temp3) <- c("COG", "Count3")
  temp <- merge(temp1, temp2, by=c("COG"), all=TRUE)
  temp <- merge(temp, temp3, by=c("COG"), all=TRUE)
  rm(temp1, temp2, temp3)
  temp$Count <- rowSums(temp[,c(2:4)], na.rm=TRUE)
  temp[,c(2:4)] <- NULL
  temp <- subset(temp, COG !="NA")
  colnames(temp) <- c("COG", "Total_Count")
  output <- temp
  rm(temp)
  output
}

rnaseq_cog_counts <- cog_counter(rnaseq_cog_counts)
rnaseq_cog_counts <- subset(rnaseq_cog_counts, COG !="")
rnaseq_cog_counts$COG <- forcats::fct_inorder(rnaseq_cog_counts$COG)
rm(temp)

rnaseq_cog_fisher_tests <- list()

## Calculate hypergeometric test ####
for(i in c(1:4)){ # For each differential abundance comparison...
  rnaseq_cog_fisher_tests[[i]] <- data.frame(matrix(nrow=21, ncol=4)) # Prepare a data frame to receive over-representation data.
  rnaseq_cog_fisher_tests[[i]][,1] <- as.character(rnaseq_cog_counts$COG) # Label the data frame with the COGs.
  
  for(p in c(1:3)){
    # Do one test containing all DEGs, regardless of which treatment they are upregulated in.
    if(p==1){{temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1 & is.na(COG)=="FALSE")}}
    
    # Do a second test containing all DEGs upregulated in the "left" treatment.
    if(p==2){temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & logFC < -1 & is.na(COG)=="FALSE")}
    
    # Do a third test containing all DEGs upregulated in the "right" treatment.
    if(p==3){temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & logFC > 1 & is.na(COG)=="FALSE")}
    
    totalM <- nrow(subset(rnaseq_annotations, COG != "")) # Identify the total number of transcripts with a COG class assignment.    
    totalDE <- nrow(temp) # Identify the total number of DEGs in this comparison.

    cogs <- as.factor(rnaseq_cog_counts$COG)
    
    for(j in 1:nlevels(cogs)){ # For each COG...
      total_sp <- subset(rnaseq_cog_counts, COG==levels(cogs)[j])$Total_Count # Identify the total number of transcripts assigned to that COG.
      
      total_sp_DE <- nrow(subset(temp, COG==levels(cogs)[j])) # Identify the total number of DEGs assigned to that COG.
      total_sp_DE[is.na(total_sp_DE)] <- 0 # If NA (COG not in DEG list), convert the value to zero.
      
      total_sp_notDE <- total_sp - total_sp_DE # Identify the total number of DEGs in that COG.
      total_DE_notsp <- totalDE - total_sp_DE # Identify the total number of DEGs *not* in that COG.
      total_notsp_notDE <- totalM - total_sp - total_DE_notsp # Identify the total number of genes that are *not* DE and *not* in that COG.
      
      # This produces data for a contingency table that looks like this:
      #             in COG           not in COG
      #     DE     total_sp_DE      total_DE_notsp  totalDE
      # Not DE  total_sp_notDE   total_notsp_notDE  total_notDE
      # Totals        total_sp         total_notsp  totalM
      
      # Perform Fisher's exact (hypergeometric) test.
      rnaseq_cog_fisher_tests[[i]][j,p+1] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp,
                                                                  total_sp_notDE, total_notsp_notDE), nrow=2),
                                                         alternative="greater")$p.value
    }
    rm(totalM, totalDE, total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
  }
  
  # Calculate FDR-adjusted p-values.
  rnaseq_cog_fisher_tests[[i]]$Global.Adj <- p.adjust(rnaseq_cog_fisher_tests[[i]][,2], method="BH")
  rnaseq_cog_fisher_tests[[i]]$Left.Adj <- p.adjust(rnaseq_cog_fisher_tests[[i]][,3], method="BH")
  rnaseq_cog_fisher_tests[[i]]$Right.Adj <- p.adjust(rnaseq_cog_fisher_tests[[i]][,4], method="BH")
}

names(rnaseq_cog_fisher_tests) <- names(rnaseq_diff_abund)[c(1:4)]

rm(temp)

# Rename columns based on the tested comparison.
for(i in c(1,2)){
  colnames(rnaseq_cog_fisher_tests[[i]]) <- c("Global", "AMS", "NMS", "Global-adj", "AMS-adj", "NMS-adj")
}
for(i in c(3,4)){
  colnames(rnaseq_cog_fisher_tests[[i]]) <- c("Global", "CH4", "CH3OH", "Global-adj", "CH4", "CH3OH")
}



## Create data frame counting the percentage of DEGs in each COG class ####
rnaseq_numbers <- list(rnaseq_cog_counts, rnaseq_cog_counts, rnaseq_cog_counts, rnaseq_cog_counts)

for(i in 1:4){
  temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)
  temp <- cog_counter(temp)
  temp <- subset(temp, COG != "")
  colnames(temp) <- c("COG", "DE")
  rnaseq_numbers[[i]] <- merge(rnaseq_numbers[[i]], temp, by="COG", all=TRUE)
  
  temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & logFC < -1)
  temp <- cog_counter(temp)
  temp <- subset(temp, COG != "")
  colnames(temp) <- c("COG", "Left")
  rnaseq_numbers[[i]] <- merge(rnaseq_numbers[[i]], temp, by="COG", all=TRUE)
  
  temp <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & logFC > 1)
  temp <- cog_counter(temp)
  temp <- subset(temp, COG != "")
  colnames(temp) <- c("COG", "Right")
  rnaseq_numbers[[i]] <- merge(rnaseq_numbers[[i]], temp, by="COG", all=TRUE)
}

names(rnaseq_numbers) <- tests[c(1:4)]

for(i in c(1:4)){
  if(i %in% c(1,2)){
    colnames(rnaseq_numbers[[i]]) <- c("COG", "Total_Count", "DE", "AMS", "NMS")
  }
  if(i %in% c(3,4)){
    colnames(rnaseq_numbers[[i]]) <- c("COG", "Total_Count", "DE", "CH4", "CH3OH")
  }
  rnaseq_numbers[[i]][,4] <- 100*rnaseq_numbers[[i]][,4] / rnaseq_numbers[[i]]$Total_Count
  rnaseq_numbers[[i]][,5] <- 100*rnaseq_numbers[[i]][,5] / rnaseq_numbers[[i]]$Total_Count
  rnaseq_numbers[[i]][is.na(rnaseq_numbers[[i]])] <- 0
  rnaseq_numbers[[i]]$Fill <- 25 - rnaseq_numbers[[i]][,4] - rnaseq_numbers[[i]][,5]
  rnaseq_numbers[[i]][,c("DE", "Total_Count")] <- NULL
  rnaseq_numbers[[i]] <- reshape2::melt(rnaseq_numbers[[i]])
  rnaseq_numbers[[i]] <- subset(rnaseq_numbers[[i]], !(COG %in% c("A", "B", "D", "S")))
  rnaseq_numbers[[i]]$COG <- factor(rnaseq_numbers[[i]]$COG)
  colnames(rnaseq_numbers[[i]]) <- c("COG", "Treatment", "Percentage")
}


######################################## III. DOWNSTREAM METABOLOME ANALYSIS ################
#### Differential abundance ####
metab_diff_abund <- list()

for(p in 1:4){ # For each pairwise comparison among treatments
  temp <- as.data.frame(t(metab_counts$imputed))
  temp$Carbon <- metab_sample_data$Carbon
  temp$Nitrogen <- metab_sample_data$Nitrogen
  temp$Treatment <- metab_sample_data$Treatment
  
  # Because metabolites were not detected in some treatments, some have a SD of zero.
  # To work around that for the sake of statistical testing, this code adds a "jitter"
  # around any metabolite abundances with a standard deviation of zero.
  sds <- temp[,c(1:(ncol(temp)-3),ncol(temp))] %>% 
    dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(sd=sd)) %>%
    as.data.frame()
  rownames(sds) <- sds$Treatment
  sds$Treatment <- NULL
  
  sds <- as.data.frame(t(sds))
  rownames(sds) <- metab_annotations$MetaboliteID
  for(i in c(1:2)){sds[,i] <- as.numeric(as.character(sds[,i]))}
  sds <- subset(sds, `CH4-AMS`==0 | `CH4-NMS`==0 | `MeOH-AMS`==0 | `MeOH-NMS`==0)
  for(i in c(1:nrow(sds))){
    temp[,rownames(sds)[i]] <- jitter(temp[,rownames(sds)[i]])
  }
  temp$Treatment <- NULL
  
  if(p==1){temp <- subset(temp, Carbon=="CH4")} # Subset to the relevant treatments.
  if(p==2){temp <- subset(temp, Carbon=="CH3OH")}
  if(p==3){temp <- subset(temp, Nitrogen=="AMS")}
  if(p==4){temp <- subset(temp, Nitrogen=="NMS")}
  
  if(p %in% c(1:2)){temp$Carbon <- NULL}
  if(p %in% c(3:4)){temp$Nitrogen <- NULL}
  
  diff.abund.results <- data.frame()
  
  # Perform differential abundance tests using an ANOVA and save the results.
  for(i in 1:(ncol(temp)-1)){
    test <- t.test(temp[,i] ~ temp[,ncol(temp)], var.equal=TRUE)
    diff.abund.results[i,1] <- test$statistic
    diff.abund.results[i,2] <- test[["parameter"]][[1]]
    diff.abund.results[i,3] <- test$p.value
  }
  colnames(diff.abund.results) <- c("F", "df", "p")
  rownames(diff.abund.results) <- rownames(diff.abund.results$MetaboliteID)
  
  # Calculate BH-corrected p-values.
  diff.abund.results$p.adj <- p.adjust(diff.abund.results$p, method="BH")
  diff.abund.results$log10q <- -1*log(diff.abund.results$p.adj, base=10)
  
  # Add average metabolite abundance for each treatment in the comparison.
  if(p==1){diff.abund.results <- cbind(diff.abund.results, metab_counts$avg[,c("CH4-AMS","CH4-NMS")])}
  if(p==2){diff.abund.results <- cbind(diff.abund.results, metab_counts$avg[,c("MeOH-AMS","MeOH-NMS")])}
  if(p==3){diff.abund.results <- cbind(diff.abund.results, metab_counts$avg[,c("CH4-AMS","MeOH-AMS")])}
  if(p==4){diff.abund.results <- cbind(diff.abund.results, metab_counts$avg[,c("CH4-NMS","MeOH-NMS")])}
  
  # Calculate log2FC for each comparison.
  diff.abund.results <- tibble::add_column(diff.abund.results,
                                           logFC=log(diff.abund.results[,7]/diff.abund.results[,6], base=2),
                                           .after="p.adj")
  
  # Add annotation information.
  diff.abund.results$MetaboliteID <- rownames(diff.abund.results)
  diff.abund.results <- merge(metab_annotations[,c("MetaboliteID","KEGG", "Superfamily", "Subfamily", "Biochemical")],
                              diff.abund.results,
                              by="MetaboliteID", all=TRUE)
  rownames(diff.abund.results) <- diff.abund.results$MetaboliteID
  
  metab_diff_abund[[p]] <- diff.abund.results
  rm(sds, temp, diff.abund.results)
}

#### Interaction effects ####
temp <- as.data.frame(t(metab_counts$imputed))
temp$Carbon <- metab_sample_data$Carbon
temp$Nitrogen <- metab_sample_data$Nitrogen

diff.abund.results <- data.frame()

# Perform differential abundance tests using an ANOVA and save the results.
for(i in 1:(ncol(temp)-2)){
  test <- anova(lm(temp[,i] ~ Carbon*Nitrogen, data=temp))
  diff.abund.results[i,1] <- test[3,4]
  diff.abund.results[i,2] <- test[3,1]
  diff.abund.results[i,3] <- test[3,5]
  rm(test)
}
colnames(diff.abund.results) <- c("F", "df", "p")
rownames(diff.abund.results) <- rownames(diff.abund.results$MetaboliteID)

# Calculate BH-corrected p-values.
diff.abund.results$p.adj <- p.adjust(diff.abund.results$p, method="BH")
diff.abund.results$log10q <- -1*log(diff.abund.results$p.adj, base=10)

# Add average metabolite abundance for each treatment in the comparison.
diff.abund.results <- cbind(diff.abund.results,
                            metab_counts$avg)

# Add annotation information.
diff.abund.results$MetaboliteID <- rownames(diff.abund.results)
diff.abund.results <- merge(metab_annotations[,c("MetaboliteID","KEGG", "Superfamily", "Subfamily", "Biochemical")],
                            diff.abund.results,
                            by="MetaboliteID", all=TRUE)
rownames(diff.abund.results) <- diff.abund.results$MetaboliteID

metab_diff_abund[[5]] <- diff.abund.results
rm(temp, diff.abund.results)

names(metab_diff_abund) <- tests[c(1:5)]

#### Super- and sub-family over-representation analysis ####
## Superfamily ####
metab_fisher_superfamily <- list()

for(i in c(1:4)){ # For each pairwise comparison among treatments
  # Create a data frame to store the results.
  metab_fisher_superfamily[[i]] <- data.frame(matrix(nrow=9, ncol=2))
  metab_fisher_superfamily[[i]][,1] <- as.character(levels(metab_annotations$Superfamily))
  
  for(p in c(1:3)){ # Test for (1) among *all* DAMs, (2) among DAMs enriched in the left, (3) among DAMs enriched on the right
    
    if(p==1){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)} # Subset to only DAMs.
    if(p==2){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC < -1)}
    if(p==3){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC > 1)}
    temp$Superfamily <- factor(temp$Superfamily)
    
    totalDE <- nrow(temp) # Total number of differentially abundant metabolites in this comparison.
    totalM <- 341 # Total number of metabolites in the experiment.
    
    for(j in 1:nlevels(metab_annotations$Superfamily)){ # For each superfamily
      total_sp <- nrow(subset(metab_annotations, Superfamily==levels(temp$Superfamily)[j])) # Total number in superfamily.
      total_sp_DE <- nrow(subset(temp, Superfamily==levels(temp$Superfamily)[j])) # Total DAMs in superfamily.
      total_sp_notDE <- total_sp - total_sp_DE # Total in superfamily that are not DAMs.
      total_DE_notsp <- totalDE - total_sp_DE # Total DAMs that are not in superfamily.
      total_notsp_notDE <- totalM - total_sp - total_DE_notsp # Total metabolites that are not in superfamily *and* not DAMs.
      
      # Calculate Fisher test.
      metab_fisher_superfamily[[i]][j,p+1] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE),
                                                          nrow=2),
                                                   alternative = "greater")$p.value
    }
    
    rm(temp, totalM, totalDE, total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
  }
}

# Clean data (add column names and list names)
for(i in 1:length(metab_fisher_superfamily)){
  if(i %in% c(1,2)){colnames(metab_fisher_superfamily[[i]]) <- c("superfamily", "Global", "AMS", "NMS")}
  if(i %in% c(3,4)){colnames(metab_fisher_superfamily[[i]]) <- c("superfamily", "Global", "Methane", "Methanol")}
}

names(metab_fisher_superfamily) <- tests[c(1:4)]

## Subfamily ####
metab_fisher_subfamily <- list()

for(i in c(1:4)){ # For each pairwise comparison among treatments
  metab_fisher_subfamily[[i]] <- data.frame(matrix(nrow=45, ncol=2))
  metab_fisher_subfamily[[i]][,1] <- as.character(levels(metab_annotations$Subfamily))
  
  for(p in c(1:3)){
    if(p==1){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)}
    if(p==2){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC < -1)}
    if(p==3){temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC > 1)}

    temp$Subfamily <- factor(temp$Subfamily)
    
    totalDE <- nrow(temp)
    totalM <- 341
    
    for(j in 1:nlevels(metab_annotations$Subfamily)){
      total_sp <- nrow(subset(metab_annotations, Subfamily==levels(metab_annotations$Subfamily)[j]))
      total_sp_DE <- nrow(subset(temp, Subfamily==levels(metab_annotations$Subfamily)[j]))
      total_sp_notDE <- total_sp - total_sp_DE
      total_DE_notsp <- totalDE - total_sp_DE
      total_notsp_notDE <- totalM - total_sp - total_DE_notsp
      metab_fisher_subfamily[[i]][j,1] <- levels(metab_annotations$Subfamily)[j]
      metab_fisher_subfamily[[i]][j,p+1] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2),
                                                            alternative = "greater")$p.value
    }
    rm(temp, totalM, totalDE, total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
  }
}

# Clean data
for(i in 1:length(metab_fisher_subfamily)){
  colnames(metab_fisher_subfamily[[i]]) <- c("Subfamily", "Global", "Left", "Right")
  metab_fisher_subfamily[[i]] <- merge(metab_fisher_subfamily[[i]],
                                           metab_annotations[,c("Subfamily", "Superfamily")],
                                           by="Subfamily", all.y=FALSE)
  metab_fisher_subfamily[[i]] <- dplyr::distinct(metab_fisher_subfamily[[i]])
  metab_fisher_subfamily[[i]] <- metab_fisher_subfamily[[i]][,c(5,1:4)]
  if(i %in% c(1,2)){
    colnames(metab_fisher_subfamily[[i]]) <- c("Superfamily", "Subfamily", "Global", "AMS", "NMS")
  }
  if(i %in% c(3,4)){
    colnames(metab_fisher_subfamily[[i]]) <- c("Superfamily", "Subfamily", "Global", "Methane", "Methanol")
  }
}
names(metab_fisher_subfamily) <- tests[c(1:4)]




## Create data frame counting the percentage of DAMs in each superfamily ####
metab_numbers <- list()

for(i in 1:4){
  total <- metab_annotations %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(total) <- c("Superfamily", "Total_Count")
  
  temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "DE")
  metab_numbers[[i]] <- merge(total, temp, by="Superfamily", all=TRUE)
  
  temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC < -1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "Left")
  metab_numbers[[i]] <- merge(metab_numbers[[i]], temp, by="Superfamily", all=TRUE)
  
  temp <- subset(metab_diff_abund[[i]], p.adj < 0.01 & logFC > 1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "Right")
  metab_numbers[[i]] <- merge(metab_numbers[[i]], temp, by="Superfamily", all=TRUE)
  
  rm(temp, total)
}

names(metab_numbers) <- tests[c(1:4)]

for(i in c(1:4)){
  if(i %in% c(1,2)){
    colnames(metab_numbers[[i]]) <- c("Superfamily", "Total_Count", "DE", "AMS", "NMS")
  }
  if(i %in% c(3,4)){
    colnames(metab_numbers[[i]]) <- c("Superfamily", "Total_Count", "DE", "CH4", "CH3OH")
  }
  metab_numbers[[i]][,4] <- 100*metab_numbers[[i]][,4] / metab_numbers[[i]]$Total_Count
  metab_numbers[[i]][,5] <- 100*metab_numbers[[i]][,5] / metab_numbers[[i]]$Total_Count
  metab_numbers[[i]][is.na(metab_numbers[[i]])] <- 0
  metab_numbers[[i]]$Fill <- 100 - metab_numbers[[i]][,4] - metab_numbers[[i]][,5]
  metab_numbers[[i]][,c("DE", "Total_Count")] <- NULL
  metab_numbers[[i]] <- reshape2::melt(metab_numbers[[i]])
  metab_numbers[[i]] <- subset(metab_numbers[[i]], !(Superfamily %in% c("Hormone", "Xenobiotics")))
  metab_numbers[[i]]$Superfamily <- factor(metab_numbers[[i]]$Superfamily,
                                           levels=c("Amino acid", "Carbohydrate", "Cofactors",
                                                    "Lipids", "Nucleotide", "Peptide", "Sec. met."))
  colnames(metab_numbers[[i]]) <- c("Superfamily", "Treatment", "Percentage")
}


######################################## IV. KEGG OVER-REPRESENTATION ANALYSIS #############
## Create pooled lists of DEGs and DAMs. ####

# The four data frames designed from this code are designed for downstream pathway over-representation analysis (see later)
# In addition, it creates a vector ("Query_Input") that is designed to be used directly with the KEGG search and color pathway tool.
# Website:  https://www.genome.jp/kegg/tool/map_pathway2.html # 

integr_kegg_search_color <- list()

for(i in 1:4){ # For each comparison
  temp_sub_rna <- subset(rnaseq_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1) # Identify DEGs
  temp_sub_met <- subset(metab_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1) # Identify DAMs
  
  temp_sub_rna_1 <- subset(temp_sub_rna, logFC < 0) # Left DEGs
  temp_sub_rna_2 <- subset(temp_sub_rna, logFC > 0) # Right DEGs
  temp_sub_met_1 <- subset(temp_sub_met, logFC < 0) # Left DAMs
  temp_sub_met_2 <- subset(temp_sub_met, logFC > 0) # Right DAMs
  
  temp_KO_rna_1 <- subset(rnaseq_mappings, Accession %in% rownames(temp_sub_rna_1)) # Only transcripts with a KEGG term
  temp_KO_rna_2 <- subset(rnaseq_mappings, Accession %in% rownames(temp_sub_rna_2))
  temp_KO_met_1 <- subset(metab_annotations, MetaboliteID %in% rownames(temp_sub_met_1) & is.na(KEGG)=="FALSE") # Only metabolites with a KEGG term
  temp_KO_met_2 <- subset(metab_annotations, MetaboliteID %in% rownames(temp_sub_met_2) & is.na(KEGG)=="FALSE")
  
  temp_KO_rna_1$color <- c(" red, red") # Define colors for left DEGs
  temp_KO_rna_2$color <- c(" blue, blue") # Colors for right DEGs
  temp_KO_met_1$color <- c(" red, red") # Colors for left DAMs
  temp_KO_met_2$color <- c(" blue, blue") # Colors for right DAMs
  
  temp_KO_rna_1$key <- paste0(temp_KO_rna_1$KEGG_Enzymes, temp_KO_rna_1$color) # Concatenate KEGG ID and color codes.
  temp_KO_rna_2$key <- paste0(temp_KO_rna_2$KEGG_Enzymes, temp_KO_rna_2$color)
  temp_KO_met_1$key <- paste0(temp_KO_met_1$KEGG, temp_KO_met_1$color)
  temp_KO_met_2$key <- paste0(temp_KO_met_2$KEGG, temp_KO_met_2$color)
  
  temp_KO_master <- data.frame(List = c(as.character(temp_KO_rna_1$KEGG_Enzymes), # Summarize into a data frame.
                                        as.character(temp_KO_rna_2$KEGG_Enzymes),
                                        as.character(temp_KO_met_1$KEGG),
                                        as.character(temp_KO_met_2$KEGG)),
                               Query_Input = c(temp_KO_rna_1$key,
                                               temp_KO_rna_2$key,
                                               temp_KO_met_1$key,
                                               temp_KO_met_2$key))
  
  integr_kegg_search_color[[i]] <- temp_KO_master
  
  rm(temp_KO_master, temp_sub_rna, temp_sub_met, temp_sub_rna_1, temp_sub_rna_2, temp_sub_met_1, temp_sub_met_2,
     temp_KO_rna_1, temp_KO_rna_2, temp_KO_met_1, temp_KO_met_2)
}

names(integr_kegg_search_color) <- tests[c(1:4)]

## KEGG pathway enrichment (transcripts) ####
# Calculate enrichment.
temp_models <- list()

for(i in 1:length(integr_kegg_search_color)){ # For each pairwise comparison
  temp <- integr_kegg_search_color[[i]]
  temp[] <- lapply(temp, function(x) replace(x, grepl("C", x), NA)) # Subset to only transcripts
  temp <- subset(temp, is.na(List)=="FALSE")
  temp_models[[i]] <- enrichKEGG(gene = as.character(temp$List), # Calculate over-representation
                                 organism = "ko", # Reference database
                                 universe = as.character(rnaseq_mappings$KEGG_Enzymes), # Only detected transcripts.
                                 pvalueCutoff = 0.05, pAdjustMethod = "BH")
  rm(temp)
}

names(temp_models) <- c(tests[1:4])

# Extract raw data tables from the results
rnaseq_kegg_enrichment <- list()
for(i in 1:length(temp_models)){
  df <- temp_models[[i]]@result
  
  # Separate data to calculate an enrichment factor.
  df <- df %>% tidyr::separate(col=GeneRatio, into=c("Hits","Queries"), sep="/")
  df <- df %>% tidyr::separate(col=BgRatio, into=c("Local_Total", "Universe"), sep="/")
  
  df$Hits <- as.numeric(df$Hits)
  df$Queries <- as.numeric(df$Queries)
  df$Local_Total <- as.numeric(df$Local_Total)
  df$Universe <- as.numeric(df$Universe)
  
  df$RichFactor <- df$Hits / df$Local_Total
  
  rnaseq_kegg_enrichment[[i]] <- df
  rm(df)
}

names(rnaseq_kegg_enrichment) <- c(tests[1:4])
rm(temp_models)

## KEGG pathway enrichment (metabolites) ####
# For this analysis, we used the 'pathway analysis' function in the online MetaboAnalyst interface.
# Link: https://www.metaboanalyst.ca/MetaboAnalyst/home.xhtml
# The metabolites were input from the 'integr_kegg_search_color' spreadsheets designed above.
######################################## V. SUPPLEMENT: INTERNAL RELATIVIZATION ANALYSIS ################
#### Prepare relativized data ####
norm_data <- read.csv("raw_data/metab_counts.csv")

rownames(norm_data) <- norm_data$MetaboliteID
norm_data$MetaboliteID <- NULL

# Impute missing data points.
norm_data <- as.data.frame(norm_data)
for(i in 1:nrow(norm_data)){
  rowmin <- 0.5*min(norm_data[i,], na.rm=TRUE)
  for(j in 1:ncol(norm_data)){
    if(is.na(norm_data[i,j])=="TRUE"){
      norm_data[i,j] <- rowmin
    }
  }
  rm(rowmin)
}

# Add annotations.
norm_data <- merge(metab_annotations[,c("MetaboliteID", "Superfamily", "Subfamily")],
                   norm_data,
                   by=0, all.y=TRUE)
rownames(norm_data) <- norm_data$MetaboliteID
norm_data[,c("Row.names", "MetaboliteID")] <- NULL

norm_data_rel <- data.frame(matrix(ncol=ncol(norm_data)))
colnames(norm_data_rel) <- colnames(norm_data)

# Calculate relativized internal abundances within each "superfamily."
# e.g., the abundance of a metabolite is expressed as a fraction of the total peak area of each superfamily.
for(i in 1:nlevels(norm_data$Superfamily)){
  temp <- subset(norm_data, Superfamily == levels(norm_data$Superfamily)[i])
  temp[,c(3:ncol(temp))] <- 100*t(decostand(t(temp[,c(3:ncol(temp))]), method = "total"))
  norm_data_rel <- rbind(norm_data_rel, temp)
  rm(temp)
}
norm_data_rel <- subset(norm_data_rel, is.na(Superfamily)=="FALSE" &
                          Superfamily !="Hormone" & Superfamily !="Xenobiotics")
norm_data_rel[,c(1:2)] <- NULL

# Calculate mean relative abundances of metabolites in each treatment.
temp <- as.data.frame(t(norm_data_rel))
temp$Treatment <-  metab_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(mean=mean))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(norm_data_rel)
temp$MetaboliteID <- rownames(temp)
norm_data_avg <- temp
rm(temp)

# Calculate standard deviations of metabolite abundances per treatment.
temp <- as.data.frame(t(norm_data_rel))
temp$Treatment <-  metab_sample_data$Treatment
temp <- temp %>% dplyr::group_by(Treatment) %>% dplyr::summarise_all(.funs=list(sd=sd))
temp$Treatment <- NULL
temp <- as.data.frame(t(temp))
colnames(temp) <- treatment_order
rownames(temp) <- rownames(norm_data_rel)
temp$MetaboliteID <- rownames(temp)
norm_data_sd <- temp
rm(temp, i, j)

#### Calculate differential abundance ####
norm_diff_abund <- list()

for(p in 1:4){ # For each pairwise comparison among treatments
  temp <- as.data.frame(t(norm_data_rel))
  temp$Carbon <- metab_sample_data$Carbon
  temp$Nitrogen <- metab_sample_data$Nitrogen

  if(p==1){temp <- subset(temp, Carbon=="CH4")} # Subset to the relevant treatments.
  if(p==2){temp <- subset(temp, Carbon=="CH3OH")}
  if(p==3){temp <- subset(temp, Nitrogen=="AMS")}
  if(p==4){temp <- subset(temp, Nitrogen=="NMS")}
  
  if(p %in% c(1:2)){temp$Carbon <- NULL}
  if(p %in% c(3:4)){temp$Nitrogen <- NULL}
  
  if(p %in% c(1,2)){temp$Carbon <- NULL}
  if(p %in% c(3,4)){temp$Nitrogen <- NULL}
  
  diff.abund.results <- data.frame()
  
  # Perform differential abundance tests using an ANOVA and save the results.
  for(i in 1:(ncol(temp)-1)){
    test <- t.test(temp[,i] ~ temp[,ncol(temp)], var.equal=TRUE)
    diff.abund.results[i,1] <- test$statistic
    diff.abund.results[i,2] <- test[["parameter"]][[1]]
    diff.abund.results[i,3] <- test$p.value
  }
  colnames(diff.abund.results) <- c("F", "df", "p")
  rownames(diff.abund.results) <- rownames(diff.abund.results$MetaboliteID)
  
  # Calculate BH-corrected p-values.
  diff.abund.results$p.adj <- p.adjust(diff.abund.results$p, method="BH")
  diff.abund.results$log10q <- -1*log(diff.abund.results$p.adj, base=10)
  
  # Add average metabolite abundance for each treatment in the comparison.
  if(p==1){diff.abund.results <- cbind(diff.abund.results, norm_data_avg[,c("CH4-AMS","CH4-NMS")])}
  if(p==2){diff.abund.results <- cbind(diff.abund.results, norm_data_avg[,c("MeOH-AMS","MeOH-NMS")])}
  if(p==3){diff.abund.results <- cbind(diff.abund.results, norm_data_avg[,c("CH4-AMS","MeOH-AMS")])}
  if(p==4){diff.abund.results <- cbind(diff.abund.results, norm_data_avg[,c("CH4-NMS","MeOH-NMS")])}
  
  # Calculate log2FC for each comparison.
  diff.abund.results <- tibble::add_column(diff.abund.results,
                                           logFC=log(diff.abund.results[,7]/diff.abund.results[,6], base=2),
                                           .after="p.adj")
  
  # Add annotation information.
  diff.abund.results$MetaboliteID <- rownames(diff.abund.results)
  diff.abund.results <- merge(metab_annotations[,c("MetaboliteID","KEGG", "Superfamily", "Subfamily", "Biochemical")],
                              diff.abund.results,
                              by="MetaboliteID", all=TRUE)
  rownames(diff.abund.results) <- diff.abund.results$MetaboliteID
  
  norm_diff_abund[[p]] <- diff.abund.results
  rm(temp, diff.abund.results)
}

names(norm_diff_abund) <- tests[c(1:4)]

#### Count the number of DAMs ####
norm_numbers <- list()

for(i in 1:4){
  total <- metab_annotations %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(total) <- c("Superfamily", "Total_Count")
  
  temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "DE")
  norm_numbers[[i]] <- merge(total, temp, by="Superfamily", all=TRUE)
  
  temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC < -1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "Left")
  norm_numbers[[i]] <- merge(norm_numbers[[i]], temp, by="Superfamily", all=TRUE)
  
  temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC > 1)
  temp <- temp %>% dplyr::group_by(Superfamily) %>% dplyr::count()
  colnames(temp) <- c("Superfamily", "Right")
  norm_numbers[[i]] <- merge(norm_numbers[[i]], temp, by="Superfamily", all=TRUE)
  
  rm(temp, total)
}

names(norm_numbers) <- tests[c(1:4)]

for(i in c(1:4)){
  if(i %in% c(1,2)){
    colnames(norm_numbers[[i]]) <- c("Superfamily", "Total_Count", "DE", "AMS", "NMS")
  }
  if(i %in% c(3,4)){
    colnames(norm_numbers[[i]]) <- c("Superfamily", "Total_Count", "DE", "CH4", "CH3OH")
  }
  norm_numbers[[i]][,4] <- 100*norm_numbers[[i]][,4] / norm_numbers[[i]]$Total_Count
  norm_numbers[[i]][,5] <- 100*norm_numbers[[i]][,5] / norm_numbers[[i]]$Total_Count
  norm_numbers[[i]][is.na(norm_numbers[[i]])] <- 0
  norm_numbers[[i]]$Fill <- 100 - norm_numbers[[i]][,4] - norm_numbers[[i]][,5]
  norm_numbers[[i]][,c("DE", "Total_Count")] <- NULL
  norm_numbers[[i]] <- reshape2::melt(norm_numbers[[i]])
  norm_numbers[[i]] <- subset(norm_numbers[[i]], !(Superfamily %in% c("Hormone", "Xenobiotics")))
  norm_numbers[[i]]$Superfamily <- factor(norm_numbers[[i]]$Superfamily,
                                           levels=c("Amino acid", "Carbohydrate", "Cofactors",
                                                    "Lipids", "Nucleotide", "Peptide", "Sec. met."))
  colnames(norm_numbers[[i]]) <- c("Superfamily", "Treatment", "Percentage")
}

#### Super- and sub-family over-representation analysis ####
## Superfamily ####
norm_fisher_superfamily <- list()

for(i in c(1:4)){ # For each pairwise comparison among treatments
  # Create a data frame to store the results.
  norm_fisher_superfamily[[i]] <- data.frame(matrix(nrow=9, ncol=2))
  norm_fisher_superfamily[[i]][,1] <- as.character(levels(metab_annotations$Superfamily))
  
  for(p in c(1:3)){ # Test for (1) among *all* DAMs, (2) among DAMs enriched in the left, (3) among DAMs enriched on the right
    
    if(p==1){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)} # Subset to only DAMs.
    if(p==2){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC < -1)}
    if(p==3){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC > 1)}
    temp$Superfamily <- factor(temp$Superfamily)
    
    totalDE <- nrow(temp) # Total number of differentially abundant metabolites in this comparison.
    totalM <- 341 # Total number of metabolites in the experiment.
    
    for(j in 1:nlevels(metab_annotations$Superfamily)){ # For each superfamily
      total_sp <- nrow(subset(metab_annotations, Superfamily==levels(temp$Superfamily)[j])) # Total number in superfamily.
      total_sp_DE <- nrow(subset(temp, Superfamily==levels(temp$Superfamily)[j])) # Total DAMs in superfamily.
      total_sp_notDE <- total_sp - total_sp_DE # Total in superfamily that are not DAMs.
      total_DE_notsp <- totalDE - total_sp_DE # Total DAMs that are not in superfamily.
      total_notsp_notDE <- totalM - total_sp - total_DE_notsp # Total metabolites that are not in superfamily *and* not DAMs.
      
      # Calculate Fisher test.
      norm_fisher_superfamily[[i]][j,p+1] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE),
                                                                 nrow=2),
                                                          alternative = "greater")$p.value
    }
    
    rm(temp, totalM, totalDE, total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
  }
}

# Clean data (add column names and list names)
for(i in 1:length(norm_fisher_superfamily)){
  if(i %in% c(1,2)){colnames(norm_fisher_superfamily[[i]]) <- c("superfamily", "Global", "AMS", "NMS")}
  if(i %in% c(3,4)){colnames(norm_fisher_superfamily[[i]]) <- c("superfamily", "Global", "Methane", "Methanol")}
}

names(norm_fisher_superfamily) <- tests[c(1:4)]

## Subfamily ####
norm_fisher_subfamily <- list()

for(i in c(1:4)){ # For each pairwise comparison among treatments
  norm_fisher_subfamily[[i]] <- data.frame(matrix(nrow=45, ncol=2))
  norm_fisher_subfamily[[i]][,1] <- as.character(levels(metab_annotations$Subfamily))
  
  for(p in c(1:3)){
    if(p==1){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & abs(logFC) > 1)}
    if(p==2){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC < -1)}
    if(p==3){temp <- subset(norm_diff_abund[[i]], p.adj < 0.01 & logFC > 1)}
    
    temp$Subfamily <- factor(temp$Subfamily)
    
    totalDE <- nrow(temp)
    totalM <- 341
    
    for(j in 1:nlevels(metab_annotations$Subfamily)){
      total_sp <- nrow(subset(metab_annotations, Subfamily==levels(metab_annotations$Subfamily)[j]))
      total_sp_DE <- nrow(subset(temp, Subfamily==levels(metab_annotations$Subfamily)[j]))
      total_sp_notDE <- total_sp - total_sp_DE
      total_DE_notsp <- totalDE - total_sp_DE
      total_notsp_notDE <- totalM - total_sp - total_DE_notsp
      norm_fisher_subfamily[[i]][j,1] <- levels(metab_annotations$Subfamily)[j]
      norm_fisher_subfamily[[i]][j,p+1] <- fisher.test(matrix(c(total_sp_DE, total_DE_notsp, total_sp_notDE, total_notsp_notDE), nrow=2),
                                                        alternative = "greater")$p.value
    }
    rm(temp, totalM, totalDE, total_sp, total_sp_DE, total_sp_notDE, total_DE_notsp, total_notsp_notDE)
  }
}

# Clean data
for(i in 1:length(norm_fisher_subfamily)){
  colnames(norm_fisher_subfamily[[i]]) <- c("Subfamily", "Global", "Left", "Right")
  norm_fisher_subfamily[[i]] <- merge(norm_fisher_subfamily[[i]],
                                       metab_annotations[,c("Subfamily", "Superfamily")],
                                       by="Subfamily", all.y=FALSE)
  norm_fisher_subfamily[[i]] <- dplyr::distinct(norm_fisher_subfamily[[i]])
  norm_fisher_subfamily[[i]] <- norm_fisher_subfamily[[i]][,c(5,1:4)]
  if(i %in% c(1,2)){
    colnames(norm_fisher_subfamily[[i]]) <- c("Superfamily", "Subfamily", "Global", "AMS", "NMS")
  }
  if(i %in% c(3,4)){
    colnames(norm_fisher_subfamily[[i]]) <- c("Superfamily", "Subfamily", "Global", "Methane", "Methanol")
  }
}
names(norm_fisher_subfamily) <- tests[c(1:4)]




