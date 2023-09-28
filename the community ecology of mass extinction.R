rm(list=ls())
setwd("~/R/MRes summer project")

{library(dplyr)
  library(tidyr)
  library(ape)
  library(vegan)
  library(pez)
  library(ggplot2)
  library(tibble)
  library(readxl)
  library(FD)
  library(readxl)
  library(segmented)
  library(PDtoolkit)
  library(zoo)
  library(ggpubr)
  library(stat)
}

############# Import species abundance, functional and taxonomy data ###########

data = read.csv('Danise_etal_2013_dataset.csv')
taxonomy = read.csv('Danise_etal_2013_taxonomy.csv')
taxonomy_changed = read.csv('Danise_etal_2013_taxonomy_changed.csv')
trait1 = read.csv('Lauren_James_functional_traits.csv',stringsAsFactors=TRUE)
size_pseudomy = read.csv('size_pseudomy.csv')
size_palaeo = read.csv('size_palaeo.csv')
colnames(size_pseudomy) <- colnames(size_palaeo)
size_data <- rbind(size_pseudomy,size_palaeo)

# Transpose the data frame and rename columns
rownames(data) <- (data[, 1])
data <- data[, -1]
df <- data.frame(t(data))
colnames(df) <- gsub("(?<=\\w)\\.(?=\\w)", " ", colnames(df), perl = TRUE)
names(df)[61] <- "Strat_Height"

# Calculate the species richness
df$species_richness <- apply(df[, 1:56], 1, function(x) sum(x != 0))

# Calculate the Simpson's diversity
df[, 1:56] <- sapply(df[, 1:56], as.numeric)
df$simpson <- diversity(df[, 1:56], index = 'simpson')


######################## Calculate phylogenetic diversity ######################

# Make the phylogeny from taxonomy
taxonomy <- as.data.frame(lapply(taxonomy, as.factor))
frm <- ~Phylum/Class/Order/Superfamily/Family/Genus/Species
phylogeny <- as.phylo(frm, data=taxonomy, collapse = FALSE)
plot.phylo(phylogeny)

# Make the revised phylogeny
taxonomy_changed <- as.data.frame(lapply(taxonomy_changed, as.factor))
phylogeny_changed <- as.phylo(frm, data=taxonomy_changed, collapse = FALSE)

# Add edge length for each phylogeny
phylogeny_el<-phylogeny
phylogeny_el$edge.length <- branching.times(phylogeny_el)
phylogeny_el$edge.length <- rep(1, length(phylogeny_el$edge.length))

phylogeny_changed_el <- phylogeny_changed
phylogeny_changed_el$edge.length <- branching.times(phylogeny_changed_el)
phylogeny_changed_el$edge.length <- rep(1, length(phylogeny_changed_el$edge.length))

# Create comparative data for original phylogeny
c.data <- comparative.comm(multi2di(phylogeny), df_matrix)
c.data_el <- comparative.comm(multi2di(phylogeny_el), df_matrix)
c.data_el[["phy"]][["edge.length"]] <-  rep(1, length(c.data_el[["phy"]][["edge.length"]]))

# Create comparative data for revised taxonomy
c.data_changed <- comparative.comm(multi2di(phylogeny_changed), df_matrix)
c.data_changed_el <- comparative.comm(multi2di(phylogeny_changed_el), df_matrix)
c.data_changed_el[["phy"]][["edge.length"]] <-  rep(1, length(c.data_changed_el[["phy"]][["edge.length"]]))

# Calculate Faith's pd
faiths <- .pd(c.data_el,abundance.weighted = TRUE) 
faiths_changed<-.pd(c.data_changed_el,abundance.weighted = TRUE)

# Combine results into a single data frame
df <- tibble::rownames_to_column(df, "Site")
site_names <- rownames(df_matrix)
pd_results <- data.frame(Site = site_names, faiths)
pd_results <- pd_results %>% left_join(MNTD,by = c("Site"))
pd_results <- pd_results %>% left_join(df[ , c('Site',"AGES", "Strat_Height")], by = c('Site'))

pd_results$AGES <- as.numeric(pd_results$AGES)
pd_results$Strat_Height <- as.numeric(pd_results$Strat_Height)
pd_results$AGES_r <- -pd_results$AGES

# Combine results into a single data frame for revised taxonomy
pd_changed <- data.frame(Site = site_names, faiths_changed)
pd_changed <- pd_changed %>% left_join(MNTD_changed,by = c("Site"))
pd_changed <- pd_changed %>% left_join(df[ , c('Site',"AGES", "Strat_Height")], by = c('Site'))

pd_changed$AGES <- as.numeric(pd_changed$AGES)
pd_changed$Strat_Height <- as.numeric(pd_changed$Strat_Height)
pd_changed$AGES_r <- -pd_changed$AGES


####################### Calculate functional trait diversity####################

# Match ecological traits data with abundance data 
df_fd <- df[, -c(57:65)]
df_matrix<-as.matrix(df_fd)

rownames(trait1) <- trait1[, 1]
trait1 <- trait1[, -1]
common_species <- intersect(rownames(trait1), colnames(df_matrix))
trait1 <- trait1[common_species, ]
df_matrix <- df_matrix[, common_species]

# Calculate ecological trait FDis
FDis_eco<-dbFD(trait1,df_matrix)

# Combine with stratigraphic height and ages
FDis_trait<- as.data.frame(FDis_eco[["FDis"]])
FDis_trait <- tibble::rownames_to_column(FDis_trait, "Site")
FDis_trait <- FDis_trait %>% left_join(df[ , c('Site',"AGES", "Strat_Height")], by = c('Site'))
colnames(FDis_trait)[2]<- 'FDis_eco'
FDis_trait$AGES <- as.numeric(FDis_trait$AGES)
FDis_trait$Strat_Height <- as.numeric(FDis_trait$Strat_Height)
FDis_trait$AGES_r <- -FDis_trait$AGES


################################### Size analysis ##############################

# Import the size data, read all sheets into a list
species_names <- excel_sheets("size_Seb Newell Body Size Data Final_2021.xlsx") 
all_sheets <- lapply(species_names, function(sheet) {
  read_excel("size_Seb Newell Body Size Data Final_2021.xlsx", sheet = sheet)
})

size_data <- bind_rows(all_sheets)
size_data$`Strat. Height` <- as.numeric(size_data$`Strat. Height`)
colnames(size_data)[5] <- 'Strat_Height'

# Modify abundance data for size analysis
df_size <- data.frame(t(data))
df_size <- df_size[,-c(57:60,62:63)]
transformed_abundance <- df_size %>% gather(Species, Abundance, -STRATIGRAPHIC.HEIGHT..m.) 
colnames(transformed_abundance)[1] <- 'Strat_Height'
transformed_abundance$Strat_Height <- as.numeric(transformed_abundance$Strat_Height)
transformed_abundance$Abundance <- as.numeric(transformed_abundance$Abundance)

######### Calculate community weighted mean for geometric mean #################

size_joined <- transformed_abundance %>%
  left_join(size_data, by = c("Species", "Strat_Height"), multiple = "all")
size_joined$Abundance <- as.numeric(size_joined$Abundance)
size_joined$Geometric.Mean <- as.numeric(size_joined$Geometric.Mean)

community_size <- size_joined %>% group_by(Strat_Height) %>%
  summarize(CWM_Size = weighted.mean(Geometric.Mean,
                                     w = Abundance / sum(Abundance),na.rm=TRUE))

# Interpolate missing data using na.approx
interpolated_community_size <- as.data.frame(na.approx(community_size))
colnames(interpolated_community_size)[2] <- 'interpolated_CWM'

#################### Estimate missing data using species mean ##################

species_mean_size <- size_data %>% group_by(Species) %>%
  summarize(species_mean_size = mean(Geometric.Mean))

# Join abundance data with species mean size data
size_estimate_data <- transformed_abundance %>%
  left_join(species_mean_size, by = "Species")

# Calculate the interpolated community weighted mean size
estimated_community_data <- size_estimate_data %>%
  group_by(Strat_Height) %>%
  summarize(estimated_CWM = weighted.mean(species_mean_size, w = Abundance / sum (Abundance), na.rm = TRUE))
combined_estimated <- community_size %>%
  left_join(estimated_community_data, by = "Strat_Height")
combined_estimated$CWM_estimated <- ifelse(is.na(combined_estimated$CWM_Size), 
                                           combined_estimated$estimated_CWM, combined_estimated$CWM_Size)
combined_estimated<-combined_estimated[, -c(2:3)]

# Use na.approx to estimate the rest missing site
full_interpolated <- as.data.frame(na.approx(combined_estimated))
colnames(full_interpolated)[4] <- 'full_interpolated'
full_interpolated<-full_interpolated[, -c(2:3)]


########################## Time series analysis ################################

#Integrate data for analysis 
{data_for_analysis <- df[ , c("AGES", "Strat_Height","species_richness", "simpson")] 
data_for_analysis$AGES <- as.numeric(data_for_analysis$AGES)
data_for_analysis$AGES_r <- -data_for_analysis$AGES
data_for_analysis$Strat_Height <- as.numeric(data_for_analysis$Strat_Height)}

{data_for_analysis <- data_for_analysis %>%
    left_join(community_size,by = c("Strat_Height")) 
  data_for_analysis <- data_for_analysis %>%
    left_join(FDis_size,by = c("Strat_Height"))
  data_for_analysis <- data_for_analysis %>%
    left_join(interpolated_community_size,by = c("Strat_Height")) 
  data_for_analysis <- data_for_analysis %>%
    left_join(community_size_without,by = c("Strat_Height")) 
  data_for_analysis <- data_for_analysis %>%
    left_join(combined_estimated,by = c("Strat_Height")) 
  data_for_analysis <- data_for_analysis %>%
    left_join(full_interpolated,by = c("Strat_Height")) }

data_for_analysis <- data_for_analysis %>%
  left_join(pd_results,by = c("Strat_Height")) 
data_for_analysis <- data_for_analysis %>%
  left_join(FDis_trait,by = c("Strat_Height")) 


################################ PCA #####################################

{metric1 <- c('species_richness','simpson','CWM_Size')
metrics <- data_for_analysis[metric1]
metrics$pd_revised <- pd_changed$pd
metrics$pd.ivs <- pd_changed$pd.ivs
metrics$FDis <- FDis_trait$FDis_eco}
metrics_cleaned<-na.omit(metrics)

metrics_pca<-prcomp(metrics_cleaned,scale=TRUE)
summary(metrics_pca)
biplot(metrics_pca,xlabs = rep("x",94),main='Biodiversity metrics PCA')

############################ Models & analysis for size ########################

# Model 1.1: CWM size; 94 sites; 2 break points
model1.1. <- segmented(lm(log(CWM_Size)~ AGES_r,data = data_for_analysis), npsi=2,seg.Z=~AGES_r)
summary(model1.1)

plot(model1.1) # Plot model result and shade the extinction interval
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

data_for_analysis$Segment1.1 <- ifelse(data_for_analysis$AGES_r < -183.690, "Before",
                                         ifelse(data_for_analysis$AGES_r < -183.466, "Between", "After"))
aov1.1 <- aov(log(CWM_Size) ~ Segment1.1, data = data_for_analysis)
summary(aov1.1)
TukeyHSD(aov1.1)

# Model 1.2: CWM size; 144 sites; 2 break points
model1.2 <- segmented(lm(log(interpolated_CWM)~ AGES_r,data = data_for_analysis), seg.Z=~AGES_r,npsi=2)
summary(model1.2)

plot(model1.2)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))
data_for_analysis$Segment1.2 <- ifelse(data_for_analysis$AGES_r < -182.868, "Before",
                                         ifelse(data_for_analysis$AGES_r < -182.725, "Between", "After"))
aov1.2 <- aov(log(interpolated_CWM) ~ Segment1.2, data = data_for_analysis)
summary(aov1.2)
TukeyHSD(aov1.2)


# Model 1.3: CWM size; 137 sites
model1.3 <- segmented(lm(log(CWM_estimated)~ AGES_r,data = data_for_analysis_complete), seg.Z=~AGES_r,npsi=2)
summary(model1.3)
plot(model1.3)
data_for_analysis$Segment1.3 <- ifelse(data_for_analysis$AGES_r < -181.615, "Before",
                                         ifelse(data_for_analysis$AGES_r < -181.243, "Between", "After"))
aov1.3 <- aov(log(CWM_estimated) ~ Segment1.3, data = data_for_analysis)
summary(aov1.3)
TukeyHSD(aov1.3)


# Model 1.4: CWM size; 153 sites; 2 break points
model1.4 <- lm(log(full_interpolated)~ AGES_r,data = data_for_analysis_complete)
summary(model1.4) # Failed to detect breakpoint, results output are simple linear model results
plot(model1.4)


# Model 2.1: Functional traits; 2 break points
model2.1 <- segmented(lm(sqrt(FDis_eco)~ AGES_r,data = FDis_trait), seg.Z=~AGES_r,npsi=2)
summary(model2.1)

plot(model2.1)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

FDis_trait$Segment2.1 <- ifelse(FDis_trait$AGES_r < -182.926, "Before",
                                  ifelse(FDis_trait$AGES_r < -182.697, "Between", "After"))
aov2.1 <- aov(sqrt(FDis_eco) ~ Segment2.1, data = FDis_trait)
summary(aov2.1)
TukeyHSD(aov2.1)

# Model 3.1: Phylogenetic diversity; original taxonomy; 2 break points
model3.1 <- segmented(lm(log(pd)~ AGES_r,data = pd_results), seg.Z=~AGES_r,npsi=2)
summary(model3.1)

plot(model3.1)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

pd_results$Segment3.1 <- ifelse(pd_results$AGES_r < -182.439, "Before",
                                  ifelse(pd_results$AGES_r < -181.554, "Between", "After"))
aov3.1 <- aov(log(pd) ~ Segment3.1, data = pd_results)
summary(aov3.1)
TukeyHSD(aov3.1)


# Model 3.2: Phylogenetic diversity; changed taxonomy; 2 break points
model3.2 <- segmented(lm(log(pd)~ AGES_r,data = pd_changed),seg.Z=~AGES_r, npsi=2)
summary(model3.2)

plot(model3.2)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

pd_changed$Segment3.2 <- ifelse(pd_changed$AGES_r < -182.427, "Before",
                                  ifelse(pd_changed$AGES_r < -181.594, "Between", "After"))
aov3.2 <- aov(log(pd) ~ Segment3.2, data = pd_changed)
summary(aov3.2)
TukeyHSD(aov3.2)


# Model 4: Simpson index and species richness
model4.1 <- segmented(lm(log(species_richness)~ AGES_r,data = data_for_analysis), seg.Z=~AGES_r,npsi=2)
summary(model4.1)

plot(model4.1)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

data_for_analysis$Segment4.1 <- ifelse(data_for_analysis$AGES_r < -182.970, "Before",
                                         ifelse(data_for_analysis$AGES_r < -182.713, "Between", "After"))
aov4.1 <- aov(log(species_richness) ~ Segment4.1, data = data_for_analysis)
summary(aov4.1)
TukeyHSD(aov4.1)

model4.2 <- segmented(lm(simpson~ AGES_r,data = data_for_analysis),seg.Z=~AGES_r,npsi=2)
summary(model4.2)
plot(model4.2)
data_for_analysis$Segment4.2 <- ifelse(data_for_analysis$AGES_r < -182.982, "Before",
                                         ifelse(data_for_analysis$AGES_r < -182.709, "Between", "After"))
aov4.2.1 <- aov(simpson ~ Segment4.2, data = data_for_analysis)
summary(aov4.2)
TukeyHSD(aov4.2)


# Model_all: Everything
{metric1 <- c('species_richness','simpson','CWM_Size','interpolated_CWM','CWM_estimated','full_interpolated','AGES_r')
  metrics_all <- data_for_analysis[metric1]
  metrics_all$pd_revised <- pd_changed$pd
  metrics_all$FDis <- FDis_trait$FDis_eco}
metrics_all_cleaned<-na.omit(metrics_all)

# Use CWM from 94 horizons
model_all1 <- segmented(lm(log(CWM_Size)+ log(species_richness)+simpson
                             + sqrt(FDis)+ log(pd_revised) ~ AGES_r, data=metrics_all),npsi=2,seg.Z=~AGES_r)
summary(model_all1)

plot(model_all1)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

metrics_all1$Segment_all1 <- ifelse(metrics_all$AGES_r < -183.533, "Before",
                                     ifelse(metrics_all$AGES_r < -182.035, "Between", "After"))
aov_all1 <- aov(log(full_interpolated)+ log(species_richness)+simpson
                  + sqrt(FDis)+ log(pd_revised) ~ Segment_all1, data = metrics_all)
summary(aov_all1)
TukeyHSD(aov_all1)

# Use CWM from 144 horizons
model_all2 <- segmented(lm(log(interpolated_CWM)+ log(species_richness)+simpson
                           + sqrt(FDis)+ log(pd_revised) ~ AGES_r, data=metrics_all),npsi=2,seg.Z=~AGES_r)
summary(model_all2)
plot(model_all2)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

metrics_all2$Segment_all2 <- ifelse(metrics_all$AGES_r < -183.533, "Before",
                                    ifelse(metrics_all$AGES_r < -182.035, "Between", "After"))
aov_all2 <- aov(log(full_interpolated)+ log(species_richness)+simpson
                + sqrt(FDis)+ log(pd_revised) ~ Segment_all2, data = metrics_all)
summary(aov_all2)
TukeyHSD(aov_all2)

# Use CWM from 137 horizons
model_all3 <- segmented(lm(log(CWM_estimated)+ log(species_richness)+simpson
                           + sqrt(FDis)+ log(pd_revised) ~ AGES_r, data=metrics_all),npsi=2,seg.Z=~AGES_r)
summary(model_all3)

plot(model_all3)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

metrics_all3$Segment_all3 <- ifelse(metrics_all$AGES_r < -183.533, "Before",
                                    ifelse(metrics_all$AGES_r < -182.035, "Between", "After"))
aov_all3 <- aov(log(full_interpolated)+ log(species_richness)+simpson
                + sqrt(FDis)+ log(pd_revised) ~ Segment_all3, data = metrics_all)
summary(aov_all3)
TukeyHSD(aov_all3)

# Use CWM from 153 horizons
model_all4 <- segmented(lm(log(full_interpolated)+ log(species_richness)+simpson
                           + sqrt(FDis)+ log(pd_revised) ~ AGES_r, data=metrics_all),npsi=2,seg.Z=~AGES_r)
summary(model_all4)

plot(model_all4)
rect(-182.7,-20, -183.3, 20,col=rgb(red = 0, green = 0.1, blue = 1, alpha = 0.2))

metrics_all4$Segment_all4 <- ifelse(metrics_all$AGES_r < -183.533, "Before",
                                    ifelse(metrics_all$AGES_r < -182.035, "Between", "After"))
aov_all4 <- aov(log(full_interpolated)+ log(species_richness)+simpson
                + sqrt(FDis)+ log(pd_revised) ~ Segment_all4, data = metrics_all)
summary(aov_all4)
TukeyHSD(aov_all4)


################################# Plot #########################################

# Create the plot for phylogenetic diversity
plot_changed <-  ggplot(pd_changed, aes(x = AGES_r, y = pd)) +
  geom_line() +
  geom_ribbon(
    aes(
      xmin = -182.7,
      xmax = -183.3,
    ),
    fill = 'brown1',
    alpha = 0.2
  )+
  geom_vline(xintercept = -182.524,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  geom_vline(xintercept = -181.505,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  # Customize labels and titles
  labs(
    x = "Time",
    y = "Phylogenetic diversity (revised)",
  ) +
  coord_flip() +
  theme_pubr()
# Print the plot
print(plot_changed)

# Create the plot for species richness
plot_sr <-  ggplot(data_for_analysis, aes(x = AGES_r, y = species_richness)) +
  geom_line() +
  geom_ribbon(
    aes(
      xmin = -182.7,
      xmax = -183.3,
    ),
    fill = 'brown1',
    alpha = 0.2
  )+
  geom_vline(xintercept = -182.970,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  geom_vline(xintercept = -182.713,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  # Customize labels and titles
  labs(
    x = "Time",
    y = "Species richness",
  ) +
  coord_flip() +
  theme_pubr()
# Print the plot
print(plot_sr)

# Create the plot for simpson index
plot_simpson <-  ggplot(data_for_analysis, aes(x = AGES_r, y = simpson)) +
  geom_line() +
  geom_ribbon(
    aes(
      xmin = -182.7,
      xmax = -183.3,
    ),
    fill = 'brown1',
    alpha = 0.2
  )+
  geom_vline(xintercept = -182.982,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  geom_vline(xintercept = -182.709,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  # Customize labels and titles
  labs(
    x = "Time",
    y = "Simpson index of diversity",
  ) +
  coord_flip() +
  theme_pubr()
# Print the plot
print(plot_simpson)


# Create the plot for functional dispersion
plot_fd <-  ggplot(FDis_trait, aes(x = AGES_r, y = FDis_eco)) +
  geom_line() +
  geom_ribbon(
    aes(
      xmin = -182.7,
      xmax = -183.3,
    ),
    fill = 'brown1',
    alpha = 0.2
  )+
  geom_vline(xintercept = -182.926,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  geom_vline(xintercept = -182.697,color='dodgerblue2', linetype = 'dashed',
             linewidth = 1)+
  # Customize labels and titles
  labs(
    x = "Time",
    y = "Functional dispersion",
  ) +
  coord_flip() +
  theme_pubr()
# Print the plot
print(plot_fd)


