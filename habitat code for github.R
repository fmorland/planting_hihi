## habitat analysis 
library(readr)
library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(reshape2)
library(ggplot2)
library(ggsignif)
library(viridis)

#### import & sort out data ####

data <- read.csv(habitat_code1_data.csv)

breedingdata1 <- read_csv("breedingdata1.csv", 
                          col_types = cols(...1 = col_skip()))

planting_data <- read_csv("Planting data6.csv", 
                          col_types = cols(Year = col_integer(), 
                                           Notes = col_skip()))

planting_data$Number_planted<-as.numeric(planting_data$Number_planted)

# summarise planting data ####

planting_summary<-planting_data %>%
  group_by(Bush,Species) %>%
  summarise(Total=sum(Number_planted))

# calculate species richness of planted trees
species_richness<-planting_data %>%
  group_by(Bush) %>%
  summarise(species_richness=n_distinct(Species))

# change species names to lower case
planting_summary$Species<-tolower(planting_summary$Species)


species_totals<-planting_data %>%
  group_by(Species) %>%
  summarise(Total=sum(na.omit(Number_planted)))


planting_summary2<-dcast(planting_summary,Bush~Species, fun.aggregate = sum)


planting_summary2<-left_join(planting_summary2, species_richness)

colnames(planting_summary2)[64]<-"planting_diversity"
colnames(planting_summary2)[65]<-"planting_richness"
colnames(planting_summary2)[1]<-"bush"

#calculate most abundant species planted in a bush area
planting_summary2$max_species<-colnames(planting_summary2)[apply(planting_summary2,1,which.max)]

# calculate number of trees planted 

planting_count<-planting_summary %>%
  group_by(Bush) %>%
  summarise(Total_count=sum(na.omit(Total)))

#calculate shannon diversity of plants planted at each bush 

library(vegan)

planting_summary2[is.na(planting_summary2)] <- 0

planting_summary2$planting_diversity<-apply(planting_summary2[,2:55], 1, FUN=diversity)

colnames(planting_summary2)[64]<-"planting_diversity"
colnames(planting_summary2)[65]<-"planting_richness"
colnames(planting_summary2)[1]<-"bush"


# put planting data into full dataframe
data<-left_join(data,planting_summary2[,c(1,63:66)],by="bush")


# calculate shannon diversity in regenerating forest ####

data$shannon_diversity<-apply(data[,5:35], 1, FUN=diversity)

# join breeding data and habitat data ####

breedingdata1<-left_join(breedingdata1,select(data,c(nest,species_richness,shannon_diversity)))

# PERMANOVA analysis ####

# aggregate data into community matrix at bush level 

planting_summary3<-na.omit(planting_summary2)

community_matrix <- data %>%
  filter(!is.na(bush)) %>%
  group_by(bush) %>%
  summarise(across(CAB:SIL, sum)) %>%
  as.data.frame()

rownames(community_matrix) <- community_matrix$bush
community_matrix$bush <- NULL

str(community_matrix)

# only keep planting data for bush sections for which there is community composition data (26)
planting_summary3 <- planting_summary3[planting_summary3$bush %in% rownames(community_matrix), ]

# PERMANOVA of species richness of tree plantings on regenerating forest tree communities ####
adonis2(
  community_matrix ~ planting_richness,
  data = planting_summary3,
  method = "bray"
)
# richness of plantings significantly affected current community composition

# PERMANOVA of species diversity of tree plantings on regenerating forest tree communities 
adonis2(
  community_matrix ~ planting_diversity,
  data = planting_summary3,
  method = "bray"
)

# which species are driving the shift in community? ####
  
  species_cor <- apply(community_matrix, 2, function(x)
    cor(x, planting_summary3$planting_richness, method = "spearman"))

species_tests <- lapply(community_matrix, function(x)
  cor.test(x,
           planting_summary3$planting_richness,
           method = "spearman"))
species_tests

species_results <- data.frame(
  species = colnames(community_matrix),
  rho = sapply(species_tests, function(x) x$estimate),
  p.value = sapply(species_tests, function(x) x$p.value)
)

species_results

sort(species_cor, decreasing = TRUE)

species_fit <- envfit(nmds, community_matrix, permutations = 999)
plot(species_fit)

# PERMANOVA comparing planted vs regenerating tree communities #

# make matrix of planted community

planting_matrix<-planting_data %>%
  group_by(Bush,species_code) %>%
  summarise(Total=sum(Number_planted))

colnames(planting_matrix)[1]<-"bush"

setdiff(unique(planting_matrix$bush), rownames(community_matrix))

planting_matrix$bush[planting_matrix$bush == "0"] <- NA

planting_matrix<-planting_matrix %>% filter(!is.na(bush))
planting_matrix<-planting_matrix %>% filter(bush != "B12")

planting_matrix<-dcast(planting_matrix,bush~species_code, fun.aggregate = sum)

rownames(planting_matrix)<-planting_matrix$bush

all(colnames(planting_matrix) == colnames(community_matrix))

# some columns/species were planted but not found in habitat surveys

# Columns present in either planted or current communities


planting_matrix$COP<-planting_matrix$COP+planting_matrix$COP_M


# Make sure both matrices have all species columns, fill missing with 0
# Step 1: all species across both matrices
all_species <- union(colnames(planting_matrix), colnames(community_matrix))

# Step 2: add missing columns to planting_matrix and community_matrix
for (sp in setdiff(all_species, colnames(planting_matrix))) {
  planting_matrix[[sp]] <- 0
}

for (sp in setdiff(all_species, colnames(community_matrix))) {
  community_matrix[[sp]] <- 0
}

# Step 3: order columns identically
planting_matrix <- planting_matrix[, all_species]
community_matrix <- community_matrix[, all_species]

#remove pohutakawa because it is known that many pohutakawa were planted but also felled at some point
planting_matrix <- planting_matrix %>%
  select(-POH)
community_matrix <- community_matrix %>%
  select(-POH)


# test if they are different 

# Stack matrices
combined <- rbind(planting_matrix, community_matrix)
time <- factor(rep(c("Planted", "Current"),
                   each = nrow(planting_matrix)))
site <- factor(rep(rownames(planting_matrix), 2))

# Paired PERMANOVA to test if planted vs current communities are significantly different  ####
adonis2(combined ~ time,
        strata = site,
        method = "bray")

# how similar are the communities

# Example assuming rows in planting_matrix and current_aligned correspond
similarity <- numeric(nrow(planting_matrix))

for (i in 1:nrow(planting_matrix)) {
  bc_dist <- vegdist(rbind(planting_matrix[i, ], community_matrix[i, ]),
                     method = "bray")
  similarity[i] <- 1 - as.numeric(bc_dist)  # Bray–Curtis similarity
}

# Check results
similarity
mean(similarity)

## NMDS Plot of planted vs regenerated forest communities (figure 2b) ####

# plot this difference 
#  Align columns (species) in both matrices
all_species <- union(colnames(planting_matrix), colnames(community_matrix))

# Add missing species columns
for(sp in setdiff(all_species, colnames(planting_matrix))) planting_matrix[[sp]] <- 0
for(sp in setdiff(all_species, colnames(community_matrix))) community_matrix[[sp]] <- 0

# Reorder columns
planting_matrix <- planting_matrix[, all_species]
community_matrix <- community_matrix[, all_species]

#  Combine matrices for NMDS
combined <- rbind(planting_matrix, community_matrix)
time <- factor(rep(c("Planted", "Contemporary"), each = nrow(planting_matrix)))
site <- factor(rep(rownames(planting_matrix), 2))  # pairing info

#  Run NMDS 
set.seed(123)  # for reproducibility
nmds <- metaMDS(combined, distance = "bray", k = 2, trymax = 100)
nmds$stress
# Extract NMDS scores for plotting
scores_df <- as.data.frame(scores(nmds, display = "sites"))
scores_df$time <- time
scores_df$site <- site

#  extract hull values ####

nmds_hull <- 
  scores_df %>% 
  group_by(time) %>%
  slice(chull(NMDS1, NMDS2))

#  Plot NMDS with points coloured by time 
ggplot(scores_df, aes(NMDS1, NMDS2, colour = time)) +
  geom_point(size = 3) +
  # Draw lines connecting planted -> current for each bush
  geom_line(aes(group = site), alpha = 0.5, colour = "grey50") +
  theme_minimal() +
  scale_color_discrete(palette="Set2", labels=c("Planted","Regenerating"))+
  labs(colour = "Community")+
  annotate("text", x = 0.25, y =2, 
           label = paste("2d stress =", round(nmds$stress, 3))) +
  geom_polygon(data = nmds_hull, 
               aes(x = NMDS1, y = NMDS2, fill = time, group = time), 
               alpha = 0.30) +
  scale_fill_discrete(palette="Set2", guide="none")+
  theme(text = element_text(size = 14),
        axis.text = element_blank(),
        axis.ticks = element_blank())

### Principal component analysis ####

# 1 # make nestbox column and create new df of top 5 species (in current) and heights, and  with data 

habitat_data_pca <- data[,c(36:45)]

# Keep only columns with non-zero variance
habitat_data_pca <- habitat_data_pca[, apply(habitat_data_pca, 2, var, na.rm = TRUE) != 0]


# vegetation density PCA
m1 <- prcomp(habitat_data_pca[,c(29:38)], scale = TRUE, center = TRUE)

par(mfrow = c(1,1))
biplot(m1)

summary(m1)

m1$rotation[,1:5]

## save PC scores and put them into data
pc_scores <- as.data.frame(m1$x[, 1:3])

data_wrk <- data %>%
  group_by(bush, nestbox) %>%
  slice_sample(n = 1) %>%
  ungroup()

#put PCA scores into habitat data

data_wrk <- cbind(data_wrk, pc_scores)

#put PCA scores into breeding data

breedingdata1<-left_join(breedingdata1,data_wrk[,c(65:68)])

#### MCMCGLMM models ####

library(MCMCglmm)

# effects of planting on regenerated forest characteristics ####

# put Total_count column into data
names(planting_count)[1]<-"bush"
data_wrk<-left_join(data_wrk, planting_count)

# calculate median year of planting and put into data
year_of_planting<-planting_data %>%
  group_by(Bush) %>%
  summarise(median(Year))
colnames(year_of_planting)<-c("bush", "year_of_planting")
data_wrk<-left_join(data_wrk,year_of_planting, by="bush")

# import bush area data and calculated total trees planted per area

bush_areas<-read_csv("bush_areas.csv")

data_wrk<-left_join(data_wrk,bush_areas)

data_wrk$planted_per_area<-(data_wrk$Total_count/data_wrk$area_of_bush)

# create dataframe for modelling 
model.data.plant<-data.frame(data_wrk$bush, data_wrk$nest, data_wrk$collector, data_wrk$shannon_diversity, data_wrk$species_richness, 
                             data_wrk$planting_diversity, data_wrk$planting_richness,  data_wrk$max_species, 
                             data_wrk$max_species_no_pohu, data_wrk$Total_count, data_wrk$planted_per_area, data_wrk$year_of_planting)
#rename columns
names(model.data.plant)<-c("bush", "nest","collector" ,"shannon_diversity","species_richness",
                           "planting_diversity","planting_richness", "max_species", "max_species_no_pohu", "Total_count", "planted_per_area", "year_of_planting")
# remove NA values
model.data.plant<-na.omit(model.data.plant)

# create priors
priors2 <- list(
  R = list(V = 1, nu = 0.002),
  G = list(
    G1 = list(V = 1, nu = 0.002),  # bush
    G2 = list(V = 1, nu = 0.002)   # year
  )
)

# species diversity models ####
modelplant1 <- MCMCglmm(fixed = shannon_diversity ~ planting_diversity, random=~bush + collector, data = model.data.plant,
                        family = "gaussian",
                        prior = priors2, verbose=FALSE,
                        nitt= 100000, thin=100, burnin=1000)

autocorr(modelplant1$Sol)
plot(modelplant1$Sol)
#good
autocorr(modelplant1$VCV)
plot(modelplant1$VCV)

summary(modelplant1)

modelplant1.a <- MCMCglmm(fixed = shannon_diversity ~ planting_richness, random=~bush + collector, data = model.data.plant,
                        family = "gaussian",
                        prior = priors2, verbose=FALSE,
                        nitt= 100000, thin=100, burnin=1000)

autocorr(modelplant1.a$Sol)
plot(modelplant1.a$Sol)
#good
autocorr(modelplant1.a$VCV)
plot(modelplant1.a$VCV)

summary(modelplant1.a)

# most abundantly planted species vs contemporary diversity
modelplant2a <- MCMCglmm(fixed = shannon_diversity ~ max_species, random=~bush, data = model.data.plant,
                         family = "gaussian",
                         prior = priors1, verbose=FALSE,
                         nitt= 100000, thin=100, burnin=1000)
autocorr(modelplant2a$Sol)
plot(modelplant2a$Sol)
#good
autocorr(modelplant2a$VCV)
plot(modelplant2a$VCV)

summary(modelplant2a)

# null model for effect of most abundant species on species richness
model_null <- MCMCglmm(
  shannon_diversity ~ 1,
  random = ~ bush,
  data = data,
  prior = priors1,
  verbose = FALSE
)
model_null$DIC -
  modelplant2c$DIC

# number of trees planted
effort_model1 <- MCMCglmm(fixed = shannon_diversity ~ planted_per_area, random=~bush, data = model.data.plant,
                          family = "gaussian",
                          prior = priors1, verbose=FALSE,
                          nitt= 100000, thin=100, burnin=1000)
autocorr(effort_model1$Sol)
plot(effort_model1$Sol)
#good
autocorr(effort_model1$VCV)
plot(effort_model1$VCV)
summary(effort_model1)

# year of planting vs contemporary diversity 
modelplant_year1 <- MCMCglmm(fixed = shannon_diversity ~ year_of_planting, random=~bush + collector, data = model.data.plant,
                             family = "gaussian",
                             prior = priors2, verbose=FALSE,
                             nitt= 100000, thin=100, burnin=1000)
autocorr(modelplant_year1$Sol)
plot(modelplant_year1$Sol)
#good
autocorr(modelplant_year1$VCV)
plot(modelplant_year1$VCV)
summary(modelplant_year1)



# species richness models ####

modelplant3a <- MCMCglmm(fixed = species_richness ~ planting_richness, random=~bush + collector, data = model.data.plant,
                          family = "gaussian",
                          prior = priors2, verbose=FALSE,
                          nitt= 100000, thin=100, burnin=1000)

autocorr(modelplant3a$Sol)
plot(modelplant3a$Sol)
#good
autocorr(modelplant3a$VCV)
plot(modelplant3a$VCV)

summary(modelplant3a)

modelplant3b <- MCMCglmm(fixed = species_richness ~ planting_diversity, random=~bush + collector, data = model.data.plant,
                          family = "gaussian",
                          prior = priors2, verbose=FALSE,
                          nitt= 100000, thin=100, burnin=1000)

autocorr(modelplant3b$Sol)
plot(modelplant3b$Sol)
#good
autocorr(modelplant3b$VCV)
plot(modelplant3b$VCV)

summary(modelplant3b)

# total trees planted
modelplant2 <- MCMCglmm(fixed = species_richness ~ Total_count, random=~bush, data = model.data.plant,
                        family = "gaussian",
                        prior = priors1, verbose=FALSE,
                        nitt= 100000, thin=100, burnin=1000)


autocorr(modelplant2$Sol)
plot(modelplant2$Sol)
#good
autocorr(modelplant2$VCV)
plot(modelplant2$VCV)

summary(modelplant2)

# most abundant trees planted
modelplant2b <- MCMCglmm(fixed = species_richness ~ max_species, random=~bush, data = model.data.plant,
                         prior = priors1, verbose=FALSE,
                         nitt= 100000, thin=100, burnin=1000)
summary(modelplant2b)

autocorr(modelplant2b$Sol)
plot(modelplant2b$Sol)
#good
autocorr(modelplant2b$VCV)
plot(modelplant2b$VCV)

# null model for effect of most abundant species on species richness
model_null <- MCMCglmm(
  species_richness ~ 1,
  random = ~ bush,
  data = data,
  prior = priors1,
  verbose = FALSE
)
model_null$DIC -
  modelplant2b$DIC

# DIC difference quite large so obvious overall effect of max_species
# which species are different to each other?

library(emmeans)
emm <- emmeans(modelplant2b, ~ max_species, data=model.data.plant)
emm
pairs(emm)

# number of trees planted

effort_model2 <- MCMCglmm(fixed = species_richness ~ planted_per_area, random=~bush, data = model.data.plant,
                          family = "gaussian",
                          prior = priors1, verbose=FALSE,
                          nitt= 100000, thin=100, burnin=1000)
summary(effort_model2)

autocorr(effort_model2$Sol)
plot(effort_model2$Sol)
#good
autocorr(effort_model2$VCV)
plot(effort_model2$VCV)

# year of planting
modelplant_year2 <- MCMCglmm(fixed = species_richness ~ year_of_planting, random=~bush + collector, data = model.data.plant,
                             family = "gaussian",
                             prior = priors2, verbose=FALSE,
                             nitt= 100000, thin=100, burnin=1000)
autocorr(modelplant_year2$Sol)
plot(modelplant_year2$Sol)
#good
autocorr(modelplant_year2$VCV)
plot(modelplant_year2$VCV)

summary(modelplant_year2)


### models of habitat density (effects of planting) ####

density_plus_planting<-density_plus_planting[,-c(36:40)]
density_plus_planting<-na.omit(density_plus_planting)

# planting diversity 

pc_model1 <- MCMCglmm(fixed = PC1 ~ planting_diversity, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
autocorr(pc_model1$Sol)
plot(pc_model1$Sol)
#good
autocorr(pc_model1$VCV)
plot(pc_model1$VCV)
summary(pc_model1) 
# significant effect of planting diversity on PC1 (negative effect)
pc_model2 <- MCMCglmm(fixed = PC2 ~ planting_diversity, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
autocorr(pc_model2$Sol)
plot(pc_model2$Sol)
#good
autocorr(pc_model2$VCV)
plot(pc_model2$VCV)
summary(pc_model2) 
# no effect on PC2

pc_model3 <- MCMCglmm(fixed = PC3 ~ planting_diversity, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
autocorr(pc_model3$Sol)
plot(pc_model3$Sol)
#good
autocorr(pc_model3$VCV)
plot(pc_model3$VCV)
summary(pc_model3) 
# or PC3

# now species richness

pc_model1a <- MCMCglmm(fixed = PC1 ~ planting_richness, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
autocorr(pc_model1a$Sol)
plot(pc_model1a$Sol)
#good
autocorr(pc_model1a$VCV)
plot(pc_model1a$VCV
summary(pc_model1a) #no effect

pc_model2a <- MCMCglmm(fixed = PC2 ~ planting_richness, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
autocorr(pc_model2a$Sol)
plot(pc_model2a$Sol)
#good
autocorr(pc_model2a$VCV)
plot(pc_model2a$VCV
summary(pc_model2a) # no effect
# no effect on PC2

pc_model3a <- MCMCglmm(fixed = PC3 ~ planting_richness, random=~bush, data = density_plus_planting,
                      family = "gaussian",
                      prior = priors1, verbose=FALSE,
                      nitt= 100000, thin=100, burnin=1000)
utocorr(pc_model3a$Sol)
plot(pc_model3a$Sol)
#good
autocorr(pc_model3a$VCV)
plot(pc_model3a$VCV
summary(pc_model3a) 
# no effect
# 1. make priors and sort out data for breeding models ####
a <- 1000
priors <- list(R = list(V = diag(1), nu = 0.002),
               G = list(G1 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G2 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a),
                        G3 = list(V = diag(1), nu = 1, alpha.mu = 0, alpha.V = diag(1)*a)))

# exclude NAs from data, make new dataframe with only variables of interest
# make model data frame, PC1,2 and 3 are from habitat density PCA
model.data<-data.frame(breedingdata1$eggs,breedingdata1$fledged, breedingdata1$hatched,breedingdata1$unhatched,breedingdata1$shannon_diversity,
                       breedingdata1$species_richness, breedingdata1$nearest_feeder, breedingdata1$female.code, breedingdata1$year, breedingdata1$female.age,
                       breedingdata1$PC1, breedingdata1$PC2, breedingdata1$PC3)
#rename columns
names(model.data)<-c("eggs","fledged", "hatched","unhatched", "shannon_diversity","species_richness",  "nearest_feeder", "female.code", "year", "female.age", "PC1", "PC2", "PC3")
# remove NA values
model.data<-na.omit(model.data)

model.data<-filter(model.data, year>=2016)

# models of clutch size vs regenerating forest tree species diversity and species richness ####

model1 <- MCMCglmm(fixed = eggs ~ shannon_diversity , random=~female.code+female.age+year, data = na.omit(breedingdata1),
                   family = "poisson",
                   prior = priors, verbose=FALSE,      
                   nitt= 100000, thin=100, burnin=1000)
)

# use these plots to check that variation around the mean is constant
autocorr(model1$Sol)
# these values show autocorrelation of the mean of successive samples from the chain
# what you want is LOW correlation - they should be independent. 
# these values are a bit high
plot(model1$Sol)
# plot also shows slight trends

# use these plots to check for autocorrelation in variance from successive samples in the chain (trends in left hand plots)
plot(model1$VCV)
# a bit of a pattern.. 
autocorr(model1$VCV)
# these values show autocorrelation between variances from successive samples in the chain, again low correlation is good
# here we have slightly high values which suggests the model needs to be run for more iterations, lets try it

summary(model1)
# with nearest feeder
model1 <- MCMCglmm(fixed = eggs ~ shannon_diversity * nearest_feeder, random=~female.code+female.age+year, data = model.data,
                   family = "poisson",
                   prior = priors, verbose=FALSE,
                   nitt= 100000, thin=100, burnin=1000)
#without nearest feeder
model1.1 <- MCMCglmm(fixed = eggs ~ shannon_diversity, random=~female.code+female.age+year, data = model.data,
                     family = "poisson",
                     prior = priors, verbose=FALSE,
                     nitt= 1000000, thin=1000, burnin=10000)
# try diagnostics again

autocorr(model1$Sol)
plot((model1$Sol))
#good
autocorr(model1$VCV)
plot(model1$VCV)
# good

summary(model1)
summary(model1.1)


# with species richness
# with nearest feeder included
model1.a <- MCMCglmm(fixed = eggs ~ species_richness * nearest_feeder, random=~female.code+female.age+year, data = model.data,
                     family = "poisson",
                     prior = priors, verbose=FALSE,
                     nitt= 100000, thin=100, burnin=1000)
#without feeder variable
model1.b <- MCMCglmm(fixed = eggs ~ species_richness, random=~female.code+female.age+year, data = model.data,
                     family = "poisson",
                     prior = priors, verbose=FALSE,
                     nitt= 100000, thin=100, burnin=1000)

summary(model1.a)
summary(model1.b)

# models of fledging numbers vs regenerating forest tree species diversity and species richness ####
#species richness
# with feeders
model2 <- MCMCglmm(fixed = fledged ~ species_richness * nearest_feeder, random=~female.code+female.age+year, data = model.data,
                   family = "poisson",
                   prior = priors, verbose=FALSE,
                   nitt= 100000, thin=100, burnin=1000)
#without feeder
model2.a <- MCMCglmm(fixed = fledged ~ species_richness, random=~female.code+female.age+year, data = model.data,
                     family = "poisson",
                     prior = priors, verbose=FALSE,
                     nitt= 100000, thin=100, burnin=1000)

autocorr(model2$Sol)
#good
autocorr(model2$VCV)
# good

autocorr(model2.a$Sol)
#good
autocorr(model2.a$VCV)
# good

summary(model2)
summary(model2.a)

# species diversity 
# with feeder distance
model2.b <- MCMCglmm(fixed = fledged ~ shannon_diversity + nearest_feeder, random=~female.code+female.age+year, data = model.data,
                     family = "poisson",
                     prior = priors, verbose=FALSE,
                     nitt= 100000, thin=100, burnin=1000)
# without feeder distance
model2c <- MCMCglmm(fixed = fledged ~ shannon_diversity, random=~female.code+female.age+year, data = model.data,
                    family = "poisson",
                    prior = priors, verbose=FALSE,
                    nitt= 100000, thin=100, burnin=1000)

autocorr(model2.b$Sol)
#good
autocorr(model2.b$VCV)
# good

autocorr(model2c$Sol)
#good
autocorr(model2c$VCV)
# good

# models of hatching success vs regenerating forest tree species diversity and species richness ####

# hatching vs diversity with feeder
model3 <- MCMCglmm(fixed = cbind(hatched,unhatched) ~ shannon_diversity * nearest_feeder, random=~female.code+female.age+year, data = model.data,
                   family = "multinomial2",
                   prior = priors, verbose=FALSE,
                   nitt= 1000000, thin=1000, burnin=10000)
# hatching vs diversity no feeder
model3a <- MCMCglmm(fixed = cbind(hatched,unhatched) ~ shannon_diversity , random=~female.code+female.age+year, data = model.data,
                    family = "multinomial2",
                    prior = priors, verbose=FALSE,
                    nitt= 100000, thin=100, burnin=1000)

autocorr(model3$Sol)
#good
autocorr(model3$VCV)
# good

summary(model3)

autocorr(model3a$Sol)
#good
autocorr(model3a$VCV)
# good

summary(model3a)
# species richness on hatching success

model3b <- MCMCglmm(fixed = cbind(hatched,unhatched) ~ species_richness * nearest_feeder, random=~female.code+female.age, data = model.data,
                    family = "multinomial2",
                    prior = priors, verbose=FALSE,
                    nitt= 1000000, thin=1000, burnin=10000)
#without feeder
model3c <- MCMCglmm(fixed = cbind(hatched,unhatched) ~ species_richness , random=~female.code+female.age+year, data = model.data,
                    family = "multinomial2",
                    prior = priors, verbose=FALSE,
                    nitt= 1000000, thin=1000, burnin=10000)

autocorr(model3b$Sol)
#good
autocorr(model3b$VCV)
# good

summary(model3b)

autocorr(model3c$Sol)
#good
autocorr(model3c$VCV)

summary(model3c)

# models of breeding success vs vegetation density PCA scores ####

# hatching success
model4 <- MCMCglmm(fixed = cbind(hatched,unhatched) ~ PC1 + PC2 + PC3 , random=~female.code+female.age+year, data = filter(model.data, year>= 2016),
                   family = "multinomial2",
                   prior = priors, verbose=FALSE,
                   nitt= 100000, thin=100, burnin=1000)

summary(model4)

autocorr(model4$Sol)
#good
autocorr(model4$VCV)
# good

# clutch size
model4a <- MCMCglmm(fixed = eggs ~ PC1 + PC2 + PC3 , random=~female.code+female.age+year, data = filter(model.data, year>= 2016),
                    family = "poisson",
                    prior = priors, verbose=FALSE,
                    nitt= 100000, thin=100, burnin=1000)

summary(model4a)

autocorr(model4$Sol)
#good
autocorr(model4$VCV)
# good

#fledgling numbers
model4b <- MCMCglmm(fixed = fledged ~ PC1 + PC2 + PC3 , random=~female.code+female.age+year, data = filter(model.data, year>= 2016),
                    family = "poisson",
                    prior = priors, verbose=FALSE,
                    nitt= 100000, thin=100, burnin=1000)

summary(model4b)

autocorr(model4a$Sol)
#good
autocorr(model4a$VCV)
# good

# figure 2a plot ####

# surveyed total ###

surveyed_total <- data[,5:35] %>%
  summarise(across(everything(), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "column",
    values_to = "total_sum")%>%
  arrange(desc(total_sum))

surveyed_total$variable<-"surveyed"
surveyed_total_plot<-surveyed_total[1:10,]


planting_total <-planting_summary2 %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(),
               names_to = "column",
               values_to = "total_sum") %>%
  arrange(desc(total_sum))

planting_total$variable<-"planted"
planting_total_plot<-planting_total[2:11,]

colnames(planting_total_plot)[1]<-"latin_name"
colnames(surveyed_total_plot)[1]<-"code"


library(readr)
tree_names <- read_csv("~/Library/CloudStorage/GoogleDrive-f.morland@sheffield.ac.uk/My Drive/habitat stuff 2026/tree names.csv")
View(tree_names)

surveyed_total_plot<-left_join(surveyed_total_plot, tree_names)

community_composition<-bind_rows(planting_total_plot,surveyed_total_plot[,-c(1,4)])

library(stringr)
community_composition$latin_name <- str_to_sentence(community_composition$latin_name)

community_composition$latin_name[community_composition$latin_name == "Pittisporum crassifolium"] <- 
  "Pittosporum crassifolium"


# colour palette 
library(RColorBrewer)
set2_13 <- colorRampPalette(
  brewer.pal(8, "Set2")
)(16)

set2_13
set.seed(123)
library(ggtext)

ggplot(data=community_composition, aes(fill=latin_name, x=variable, y=total_sum))+
  geom_bar(position="fill", stat="identity", width=0.4)+
  theme_minimal()+
  labs(
    x = "Community",
    y = "Proportion",
    fill = "Species"
  ) +
  scale_fill_manual(
    values = sample(set2_13),
    labels = function(x) paste0("<i>", x, "</i>")
  ) +
  theme(
    legend.text = element_markdown()
  )+
  scale_x_discrete(labels = c(
    "planted" = "Planted",
    "surveyed" = "Regenerating"))

## figure 3 plot ####

library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# Helper function for MCMCglmm objects
extract_mcmc_post <- function(model, response_name){
  as.data.frame(model$Sol) %>%         # fixed effect posteriors
    select(contains("PC")) %>%         # keep only PC1, PC2, PC3
    pivot_longer(cols = everything(), names_to = "term", values_to = "posterior") %>%
    group_by(term) %>%
    summarise(
      post.mean = mean(posterior),
      l95 = quantile(posterior, 0.025),
      u95 = quantile(posterior, 0.975),
      .groups = "drop"
    ) %>%
    mutate(response = response_name)
}

# Apply to each model
df1 <- extract_mcmc_post(model4a, "clutch_size")
df2 <- extract_mcmc_post(model4, "hatching_success")
df3 <- extract_mcmc_post(model4b, "fledglings")

# Combine
forest_data <- bind_rows(df1, df2, df3)

forest_data <- forest_data %>%
  mutate(term = factor(term, levels = c("PC2", "PC1", "PC3")))

ggplot(forest_data, aes(x = post.mean, y = term, color = response)) +
  geom_point(position = position_dodge(width = 0.7), size = 3) +
  geom_errorbarh(aes(xmin = l95, xmax = u95),
                 position = position_dodge(width = 0.7), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(x = "Posterior mean (log scale)", y = "PC predictor", color = "Response") +
  scale_color_discrete(palette="Set2", labels=c("Clutch Size", "Fledglings", "Hatching Success"))+
  theme_minimal() +
  theme(text = element_text(size = 12))
