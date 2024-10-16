## ---- 1. IMPORT PACKAGE ----

#{{multivariate analysis}}
library(FactoMineR) 
library(missMDA) #permits to handle missing values in principal component methods (PCA, CA, MCA, MFA, FAMD)
library(ade4)
library(factoextra) #extract information from multivariate analysis
library(mixOmics) #discriminant analysis PLS-DA and sPLS-DA
####library(uwot) #Uniform Manifold Approximation and Projection (UMAP)
####library(phateR) #Markov Affinity-Based Graph Imputation

#{{hierarchical clustering with bootstrap values}}
####library(pvclust)

#{{multivariate compositional analysis}}
####library (ALDEx2)
#{{multivariate Count-Based Differential Abundance Analysis}}
library(edgeR)

#{{correlation figure}}
library(corrplot)
library(ggcorrplot)

#{{dendogram graph}}
library(ggdendro)
library(dendextend)

#{{Circular and network graph}}
####library(ggraph)
####library(igraph)

##{{Heatmap graph}}
####library(ComplexHeatmap)

#{{vendiagram analysis}}
library(ggVennDiagram)
library(UpSetR)

#{{ecological analysis; diversity indices, permanova on distance, multivariate analysis}}
library(vegan)
library(pairwiseAdonis)
library(hilldiv)

#{{normality test}}
library(nortest)

#{{data distribution}}
####library(performance)

#{{Indicator species analysis}}
library(indicspecies)

#{{Import functional database}}
library(metagMisc) #use to import FUNGuild, but many other tips as
#pairwise dissimilarity boxplots, prevalence plots, diversity profiles based on Hill numbers
library(fungaltraits) #use to import FunFun database

#{{bi-network analysis}}
library(bipartite) 

#{{all inclusive packages}}
####library(phyloseq)
####library(microeco)
####library(microbiomeExplorer)

#{{color palette}}
library(RColorBrewer)
library(viridis)
library(ggsci)

#{{figure and text manipulation extra}}
library(cowplot) 
library(ggrepel)
library(scales)
library(ggforce)
library(ggpubr)
library(patchwork)
library(ggtext) #allow text modification, for example name in italic with *name* 
#{{Marginal distribution}}
library(ggside)
#{{convert plot in ggplot object}}
library(ggplotify)

#{{table manipulation and graph, including ggplot2}}
library(tidyverse)

## ---- 2. VARIABLE SELECTION ----
#{{.table files correspond to several projects, need to select ORACLE samples}}

input1 <- "data/metadata.txt"
input3 <- "data/Oracle_soils_16S_341F_785R_87_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"
input4b <- "data/Oracle_run_20200106_ITS2_ITS86F_ITS4_169_samples.OTU.filtered.cleaved.nosubstringOTUs.mumu.table2"

#{{set up parameter}}
bact_taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family", 
                           "genus", "species")
its_taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family",
                          "genus", "species")
seed <- 1 #set up the seed for statistics
percentage_threshold <- 1

## ---- 3. IMPORT AND FORMAT METADATA ----
metadata <- input1 %>%
  read_tsv() %>%
  rename_with(~gsub(".", "_", ., fixed = TRUE))

## ---- 4. IMPORT AND FORMAT OTU ---- 
## -------- 4.1. bacteria (soils) --------

bact_OTUs <- input3 %>% 
  #{{import data}}
  read_tsv() %>%
  #{{format column names}}
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~gsub("OR", "BU",.)) %>%
  rename(CONT_BU_ext = BU_T_extr) %>%
  rename(CONT_BU_PCR1 = T_1) %>%
  rename(CONT_BU_PCR2 = T_2) %>%
  #{{select samples from ORACLE project}}
  select(OTU, taxonomy, starts_with("CONT", ignore.case = FALSE), 
         starts_with("BU", ignore.case = FALSE)) %>%
  #{{keep only resequenced samples when available}}
  select(-BU_SER_10, -BU_ST_10, -BU_YL_19, -BU_TS_3) %>%
  #{{homogeneize resequenced sample names}}
  rename(BU_SER_10 = BU_SER_10_S) %>%
  rename(BU_ST_10 = BU_ST_10_S) %>%
  rename(BU_YL_19 = BU_YL_19_S) %>%
  rename(BU_TS_3 = BU_TS_3_S) %>%
  rename_with(~sub("\\_1$", "_01", .)) %>%
  rename_with(~sub("\\_2$", "_02", .)) %>%
  rename_with(~sub("\\_3$", "_03", .)) %>%
  rename_with(~sub("\\_4$", "_04", .)) %>%
  rename_with(~sub("\\_5$", "_05", .)) %>%
  rename_with(~sub("\\_6$", "_06", .)) %>%
  rename_with(~sub("\\_7$", "_07", .)) %>%
  rename_with(~sub("\\_8$", "_08", .)) %>%
  rename_with(~sub("\\_9$", "_09", .)) %>%
  #{{sequences with "No_hit" has to be deleted}}
  filter(grepl("^Bacteria|^Archaea",taxonomy)) %>%
  separate(taxonomy, bact_taxonomic_levels, sep = "[|]", extra = "drop") %>%
  #{{sequences affiliated to chloroplast and mitonchondria are excluded}}
  filter(order != "Chloroplast", family != "Mitochondria") %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE), 
         starts_with("BU", ignore.case = FALSE)) %>%
  droplevels()

## -------- 4.3. fungi (soils) --------

ITS2_OTUs <- input4b %>%
  #{{import data}}
  read_tsv() %>%
  #{{format column names}}
  rename_with(~sub("\\_.*", "", .)) %>%
  rename_with(~gsub("-", "_", .)) %>%
  rename_with(~gsub("OR", "BU",.)) %>%
  rename(CONT_BU_ext = BU_T_extr) %>%
  rename(CONT_BU_PCR1 = T_1) %>%
  rename(CONT_BU_PCR2a = T_2) %>%
  rename(CONT_BU_PCR2b = T_3) %>%
  #{{select samples from ORACLE project}}
  select(OTU, taxonomy, starts_with("CONT", ignore.case = FALSE),
         starts_with("BU_RT_", ignore.case = FALSE),
         starts_with("BU_SERL_", ignore.case = FALSE),
         starts_with("BU_ST_", ignore.case = FALSE),
         starts_with("BU_TS_", ignore.case = FALSE),
         starts_with("BU_YL_", ignore.case = FALSE)) %>%
  #{{homogeneize resequenced sample names}}
  rename_with(~sub("\\_1$", "_01", .)) %>%
  rename_with(~sub("\\_2$", "_02", .)) %>%
  rename_with(~sub("\\_3$", "_03", .)) %>%
  rename_with(~sub("\\_4$", "_04", .)) %>%
  rename_with(~sub("\\_5$", "_05", .)) %>%
  rename_with(~sub("\\_6$", "_06", .)) %>%
  rename_with(~sub("\\_7$", "_07", .)) %>%
  rename_with(~sub("\\_8$", "_08", .)) %>%
  rename_with(~sub("\\_9$", "_09", .)) %>%
  #{{"Fungi" are selected}}
  filter(grepl("Fungi", taxonomy)) %>%
  #{{select only OTUs and controls}}
  select(OTU, starts_with("CONT", ignore.case = FALSE),
         starts_with("BU", ignore.case = FALSE)) %>%
  droplevels()

## ---- 5. IMPORT AND FORMAT TAXONOMIC DATA ----
## -------- 5.1. bacteria (soils) --------

bact_tax <- bact_OTUs %>%
  #{{merge taxonomy}}
  left_join(select(read_tsv(input3), OTU, taxonomy), by = "OTU") %>%
  #{{select only OTU number and taxonomy}}
  select(OTU,taxonomy) %>%  
  separate(taxonomy, bact_taxonomic_levels, sep = "[|]", extra = "drop") %>%
  #{{replace missing taxonomic information and homogeneize the terms}}
  mutate_if(is.character, ~sub("uncultured", "unidentified", .)) %>% 
  mutate_if(is.character, ~sub("uncultured_bacterium", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("unidentified_bacterium", "unidentified", .)) %>%
  mutate_if(is.character, ~sub("[*]", "unidentified", .)) %>%
  mutate_if(is.character , replace_na, replace = "unidentified") %>%
  #{{replace unidentified by the higher taxonomic rank}}
  mutate(phylum = if_else(str_detect(phylum, "unidentified"), 
                          str_c("unidentified_", kingdom), phylum)) %>%
  mutate(class = if_else(str_detect(class, "unidentified"), 
                         str_c("unidentified_", phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "unidentified"), 
                         str_c("unidentified_", class), order)) %>%
  mutate(family = if_else(str_detect(family, "unidentified"), 
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "unidentified"), 
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "unidentified"), 
                           str_c("unidentified_", genus), species)) %>%
  mutate(phylum = if_else(str_detect(phylum, "metagenome"), 
                          str_c("unidentified_", kingdom), phylum)) %>%
  mutate(class = if_else(str_detect(class, "metagenome"), 
                         str_c("unidentified_", phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "metagenome"), 
                         str_c("unidentified_", class), order)) %>%
  mutate(family = if_else(str_detect(family, "metagenome"), 
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "metagenome"), 
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "metagenome"),
                           str_c("unidentified_", genus), species)) %>%
  #{{homogenize unidentified terms for taxonomic rank}}
  mutate(class = str_replace(class,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{homogenize terms for unidentified or uncharacterized species (sp)}}
  mutate(species = if_else(str_detect(species, "_sp"), 
                           str_c("unidentified_", genus), species)) %>% 
  mutate(species = str_replace(species, "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{convert OTU column in character from numeric}}
  mutate_if(is.numeric, as.character)

## -------- 5.3. fungi (soils) --------

ITS2_tax <- ITS2_OTUs %>% 
  #{{merge taxonomy}}
  left_join(select(read_tsv(input4b), OTU, taxonomy), by = "OTU") %>%
  #{{select only OTU number and taxonomy}}
  select(OTU, taxonomy) %>%
  #{{"Fungi" are selected}}
  separate(taxonomy, its_taxonomic_levels, sep = "[|]", extra = "drop") %>%
  #{{remove "k__" bits, etc in taxonomic names and replace missing taxonomic information}}
  mutate_if(is.character, ~sub("^[kpcofgs]__", "", .)) %>%
  mutate_if(is.character, ~sub("[*]", "unidentified", .)) %>%
  #{{replace unidentified by the higher taxonomic rank}}
  mutate(phylum = if_else(str_detect(phylum, "unidentified"),
                          str_c("unidentified_", kingdom), phylum)) %>%
  mutate(class = if_else(str_detect(class, "unidentified"),
                         str_c("unidentified_",phylum), class)) %>%
  mutate(order = if_else(str_detect(order, "unidentified"),
                         str_c("unidentified_",class), order)) %>%
  mutate(family = if_else(str_detect(family, "unidentified"),
                          str_c("unidentified_", order), family)) %>%
  mutate(genus = if_else(str_detect(genus, "unidentified"),
                         str_c("unidentified_", family), genus)) %>%
  mutate(species = if_else(str_detect(species, "unidentified"),
                           str_c("unidentified_",genus), species)) %>%
  #{{homogneize unidentified terms for taxonomic rank}}
  mutate(class = str_replace(class,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(order = str_replace(order,
                             "unidentified_unidentified_",
                             "unidentified_")) %>%
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(family = str_replace(family,
                              "unidentified_unidentified_",
                              "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(genus = str_replace(genus,
                             "unidentified_unidentified_",
                             "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>% 
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_unidentified_",
                               "unidentified_")) %>%
  mutate(species = str_replace(species,
                               "unidentified_unidentified_",
                               "unidentified_")) %>% 
  #{{homogenize terms for unidentified or uncharacterized species (sp)}}
  mutate(species = if_else(str_detect(species, "_sp"), 
                           str_c("unidentified_", genus), species)) %>% 
  mutate(species = str_replace(species, "unidentified_unidentified_",
                               "unidentified_")) %>%
  #{{convert OTU column in character from numeric}}
  mutate_if(is.numeric,as.character)

## ---- 6. DECONTAMINATION OF DATA ----

#{{For now, I decide to compute the total number of reads for each OTU
# in control samples. That sum is then substracted from occurrences
# of that OTU in true samples. OTUs present in control samples. The
# rationale is as follows:
# - strong in control samples, weak in true samples (contamination
#   specific to the control samples, will be eliminated by the
#   substraction, i.e, final abundance is 0),
# - present in control samples, present in true samples (systematic
#   contamination, will be mitigated by the substraction)
# - weak in control samples, strong in true samples (cross-talk, will
#   be eliminated/mitigated by the substraction)
# Finally, control samples should be eliminated from the statistical
# analysis now that the substraction is done}}

## -------- 6.1. bacteria (soils) --------

#{{Extract control samples}}
d <- replace(bact_OTUs, bact_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  gather("samples", "reads", -OTU) %>%
  filter(!is.na(reads))

#{{Extract samples}}

bact_OTUs_decont <- replace(bact_OTUs, bact_OTUs == 0, NA) %>%
  gather("samples", "reads", -c("OTU", starts_with("CONT"))) %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  spread(samples, reads, fill = 0) %>%
  select(-starts_with("CONT")) 
bact_OTUs_decont
#{{delete previous data}}
rm(d)
rm(bact_OTUs)

## -------- 6.3. fungi (soils) --------

#{{Extract control samples}}
d <- replace(ITS2_OTUs, ITS2_OTUs == 0, NA) %>%
  select(OTU, starts_with("CONT")) %>%
  gather("samples", "reads", -OTU) %>%
  filter(!is.na(reads))

ITS2_OTUs_decont <- replace(ITS2_OTUs, ITS2_OTUs == 0, NA) %>%
  gather("samples", "reads", -c("OTU", starts_with("CONT"))) %>%
  filter(!is.na(reads)) %>%
  #{{merge with control samples}}
  left_join(count(d, OTU, wt = reads), by = "OTU") %>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  spread(samples, reads, fill = 0) %>%
  select(-starts_with("CONT"))

#{{delete previous data}}
rm(d)
rm(ITS2_OTUs)

## ---- 7. FORMAT OTU DATA FOR RAREFACTION ----
## -------- 7.1. bacteria (soils) --------

bact_OTUs_t <- bact_OTUs_decont %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

#{{delete previous data}}
rm(bact_OTUs_decont)

## -------- 7.3. fungi (soils) --------

ITS2_OTUs_t <- ITS2_OTUs_decont %>%
  column_to_rownames(var = "OTU") %>%
  #{{transpose data}}
  t() %>% 
  as_tibble(rownames = NA)

#{{delete previous data}}
rm(ITS2_OTUs_decont)


## ---- 8. SELECT DATA RAREFACTION THRESHOLD ----
## -------- 8.1. bacteria (soils) -------- 

#{{Identification of Low samples and outliers + rarefaction threshold with smallest}}
quantile(rowSums(bact_OTUs_t))
min(rowSums(bact_OTUs_t))
sort(rowSums(bact_OTUs_t))
head(sort(rowSums(bact_OTUs_t)))
smallest_bact <- min(rowSums(bact_OTUs_t))

## -------- 8.3. fungi (soils) --------

#{{Identification of Low samples and outliers + rarefaction threshold with smallest}}
quantile(rowSums(ITS2_OTUs_t))
min(rowSums(ITS2_OTUs_t))
sort(rowSums(ITS2_OTUs_t))
head(sort(rowSums(ITS2_OTUs_t)))
smallest_its2 <- min(rowSums(ITS2_OTUs_t))

## ---- 9. DATA RAREFACTION (fixed seed) ------
## -------- 9.1. bacteria (soils) --------

set.seed(seed)
#{{rarefaction multiple times}}
bact_OTUs_rarefied_t <-  bact_OTUs_t %>% 
  rrarefy(smallest_bact) %>% 
  #{{create a new column with sum of reads for each OTU in all samples}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{delete all sample with = 0 in rarefied samples}}
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  filter(total != 0) %>%
  select(-total) %>%
  column_to_rownames(var = "OTU") %>% 
  t() %>%
  as_tibble(rownames = NA)


#{{export data}}
output <- "bact_OTUs_soils_final.txt"

tmp <- bact_OTUs_rarefied_t %>%
  rownames_to_column(var = "soil_code")

write_delim(tmp, file = output, delim = "\t", col_names = TRUE)

#{{delete previous data}}
rm(bact_OTUs_t) 

## -------- 9.3. fungi (soils) --------

set.seed(seed)
ITS2_OTUs_rarefied_t <- ITS2_OTUs_t %>% 
  rrarefy(smallest_its2) %>%
  #{{create a new column with sum of reads for each OTU in all samples}}
  t() %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "OTU") %>%
  #{{delete all sample with = 0 in rarefied samples}}
  mutate(total = rowSums(across(where(is.numeric)))) %>% 
  filter(total != 0) %>%
  select(-total) %>%
  column_to_rownames(var = "OTU") %>% 
  t() %>%
  as_tibble(rownames = NA)

#{{export data}}
output <- "its2_OTUs_soils_final.txt"

tmp <- ITS2_OTUs_rarefied_t %>%
  rownames_to_column(var = "soil_code")

write_delim(tmp, file = output, delim = "\t", col_names = TRUE)

#{{delete table data}}
rm(ITS2_OTUs_t)

## ---- 10. ALPHA-DIVERSITY ----
## -------- 10.1. bacteria (soils) --------
## ------------ 10.1.1. estimate usual diversity index -----------

#{{richness}}
R <- bact_OTUs_rarefied_t %>%
  estimateR() %>% 
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "soil_code") %>% 
  select(soil_code, S.obs)

#{{Shannon index H; richness + evenness}}
H <- bact_OTUs_rarefied_t %>% 
  diversity(index = "shannon", MARGIN = 1, base =
              exp(1)) %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{Pielou’s index of evenness; 0-1, 1 = max. evenness}}
S <- bact_OTUs_rarefied_t %>% 
  specnumber()
J <- H %>% 
  mutate(value = value/log(S))

#{{Simpson's D index; richness + evenness; 0-1; 1 - D rises as evenness increases}}
D <- bact_OTUs_rarefied_t %>% 
  diversity("simpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")
inv_D <- bact_OTUs_rarefied_t %>% 
  diversity("invsimpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 10.1.2. estimate effective number of species - Hill numbers -----------

#This equation has a parameter q that defines its sensitivity to rare species:
#low values of q favor rare species, high values of q favor abundant species. 
#For example, Shannon diversity is of order q = 1, and for Simpson diversity q = 2. 
#When q = 0, diversity = S (richness), because rare species are treated the same as abundant ones.

##{{richness as hill numbers}}
dR <- bact_OTUs_rarefied_t %>% 
  t() %>% 
  hill_div(qvalue = 0) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{richness + evenness (= shannon diversity) as hill numbers; rare species}}
dREr <- bact_OTUs_rarefied_t %>%
  t() %>%
  hill_div(qvalue = 1) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

##{{richness + evenness (= inverse Simpson = simpson diversity) as hill numbers; abundant species}}
dREa <- bact_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 2) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

#{{evenness (= shannon evenness) as hill ratio}}
eDRr <- dREr %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)

#{{evenness (= ssimpson evenness) as hill ratio}}
eDRa  <- dREa %>%
  left_join(dR, by = "soil_code") %>%
  mutate(value = value.x / value.y) %>%
  select(soil_code,value)

## ------------ 10.1.3. format diversity table (usual index + Hill numbers) -----------

#{{merge data}}
Indices_bact <- R %>%
  left_join(H, by = "soil_code") %>% 
  left_join(inv_D, by = "soil_code") %>%
  left_join(J, by = "soil_code") %>%
  left_join(dR, by = "soil_code") %>%
  left_join(dREr, by = "soil_code") %>%
  left_join(dREa, by = "soil_code") %>%
  left_join(eDRr, by = "soil_code") %>%
  left_join(eDRa, by = "soil_code")
names(Indices_bact) <- c("soil_code", "Richness", "Shannon", "Inv_Simpson", "Pielou", 
                           "Hill_Richness","Hill_Shannon","Hill_Inv_Simpson", 
                           "Hill_Shannon_evenness", "Hill_Simpson_evenness")

write_delim(Indices_bact, file = "data/divesities_bact.tsv", delim = "\t", col_names = TRUE)

## -------- 10.3. fungi (soils) --------
## ------------ 10.3.1. estimate usual diversity index -----------

#{{richness}}
R <- ITS2_OTUs_rarefied_t %>%
  estimateR() %>% 
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "plant_code") %>% 
  select(plant_code, S.obs)

#{{Shannon index H; richness + evenness}}
H <- ITS2_OTUs_rarefied_t %>% 
  diversity(index = "shannon", MARGIN = 1, base =
              exp(1)) %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")

#{{Pielou's index of evenness; 0-1, 1 = max. evenness}}
S <- ITS2_OTUs_rarefied_t %>% 
  specnumber()
J <- H %>% 
  mutate(value = value/log(S))

#{{Simpson's D index; richness + evenness; 0-1; 1 - D rises as evenness increases}}
D <- ITS2_OTUs_rarefied_t %>% 
  diversity("simpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")
inv_D <- ITS2_OTUs_rarefied_t %>% 
  diversity("invsimpson") %>% 
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")

## ------------ 10.3.2. estimate effective number of species - Hill numbers -----------

#This equation has a parameter q that defines its sensitivity to rare species:
#low values of q favor rare species, high values of q favor abundant species. 
#For example, Shannon diversity is of order q = 1, and for Simpson diversity q = 2. 
#When q = 0, diversity = S (richness), because rare species are treated the same as abundant ones.

##{{richness as hill numbers}}
dR <- ITS2_OTUs_rarefied_t %>% 
  t() %>% 
  hill_div(qvalue = 0) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")

#{{richness + evenness (= Shannon entropy) as hill numbers; rare species}}
dREr <- ITS2_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 1) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")

##r{{richness + evenness (= inverse Simpson) as hill numbers; abundant species}}
dREa <- ITS2_OTUs_rarefied_t %>% 
  t() %>%
  hill_div(qvalue = 2) %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "plant_code")

#{{evenness (= shannon evenness) as hill ratio}}
eDRr <- dREr %>%
  left_join(dR, by = "plant_code") %>%
  mutate(value = value.x / value.y) %>%
  select(plant_code,value)

#{{evenness (= ssimpson evenness) as hill ratio}}
eDRa  <- dREa %>%
  left_join(dR, by = "plant_code") %>%
  mutate(value = value.x / value.y) %>%
  select(plant_code,value)

## ------------ 10.3.3. format diversity table (usual index + Hill numbers) -----------

#{{merge data}}
Indices_its2 <- R %>%
  left_join(H, by = "plant_code") %>% 
  left_join(inv_D, by = "plant_code") %>%
  left_join(J, by = "plant_code") %>%
  left_join(dR, by = "plant_code") %>%
  left_join(dREr, by = "plant_code") %>%
  left_join(dREa, by = "plant_code") %>%
  left_join(eDRr, by = "plant_code") %>%
  left_join(eDRa, by = "plant_code")
names(Indices_its2) <- c("plant_code", "Richness", "Shannon", "Inv_Simpson", "Pielou", 
                           "Hill_Richness","Hill_Shannon","Hill_Inv_Simpson", 
                           "Hill_Shannon_evenness", "Hill_Simpson_evenness")

write_delim(Indices_its2, file = "data/divesities_fungi.tsv", delim = "\t", col_names = TRUE)
## A FAIRE ---- 11. FORMAT FUNCTIONAL MICROBIAL DATA ----
## -------- 11.1. bacteria (soils) ----
## ------------ 11.1.1. major legume nodulating N-fixers (Shamseldin et al 2017, Symbiosis) ----

bact_leg_Nfixers_rarefied_t <- bact_OTUs_rarefied_t %>%
  t() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  left_join(select(bact_tax, OTU, genus), by = "OTU") %>%
  #{{select major nodulating N-fixer taxonomic groups}}
  filter(genus == "Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium" | genus == "Mesorhizobium" | genus == "Sinorhizobium" |
           genus == "Ensifer" | genus == "Bradyrhizobium" | genus == "Microvirga") %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 11.1.2. major actinorhizal nodulating N-fixers (???, Symbiosis) ----

bact_act_Nfixers_rarefied_t <- bact_OTUs_rarefied_t %>%
  t() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  left_join(select(bact_tax, OTU, genus), by = "OTU") %>%
  #{{select major nodulating N-fixer taxonomic groups}}
  filter(genus == "Frankia") %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 11.1.3. AOB Nitrifiers ----

bact_AOB_rarefied_t <- bact_OTUs_rarefied_t %>%
  t() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  left_join(select(bact_tax, OTU, genus), by = "OTU") %>%
  #{{select AOB taxonomic groups}}
  filter(genus == "Nitrosomonas" | genus == "Nitrosospira" | genus == "Nitrosococcus" ) %>% 
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## ------------ 11.1.4. NOB Nitrifiers ----
bact_NOB_rarefied_t <- bact_OTUs_rarefied_t %>%
  t() %>% 
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>% 
  left_join(select(bact_tax, OTU, genus), by = "OTU") %>%
  #{{select NOB taxonomic groups}}
  filter(genus == "Nitrobacter" | genus == "Nitrospira" | genus == "Nitrococcus" ) %>%
  select(-OTU) %>%
  group_by(genus) %>%
  summarise_all(sum) %>%
  column_to_rownames(var = "genus") %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  rownames_to_column(var = "soil_code")

## -------- 11.2. fungi (soils) ----
## ---- 12. TAXONOMIC COUNT ----
## ---- 13. STATISTICS ON ENVIRONMENTAL DATA ----
## -------- 13.1. soils ----
## ------------ 13.1.1. normality test {normtest} ------

#{{soil data}}
tmp <- soils %>%
  #{{merge with metadata}}
  left_join(metadata, by = "spatial_code") %>%
  #{{delete duplicates}}
  distinct(spatial_code, .keep_all = TRUE) %>% 
  select(-soil_code) %>% 
  #{{subdata to select}}
  filter(land_use == "conventional")

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$pH)
hist(tmp$pH)
ggqqplot(tmp$Myc_rate, ylab = "pH")


#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(logit_trans(tmp$pH))
hist(sqrt(soils_factor$pH))
ggqqplot(sqrt(tmp$pH), ylab = "pH")

## ------------ 13.1.2. Correlation analyis {ggcorrplot} ----
## ---------------- 13.1.2.1. land use & orchard position ----

#{{select file name}}
output <- "Correlation_soil_landuse_position_spearman.pdf"

#{{Select conventional R data}}
tmp <- soils %>%
  #{{merge with metadata}}
  left_join(metadata, by = "spatial_code") %>%
  #{{delete duplicates}}
  distinct(spatial_code, .keep_all = TRUE) %>% 
  select(-soil_code) %>%
  #{{select the land_use}}
  filter(land_use == "conventional" & orchard_position == "row") %>%
  column_to_rownames("spatial_code") %>% 
  #{{#select only soil parameters}}
  select(where(is.numeric))

#{{Compute a correlation matrix, select the method}}
cor_soil <- cor(tmp, method = "spearman")
corr <- round(cor_soil, 1)

#{{Compute a correlation pvalues matrix}}
p.value <- cor_pmat(tmp)

#{{plot}}
corr1 <- ggcorrplot(corr, p.mat = p.value, sig.level = 0.05, 
                    hc.order = FALSE, hc.method = "ward2",
                    type = "upper", insig = "blank", lab = TRUE, outline.color = "grey",
                    legend.title = "Spearman",
                    ggtheme = ggplot2::theme_gray) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.text.y = element_text(vjust = 0.25))


#{{Select conventional R data}}
tmp <- soils %>%
  #{{merge with metadata}}
  left_join(metadata, by = "spatial_code") %>%
  #{{delete duplicates}}
  distinct(spatial_code, .keep_all = TRUE) %>% 
  select(-soil_code) %>%
  #{{select the land_use}}
  filter(land_use == "conventional" & orchard_position == "inter-row") %>%
  column_to_rownames("spatial_code") %>% 
  #{{#select only soil parameters}}
  select(where(is.numeric))

#{{Compute a correlation matrix, select the method}}
cor_soil <- cor(tmp, method = "spearman")
corr <- round(cor_soil, 1)

#{{Compute a correlation pvalues matrix}}
p.value <- cor_pmat(tmp)

#{{plot}}
corr2 <- ggcorrplot(corr, p.mat = p.value, sig.level = 0.05, 
                    hc.order = FALSE, hc.method = "ward2",
                    type = "upper", insig = "blank", lab = TRUE, outline.color = "grey",
                    legend.title = "Spearman",
                    ggtheme = ggplot2::theme_gray) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.text.y = element_text(vjust = 0.25))

#{{Select organic R data}}
tmp <- soils %>%
  #{{merge with metadata}}
  left_join(metadata, by = "spatial_code") %>%
  #{{delete duplicates}}
  distinct(spatial_code, .keep_all = TRUE) %>% 
  select(-soil_code) %>%
  #{{select the land_use}}
  filter(land_use == "organic" & orchard_position == "row") %>%
  column_to_rownames("spatial_code") %>% 
  #{{#select only soil parameters}}
  select(where(is.numeric))

#{{Compute a correlation matrix, select the method}}
cor_soil <- cor(tmp, method = "spearman")
corr <- round(cor_soil, 1)

#{{Compute a correlation pvalues matrix}}
p.value <- cor_pmat(tmp)

#{{plot}}
corr3 <- ggcorrplot(corr, p.mat = p.value, sig.level = 0.05, 
                    hc.order = FALSE, hc.method = "ward2",
                    type = "upper", insig = "blank", lab = TRUE, outline.color = "grey",
                    legend.title = "Spearman",
                    ggtheme = ggplot2::theme_gray) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.text.y = element_text(vjust = 0.25))

#{{Select organic R data}}
tmp <- soils %>%
  #{{merge with metadata}}
  left_join(metadata, by = "spatial_code") %>%
  #{{delete duplicates}}
  distinct(spatial_code, .keep_all = TRUE) %>% 
  select(-soil_code) %>%
  #{{select the land_use}}
  filter(land_use == "organic" & orchard_position == "inter-row") %>%
  column_to_rownames("spatial_code") %>% 
  #{{#select only soil parameters}}
  select(where(is.numeric))

#{{Compute a correlation matrix, select the method}}
cor_soil <- cor(tmp, method = "spearman")
corr <- round(cor_soil, 1)

#{{Compute a correlation pvalues matrix}}
p.value <- cor_pmat(tmp)

#{{plot}}
corr4 <- ggcorrplot(corr, p.mat = p.value, sig.level = 0.05, 
                    hc.order = FALSE, hc.method = "ward2",
                    type = "upper", insig = "blank", lab = TRUE, outline.color = "grey",
                    legend.title = "Spearman",
                    ggtheme = ggplot2::theme_gray) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
        axis.text.y = element_text(vjust = 0.25))

#{{assembly plots}}

((corr1 | corr2) / (corr3 | corr4)) +
  plot_annotation(tag_levels = 'A') +
  plot_layout(guides = "collect") & theme(legend.position = "right")


width <- 15
height <- 15

ggsave(file = output, width = width , height = height)

## ------------ 13.1.3. factor analysis of mixed data {FactoMineR and factoextra} ----
## ---------------- 13.1.3.1. format data ----------------

#{{select data}}
tmp <- metadata %>%
  select(-station_code) %>%
  
  drop_na() %>%
  column_to_rownames("soil_code")

## ---------------- 13.1.3.2. famd ----
set.seed(seed)
#{{Generate FAMD matrix}}
#{{possibilty to manage NA data before FAMD with missMDA}}
res_famd <- tmp %>%
  select(-nivins, -perschar, -nbract, 
         -petrum, -outtrans, -pet.out,
         -suptot, -dismai,-sexe, -langues) %>% #socio-eco parameters are excluded, only soil management + characteristic are kept
  FAMD(graph = FALSE) %>%

#{{visualize inertia af each axis}}
fviz_eig(res_famd, addlabels = TRUE, ylim = c(0, 50))

## -------------------- 13.3.1.2.1. plot famd (samples) ----

#{{select file name}}
output <- "FAMD_quali_soils_legacy_.pdf"

##{{extract and format data}}
eigenvalues <- res_famd %>%
  get_eigenvalue() %>% #extract inertia (eigen) values
  as_tibble()

set.seed(seed)
##{{select famd with qualitative variable axis to plot}}
quali_var <- fviz_famd_var(res_famd, "quali.var", col.var = "contrib", 
                            gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE,
                           legend.title = "Contribution (%)") +
  ggtitle("") +
  xlab(paste("Axis 1\nVariance = ", round(eigenvalues$variance.percent[1], digits = 1), " %", sep = "")) +
  ylab(paste("Axis 2\nVariance = ", round(eigenvalues$variance.percent[2], digits = 1), " %", sep = "")) 

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

#{{select file name}}
output <- "FAMD_quali_soils_legacy_.pdf"

set.seed(seed)
##{{select famd with quantitative variable axis}}
quanti_var <- fviz_famd_var(res_famd, "quanti.var", col.var = "contrib", 
                           gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                           repel = TRUE,
                           legend.title = "Contribution (%)") +
  ggtitle("") +
  xlab(paste("Axis 1\nVariance = ", round(eigenvalues$variance.percent[1], digits = 1), " %", sep = "")) +
  ylab(paste("Axis 2\nVariance = ", round(eigenvalues$variance.percent[2], digits = 1), " %", sep = "")) 

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)


## ---- 14. A TAXONOMIC OR FUNCTIONAL DISTRIBUTION ----
## -------- 14.1. bacteria (soils) ----
## ------------ 14.1.1. taxonomy ----
## ---------------- 14.1.1.1. locality ----
## -------------------- 14.1.1.1.1. barchart ----
#{{select file name and corresponding conditions}}
output <- "Barchart_bact_phylum_locality.pdf"

#{{simpson palette visualisation}}
show_col(pal_simpsons("springfield")(16))

tax <- "phylum"
fact_l <- "locality"

#{{select and re-order the factors}}
flevels_locality <- c("Raguitenga", "Tingressene", "Boussouma", "Sera", "Yilou")

bar1 <- bact_OTUs_rarefied_t %>%
  t() %>%
  as_tibble(rownames = NA) %>% 
  rownames_to_column(var = "OTU") %>%
  #{{select the taxonomic rank to analyse}}
  left_join(select(bact_tax, OTU, all_of(tax)), by = "OTU") %>%
  select(-OTU) %>%
  #{{group by taxonomic rank}}
  group_by_at(tax) %>% #group_by_at replace group_by to pass a variable as column name  
  summarise_all(sum) %>%
  column_to_rownames(var = tax) %>%
  t() %>%
  as_tibble(rownames = NA) %>%
  #{{select condition to plot}}
  rownames_to_column(var = "soil_code") %>%
  left_join(select(metadata, soil_code, all_of(fact_l)), 
            by = "soil_code") %>%
  select(-soil_code) %>%
  #{{transform data for plotting}}
  gather(key = "taxonomy", value = "reads", -all_of(fact_l)) %>%
  #{{determine relative abundance of each taxonomic rank per conditionl}}
  group_by_at(fact_l) %>%  #group_by_at replace group_by to import column name from vector
  mutate(Rel_abund = reads/sum(reads)) %>%
  ungroup() %>%
  #{{all taxa less than x% in abundance for each plant species are renamed}}
  mutate(taxonomy = if_else(Rel_abund < 0.001, "< 0.1%", taxonomy)) %>%
  #{{barchart plot}}
  ggplot(aes(x = locality, y = reads)) +
  #{{proportional stacked visualisation}}
  geom_col(position = "fill", ####color = "black",
           aes(fill = taxonomy)) +
  ggtitle("") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        legend.title = element_blank(),
        legend.position = "bottom",
        legend.text = element_text(hjust = 0, size = 12),
        axis.text.x = element_text(angle = 0, vjust = 0.9, hjust = 0.5, size = 12, colour = "black"),
        axis.title = element_text(size = 16)) +
  scale_x_discrete(limits = flevels_locality,
                   labels = c("Raguitenga\n(Korsimoro)", 
                              "Tingressene\n(Korsimoro)",
                              "Boussouma\n(Boussouma)", 
                              "Sera\n(Boussouma)",
                              "Yilou\n(Guibare)")) +
  scale_y_continuous(labels = c("0", "25", "50", "75", "100")) +
  #{{palette automatic}}
  ####scale_fill_nejm() +
  ####scale_fill_npg() +
  ####scale_fill_startrek() +
  ####scale_fill_viridis_d() +
  ####scale_fill_brewer() +
  scale_fill_simpsons(labels = c("<0.1%",
                            expression(italic("Acidobacteria")),
                            expression(italic("Actinobacteriota")),
                            expression(italic("Bacteroidota")),
                            expression(italic("Chloroflexi")),
                            expression(italic("Cyanobacteria")),
                            expression(italic("Entotheonellaeota")),
                            expression(italic("Firmicutes")),
                            expression(italic("Gemmatimonadota")),
                            expression(italic("Methylomirabilota")),
                            expression(italic("Myxococcota")),
                            expression(italic("Planctomycetota")),
                            expression(italic("Proteobacteria")),
                            "un. Bacteria",
                            expression(italic("Verrucomicrobiota")))) +
  ylab("Relative abundance (%)") +
  xlab("")

width <- 9
height <- 9

ggsave(file = output, width = width , height = height)

## ------------ 14.1.2. N-cycling ----
## ---------------- 14.1.2.1. locality ----
## -------------------- 14.1.2.1.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}

tmp <- bact_leg_Nfixers_rarefied_t %>%
  left_join(bact_act_Nfixers_rarefied_t , by = "soil_code") %>%
  left_join(bact_AOB_rarefied_t , by = "soil_code") %>%
  left_join(bact_NOB_rarefied_t , by = "soil_code") %>%
  mutate(Leg_N_fixers = `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
           Bradyrhizobium +
           #Ensifer +
           Mesorhizobium +
           Microvirga) %>%
  mutate(Act_N_fixers = Frankia) %>%
  mutate(AOB_nitrifiers = #Nitrosococcus + 
           Nitrosomonas + Nitrosospira) %>%
  mutate(NOB_nitrifiers = #Nitrobacter + Nitrococcus
           Nitrospira) %>%
  select(soil_code, Leg_N_fixers, Act_N_fixers, AOB_nitrifiers, NOB_nitrifiers) %>%
  #{{merge with metadata}}
  left_join(metadata, by = "soil_code") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Leg_N_fixers)
hist(tmp$Leg_N_fixers)
ggqqplot(tmp$Leg_N_fixers, ylab = "abundance")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(sqrt(tmp$Leg_N_fixers))
hist(sqrt(tmp$Leg_N_fixers))
ggqqplot(sqrt(tmp$Leg_N_fixers), ylab = "abundance")

## ------------------------ 14.1.2.1.1.1. anova and tukey tests {stats}---- 
## ------------------------ 14.1.2.1.1.2. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the global data}}
{
  
  tmp <- bact_leg_Nfixers_rarefied_t %>%
    left_join(bact_act_Nfixers_rarefied_t , by = "soil_code") %>%
    left_join(bact_AOB_rarefied_t , by = "soil_code") %>%
    left_join(bact_NOB_rarefied_t , by = "soil_code") %>%
    mutate(Leg_N_fixers = `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
             Bradyrhizobium +
             #Ensifer +
             Mesorhizobium +
             Microvirga) %>%
    mutate(Act_N_fixers = Frankia) %>%
    mutate(AOB_nitrifiers = #Nitrosococcus + 
             Nitrosomonas + Nitrosospira) %>%
    mutate(NOB_nitrifiers = #Nitrobacter + Nitrococcus
             Nitrospira) %>%
    select(soil_code, Leg_N_fixers, Act_N_fixers, AOB_nitrifiers, NOB_nitrifiers) %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    droplevels()
}

#{{locality effect}}
{
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  krus_loc_leg <- kruskal.test(Leg_N_fixers~locality,
                               data = tmp, na.action = na.omit)
  
  krus_loc_act <- kruskal.test(Act_N_fixers~locality,
                               data = tmp, na.action = na.omit)
  
  krus_loc_aob <- kruskal.test(AOB_nitrifiers~locality,
                               data = tmp, na.action = na.omit)
  
  krus_loc_nob <- kruskal.test(NOB_nitrifiers~locality,
                               data = tmp, na.action = na.omit)
  
  #{{multi groups}}
  dunn_loc_leg <- dunn.test::dunn.test(tmp$Leg_N_fixers, tmp$locality,
                                        method = "bh") 
  
  dunn_loc_act <- dunn.test::dunn.test(tmp$Act_N_fixers, tmp$locality,
                                       method = "bh")
  
  dunn_loc_aob <- dunn.test::dunn.test(tmp$AOB_nitrifiers, tmp$locality,
                                       method = "bh")
  
  dunn_loc_nob <- dunn.test::dunn.test(tmp$NOB_nitrifiers, tmp$locality,
                                       method = "bh")
  }
  
## -------------------- 14.1.2.1.2. barchart ----
## -------------------- 14.1.2.1.3. boxplot ---- 

#{{select file name}}
output <- "Boxplot_bact_Ncycling_locality.pdf"

#{{select and re-order the factors}}
flevels_locality <- c("Raguitenga", "Tingressene", "Boussouma", "Sera", "Yilou")

#{{jco palette visualisation}}
show_col(pal_jco("default")(10))
#{{select color palette}}
colpalette_locality <- c("#0073C2FF","#7AA6DCFF", "#A73030FF", "#CD534CFF", "#EFC000FF")


#{{igv palette visualisation}}
show_col(pal_igv("default")(16))
#{{select color palette}}
colpalette_Ncycling <- c("#924822FF", "#BA6338FF", "#D58F5CFF", "#E4AF69FF")

#{{format kruskall statistic}}
krus_loc_leg_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_leg$p.value, digits = 4),
                                 sep = "")

krus_loc_act_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_act$p.value, digits = 4),
                                 sep = "")

krus_loc_aob_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_aob$p.value, digits = 4),
                                 sep = "")

krus_loc_nob_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_nob$p.value, digits = 4),
                                 sep = "")

#{{format dunn statistic}}
dunn_loc_leg_BR_pvalue <- dunn_loc_leg$P.adjusted[1]
dunn_loc_leg_BS_pvalue <- dunn_loc_leg$P.adjusted[2]
dunn_loc_leg_RS_pvalue <- dunn_loc_leg$P.adjusted[3]
dunn_loc_leg_BT_pvalue <- dunn_loc_leg$P.adjusted[4]
dunn_loc_leg_RT_pvalue <- dunn_loc_leg$P.adjusted[5]
dunn_loc_leg_ST_pvalue <- dunn_loc_leg$P.adjusted[6]
dunn_loc_leg_BY_pvalue <- dunn_loc_leg$P.adjusted[7]
dunn_loc_leg_RY_pvalue <- dunn_loc_leg$P.adjusted[8]
dunn_loc_leg_SY_pvalue <- dunn_loc_leg$P.adjusted[9]
dunn_loc_leg_TY_pvalue <- dunn_loc_leg$P.adjusted[10]

dunn_loc_act_BR_pvalue <- dunn_loc_act$P.adjusted[1]
dunn_loc_act_BS_pvalue <- dunn_loc_act$P.adjusted[2]
dunn_loc_act_RS_pvalue <- dunn_loc_act$P.adjusted[3]
dunn_loc_act_BT_pvalue <- dunn_loc_act$P.adjusted[4]
dunn_loc_act_RT_pvalue <- dunn_loc_act$P.adjusted[5]
dunn_loc_act_ST_pvalue <- dunn_loc_act$P.adjusted[6]
dunn_loc_act_BY_pvalue <- dunn_loc_act$P.adjusted[7]
dunn_loc_act_RY_pvalue <- dunn_loc_act$P.adjusted[8]
dunn_loc_act_SY_pvalue <- dunn_loc_act$P.adjusted[9]
dunn_loc_act_TY_pvalue <- dunn_loc_act$P.adjusted[10]

dunn_loc_aob_BR_pvalue <- dunn_loc_aob$P.adjusted[1]
dunn_loc_aob_BS_pvalue <- dunn_loc_aob$P.adjusted[2]
dunn_loc_aob_RS_pvalue <- dunn_loc_aob$P.adjusted[3]
dunn_loc_aob_BT_pvalue <- dunn_loc_aob$P.adjusted[4]
dunn_loc_aob_RT_pvalue <- dunn_loc_aob$P.adjusted[5]
dunn_loc_aob_ST_pvalue <- dunn_loc_aob$P.adjusted[6]
dunn_loc_aob_BY_pvalue <- dunn_loc_aob$P.adjusted[7]
dunn_loc_aob_RY_pvalue <- dunn_loc_aob$P.adjusted[8]
dunn_loc_aob_SY_pvalue <- dunn_loc_aob$P.adjusted[9]
dunn_loc_aob_TY_pvalue <- dunn_loc_aob$P.adjusted[10]

dunn_loc_nob_BR_pvalue <- dunn_loc_nob$P.adjusted[1]
dunn_loc_nob_BS_pvalue <- dunn_loc_nob$P.adjusted[2]
dunn_loc_nob_RS_pvalue <- dunn_loc_nob$P.adjusted[3]
dunn_loc_nob_BT_pvalue <- dunn_loc_nob$P.adjusted[4]
dunn_loc_nob_RT_pvalue <- dunn_loc_nob$P.adjusted[5]
dunn_loc_nob_ST_pvalue <- dunn_loc_nob$P.adjusted[6]
dunn_loc_nob_BY_pvalue <- dunn_loc_nob$P.adjusted[7]
dunn_loc_nob_RY_pvalue <- dunn_loc_nob$P.adjusted[8]
dunn_loc_nob_SY_pvalue <- dunn_loc_nob$P.adjusted[9]
dunn_loc_nob_TY_pvalue <- dunn_loc_nob$P.adjusted[10]

#{{boxplot for functional groups for locality}}
{
  fact_l <- "locality"
  
  box_leg_nfix_locality <- bact_leg_Nfixers_rarefied_t %>%
    mutate(Leg_N_fixers = `Allorhizobium-Neorhizobium-Pararhizobium-Rhizobium` + 
             Bradyrhizobium +
             #Ensifer +
             Mesorhizobium +
             Microvirga) %>%
    select(soil_code, Leg_N_fixers) %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{plot}}
    ggplot(aes(x = get(fact_l), y = Leg_N_fixers)) +
    geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.1),
                alpha = 1/5, cex = 2.5) +
    theme_grey(base_size = 12) + 
    scale_x_discrete(limits = flevels_locality,
                     labels = c("Raguitenga\n(Korsimoro)",
                                "Tingressene\n(Korsimoro)",
                                "Boussouma\n(Boussouma)", 
                                "Sera\n(Boussouma)",
                                "Yilou\n(Guibare)")) +
    scale_fill_manual(limits = flevels_locality,
                      values = colpalette_locality) +
    scale_y_log10(limits = c(1, 1000)) +
    xlab("") + 
    ylab("Relative adundance (log10[reads])") +
    labs(fill = "Locality") +
    labs(title =  "Legume N-fixers", subtitle = krus_loc_leg_annotation) +
    theme(plot.title = element_text(color = colpalette_Ncycling[1], hjust = 0, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    annotate("text", x = 3.5, y = 20,
             label = paste("Dunn test P = ", round(dunn_loc_leg_BS_pvalue, digits = 5), sep = ""),
                     fontface = "italic", size = 3) +
    geom_segment(aes(x = 3.1, xend = 3.9, y = 30, yend = 30))
    
    
  box_act_nfix_locality <- bact_act_Nfixers_rarefied_t %>%
    mutate(Act_N_fixers = Frankia) %>%
    select(soil_code, Act_N_fixers) %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{plot}}
    ggplot(aes(x = get(fact_l), y = Act_N_fixers)) +
    geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.1),
                alpha = 1/5, cex = 2.5) +
    theme_grey(base_size = 12) + 
    scale_x_discrete(limits = flevels_locality,
                     labels = c("Raguitenga\n(Korsimoro)", 
                                "Tingressene\n(Korsimoro)",
                                "Boussouma\n(Boussouma)", 
                                "Sera\n(Boussouma)",
                                "Yilou\n(Guibare)")) +
    scale_fill_manual(limits = flevels_locality,
                      values = colpalette_locality) +
    scale_y_log10(limits = c(1, 1000)) +
    xlab("") + 
    ylab("Relative adundance (log10[reads])") +
    labs(fill = "Locality") +
    labs(title =  "Actinorhizal N-fixers", subtitle = krus_loc_act_annotation) +
    theme(plot.title = element_text(color = colpalette_Ncycling[2], hjust = 0, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0))
    #annotate("text", x = 3.5, y = 20,
    #         label = paste("Dunn test P = ", round(dunn_loc_act_BS_pvalue, digits = 5), sep = ""),
    #         fontface = "italic", size = 3) +
    #geom_segment(aes(x = 3.1, xend = 3.9, y = 30, yend = 30))
  
  box_aob_locality <- bact_AOB_rarefied_t %>%
    mutate(AOB_nitrifiers = #Nitrosococcus + 
             Nitrosomonas + Nitrosospira) %>%
    select(soil_code, AOB_nitrifiers) %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{plot}}
    ggplot(aes(x = get(fact_l), y = AOB_nitrifiers)) +
    geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.1),
                alpha = 1/5, cex = 2.5) +
    theme_grey(base_size = 12) + 
    scale_x_discrete(limits = flevels_locality,
                     labels = c("Raguitenga\n(Korsimoro)", 
                                "Tingressene\n(Korsimoro)",
                                "Boussouma\n(Boussouma)", 
                                "Sera\n(Boussouma)",
                                "Yilou\n(Guibare)")) +
    scale_fill_manual(limits = flevels_locality,
                      values = colpalette_locality) +
    scale_y_log10(limits = c(1, 1000)) +
    xlab("") + 
    ylab("Relative adundance (log10[reads])") +
    labs(fill = "Locality") +
    labs(title =  "AOB nitrifiers", subtitle = krus_loc_aob_annotation) +
    theme(plot.title = element_text(color = colpalette_Ncycling[3], hjust = 0, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
    annotate("text", x = 3.5, y = 30,
           label = paste("Dunn test P = ", round(dunn_loc_aob_BS_pvalue, digits = 5), sep = ""),
           fontface = "italic", size = 3) +
    geom_segment(aes(x = 3.1, xend = 3.9, y = 20, yend = 20))
  
  box_nob_locality <- bact_NOB_rarefied_t %>%
    mutate(NOB_nitrifiers = #Nitrobacter + Nitrococcus
             Nitrospira) %>%
    select(soil_code, NOB_nitrifiers) %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    #{{plot}}
    ggplot(aes(x = get(fact_l), y = NOB_nitrifiers)) +
    geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
    geom_jitter(position = position_jitter(0.1),
                alpha = 1/5, cex = 2.5) +
    theme_grey(base_size = 12) + 
    scale_x_discrete(limits = flevels_locality,
                     labels = c("Raguitenga\n(Korsimoro)", 
                                "Tingressene\n(Korsimoro)",
                                "Boussouma\n(Boussouma)", 
                                "Sera\n(Boussouma)",
                                "Yilou\n(Guibare)")) +
    scale_fill_manual(limits = flevels_locality,
                      values = colpalette_locality) +
    scale_y_log10(limits = c(1, 1000)) +
    xlab("") + 
    ylab("Relative adundance (log10[reads])") +
    labs(fill = "Locality") +
    labs(title =  "NOB nitrifiers", subtitle = krus_loc_nob_annotation) +
    theme(plot.title = element_text(color = colpalette_Ncycling[4], hjust = 0, size = 16, face = "bold"),
          plot.subtitle = element_text(size = 12),
          legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
          axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) 
    #annotate("text", x = 3.5, y = 20,
    #       label = paste("Dunn test P = ", round(dunn_loc_aob_BS_pvalue, digits = 5), sep = ""),
    #       fontface = "italic", size = 3) +
    #geom_segment(aes(x = 3.1, xend = 3.9, y = 30, yend = 30))
  }

#{{assembly of boxplots}}
{
  (box_leg_nfix_locality / box_act_nfix_locality / box_aob_locality / box_nob_locality) + 
    #plot_annotation(tag_levels = 'A') +
    plot_layout(guides = "collect")
  
  width <- 9
  height <- 12
  
  ggsave(file = output, width = width , height = height)
}
  
## -------------------- 14.1.2.1.4. Multiple plots ---- 
## -------- 14.2. fungi (soils) ----
## ---- 17. ALPHADIVERSITY + VARIANCE ANALYSIS ----
# p-values and fdr (http://www.nonlinear.com/support/progenesis/comet/faq/v2.0/pq-values.aspx)
# p-values and fdr (https://stats.stackexchange.com/questions/198073/correction-for-multiple-testing-on-a-modest-number-of-tests-10-20-with-fdr)
# p-values and fdr http://www.biostathandbook.com/multiplecomparisons.html
# p-values and fdr http://rcompanion.org/rcompanion/f_01.html

## -------- 17.1. bacteria (soils) --------
## ------------ 17.1.1. normality test {normtest} ------------
#{{Shapiro-Francia test (shapiro.test); Anderson–Darling test (ad.test); Cramer–von Mises test (cvm.test); Lilliefors test (lillie.test); Pearson chi-squared test for the composite hypothesis of normality (pearson.test)}}

#{{select the data}}
tmp <- Indices_bact %>%
  #{{merge with metadata}}
  left_join(metadata, by = "soil_code") %>%
  droplevels()

#{{if Shapiro-Francia test with p<0.05, data need to be transform to reach normality}}
set.seed(seed)
shapiro.test(tmp$Hill_Richness)
hist(tmp$Hill_Richness)
ggqqplot(tmp$Shannon, ylab = "Richness")

#{{select tha data transformation, sqrt(temp),sqrt(sqrt(tmp)), log1p(tmp),log10(tmp), asin(sqrt(tmp/100))}}
set.seed(seed)
shapiro.test(log1p(tmp$Hill_Richness))
hist(log1p(tmp$Hill_Richness))
ggqqplot(log1p(tmp$Hill_Richness), ylab = "Number of OTUs")

## ------------ 17.1.2. anova and tukey tests {stats}---- 
## ------------ 17.1.3. kruskall Wallis {RVAideMemoire} and Post-hoc Dunn tests {FSA} ---- 

#{{select the global data}}
{
   tmp <- Indices_bact %>%
    #{{merge with metadata}}
    left_join(metadata, by = "soil_code") %>%
    droplevels()
  }
  
#{{locality effect}}
{ 
  #{{perform kruskal test and dunn posthoc test for unequal numbers of observations}}
  #{{two groups}}
  
  krus_loc_dR <- kruskal.test(Hill_Richness~locality,
                               data = tmp, na.action = na.omit)
  
  krus_loc_dREr <- kruskal.test(Hill_Shannon~locality,
                                 data = tmp, na.action = na.omit)
  
  krus_loc_dREa <- kruskal.test(Hill_Inv_Simpson~locality,
                                 data = tmp, na.action = na.omit)
  
  krus_loc_eDRr <- kruskal.test(Hill_Shannon_evenness~locality,
                                 data = tmp, na.action = na.omit)
  
  krus_loc_eDRa <- kruskal.test(Hill_Simpson_evenness~locality,
                                 data = tmp, na.action = na.omit)
  
  #{{multi groups}}
  dunn_loc_dR <- dunn.test::dunn.test(tmp$Hill_Richness, tmp$locality,
                                       method = "bh") 
  
  dunn_loc_dREr <- dunn.test::dunn.test(tmp$Hill_Shannon, tmp$locality,
                                       method = "bh")
  
  dunn_loc_dREa <- dunn.test::dunn.test(tmp$Hill_Inv_Simpson, tmp$locality,
                                       method = "bh")
  
  dunn_loc_eDRr <- dunn.test::dunn.test(tmp$Hill_Shannon_evenness, tmp$locality,
                                       method = "bh")
  
  dunn_loc_eDRa <- dunn.test::dunn.test(tmp$Hill_Simpson_evenness, tmp$locality,
                                       method = "bh")
  
}  

## ------------ 17.1.4. boxplot ------------

#{{select and re-order the factors}}
flevels_locality <- c("Raguitenga", "Tingressene", "Boussouma", "Sera", "Yilou")

#{{jco palette visualisation}}
show_col(pal_jco("default")(10))
#{{select color palette}}
colpalette_locality <- c("#0073C2FF","#7AA6DCFF", "#A73030FF", "#CD534CFF", "#EFC000FF")


#{{format kruskall statistic}}
krus_loc_dR_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_dR$p.value, digits = 4),
                                 sep = "")

krus_loc_dREr_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_dREr$p.value, digits = 4),
                                 sep = "")

krus_loc_dREa_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_dREa$p.value, digits = 4),
                                 sep = "")

krus_loc_eDRr_annotation <- paste("Kruskal-Wallis rank sum test",
                                 "\nLocality P = ", round(krus_loc_eDRr$p.value, digits = 4),
                                 sep = "")
krus_loc_eDRa_annotation <- paste("Kruskal-Wallis rank sum test",
                                  "\nLocality P = ", round(krus_loc_eDRa$p.value, digits = 4),
                                  sep = "")

#{{select file name}}
output <- "Boxplot_bact_diversity_locality.pdf"
{
  #{{boxplot for locality effect}}
  {
    fact_l <- "locality"
    
    box_loc_dR <- Indices_bact %>%
      #{{merge with metadata}}
      left_join(metadata, by = "soil_code") %>%
      ggplot(aes(x = get(fact_l), y = Hill_Richness)) +
      geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.1),
                  alpha = 1/5, cex = 2.5) +
      theme_grey(base_size = 12) +
      ylim(c(0,6000)) +
      scale_x_discrete(limits = flevels_locality,
                       labels = c("Raguitenga\n(Korsimoro)",
                                  "Tingressene\n(Korsimoro)",
                                  "Boussouma\n(Boussouma)", 
                                  "Sera\n(Boussouma)",
                                  "Yilou\n(Guibare)")) +
      scale_fill_manual(limits = flevels_locality,
                        values = colpalette_locality) +
      labs(fill = "Locality") +
      labs(title =  "", subtitle = krus_loc_dR_annotation) +
      theme(plot.title = element_text(hjust = 0, size = 16),
            plot.subtitle = element_text(hjust = 0, size = 10),
            legend.position = "right",
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
      xlab("") + 
      ylab("Richness\n(q = 0)")
    
    box_loc_dREr <- Indices_bact %>%
      #{{merge with metadata}}
      left_join(metadata, by = "soil_code") %>%
      ggplot(aes(x = get(fact_l), y = Hill_Shannon)) +
      geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.1),
                  alpha = 1/5, cex = 2.5) +
      theme_grey(base_size = 12) +
      ylim(c(0,3000)) +
      scale_x_discrete(limits = flevels_locality,
                       labels = c("Raguitenga\n(Korsimoro)",
                                  "Tingressene\n(Korsimoro)",
                                  "Boussouma\n(Boussouma)", 
                                  "Sera\n(Boussouma)",
                                  "Yilou\n(Guibare)")) +
      scale_fill_manual(limits = flevels_locality,
                        values = colpalette_locality) +
      labs(fill = "Locality") +
      labs(title =  "", subtitle = krus_loc_dREr_annotation) +
      theme(plot.title = element_text(hjust = 0, size = 16),
            plot.subtitle = element_text(hjust = 0, size = 10),
            legend.position = "right",
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
      xlab("") +
      ylab("Shannon diversity\n(q = 1)")
    
    box_loc_dREa <- Indices_bact %>%
      #{{merge with metadata}}
      left_join(metadata, by = "soil_code") %>%
      ggplot(aes(x = get(fact_l), y = Hill_Inv_Simpson)) +
      geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.1),
                  alpha = 1/5, cex = 2.5) +
      theme_grey(base_size = 12) +
      ylim(c(0,3000)) +
      scale_x_discrete(limits = flevels_locality,
                       labels = c("Raguitenga\n(Korsimoro)",
                                  "Tingressene\n(Korsimoro)",
                                  "Boussouma\n(Boussouma)", 
                                  "Sera\n(Boussouma)",
                                  "Yilou\n(Guibare)")) +
      scale_fill_manual(limits = flevels_locality,
                        values = colpalette_locality) +
      labs(fill = "Locality") +
      labs(title =  "", subtitle = krus_loc_dREa_annotation) +
      theme(plot.title = element_text(hjust = 0, size = 16),
            plot.subtitle = element_text(hjust = 0, size = 10),
            legend.position = "right",
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
      xlab("") + 
      ylab("Simpson diversity\n(q = 2)")
    
    box_loc_eDRr <- Indices_bact %>%
      #{{merge with metadata}}
      left_join(metadata, by = "soil_code") %>%
      ggplot(aes(x = get(fact_l), y = Hill_Shannon_evenness)) +
      geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.1),
                  alpha = 1/5, cex = 2.5) +
      theme_grey(base_size = 12) +
      ylim(c(0,0.6)) +
      scale_x_discrete(limits = flevels_locality,
                       labels = c("Raguitenga\n(Korsimoro)",
                                  "Tingressene\n(Korsimoro)",
                                  "Boussouma\n(Boussouma)", 
                                  "Sera\n(Boussouma)",
                                  "Yilou\n(Guibare)")) +
      scale_fill_manual(limits = flevels_locality,
                        values = colpalette_locality) +
      labs(fill = "Locality") +
      labs(title =  "", subtitle = krus_loc_eDRr_annotation) +
      theme(plot.title = element_text(hjust = 0, size = 16),
            plot.subtitle = element_text(hjust = 0, size = 10),
            legend.position = "right",
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
      xlab("") + 
      ylab("Evenness\nShannon-based")
    
    box_loc_eDRa <- Indices_bact %>%
      #{{merge with metadata}}
      left_join(metadata, by = "soil_code") %>%
      ggplot(aes(x = get(fact_l), y = Hill_Simpson_evenness)) +
      geom_boxplot(aes(fill = get(fact_l)), outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.1),
                  alpha = 1/5, cex = 2.5) +
      theme_grey(base_size = 12) +
      ylim(c(0,0.6)) +
      scale_x_discrete(limits = flevels_locality,
                       labels = c("Raguitenga\n(Korsimoro)",
                                  "Tingressene\n(Korsimoro)",
                                  "Boussouma\n(Boussouma)", 
                                  "Sera\n(Boussouma)",
                                  "Yilou\n(Guibare)")) +
      scale_fill_manual(limits = flevels_locality,
                        values = colpalette_locality) +
      labs(fill = "Locality") +
      labs(title =  "", subtitle = krus_loc_eDRa_annotation) +
      theme(plot.title = element_text(hjust = 0, size = 16),
            plot.subtitle = element_text(hjust = 0, size = 10),
            legend.position = "right",
            axis.text.x = element_text(size = 12, angle = 0, hjust = 0.5, vjust = 0.5),
            axis.text.y = element_text(size = 12, angle = 0, hjust = 1, vjust = 0)) +
      xlab("") + 
      ylab("Evenness\nSimpson-based")
    }
  
  #{{assembly of boxplots}}
  
  (box_loc_dR / box_loc_dREr / box_loc_dREa / box_loc_eDRr / box_loc_eDRa) + 
    plot_layout(guides = "auto") & theme(legend.position = "none")
  
  width <- 9
  height <- 15
  
  ggsave(file = output, width = width , height = height)
  
}
