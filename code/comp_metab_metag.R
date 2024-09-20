# ---- Packages ----
library(tidyverse)

main_theme = theme_minimal()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22, face="italic", angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(colour = "black", size=22, face="italic"),
        legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),
        
        legend.text = element_text(colour = "black", size=18),
        axis.title= element_text(size=28),
        strip.text = element_text(colour = "black", size=15, face ="italic"))

# ---- Load data ----

metab1 = "data/metabarcoding_vs_metagenomics/oracle_metaBt0_16S_samples_run_20190715_kraken2_assignment_genuses.tsv"
metab2 = "data/metabarcoding_vs_metagenomics/oracle_metaBt0_16S_samples_run_20200106_kraken2_assignment_genuses.tsv"
metag1 = "data/metabarcoding_vs_metagenomics/oracle_metaGt0_SAMA_12_samples_kraken2_SSU_assignment_genuses.tsv"
metag2 = "data/metabarcoding_vs_metagenomics/oracle_metaGt0_SAMA_21_1_samples_kraken2_SSU_assignment_genuses.tsv"
metag3 = "data/metabarcoding_vs_metagenomics/oracle_metaGt0_SAMA_21_2_samples_kraken2_SSU_assignment_genuses.tsv"

taxonomic_levels <- c("kingdom", "phylum", "class", "order", "family", 
                      "genus", "species")
## ---- Arrange metaB ----

# Select dat from bracken method, and rename columns 
metab1%>%
  read_tsv(col_names = T, skip = 1)%>%
  select(which(str_detect(colnames(.),"bracken_genuses")), "#OTU ID", "taxonomy")%>%
  rename_with(
    .%>%
    str_remove( "_bracken_genuses")%>%
    str_replace_all("^OR-", "BU_") %>% # Replace "OR-" with "BU_"
    str_replace_all("-", "_")%>% 
    str_remove("_S\\d+")%>%
    str_replace("((?<!.)T_)(\\d)", "CONT_BU_PCR_\\2")%>%
    str_replace("(BU_T_extr)", "CONT_BU_ext")%>%
    str_replace("#OTU ID", "OTU"))%>%
  separate_wider_delim(taxonomy, names = taxonomic_levels, delim = "; ", cols_remove = F)%>%
  filter(kingdom == "k__Bacteria")-> metab1_table
metab1_table$taxonomy
# Disjoint R1 and R2 and decontaminate

metab1_table%>%
  select(!ends_with("R2"))%>%
  rename_with(
   .%>% str_remove("_R1")
  )%>%
  replace(. == 0, NA) %>%
  pivot_longer(starts_with("BU"), names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads)) %>%
  #{{merge control samples}}
  mutate( n = rowSums(across(starts_with("CONT")), na.rm = T))%>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  pivot_wider(values_from = reads, names_from = samples, values_fill = 0) %>%
  select(-starts_with("CONT"))%>%
  select(order(colnames(.)))-> metab1_table_decont_R1

metab1_table%>%
  select(!ends_with("R1"))%>%
  rename_with(
    .%>% str_remove("_R2")
  )%>%
  replace(. == 0, NA) %>%
  pivot_longer(starts_with("BU"), names_to = "samples", values_to = "reads") %>%
  filter(!is.na(reads)) %>%
  #{{merge control samples}}
  mutate( n = rowSums(across(starts_with("CONT")), na.rm = T))%>%
  #{{substract abundance of control samples}}
  mutate(reads = case_when(
    is.na(n)  ~ reads,
    n > reads ~ 0,
    TRUE      ~ reads - n)) %>%
  select(-n) %>%
  pivot_wider(values_from = reads, names_from = samples, values_fill = 0) %>%
  select(-starts_with("CONT"))%>%
  select(order(colnames(.)))-> metab1_table_decont_R2
metab1_table_decont_R2%>%
  rename_with(
    .%>% str_c("_R2"),
    starts_with("BU")
  )%>%
  full_join(metab1_table_decont_R1%>%
              rename_with(
                .%>% str_c("_R1"),
                starts_with("BU")
              )%>%
              select(starts_with("BU"), OTU), by = "OTU")->metab1_table_decont
metab1_table_decont%>%colnames
## ---- analysis ----

### ---- coinertie ----
#Compaire 
library(ade4)
data(doubs)
colnames(metab1_table_decont_R1)
all(colnames(metab1_table_decont_R2) == colnames(metab1_table_decont_R1)) # need to be true for the next step

metab1_table_decont_R1%>%
 bind_rows(metab1_table_decont_R2%>%
         filter(!(OTU %in% metab1_table_decont_R1$OTU)))%>%
  arrange(OTU)-> metab1_table_decont_R1

metab1_table_decont_R2%>%
  bind_rows(metab1_table_decont_R1%>%
              filter(!(OTU %in% metab1_table_decont_R2$OTU)))%>%
  arrange(OTU) -> metab1_table_decont_R2

metab1_table_decont_R1%>%
  column_to_rownames("OTU")%>%
  select(starts_with("BU"))%>%
  decostand(method = "hellinger")%>%
  t()%>%
  dudi.pca( scale = TRUE, scan = FALSE, nf = 3) -> dudi1

metab1_table_decont_R2%>%
  column_to_rownames("OTU")%>%
  select(starts_with("BU"))%>%
  decostand(method = "hellinger")%>%
  t()%>%
  dudi.pca(scale = FALSE, scan = FALSE, nf = 2) -> dudi2

coin1 <- coinertia(dudi1,dudi2, scan = FALSE, nf = 2)

summary(coin1)
plot(coin1)
metab1_table_decont_R2%>%
  filter(OTU == 1688)
metab1_table_decont_R1%>%
  filter(OTU == 1688)

metab1_table_decont_R1%>%
  column_to_rownames("OTU")%>%
  select(starts_with("BU"))%>%
  decostand(method = "hellinger")%>%
  dudi.pca( scale = TRUE, scan = FALSE, nf = 3) -> dudi3

metab1_table_decont_R2%>%
  column_to_rownames("OTU")%>%
  select(starts_with("BU"))%>%
  decostand(method = "hellinger")%>%
  dudi.pca(scale = FALSE, scan = FALSE, nf = 2) -> dudi4

coin2 <- coinertia(dudi3,dudi4, scan = FALSE, nf = 2)
summary(coin2)
coin2$mX
coin2$mY
library(ggrepel)
compute.dist = function(df1, df2){
  data.frame(ID = as.numeric(rownames(df1)) , dist = sqrt(rowSums((df1-df2)^2)))
}
compute.dist(coin2$mX, coin2$mY) %>%
  arrange(desc(dist)) -> dist_sp_df
plot(coin2)
cbind(coin2$mX,coin2$mY) ->  df_coin_pos
colnames(df_coin_pos)=c("R1x", "R1y", "R2x", "R2y")

df_coin_pos%>%
  rownames_to_column("OTU")%>%
  mutate(OTU = as.numeric(OTU))%>%
  left_join(metab1_table_decont_R2%>%select(OTU, genus), by = "OTU")%>%
  left_join(dist_sp_df, by = join_by(OTU == ID))%>%
  arrange(desc(dist)) %>%
  slice(c(1:10, (n() - 9):n()))%>%
  mutate(diff = c(rep("Strong",10 ), rep("Weak", 10)))%>%
  
  ggplot()+
  facet_wrap(~diff)+
  geom_label_repel(aes(x= R1x, y = R1y ,label = genus, col =dist))+
  geom_segment(aes(x = R1x, y = R1y, xend = R2x, yend = R2y, col =dist),
               arrow = arrow(length = unit(0.1, "cm"), type = "closed"), cex =1)+
  
  labs(x="X",y = "Y", title = "Genera estimations differences between R1 (Tail) and R2 (Head) across 80 samples (for metaB/bracken)")+
  main_theme

# Exemple of g__Gordonia for strong
metab1_table_decont_R2%>%
  filter(genus == "g__Gordonia")%>%pull(taxonomy)
metab1_table_decont_R1%>%
  filter(genus == "g__Gordonia")

# Exemple of g__Asidisphaera for weak
metab1_table_decont_R2%>%
  filter(genus == "g__Acidisphaera")
metab1_table_decont_R1%>%
  filter(genus == "g__Acidisphaera")

### ---- metacoder ----
library(metacoder)

  

hmp_otus%>%colnames
hmp_otus$otu_id
metab1_table_decont_R1%>%
  mutate(
    taxonomy = str_remove(taxonomy, "; s__"),
    taxonomy = str_replace_all(taxonomy,"(k__|p__|c__|o__|f__|g__)(;|$)","\\1none\\2"))%>%
  select(OTU, taxonomy,starts_with("BU"))%>%
  mutate(across(where(is.numeric), ~as.integer(.x)))%>%
  parse_tax_data(class_cols = "taxonomy", # the column that contains taxonomic information
                      class_sep = "; ", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name")) -> obj1

metab1_table_decont = metab1_table_decont_R1%>%
  mutate(type = "R1")%>%
  bind_rows(metab1_table_decont_R1%>%
              mutate(type = "R2"))
  
  
metab1_table$R1_sum <- rowSums(metab1_table[, grep("R1", colnames(metab1_table))])
metab1_table$R2_sum <- rowSums(metab1_table[, grep("R2", colnames(metab1_table))])
obj1 <- parse_tax_data(metab1_table_decont,
                      class_cols = "taxonomy", # The column in the input table
                      class_sep = "; ",
                      class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                      class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))

obj2 <- parse_tax_data(metab1_table_decont_R2,
                       class_cols = "taxonomy", # The column in the input table
                       class_sep = "; ",
                       class_regex = "^([a-z]{0,1})_{0,2}(.*)$",
                       class_key = c("tax_rank" = "taxon_rank", "name" = "taxon_name"))
obj1%>%
  filter_taxa(taxon_ranks == "", supertaxa = TRUE) -> obj1_temp 
obj2%>%
  filter_taxa(taxon_ranks == "g", supertaxa = TRUE) -> obj2_temp 
obj1$data$tax_data$ty

obj1$data$tax_abund <- calc_taxon_abund(obj1, "tax_data",
                                       cols = startsWith(colnames(obj1$data$tax_data),"BU"),
                                                         groups = str_extract(str_subset(colnames(obj1$data$tax_data), "_R\\d"),"R\\d"))

obj2_temp$n_obs()%>%length()
obj1_temp$n_obs()%>%length()
obj1%>%
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE)%>%
  heat_tree(node_label = gsub(pattern = "\\[|\\]", replacement = "", taxon_names),
            node_size = n_obs,
            node_color = R1-R2,
            node_color_range = c("blue", "gray", "red"),
            node_color_axis_label = "OTU count",
            layout = "davidson-harel", initial_layout = "reingold-tilford")

metab1_table%>%
 
  select(OTU, taxonomy,starts_with("BU"))%>%
  parse_tax_data(class_cols = "taxonomy", # the column that contains taxonomic information
                 class_sep = "; ", # The character used to separate taxa in the classification
                 class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                 class_key = c(tax_rank = "info", # A key describing each regex capture group
                               tax_name = "taxon_name"),
                 datasets = list(R1 = metab1_table$R1_sum, 
                                 R2 = metab1_table$R2_sum),
                 mapping = c("{{index}}" = "{{index}}",
                             "{{index}}" = "{{index}}")) -> obj2
obj2$data$differences <- abs(metab1_table$R1_sum - metab1_table$R2_sum)
obj2$data$tax_data$taxon_id%>%length()
obj2$data$differences %>%length()
heat_tree(obj2,
          node_label = taxon_names,
          
          node_color_range = c("blue", "gray", "red"),  # blue for R2, red for R1
          node_color_trans = "linear",
          node_color_interval = c(-max(abs(obj2$data$differences), na.rm = TRUE), 
                                  max(abs(obj2$data$differences), na.rm = TRUE)),
          edge_color_interval = c(-max(abs(obj2$data$differences), na.rm = TRUE), 
                                  max(abs(obj2$data$differences), na.rm = TRUE)),
          edge_color_trans = "linear",
          title = "Differences Between R1 and R2 Methods",
          node_color_axis_label = "Difference (R1 - R2)")

set.seed(1) # This makes the plot appear the same each time it is run 
obj %>%
  filter_taxa(grepl(pattern = "^[a-zA-Z]+$", taxon_names))%>% # remove "odd" taxa
  filter_taxa(taxon_ranks == "o", supertaxa = TRUE)
mutate_obs(obj2, "info",
           new_col = "Im new",
           newer_col = paste0(new_col, "er!"))
filter_obs(obj2, target = ends_with("R1"), ! no_reads, drop_taxa = TRUE)

calc_obs_props(obj, "tax_data")
calc_n_samples(obj1, "tax_abund", groups = hmp_samples$body_site, cols = hmp_samples$sample_id)
heat_tree(obj1, 
          node_label = taxon_names,
          node_size = n_obs,
          node_color = R1 -R2, 
          node_size_axis_label = "OTU count",
          node_color_axis_label = "Samples with reads",
          layout = "davidson-harel", initial_layout = "reingold-tilford")

  