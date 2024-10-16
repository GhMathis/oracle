## ---- IMPORT PACKAGES ----
library(readxl)
library(tidyverse)
library(iNEXT)
library(GGally)
library(sbm)
library(FactoMineR)
library(vegan)
library(alluvial) #aluvial plot
library(MuMIn) # model selection
library(aricode) # NMI

library(factoextra)
library(corrplot)
library(RColorBrewer)

## ---- usefull function ----
main_theme = theme_minimal()+
  theme(line = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.border = element_blank(),
        axis.ticks =  element_line(colour = "black"),
        axis.text.x = element_text(colour = "black", size=22, face="italic", angle = 45,),
        axis.text.y = element_text(colour = "black", size=22, face="italic"),
        legend.title = element_text(colour = "black", size=20,
                                    hjust =0.5),

        legend.text = element_text(colour = "black", size=18),
        axis.title= element_text(size=28),
        strip.text = element_text(colour = "black", size=15, face ="italic"))

zscore_normalization <- function(data) {
  # Step 1: Rank the data
  ranked_data <- rank(data)
  n <- length(data)
  
  # Step 2: Calculate quantile positions
  quantile_positions <- (ranked_data) / (n + 1)
  
  # Step 3: Transform quantiles to z-scores
  mean_std_normal <- qnorm(quantile_positions)
  
  return(mean_std_normal)
}
## ---- IMPORT DATA ----

read_tsv("data/metadata.txt") %>%
  rename_with(~gsub(".", "_", ., fixed = TRUE)) -> metadata
read_tsv("data/sample_coordinates.tsv")-> coord
read_tsv("data/divesities_bact.tsv")%>%
  dplyr::select(starts_with("Hill_"),soil_code)-> div_bact

read_xlsx("data/Historique_Parcelles-Oracle.xlsx")%>%
  dplyr::select(`N° Parcelle`, ...14, ...17,...20,) %>%
  rename(soil_code = "N° Parcelle", cultu_2014 = "...14", cultu_2015 = "...17",cultu_2016 = "...20")%>%
  mutate(soil_code = paste0("BU_",str_replace(soil_code, "([A-Z]+)(\\d+)", "\\1_\\2"))) -> preced_cultu

read_xlsx("data/Soil_Parameters_oracle.xlsx")%>%
  dplyr::select(`Var niébé 2017`, `Var Sorgho 2017`)%>%
  rename(varnieb = "Var niébé 2017", varsorgh = "Var Sorgho 2017")-> variety

metadata%>%
  dplyr::select(pracul,modsem,FO_Qte ,denpoq, freqsar, NPK_qte, uree_qte)%>%
  mutate(pracul = as.factor(pracul),
                 modsem = as.factor(modsem))->agri_practice

agri_practice%>%  
  mutate(across(where(is.numeric),~decostand(.x,"standardize")[,1]))->agri_practice_standa

read_xlsx("data/Soil_Parameters_oracle.xlsx")%>%
  dplyr::select(`C (g kg-1)`,`C/N`, `CEC (cmolc/kg)`, `Ca (cmol/kg)`,`Mg (cmolc/kg)`, `K(cmolc/kg)`, `SBE (cmolc/kg)`,           
    `TS (%)`, Site...3, `N° site...4`)%>%
  rename(Site = "Site...3",N_site = "N° site...4", C = "C (g kg-1)",CN = "C/N", CEC = "CEC (cmolc/kg)", Ca = "Ca (cmol/kg)",Mg = "Mg (cmolc/kg)",K = "K(cmolc/kg)",SBE = "SBE (cmolc/kg)",           
         TS = "TS (%)" )%>%
  mutate( N_site = str_replace(as.character(N_site), "^[1-9]$", "0\\0"),
          soil_code = str_c(Site,N_site))%>%
  mutate(soil_code = str_c("BU_",Site,"_",N_site))%>%
  left_join(metadata%>%dplyr::select(pHeau,Pto,Kto,Nto,MO,soil_code), by = "soil_code")-> soil


read.table("data/bact_OTUs_soils_final.txt",header = T)%>%as_tibble -> otu_bact



## ---- CLUSTERING ----
### ---- Crops history ----


### rm NA values
preced_cultu%>%
  filter(is.na(cultu_2014) | is.na(cultu_2015) | is.na(cultu_2016))
preced_cultu%>%
  filter(!(is.na(cultu_2014) | is.na(cultu_2015) | is.na(cultu_2016)))%>%
  filter(cultu_2014 != "Ass_Ses+N") -> preced_cultu # rm NA and BU_SE_R12 because only one replicate of sesame culture

### disjonctif table usable by sbm
tab.disjonctif(preced_cultu$cultu_2014) -> cult2014
tab.disjonctif(preced_cultu$cultu_2015) -> cult2015
tab.disjonctif(preced_cultu$cultu_2016) -> cult2016

### check if colnames concord with eachother
colnames(cult2014)
colnames(cult2015)
colnames(cult2016)

### visualisation of one tab
tab.disjonctif(preced_cultu$cultu_2014)%>%
  t()%>%
  plotMyMatrix(dimLabels = c("cult2014", "soil code"), plotOptions = list(rowNames = T,colNames = F))

### create networks objects
net2014 <- defineSBM(cult2014, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))
net2015 <- defineSBM(cult2015, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))
net2016 <- defineSBM(cult2016, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))

### compute smb
estimateMultiplexSBM(list(net2014, net2015, net2016), dependent = FALSE) -> sbm_preced_cult 

plot(sbm_preced_cult)



cluster_df <- data.frame(clust_cult_history = sbm_preced_cult$memberships$cult_cult_history, soil_code = preced_cultu$soil_code)

#### Or #
read_xlsx("data/Historique_Parcelles-Oracle.xlsx")%>%
  dplyr::select(`N° Parcelle`, ...14, ...17,...20,) %>%
  rename(soil_code = "N° Parcelle", cultu_2014 = "...14", cultu_2015 = "...17",cultu_2016 = "...20")%>%
  mutate(soil_code = paste0("BU_",str_replace(soil_code, "([A-Z]+)(\\d+)", "\\1_\\2"))) -> preced_cultu

preced_cultu%>%
  filter(!(is.na(cultu_2014) | is.na(cultu_2015) | is.na(cultu_2016))) -> preced_cultu
preced_cultu%>%
  mutate(across(starts_with("cultu"), ~case_when(
    .x == "Ass_Ses+N" ~ "c+l",
    .x == "Ass_M+N" ~ "c+l",
    .x == "Ass_S+N" ~ "c+l",
    .x == "Ass_S+M+N" ~ "c+l",
    .x == "Rot_Avec_Leg" ~ "l",
    .x == "Rot_Sans_Leg" ~ "c",
    .x == "Jachère" ~ "jachere"
  ))) -> preced_cultu_cere_legu

### disjonctif table usable by sbm
tab.disjonctif(preced_cultu_cere_legu$cultu_2014) -> cult2014_bis
tab.disjonctif(preced_cultu_cere_legu$cultu_2015) -> cult2015_bis
tab.disjonctif(preced_cultu_cere_legu$cultu_2016) -> cult2016_bis

### check if colnames concord with eachother
colnames(cult2014_bis)
colnames(cult2015_bis)
colnames(cult2016_bis)
str(cult2014_bis)
cult2014_bis = cbind(cult2014_bis,jachere = 0)

### visualisation of one tab
cult2014_bis%>%
  t()%>%
  plotMyMatrix(dimLabels = c("cult2014", "soil_code"), plotOptions = list(rowNames = T,colNames = F))

### create networks objects
net2014_bis <- defineSBM(cult2014_bis, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))
net2015_bis <- defineSBM(cult2015_bis, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))
net2016_bis <- defineSBM(cult2016_bis, model = "bernoulli", dimLabels = c(col = "cult_cult_history", row = "soil_code"))

### compute smb
estimateMultiplexSBM(list(net2014_bis, net2015_bis, net2016_bis), dependent = FALSE) -> sbm_preced_cult_bis

plot(sbm_preced_cult_bis)



cluster_df <- data.frame(clust_cult_history = sbm_preced_cult_bis$memberships$cult_cult_history, soil_code = preced_cultu$soil_code)

### ---- Soil ----


zscore_normalization

soil%>%
  select(-c(soil_code, N_site, Site))%>%
  mutate_all(~zscore_normalization(.))%>%
  as.matrix()-> soil_standard_rank

smb_soil <- estimateBipartiteSBM(soil_standard_rank, model = "gaussian")

grid_id = smb_soil$memberships[2]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(col)

soil_id = smb_soil$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(row)


soil_standard_rank[soil_id$id,grid_id$id]%>%
  as.matrix()%>%
  plotMyMatrix(dimLabels = c("Soil", "soil_code"), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 6, face="italic"),
    axis.text.y = element_text(size = 5,face="italic"))+
  geom_vline(xintercept = which(!duplicated(grid_id$col))[-1]-0.5, col ="red")+
  geom_hline(yintercept = abs(which(!duplicated(soil_id$row))[-1]-nrow(soil_id)-1)+0.5, col = "red", linetype = 5)+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))

cluster_df%>%
  left_join( data.frame(clust_soil = smb_soil$memberships$row,
                       soil_code = soil$soil_code ), by = "soil_code") -> cluster_df

### ---- Variety ----



tab.disjonctif(variety$varnieb) -> varnieb
tab.disjonctif(variety$varsorgh) -> varsorgh


### visualization of one tab
varnieb%>%
  t()%>%
  plotMyMatrix(dimLabels = c("varnieb", "soil code"), plotOptions = list(rowNames = T,colNames = F))

varsorgh%>%
  t()%>%
  plotMyMatrix(dimLabels = c("varsorgh", "soil code"), plotOptions = list(rowNames = T,colNames = F))

### create networks objects
type = "bipartite"
model = "bernoulli"
netvarnieb<- defineSBM(varnieb, model, type, directed, dimLabels = c("soil_code",
                                                                                    "varnieb"))
netvarsorgh <- defineSBM(varsorgh, model, type, directed, dimLabels = c("soil_code",
                                                                              "varsorgh"))

plotMyMultipartiteMatrix(list(netvarnieb, netvarnieb), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 6, face="italic"),
    axis.text.y = element_text(size = 5,face="italic"))

### compute smb
estimateMultipartiteSBM(list(netvarnieb, netvarsorgh), estimOptions = list(initBM = FALSE)) -> sbm_variety 

plot(sbm_variety)

cluster_df%>%
  left_join( data.frame(clust_variety = sbm_variety$memberships$soil_code,
                        soil_code = soil$soil_code ), by = "soil_code") -> cluster_df

### ---- Agricultural practice ----
library(factoextra)
agri_practice_standa%>%
  select(-c(pracul, modsem))%>%
  PCA%>%
  HCPC(-1) -> test_hpc_pca

agri_practice_standa%>%
  select(-c(pracul, modsem))%>%
  PCA%>%
  fviz_pca_biplot(col.ind = test_hpc_pca$data.clust$clust)

agri_practice_standa$denpoq
corrplot(cor(agri_practice_standa[-11,c(3,4,6,7)]))
MFA_agri_practice = MFA(agri_practice_standa, group = c(2,5), c("n","s"),
    name.group=c("technical itin","chemical itin"))
summary(MFA_agri_practice)
HCPC_practice = MFA_agri_practice%>%
  HCPC(-1)
str(as.factor(HCPC_practice$data.clust$clust))
fviz_screeplot(MFA_agri_practice)
fviz_mfa_var(MFA_agri_practice, "group")
fviz_mfa_quali_biplot(MFA_agri_practice,repel = FALSE, col.var = "#E7B800",
                      habillage = HCPC_practice$data.clust$clust, addEllipses = TRUE, ellipse.level = 0.95)+
  scale_colour_manual(values = c(brewer.pal(6, "Set1")))
fviz_mfa_ind(MFA_agri_practice,col.ind = HCPC_practice$data.clust$clust) 


grid_id2 = data.frame(col = test_hpc_pca$data.clust$clust)%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(col)

soil_id = smb_soil$memberships[1]%>%
  as.data.frame()%>%
  rownames_to_column("id")%>%
  mutate(id = as.numeric(id))%>%
  arrange(row)


agri_practice[grid_id2$id,-c(1,2)]%>%
  decostand("standardize")%>%
  mutate(across(where(is.factor), ~as.numeric(.x)))%>%
  as.matrix()%>%
  t()%>%
  plotMyMatrix(dimLabels = c("Practices", "soil_code"), plotOptions = list(rowNames = T,colNames = T))+
  theme(
    legend.position="none",
    axis.text.x =element_text(angle = 45, vjust = 1, hjust = 1, size = 6, face="italic"),
    axis.text.y = element_text(size = 5,face="italic"))+
  geom_vline(xintercept = which(!duplicated(grid_id2$col))[-1]-0.5, col ="red")+
  #geom_hline(yintercept = abs(which(!duplicated(soil_id$row))[-1]-nrow(soil_id)-1)+0.5, col = "red", linetype = 5)+
  theme(strip.text.x = element_text(colour = "gray", size=15, face ="italic"),
        strip.text.y  = element_text(colour = "gray", size=15, face ="italic"),
        strip.background = element_rect(fill="white"),
        panel.border = element_rect(color = "black",  fill = NA))
#### ---- a try with sbm----
agri_practice%>%
  select(where(is.numeric))%>%
  mutate_all(~zscore_normalization(.)) -> agri_practice_normrank

tab.disjonctif(agri_practice$pracul) -> pracul_dummy
tab.disjonctif(agri_practice$modsem) -> modsem_dummy

netpracul<- defineSBM(pracul_dummy, model, type, dimLabels = c("soil_code",
                                                                     "pracul"))
netmodsem <- defineSBM(modsem_dummy, model, type, dimLabels = c("soil_code",
                                                                        "modsem"))
netamendement <- defineSBM(as.matrix(agri_practice_normrank), "gaussian", type, dimLabels = c("soil_code",
                                                                "amendement"))
estimateBipartiteSBM(as.matrix(agri_practice_normrank), model = "gaussian")-> sbm_practice
  
estimateMultipartiteSBM(as.matrix(agri_practice_normrank), estimOptions = list(initBM = FALSE)) -> sbm_practice 

plot(sbm_practice)
### ---- Bacterial community ----

pca_bact = otu_bact%>%
  column_to_rownames("soil_code")%>%as_tibble%>%
  mutate(across(everything(), ~decostand(., method = "hellinger")[,1]))%>%
  as.matrix%>%
  PCA

otu_bact[,1:10]%>%
  
  left_join(soil, by = "soil_code")%>%
  pivot_longer(starts_with("X"), names_to = "bact_names", values_to = "abondances")%>%
  pivot_longer(-c(soil_code, bact_names, abondances,Site, N_site), names_to = "soil_names", values_to = "soil_values")%>%
  ggplot()+
  #geom_point(aes(soil_values, abondances, col = bact_names))+
  geom_smooth(aes(soil_values, log1p(abondances), col = bact_names), se = F)+
  facet_wrap(~soil_names, scales = "free")+
  main_theme
  
  
  
otu_bact%>%
  column_to_rownames("soil_code")%>%
  t()%>%
  cor()-> corr_site

otu_bact%>%
  select(soil_code)%>%
  mutate(locality = str_extract(soil_code, "(?<=_)[A-Z]+(?=_)"))%>%
  pull(locality) -> locality

corrplot::corrplot(corr_site, method = "square",i)

fviz_pca_ind(pca_bact,col.ind = locality)+
  scale_colour_manual(values = c(brewer.pal(5, "Set1")))

hist(log(corr_site))

HCPC_bact = pca_bact%>%
  HCPC(-1)
str(HCPC_bact)

fviz_pca_ind(pca_bact,col.ind = HCPC_bact$data.clust$clust)+
  scale_colour_manual(values = c(brewer.pal(4, "Set1")))+
  main_theme
hist(corr_site)
corr_site%>%
  #log%>%
  estimateSimpleSBM(
    model = 'gaussian', 
    dimLabels = "soil_code",
    estimOptions = list(nbCores = 8,plot = T)) -> sbm_bact
plot(sbm_bact)
sbm_bact$storedModels %>% knitr::kable()

sbm_bact$memberships

cluster_df%>%
  left_join( data.frame(clust_bact = sbm_bact$memberships,
                        clust_bact_HCPC = HCPC_bact$data.clust$clust,
                        soil_code = colnames(corr_site) ), by = "soil_code") -> cluster_df
### ---- test PLN ----

library(PLNmodels)
library(svd)
table(colSums(otu_bact[,-1]!=0))
# remove bactria that apear only one or 2 time (not usable in a statistical point of view)
otu_bact%>%
  column_to_rownames("soil_code")%>%
  select(where(~sum(.x)!=1 && sum(.x)!=2 ))-> otu_bact_clean

SVD = propack.svd(otu_bact[,-1]%>%as.matrix, neig= 7)
str(SVD)
trunc_d =  diag(sqrt(SVD$d))


str(L)
L =(SVD$u %*% trunc_d)
egein_V = SVD$d

R = trunc_d %*% t(SVD$v)
t_R = t(R)#[,1:12]
otu_bact[,-1]%>%tibble
L%*%R%>%tibble()
temp = list(abundance = otu_bact_clean,
            covariate = soil%>%
              select(-c(Site,  N_site))%>%
              arrange(soil_code)%>%
              column_to_rownames("soil_code"))
data_bact_pln = prepare_data(temp$abundance,temp$covariate)
pln_bact_null <- PLN(Abundance ~ 1,data_bact_pln)

## ---- Diversity model ----

cluster_df%>%
  mutate(locality = as.numeric(factor((str_extract(soil_code,"((?<=^.{3})[A-Z]+)")))))-> cluster_df

### visualization of cluster

B <-
  as.data.frame(
    table(
      
      cluster_df$soil_code,
      cluster_df$clust_cult_history ,
      cluster_df$clust_soil,
      cluster_df$clust_variety,
      cluster_df$clust_bact,
      cluster_df$locality,
      cluster_df$clust_bact_HCPC)
  )

colnames(B) =  c( "soil_code","clust_cult_history", "clust_soil",  "clust_variety", "clust_bact", "locality", "clust_bact_HCPC","Freq")

w   <- which(B$Freq != 0)
B <- B[w, ]

alluvial(B[, c(2,4,3)],
         freq = B$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)
alluvial(B[, c(7,2)],
         freq = B$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)
alluvial(B[, c(7,3)],
         freq = B$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)

alluvial(B[, c(7,4)],
         freq = B$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)

alluvial(B[, c(7,6)],
         freq = B$Freq,
         alpha = 0.7,
         cex = 1.5,
         cex.axis = 1.5,
         xw = 0.2)
NMI(c1 = cluster_df$clust_bact, c2 = cluster_df$clust_cult_history)
NMI(c1 = cluster_df$clust_bact, c2 = cluster_df$clust_soil)
NMI(c1 = cluster_df$clust_bact, c2 = cluster_df$clust_variety)
NMI(c1 = cluster_df$clust_bact, c2 = cluster_df$locality)
cluster_df%>%
  mutate(across(where(is.numeric), as.factor)) -> cluster_df
Y = tab.disjonctif(cluster_df$clust_bact_HCPC)
cca(Y~clust_variety+ clust_cult_history+clust_soil, data = cluster_df)
mod_varpart_bact = varpart(Y,
                                  ~clust_variety,
                                  ~clust_cult_history, 
                                  ~clust_soil, data= cluster_df, chisquare=T,
                                  permutations = how(nperm = 10000))
plot(mod_varpart_bact,Xnames = c(list("variety"), list("history"), list("soil")),
     bg = c( "forestgreen","#CB6CE6","#C0344A"), lwd = 2, cex = 1.5, alpha = 170)

### ---- RDA ----
otu_bact%>%
  column_to_rownames("soil_code")%>%
  select(where(~sum(.x)!=1 && sum(.x)!=2 ))-> otu_bact_clean

soil%>%
  select(-c(Site,N_site,Pto))%>%
  column_to_rownames("soil_code")%>%
  decostand(., "standardize") -> soil_stand
str(soil_stand)
otu_bact_clean%>%
  decostand(., method = "hellinger")%>%
  as.matrix -> otu_bact_stand
str(otu_bact_stand)
rownames(soil_stand)
mod0_rda <- rda(otu_bact_stand ~ 1, data= soil_stand)
mod1_rda = rda(otu_bact_stand~., data =soil_stand)

ordiplot(mod1_rda, scaling = 1)
ordiplot(mod1_rda, scaling = 2)

corrplot(cor(as.matrix(soil_stand)),  method = 'square',order = 'hclust' )
vif.cca(mod1_rda)
mod2_rda = rda(otu_bact_stand~CN+ TS +K + pHeau + CEC + Kto + MO, data =soil_stand)
vif.cca(mod2_rda)
ordiplot(mod2_rda, scaling = 1, type = "text")
ordiplot(mod2_rda, scaling = 2)

rda_step_select = ordistep(object = mod0_rda, scope = formula(mod2_rda),direction = "both")
ordiplot(rda_step_select, scaling = 2)

### ---- Richness ----
#### ---- With cluster ----
cluster_df%>%
  left_join(div_bact, by = "soil_code") -> model_df

model_df%>%
  pivot_longer(starts_with(c("clust_", "locality")), names_to = "clust_type", values_to = "clust")%>%
  mutate(clust = as.factor(clust))%>%
  ggplot()+
  geom_boxplot(aes(clust, Hill_Richness))+
  facet_wrap(~clust_type)

##### ---- Choose model residual type ----
hist(model_df$Hill_Richness,breaks = 20)
bact.mod1 = glm(Hill_Richness~clust_cult_history + clust_soil + clust_variety, data = model_df )
par(mfrow = c(2,2))
plot(bact.mod1)
par(mfrow = c(1,1))
hist(bact.mod1$residuals, breaks = 20)

##### ---- Model selection ----

bact.fullmod <- glm(Hill_Richness~clust_cult_history * clust_soil * clust_variety*locality, data = model_df,
                   na.action = "na.fail")
### model selection
dd_selection = dredge(bact.fullmod)

### n models tested
length(attr(dd_selection, "model.calls"))

### Best models
subset(dd_selection, delta < 2)%>%
  as.data.frame()

#The best one
summary(get.models(dd_selection, 1)[[1]])

lapply(get.models(dd_selection, c(1,2,3,4,5)), summary)

best_model = get.models(dd_selection, 1)[[1]]
best_model
summary(best_model)

best_model$deviance/best_model$df.residual


dev_model <- best_model$deviance
dev_nul <- best_model$null.deviance
(dev_nul-dev_model)/dev_nul #R² approx

par(mfrow = c(2,2))
plot(best_model)
par(mfrow = c(1,1))
hist(best_model$residuals, breaks = 30)

#### ---- With gradient ----

soil%>%
  left_join(div_bact, by = "soil_code") -> model_grad_df

model_grad_df%>%
  pivot_longer(-starts_with(c("Hill_", "N_site","Site", "soil_code")), names_to = "soil_type", values_to = "soil_values")%>%
  ggplot()+
  geom_point(aes(soil_values, Hill_Shannon))+
  geom_smooth(aes(soil_values, Hill_Shannon))+
  facet_wrap(~soil_type,scale = "free")+
  main_theme

##### ---- Choose model residual type ----
hist(model_grad_df$Hill_Richness,breaks = 20)
bact.grad.mod1 = glm(log(Hill_Richness)~ CN + CEC + K  + TS + pHeau + Pto + Kto + MO, data = model_grad_df)
summary(bact.grad.mod1)
par(mfrow = c(2,2))
plot(bact.grad.mod1)
par(mfrow = c(1,1))
hist(bact.grad.mod1$residuals, breaks = 20)

##### ---- Model selection ----

bact.fullmod <- glm(Hill_Richness~clust_cult_history * clust_soil * clust_variety*locality, data = model_df,
                    na.action = "na.fail")
### model selection
dd_selection = dredge(bact.fullmod)

### n models tested
length(attr(dd_selection, "model.calls"))

### Best models
subset(dd_selection, delta < 2)%>%
  as.data.frame()

#The best one
summary(get.models(dd_selection, 1)[[1]])

lapply(get.models(dd_selection, c(1,2,3,4,5)), summary)

best_model = get.models(dd_selection, 1)[[1]]
best_model
summary(best_model)

best_model$deviance/best_model$df.residual


dev_model <- best_model$deviance
dev_nul <- best_model$null.deviance
(dev_nul-dev_model)/dev_nul #R² approx

par(mfrow = c(2,2))
plot(best_model)
par(mfrow = c(1,1))
hist(best_model$residuals, breaks = 30)
