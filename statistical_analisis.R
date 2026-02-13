
#               STATISTICAL ANALISIS ASAD AND TABENI (2026) 

# Paper Title: Bridging human drivers on riparian plant metacommunity structure along a dryland natural-rural-urban gradient

#Journal: Drylands, special issue Drylands of South America: Ecology without borders that integrates environment and society

# Instituto Argentino de Investigación de las Zonas Áridas (IADIZA-CONICET), Mendoza, Argentina


#===============================================================================
rm(list = ls())
setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Project_Drylands/data_analysis")

# Packages
packages = c(
  "ggplot2", "ggmcmc", "agricolae", "fitdistrplus", "MASS", "cowplot", "qqplotr", 
  "ggeffects", "GGally", "broom", "doBy", "corrplot", "DHARMa", "pROC", "multcomp", 
  "multcompView", "car", "broom.mixed", "glmmTMB", "gamlss.dist", "bayesplot", 
  "reshape2", "gridExtra", "brms", "emmeans", "DirichletReg", "readxl", "tidyr", 
  "writexl", "stats", "ggrip_modif", "data.table", "jsonlite", "curl", 
  "tidyverse", "vegan", "FD", "AICcmodavg", "nlme", "GA", "gawdis", "lme4", 
  "pbkrtest", "lmerTest", "lavaan", "readr", "betapart", "dplyr", "tidyverse", "mgcv", 
  "geosphere", "performance", "ggpubr", "patchwork", "tidymv", "sjPlot")

invisible(lapply(packages, require, character.only = TRUE))


#================ TOTAL BETA DIVERSITY CALCULATION =============================
#SUPPLEMENTARY TABLE 1

sp_cover_plot = read.csv("STable1.csv", header = TRUE, sep = ",", stringsAsFactors = TRUE)
View(sp_cover_plot)

#===============================================================================
# BRAY-CURTIS Communities Dissimilarities indices accounting for species abundance
# Partition into Balanced Variation (BC_bal → Turnover) and Abundance gradient (BC_gra → Nestedness)
#===============================================================================
df = sp_cover_plot
colnames(df)

# Delete non Species Colummns
species_cols = setdiff(names(df), c("Mosaic", "Transect", "Plot", "SiteID", "lon", "lat", "Frag_LDI"))

# Only SiteID, Mosaic y especies
sp_abund_site = df %>% select(SiteID, Mosaic, all_of(species_cols))
sp_abund_site

#Bray-Curtis per siteID
beta_bray_results = list()

for (s in unique(sp_abund_site$SiteID)) {
  
  # Extraer data del sitio
  df_s = sp_abund_site %>% filter(SiteID == s)
  abund_s = df_s %>% select(-SiteID, -Mosaic)
  
  # empy plots
  empty_rows = rowSums(abund_s) == 0
  
  if (all(empty_rows)) {
    warning(paste("SiteID", s, "no tiene ninguna especie en ningún plot → β no calculable"))
    beta_bray_results[[s]] = tibble(SiteID = s, Mosaic = unique(df_s$Mosaic),
                                     BC_total = NA, BC_bal = NA, BC_gra = NA)
    next
  }
  
  # Exclude emty plots
  abund_s2 = abund_s[!empty_rows, ]
  
  if (nrow(abund_s2) < 2) {
    warning(paste("SiteID", s, "tiene solo 1 plot con especies → β no calculable"))
    beta_bray_results[[s]] = tibble(SiteID = s, Mosaic = unique(df_s$Mosaic),
                                     BC_total = NA, BC_bal = NA, BC_gra = NA)
    next
  }
  
  # Beta diversity Bray-Curtis and its partiton
  beta_core = betapart.core.abund(abund_s2)
  beta_pair = beta.pair.abund(beta_core, index.family = "bray")
  
  beta_bray_results[[s]] = tibble(
    SiteID  = s,
    Mosaic  = unique(df_s$Mosaic),
    BC_total = mean(as.dist(beta_pair$beta.bray)),
    BC_bal   = mean(as.dist(beta_pair$beta.bray.bal)),
    BC_gra   = mean(as.dist(beta_pair$beta.bray.gra))
  )
}

# Convertir lista en data frame final
beta_bray_results = bind_rows(beta_bray_results)

View(beta_bray_results)

#BRAY-CURTIS ONLY NATIVES AND EXOTIC
#SUPPLEMENTARY TABLE 2

sp_traits = read.csv("STable2.csv", row.names = 1)
View(sp_traits)
View(df)

native_species = rownames(sp_traits)[sp_traits$Origin == "native"]
exotic_species = rownames(sp_traits)[sp_traits$Origin == "exotic"]

species_cols = setdiff(names(df), c("Mosaic", "Transect", "Plot", "SiteID", "lat", "lon", "Frag_LDI"))

#Select Native and Exotic
native_species = intersect(native_species, species_cols)
native_species
exotic_species = intersect(exotic_species, species_cols)
exotic_species

# Dataframes per pool
df_native = df %>% select(SiteID, Mosaic, all_of(native_species))
View(df_native)
df_exotic = df %>% select(SiteID, Mosaic, all_of(exotic_species))
View(df_exotic)

# Bray–Curtis per pool
calcular_beta_bray = function(df_pool) {
  resultados = list()
  
  for (s in unique(df_pool$SiteID)) {
    df_s = df_pool %>% filter(SiteID == s)
    abund_s = df_s %>% select(-SiteID, -Mosaic)
    empty_rows = rowSums(abund_s) == 0
    
    if (all(empty_rows) || nrow(abund_s[!empty_rows, ]) < 2) {
      resultados[[s]] = tibble(SiteID = s,
                                Mosaic = unique(df_s$Mosaic),
                                BC_total = NA,
                                BC_bal = NA,
                                BC_gra = NA)
    } else {
      abund_s2 = abund_s[!empty_rows, ]
      beta_core = betapart.core.abund(abund_s2)
      beta_pair = beta.pair.abund(beta_core, index.family = "bray")
      
      resultados[[s]] = tibble(SiteID = s,
                                Mosaic = unique(df_s$Mosaic),
                                BC_total = mean(as.dist(beta_pair$beta.bray)),
                                BC_bal = mean(as.dist(beta_pair$beta.bray.bal)),
                                BC_gra = mean(as.dist(beta_pair$beta.bray.gra)))
    }
  }
  bind_rows(resultados)
}

# Aply each pool
beta_bray_native = calcular_beta_bray(df_native)
beta_bray_exotic = calcular_beta_bray(df_exotic)

beta_bray_native[is.na(beta_bray_native)] = 0
beta_bray_exotic[is.na(beta_bray_exotic)] = 0

#Re-name columns of the data frames natives and exotics
names(beta_bray_native)[names(beta_bray_native) %in% c("BC_total", "BC_bal", "BC_gra")] =
  c("BC_total_nat", "BC_bal_nat", "BC_gra_nat")

names(beta_bray_exotic)[names(beta_bray_exotic) %in% c("BC_total", "BC_bal", "BC_gra")] =
  c("BC_total_exo", "BC_bal_exo", "BC_gra_exo")

View(beta_bray_results)
View(beta_bray_native)
View(beta_bray_exotic)

#Save
write.csv(beta_bray_results, file = "Beta_Bray_Total.csv", row.names = FALSE)
write.csv(beta_bray_native, file = "Beta_Bray_Native.csv", row.names = FALSE)
write.csv(beta_bray_exotic, file = "Beta_Bray_Exotic.csv", row.names = FALSE)


#Table 1 - Estimate Proportions of Turnover and Nestedness for Total, Native and Exotic Pool

#Total Pool
beta_bray_results_prop = beta_bray_results %>%
  mutate(
    prop_turnover = BC_bal / BC_total,
    prop_nestedness = BC_gra / BC_total
  )
summary_prop_beta = beta_bray_results_prop %>%
  group_by(Mosaic) %>%
  summarise(
    Turnover = mean(prop_turnover, na.rm = TRUE),
    Nestedness = mean(prop_nestedness, na.rm = TRUE)
  )
print(summary_prop_beta)

#Natives
beta_bray_native_prop = beta_bray_native %>%
  mutate(
    prop_turnover = BC_bal_nat / BC_total_nat,
    prop_nestedness = BC_gra_nat / BC_total_nat
  )
summary_prop_beta = beta_bray_native_prop %>%
  group_by(Mosaic) %>%
  summarise(
    Turnover = mean(prop_turnover, na.rm = TRUE),
    Nestedness = mean(prop_nestedness, na.rm = TRUE)
  )
print(summary_prop_beta)

#Exotic
beta_bray_exotic_prop = beta_bray_exotic %>%
  mutate(
    prop_turnover = BC_bal_exo / BC_total_exo,
    prop_nestedness = BC_gra_exo / BC_total_exo
  )
summary_prop_beta = beta_bray_exotic_prop %>%
  group_by(Mosaic) %>%
  summarise(
    Turnover = mean(prop_turnover, na.rm = TRUE),
    Nestedness = mean(prop_nestedness, na.rm = TRUE)
  )
print(summary_prop_beta)

#===============================================================================
#       STATISTICAL ANALYSIS BETA DIVERSITY PER LANDSCAPE AND FRAGMENTATION
#===============================================================================
#SUPPLEMENTARY TABLE 3

data = read.csv("STable3.csv", header=T,  sep = ",", stringsAsFactors = T)
colnames(data)
View(data)

#Figure 2 → Correlation Between Landscape and Fragmentation (explanatory variables)

# Correct Landscape order (from west to east)
data$Landscape = factor(data$Landscape,
                       levels = c("Foothill", "Peri-Urban", "Urban", "Rural", "Wetland"))

mean_ldis = data %>%
  group_by(Landscape) %>%
  summarise(mean_LDI = mean(Frag_LDI, na.rm = TRUE))

landscape_colors = c(
  "Foothill"    = "#D9c100", 
  "Peri-Urban"  = "orange",  
  "Urban"       = "#BEBADA",  
  "Rural"       = "#A6D854",  
  "Wetland"     = "#80B1D3"   
)

Figure_2=ggplot(data, aes(x = Landscape, y = Frag_LDI, fill = Landscape)) +
  geom_boxplot(width = 0.5, color = "black", outlier.shape = NA, alpha = 0.6) +
  geom_point(alpha = 0.5, size = 2, position = position_nudge(x = 0)) +
  geom_line(data = mean_ldis, aes(x = Landscape, y = mean_LDI, group = 1),
            color = "black", linewidth = 1.2, inherit.aes = FALSE) +
  scale_fill_manual(values = landscape_colors) +
  labs(x = "Landscape", y = "Fragmentation (LDI)") +
  theme_classic(base_size = 17) +
  theme(
    axis.text = element_text(color = "black"),
    legend.position = "none",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )

print(Figure_2)

#Response variables Analysis (change response variable in script)
colnames(data)
hist(data$BC_total)

normal=fitdist(data$BC_total,"norm")
lognormal=fitdist(data$BC_total,"lnorm") 
gamma=fitdist(data$BC_total, "gamma")

cdf.sp=cdfcomp(list(normal, lognormal,gamma), main="", legendtext =c("Normal", "lognormal", "Gamma"),fitcol = c("orange", "blue", "red"), plotstyle ="ggplot")+
  geom_line(size=1.2)+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16), 
        legend.position = c(0.70,0.25),
        legend.text=element_text(size=16))
qq.sp=qqcomp(list(normal,lognormal, gamma), main="",fitcol = c("orange","blue", "red"), plotstyle 	="ggplot")+
  geom_line(size=1.2)+
  theme_bw()+
  theme(axis.title=element_text(size=18), 
        axis.text = element_text(size=16),
        legend.position ="none")
#Distribution plots
plot_grid(cdf.sp, qq.sp)
gofstat(list(normal,lognormal,gamma))$aic

#Best Distribution for all Beta Indices: GAMMA (normal also possible)

#===============================================================================
#                            GLM AND GAM MODELS 
#===============================================================================

# Explanatory figure
ggplot(data, aes(x=Frag_LDI, y=BC_bal, color=Landscape)) + 
  geom_point(size=2)+
  theme_bw()+ 
  geom_smooth(method="glm",  method.args = list(family = Gamma(link = "log")),se=F,size=1.5)+
  ylab("BC")+
  theme(legend.position="right",
        legend.text=element_text(size = 14),
        legend.title = element_text(size=16, face="bold"),
        axis.text = element_text(size=16),
        axis.title=element_text(size=18))

#Explicative variables relations
colnames(data)
ggpairs(data[,c(2,5)], mapping = aes(color = Landscape))+
theme_bw()


#GLM Total Beta Diversity ~ Landscape (Total Pool, Natives and Exotics)

###TOTAL

glmBC_total = glm(BC_total~ Landscape, data = data, family = Gamma(link = "log"))

summary(glmBC_total)
anova(glmBC_total, test="Chisq")
tidy(glmBC_total, conf.int=T)   
glance(glmBC_total)
pred.BC_total_Land = ggpredict(glmBC_total, terms = c("Landscape"))
pred = ggpredict(glmBC_total, terms = "Landscape")
tab_model(glmBC_total,show.r2=T, show.ci = 0.95)   

pred_glmBC_total=ggpredict(glmBC_total, terms = c("Landscape"))
plot(pred_glmBC_total)+ 
  theme_bw()+
  labs(x="Landscape", y="BC_total")+
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
simres = simulateResiduals(glmBC_total)
plot(simres)
plotResiduals(simres, data$Landscape)
#Dispersion
testDispersion(simres)

#There is over-disperson and variance heterogenety, possible due to high variance in Beta Diversity of Peri-Urban sites
#Overall is acceptable

# Post-hoc Analysis
emm = emmeans(glmBC_total, pairwise ~ Landscape, adjust = "tukey")
cld = cld(emm$emmeans, alpha = 0.05, Letters = letters) #p-value < 0.05
data$Landscape = factor(data$Landscape, levels = c("Foothill", "Peri-Urban", "Urban", "Rural", "Wetland"))
cld$Landscape  = factor(cld$Landscape,  levels = c("Foothill", "Peri-Urban", "Urban", "Rural", "Wetland"))

palette_landscape = c(
  "Foothill"    = "#D9c100", 
  "Peri-Urban"  = "orange",  
  "Urban"       = "#BEBADA",  
  "Rural"       = "#A6D854",  
  "Wetland"     = "#80B1D3"
)

#Significance Letters
cld$Landscape = factor(cld$Landscape, levels = levels(data$Landscape))
med_vals = aggregate(BC_total ~ Landscape, data = data, median)
cld$y_pos = 0.7 #Letter Possition

#Plot
BC_total = ggplot(data, aes(x = Landscape, y = BC_total, fill = Landscape)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, color = "black", aes(fill = Landscape)) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 2) +
  geom_text(data = cld, aes(x = Landscape, y = y_pos, label = .group),
            inherit.aes = FALSE, size = 6, color = "black", nudge_x = 0.1) +
  scale_fill_manual(values = palette_landscape) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Landscape", 
       y = "dBC",
       title = "Total Species Pool") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none")

BC_total

### NATIVES

glmBC_total_nat = glm(BC_total_nat ~ Landscape, data = data)

summary(glmBC_total_nat)
anova(glmBC_total_nat, test="Chisq")
tidy(glmBC_total_nat, conf.int=T)   
glance(glmBC_total_nat)
pred.BC_total_nat_Land = ggpredict(glmBC_total_nat, terms = c("Landscape"))
pred = ggpredict(glmBC_total_nat, terms = "Landscape")
tab_model(glmBC_total_nat,show.r2=T, show.ci = 0.95)   

pred_glmBC_total_nat=ggpredict(glmBC_total_nat, terms = c("Landscape"))
plot(pred_glmBC_total_nat)+ 
  theme_bw()+
  labs(x="Landscape", y="BC_total_nat")+
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
simres = simulateResiduals(glmBC_total_nat)
plot(simres)
plotResiduals(simres, data$Landscape)
#Dispersion
testDispersion(simres)

#Overall is acceptable

# Post-hoc Analysis
emm = emmeans(glmBC_total_nat, pairwise ~ Landscape, adjust = "tukey")
cld = cld(emm$emmeans, alpha = 0.05, Letters = letters) #p-value < 0.05
cld$Landscape = factor(cld$Landscape, levels = c("Foothill", "Peri-Urban", "Urban", "Rural", "Wetland"))

palette_landscape = c(
  "Foothill"    = "#D9c100", 
  "Peri-Urban"  = "orange",  
  "Urban"       = "#BEBADA",  
  "Rural"       = "#A6D854",  
  "Wetland"     = "#80B1D3"
)

#Significance Letters
cld$Landscape = factor(cld$Landscape, levels = levels(data$Landscape))
med_vals = aggregate(BC_total_nat ~ Landscape, data = data, median)
cld$y_pos = 0.7 #Letter Possition

#Plot
BC_total_nat = ggplot(data, aes(x = Landscape, y = BC_total_nat, fill = Landscape)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, color = "black", aes(fill = Landscape)) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 2) +
  geom_text(data = cld, aes(x = Landscape, y = y_pos, label = .group),
            inherit.aes = FALSE, size = 6, color = "black", nudge_x = 0.1) +
  scale_fill_manual(values = palette_landscape) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Landscape", 
       y = "dBC",
       title = "Native Species") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none")


BC_total_nat


### EXOTICS

#Exclude Foothill from the analysis
data_noFoothill = subset(data, Landscape != "Foothill")
data_noFoothill$Landscape = droplevels(data_noFoothill$Landscape)
table(data_noFoothill$Landscape)

glmBC_total_exo = glm(BC_total_exo ~ Landscape, data = data_noFoothill)

summary(glmBC_total_exo)
anova(glmBC_total_exo, test="Chisq")
tidy(glmBC_total_exo, conf.int=T)   
glance(glmBC_total_exo)
pred.BC_total_exo_Land = ggpredict(glmBC_total_exo, terms = c("Landscape"))
pred = ggpredict(glmBC_total_exo, terms = "Landscape")
tab_model(glmBC_total_exo,show.r2=T, show.ci = 0.95)   

pred_glmBC_total_exo=ggpredict(glmBC_total_exo, terms = c("Landscape"))
plot(pred_glmBC_total_exo)+ 
  theme_bw()+
  labs(x="Landscape", y="BC_total_exo")+
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
simres = simulateResiduals(glmBC_total_exo)
plot(simres)
plotResiduals(simres, data$Landscape)
#Dispersion
testDispersion(simres)

#Overall is acceptable

# Post-hoc Analysis
emm = emmeans(glmBC_total_exo, pairwise ~ Landscape, adjust = "tukey")
cld = cld(emm$emmeans, alpha = 0.05, Letters = letters) #p-value < 0.05
cld$Landscape = factor(cld$Landscape, levels = c("Peri-Urban", "Urban", "Rural", "Wetland"))

palette_landscape_noFoothill = c(
  "Peri-Urban"  = "orange",  
  "Urban"       = "#BEBADA",  
  "Rural"       = "#A6D854",  
  "Wetland"     = "#80B1D3"
)

#Significance Letters
cld$Landscape = factor(cld$Landscape, levels = levels(data$Landscape))
med_vals = aggregate(BC_total_exo ~ Landscape, data = data, median)
cld$y_pos = 0.7 #Letter Possition

#Plot
BC_total_exo = ggplot(data_noFoothill, aes(x = Landscape, y = BC_total_exo, fill = Landscape)) +
  geom_violin(trim = FALSE, alpha = 0.6, color = NA) +
  geom_boxplot(width = 0.15, color = "black", aes(fill = Landscape)) +
  geom_jitter(width = 0.1, alpha = 0.4, size = 2) +
  geom_text(data = cld, aes(x = Landscape, y = y_pos, label = .group),
            inherit.aes = FALSE, size = 6, color = "black", nudge_x = 0.1) +
  scale_fill_manual(values = palette_landscape_noFoothill) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.05))) +
  theme_bw() +
  labs(x = "Landscape", 
       y = "dBC",
       title = "Exotic Species") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    legend.position = "none")

BC_total_exo

#Figure 3 A-B-C
print(BC_total)
print(BC_total_nat)
print(BC_total_exo)


#GAM Total Beta Diversity ~ Fragmentation (LDI) (Total Pool, Natives and Exotics)

#Total Pool

colnames(data)

gamBC_total = gam(BC_total ~ s(Frag_LDI), data = data, family = Gamma(link = "log"))
summary(gamBC_total)
anova(gamBC_total, test="Chisq")
tidy(gamBC_total, conf.int=T)   
glance(gamBC_total)

#Predictions
pred_gamBC_total= ggpredict(gamBC_total, terms = c("Frag_LDI"))
pred = ggpredict(gamBC_total, terms = "Frag_LDI")
tab_model(gamBC_total,show.r2=T, show.ci = 0.95)   

plot(pred_gamBC_total)+ 
  theme_bw()+
  labs(x="Frag_LDI", y="BC_total")+
  scale_y_continuous(limits = c(0, 1))+ 
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
#mgcv general check
gam.check(gamBC_total)
par(mfrow=c(2,2))
plot(gamBC_total, residuals = TRUE, pch = 19, cex = 0.6)

# Residuals vs Adjusted
plot(fitted(gamBC_total), residuals(gamBC_total),
     pch = 19,
     xlab = "Valores ajustados",
     ylab = "Residuales",
     main = "Residuales vs Ajustados")
abline(h = 0, col = "red")

#ACF
acf(residuals(gamBC_total), main = "ACF de residuales")

#Simres
plot(simres)

# Suavizant
k.check(gamBC_total)

#No Residuals Problems


#PLOT
plot_gam_total = ggplot(pred_gamBC_total, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#B3CDE3", alpha = 0.8) +
  geom_line(color = "#045a8d", size = 2) +
  geom_point(data = data, aes(x = Frag_LDI, y = BC_total), inherit.aes = FALSE, 
             alpha = 0.2, size = 1, color = "black") +
  annotate("text", x = 0.02, y = 0.95, hjust = 0,
           label = paste0("p < 0.001\nDev. explained = 30.1%"),
           size = 4.5, color = "black") +
  theme_bw() +
  labs(x = "Fragmentation (LDI)", 
       y = "dBC",
       title = "Total Species Pool") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 1))  

plot_gam_total


#NATIVES
colnames(data)


gamBC_total_nat = gam(BC_total_nat ~ s(Frag_LDI), data = data)
summary(gamBC_total_nat)
anova(gamBC_total_nat, test="Chisq")
tidy(gamBC_total_nat, conf.int=T)   
glance(gamBC_total_nat)

#Predictions
pred_gamBC_total_nat= ggpredict(gamBC_total_nat, terms = c("Frag_LDI"))
pred = ggpredict(gamBC_total_nat, terms = "Frag_LDI")
tab_model(gamBC_total_nat,show.r2=T, show.ci = 0.95)   

plot(pred_gamBC_total_nat)+ 
  theme_bw()+
  labs(x="Frag_LDI", y="BC_total_nat")+
  scale_y_continuous(limits = c(0, 1))+ 
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
#mgcv general check
gam.check(gamBC_total_nat)

par(mfrow=c(2,2))
plot(gamBC_total_nat, residuals = TRUE, pch = 19, cex = 0.6)

# Residuals vs Adjusted
plot(fitted(gamBC_total_nat), residuals(gamBC_total_nat),
     pch = 19,
     xlab = "Valores ajustados",
     ylab = "Residuales",
     main = "Residuales vs Ajustados")
abline(h = 0, col = "red")

#ACF
acf(residuals(gamBC_total_nat), main = "ACF de residuales")

#Simres
plot(simres)

# Suavizant
k.check(gamBC_total_nat)

#Not major Residuals Problems


#PLOT
plot_gam_total_nat = ggplot(pred_gamBC_total_nat, aes(x = x, y = predicted)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#B3CDE3", alpha = 0.8) +
  geom_line(color = "#045a8d", size = 2) +
  geom_point(data = data, aes(x = Frag_LDI, y = BC_total_nat), inherit.aes = FALSE, 
             alpha = 0.2, size = 1, color = "black") +
  annotate("text", x = 0.02, y = 0.95, hjust = 0,
           label = paste0("p < 0.001\nDev. explained = 21.4%"),
           size = 4.5, color = "black") +
  theme_bw() +
  labs(x = "Fragmentation (LDI)", 
       y = "dBC",
       title = "Native Species") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 1))  

plot_gam_total_nat

#EXOTIC

colnames(data_noFoothill)

gamBC_total_exo = gam(BC_total_exo ~ s(Frag_LDI), data = data_noFoothill)
summary(gamBC_total_exo)
anova(gamBC_total_exo, test="Chisq")
tidy(gamBC_total_exo, conf.int=T)   
glance(gamBC_total_exo)

#Predictions
pred_gamBC_total_exo= ggpredict(gamBC_total_exo, terms = c("Frag_LDI"))
pred = ggpredict(gamBC_total_exo, terms = "Frag_LDI")
tab_model(gamBC_total_exo,show.r2=T, show.ci = 0.95)   

plot(pred_gamBC_total_exo)+ 
  theme_bw()+
  labs(x="Frag_LDI", y="BC_total_exo")+
  scale_y_continuous(limits = c(0, 1))+ 
  theme(legend.position = c(0.25,0.20), 
        plot.title=element_blank(),
        axis.title=element_text(size=18), 
        axis.text = element_text(size=16))

#Residuals
#mgcv general check
gam.check(gamBC_total_exo)

par(mfrow=c(2,2))
plot(gamBC_total_exo, residuals = TRUE, pch = 19, cex = 0.6)

# Residuals vs Adjusted
plot(fitted(gamBC_total_exo), residuals(gamBC_total_exo),
     pch = 19,
     xlab = "Valores ajustados",
     ylab = "Residuales",
     main = "Residuales vs Ajustados")
abline(h = 0, col = "red")

#ACF
acf(residuals(gamBC_total_exo), main = "ACF de residuales")

#Simres
plot(simres)

# Suavizant
k.check(gamBC_total_exo)

#Not major Residuals Problems


#PLOT
plot_gam_total_exo = ggplot(pred_gamBC_total_exo, aes(x = x, y = predicted)) +
  #geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#B3CDE3", alpha = 0.8) +
  #geom_line(color = "#045a8d", size = 2) +
  geom_point(data = data, aes(x = Frag_LDI, y = BC_total_exo), inherit.aes = FALSE, 
             alpha = 0.2, size = 1, color = "black") +
  annotate("text", x = 0.02, y = 0.95, hjust = 0,
           label = paste0("p < 0.1\nDev. explained = 4.8%"),
           size = 4.5, color = "black") +
  theme_bw() +
  labs(x = "Fragmentation (LDI)", 
       y = "dBC",
       title = "Exotic Species") +
  theme(
    panel.grid = element_blank(),             
    panel.border = element_rect(color = "black", fill = NA, size = 1),  
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12, color = "black"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  ) +
  scale_y_continuous(limits = c(0, 1))  

plot_gam_total_exo

#Figure 3 D-E-F
print(plot_gam_total)
print(plot_gam_total_nat)
print(plot_gam_total_exo)

#FINAL PLOT

#Delete axis tittles and combine
BC_total_no_x = BC_total + theme(axis.title.x = element_blank())
BC_total_nat_no_x = BC_total_nat +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())
BC_total_exo_no_x = BC_total_exo +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())

plot_gam_total_no_x = plot_gam_total +
  theme(axis.title.x = element_blank(), plot.title = element_blank())
plot_gam_total_nat_no_x = plot_gam_total_nat +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title = element_blank())
plot_gam_total_exo_no_x = plot_gam_total_exo +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), plot.title = element_blank())

#Letters above
row1 = cowplot::plot_grid(BC_total_no_x, BC_total_nat_no_x, BC_total_exo_no_x,
  labels = c("A", "B", "C"),
  label_size = 12, label_fontface = "bold",
  label_x = 0.98, label_y = 0.95,
  hjust = 1, vjust = 1,
  ncol = 3, align = "v", axis = "l")
row2 = cowplot::plot_grid(plot_gam_total_no_x, plot_gam_total_nat_no_x, plot_gam_total_exo_no_x,
  labels = c("D", "E", "F"),
  label_size = 12, label_fontface = "bold",
  label_x = 0.98, label_y = 1.05,
  hjust = 1, vjust = 1,
  ncol = 3, align = "v", axis = "l")

# New Axis titles
xlabel1 = ggdraw() + draw_label("", size = 16)
xlabel2 = ggdraw() + draw_label("Fragmentation (LDI)", size = 16)

# Combine
top = cowplot::plot_grid(row1, xlabel1, ncol = 1, rel_heights = c(1, 0.08))
bottom = cowplot::plot_grid(row2, xlabel2, ncol = 1, rel_heights = c(1, 0.12)) 
combined = cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1.05))


# Final Figure 3
Figure_3 = cowplot::plot_grid(top, bottom, ncol = 1, rel_heights = c(1, 1.05))
print(Figure_3)

#===============================================================================
# BETA PAIRWISE MANTEL TESTS RELATIION WITH SPATIAL AND FRAGMENTATION DISTANCE 
#===============================================================================
#Supplementary Table 1

#GENERAL SPECIES POOL

df = read.csv("STable1.csv", header = TRUE, sep = ",")
View(df)
colnames(df)

#Build a community x species matrix according to SiteID
#Automatically detect al the species columns

# No-Species columns
non_species_cols = c("Mosaic", "Transect", "Plot", "SiteID", "lon", "lat", "Frag_LDI")

#Species columns
species_cols = setdiff(names(df), non_species_cols)

# Matrix per plot SiteID + especies
sp_plot = df[, c("SiteID", species_cols)]

# Add 10 plots per SiteID (metacommunity)
sp_site = sp_plot %>%
  group_by(SiteID) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Species Matrix
sp_mat = as.data.frame(sp_site)
rownames(sp_mat) = sp_mat$SiteID
sp_mat = sp_mat[, -1]

sp_mat = sp_mat[, -1]

View(sp_mat)

# Extract Fragmentation LDI values and Coordinates 
env_site = df %>%
  select(SiteID, lon, lat, Frag_LDI) %>%
  distinct()   
env_site

#Check coincidation
env_site = env_site[match(rownames(sp_mat), env_site$SiteID), ]
all(env_site$SiteID == rownames(sp_mat))

# Spatial Distance (Haversine in km)
coords = env_site[, c("lon", "lat")]
dist_geo_m = distm(coords, coords, fun = distHaversine)
dist_geo = as.dist(dist_geo_m / 1000)   #km

#Fragmentation distance (Frag_LDI)
frag = env_site$Frag_LDI
dist_frag = dist(frag, method = "euclidean")   # |Δ Frag_LDI|

# Pair-Wise Beta diversity calculation

# Total Bray–Curtis
bc_total = vegdist(sp_mat, method = "bray")

# Partition Baselga (2013)
beta_bc = beta.pair.abund(sp_mat)

bc_bal = beta_bc$beta.bray.bal   # turnover
bc_gra = beta_bc$beta.bray.gra   # nestedness

# Mantel Tests (Beta Pairwise and Environmental and Spatial (geographic) Distance)

mantel_total_geo = mantel(bc_total, dist_geo)
mantel_bal_geo   = mantel(bc_bal, dist_geo)
mantel_gra_geo   = mantel(bc_gra, dist_geo)

mantel_total_frag = mantel(bc_total, dist_frag)
mantel_bal_frag   = mantel(bc_bal, dist_frag)
mantel_gra_frag   = mantel(bc_gra, dist_frag)

mantel_total_geo
mantel_bal_geo
mantel_gra_geo

mantel_total_frag
mantel_bal_frag
mantel_gra_frag

# Order Data for Plots
df_pairs = data.frame(
  BC_total = as.vector(bc_total),
  BC_bal   = as.vector(bc_bal),
  BC_gra   = as.vector(bc_gra),
  Spatial  = as.vector(dist_geo),
  Frag     = as.vector(dist_frag)
)

# Plots

plot_dd = function(x, y, xlab, ylab, mantel_res, index_type){
  
  # Mantel Regresión line p<0.05
  pval = mantel_res$signif
  add_lm = pval < 0.05
  
  # Text
  if(index_type == "gra"){
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 1
    hjust_value = 1
    vjust_value = 1
  } else {
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 0.02
    hjust_value = 1
    vjust_value = 0
  }
  
  g = ggplot(df_pairs, aes(x = {{x}}, y = {{y}})) +
    geom_point(alpha = 0.7, size = 1) +
    
    # Line only if p < 0.05
    { if(add_lm) geom_smooth(method = "lm", formula = y ~ log(x + 1),
                             se = FALSE, color="#045a8d", linewidth=2) } +
    
    ylim(0, 1) +
    labs(x = NULL, y = ylab) +
    
    annotate("text",
             x = xpos,
             y = ypos,
             hjust = hjust_value,
             vjust = vjust_value,
             label = paste0("r² = ", round(mantel_res$statistic, 2),
                            "\np < ", signif(pval, 2))) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
      axis.title  = element_text(size=12, face="plain"),
      axis.text   = element_text(size=11, color="black"),
      plot.margin = margin(15,5,5,5)
    )
  
  return(g)
}

pA = plot_dd(Spatial, BC_total, NULL, "Pairwise dBC", mantel_total_geo, index_type="total")
pB = plot_dd(Spatial, BC_bal,   NULL, "Pairwise dBC-bal (turnover)", mantel_bal_geo, index_type="bal")
pC = plot_dd(Spatial, BC_gra,   NULL, "PairWise dBC-gra (nestedness)", mantel_gra_geo, index_type="gra")
pD = plot_dd(Frag, BC_total, NULL, "Pairwise dBC", mantel_total_frag, index_type="total")
pE = plot_dd(Frag, BC_bal,   NULL, "Pairwise dBC-bal (turnover)", mantel_bal_frag, index_type="bal")
pF = plot_dd(Frag, BC_gra,   NULL, "Pairwise dBC-gra (nestedness)", mantel_gra_frag, index_type="gra")


# Supperior Panel Letters
panel_top = ggarrange(pA, pB, pC,ncol = 3,
                       labels = c("A", "B", "C"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.99,
                       hjust = 1, vjust = 1)

# Inferior Panel Letters
panel_bot = ggarrange(pD, pE, pF,ncol = 3,
                       labels = c("D", "E", "F"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.99,
                       hjust = 1, vjust = 1)


Figure_4 = ggarrange(
  panel_top,
  text_grob("Spatial distance (Δ km)", size=12),
  panel_bot,
  text_grob("Fragmentation distance (Δ LDI)", size=12),
  ncol=1,
  heights=c(1, 0.08, 1, 0.08)
)

Figure_4

#===============================================================================
#               SUPPLEMENTARY MATERIAL ANALYSIS AND FIGURES
#===============================================================================

#NATIVE SPECIES POOL

View(df_native)
colnames(df_native)

#Build a community x species matrix according to SiteID
#Automatically detect al the species columns

# No-Species columns
non_species_cols = c("Mosaic", "SiteID")

#Species columns
species_cols = setdiff(names(df_native), non_species_cols)

species_cols

# Matrix per plot SiteID + especies
sp_plot = df_native[, c("SiteID", species_cols)]

# Add 10 plots per SiteID (metacommunity)
sp_site = sp_plot %>%
  group_by(SiteID) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Species Matrix
sp_mat_native = as.data.frame(sp_site)
rownames(sp_mat_native) = sp_mat_native$SiteID
sp_mat_native = sp_mat_native[, -1]

sp_mat_native = sp_mat_native[, -1]

View(sp_mat_native)

# Extract Fragmentation LDI values and Coordinates 
env_site = df %>%
  select(SiteID, lon, lat, Frag_LDI) %>%
  distinct()   
env_site

#Check coincidation
env_site = env_site[match(rownames(sp_mat_native), env_site$SiteID), ]
all(env_site$SiteID == rownames(sp_mat_native))


# Spatial Distance (Haversine en km)
coords = env_site[, c("lon", "lat")]
dist_geo_m = distm(coords, coords, fun = distHaversine)
dist_geo = as.dist(dist_geo_m / 1000)   #km

#Fragmentation distance (Frag_LDI)
frag = env_site$Frag_LDI
dist_frag = dist(frag, method = "euclidean")   # |Δ Frag_LDI|

# Pair-Wise Beta diversity calculation

# Total Bray–Curtis for natives
bc_total_native = vegdist(sp_mat_native, method = "bray")

# Partition Baselga (2013)
beta_bc = beta.pair.abund(sp_mat_native)
bc_bal_native = beta_bc$beta.bray.bal   # turnover
bc_gra_native = beta_bc$beta.bray.gra   # nestedness

# Mantel Tests (Beta Pairwise and Environmental and Spatial (geographic) Distance)

mantel_total_geo = mantel(bc_total_native, dist_geo)
mantel_bal_geo   = mantel(bc_bal_native, dist_geo)
mantel_gra_geo   = mantel(bc_gra_native, dist_geo)

mantel_total_frag = mantel(bc_total_native, dist_frag)
mantel_bal_frag   = mantel(bc_bal_native, dist_frag)
mantel_gra_frag   = mantel(bc_gra_native, dist_frag)

mantel_total_geo
mantel_bal_geo
mantel_gra_geo

mantel_total_frag
mantel_bal_frag
mantel_gra_frag

# Order Data for Plots
df_pairs = data.frame(
  BC_total_native = as.vector(bc_total_native),
  BC_bal_native   = as.vector(bc_bal_native),
  BC_gra_native  = as.vector(bc_gra_native),
  Spatial  = as.vector(dist_geo),
  Frag     = as.vector(dist_frag)
)

# Plots

plot_dd = function(x, y, xlab, ylab, mantel_res, index_type){
  
  # Mantel Regresión line p<0.05
  pval = mantel_res$signif
  add_lm = pval < 0.05
  
  # Text
  if(index_type == "gra"){
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 1
    hjust_value = 1
    vjust_value = 1
  } else {
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 0.02
    hjust_value = 1
    vjust_value = 0
  }
  
  g = ggplot(df_pairs, aes(x = {{x}}, y = {{y}})) +
    geom_point(alpha = 0.3, size = 1.5) +
    
    # Line only if p < 0.05
    { if(add_lm) geom_smooth(method = "lm", formula = y ~ log(x + 1),
                             se = FALSE, color="#045a8d") } +
    
    ylim(0, 1) +
    labs(x = NULL, y = ylab) +
    
    annotate("text",
             x = xpos,
             y = ypos,
             hjust = hjust_value,
             vjust = vjust_value,
             label = paste0("r² = ", round(mantel_res$statistic, 2),
                            "\np < ", signif(pval, 2))) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
      axis.title  = element_text(size=12, face="plain"),
      axis.text   = element_text(size=11, color="black"),
      plot.margin = margin(15,5,5,5)
    )
  
  return(g)
}

pG = plot_dd(Spatial, BC_total_native, NULL, "dBC total", mantel_total_geo, index_type="total")
pH = plot_dd(Spatial, BC_bal_native,   NULL, "dBC bal (turnover)", mantel_bal_geo, index_type="bal")
pI = plot_dd(Spatial, BC_gra_native,   NULL, "dBC gra (nestedness)", mantel_gra_geo, index_type="gra")
pJ = plot_dd(Frag, BC_total_native, NULL, "dBC total", mantel_total_frag, index_type="total")
pK = plot_dd(Frag, BC_bal_native,   NULL, "dBC bal (turnover)", mantel_bal_frag, index_type="bal")
pL = plot_dd(Frag, BC_gra_native,   NULL, "dBC gra (nestedness)", mantel_gra_frag, index_type="gra")


# Supperior Panel Letters
panel_top = ggarrange(pG, pH, pI,ncol = 3,
                       labels = c("G", "H", "I"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.98,
                       hjust = 1, vjust = 1)

# Inferior Panel Letters
panel_bot = ggarrange(pJ, pK, pL,ncol = 3,
                       labels = c("J", "K", "L"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.98,
                       hjust = 1, vjust = 1)


final_plot_native = ggarrange(
  panel_top,
  text_grob("Spatial distance (Δ km)", size=12),
  panel_bot,
  text_grob("Fragmentation distance (Δ LDI)", size=12),
  ncol=1,
  heights=c(1, 0.08, 1, 0.08)
)
final_plot_native

#EXOTIC SPECIES POOL

View(df_exotic)
colnames(df_exotic)

#Build a community x species matrix according to SiteID
#Automatically detect al the species columns

# No-Species columns
non_species_cols = c("Mosaic", "SiteID")

#Species columns
species_cols = setdiff(names(df_exotic), non_species_cols)
species_cols

# Matrix per plot SiteID + especies
sp_plot = df_exotic[, c("SiteID", species_cols)]

# Add 10 plots per SiteID (metacommunity)
sp_site = sp_plot %>%
  group_by(SiteID) %>%
  summarise(across(everything(), sum), .groups = "drop")

# Species Matrix
sp_mat_exotic = as.data.frame(sp_site)
rownames(sp_mat_exotic) = sp_mat_exotic$SiteID
sp_mat_exotic = sp_mat_exotic[, -1]

sp_mat_exotic = sp_mat_exotic[, -1]
sp_mat_exotic = sp_mat_exotic[rowSums(sp_mat_exotic) > 0, ]

View(sp_mat_exotic)

# Extract Fragmentation LDI values and Coordinates 
env_site = df %>%
  select(SiteID, lon, lat, Frag_LDI) %>%
  distinct()   
env_site

#Check coincidation
env_site = env_site[match(rownames(sp_mat_exotic), env_site$SiteID), ]
all(env_site$SiteID == rownames(sp_mat_exotic))


# Spatial Distance (Haversine en km)
coords = env_site[, c("lon", "lat")]
dist_geo_m = distm(coords, coords, fun = distHaversine)
dist_geo = as.dist(dist_geo_m / 1000)   #km

#Fragmentation distance (Frag_LDI)
frag = env_site$Frag_LDI
dist_frag = dist(frag, method = "euclidean")   # |Δ Frag_LDI|

# Pair-Wise Beta diversity calculation


# Total Bray–Curtis for exotics
bc_total_exotic = vegdist(sp_mat_exotic, method = "bray")

# Partition Baselga (2013)
beta_bc = beta.pair.abund(sp_mat_exotic)
bc_bal_exotic = beta_bc$beta.bray.bal   # turnover
bc_gra_exotic = beta_bc$beta.bray.gra   # nestedness

# Mantel Tests (Beta Pairwise and Environmental and Spatial (geographic) Distance)

mantel_total_geo = mantel(bc_total_exotic, dist_geo)
mantel_bal_geo   = mantel(bc_bal_exotic, dist_geo)
mantel_gra_geo   = mantel(bc_gra_exotic, dist_geo)

mantel_total_frag = mantel(bc_total_exotic, dist_frag)
mantel_bal_frag   = mantel(bc_bal_exotic, dist_frag)
mantel_gra_frag   = mantel(bc_gra_exotic, dist_frag)

mantel_total_geo
mantel_bal_geo
mantel_gra_geo

mantel_total_frag
mantel_bal_frag
mantel_gra_frag

# Order Data for Plots
df_pairs = data.frame(
  BC_total_exotic = as.vector(bc_total_exotic),
  BC_bal_exotic   = as.vector(bc_bal_exotic),
  BC_gra_exotic  = as.vector(bc_gra_exotic),
  Spatial  = as.vector(dist_geo),
  Frag     = as.vector(dist_frag)
)

# Plots

plot_dd = function(x, y, xlab, ylab, mantel_res, index_type){
  
  # Mantel Regresión line p<0.05
  pval = mantel_res$signif
  add_lm = pval < 0.05
  
  # Text
  if(index_type == "gra"){
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 1
    hjust_value = 1
    vjust_value = 1
  } else {
    xpos = max(df_pairs[[deparse(substitute(x))]], na.rm = TRUE)
    ypos = 0.02
    hjust_value = 1
    vjust_value = 0
  }
  
  g = ggplot(df_pairs, aes(x = {{x}}, y = {{y}})) +
    geom_point(alpha = 0.3, size = 1.5) +
    
    # Line only if p < 0.05
    { if(add_lm) geom_smooth(method = "lm", formula = y ~ log(x + 1),
                             se = FALSE, color="#045a8d") } +
    
    ylim(0, 1) +
    labs(x = NULL, y = ylab) +
    
    annotate("text",
             x = xpos,
             y = ypos,
             hjust = hjust_value,
             vjust = vjust_value,
             label = paste0("r² = ", round(mantel_res$statistic, 2),
                            "\np < ", signif(pval, 2))) +
    
    theme_classic(base_size = 14) +
    theme(
      panel.border = element_rect(color="black", fill=NA, linewidth=0.5),
      axis.title  = element_text(size=12, face="plain"),
      axis.text   = element_text(size=11, color="black"),
      plot.margin = margin(15,5,5,5)
    )
  
  return(g)
}

pM = plot_dd(Spatial, BC_total_exotic, NULL, "dBC total", mantel_total_geo, index_type="total")
pN = plot_dd(Spatial, BC_bal_exotic,   NULL, "dBC bal (turnover)", mantel_bal_geo, index_type="bal")
pO = plot_dd(Spatial, BC_gra_exotic,   NULL, "dBC gra (nestedness)", mantel_gra_geo, index_type="gra")
pP = plot_dd(Frag, BC_total_exotic, NULL, "dBC total", mantel_total_frag, index_type="total")
pQ = plot_dd(Frag, BC_bal_exotic,   NULL, "dBC bal (turnover)", mantel_bal_frag, index_type="bal")
pR = plot_dd(Frag, BC_gra_exotic,   NULL, "dBC gra (nestedness)", mantel_gra_frag, index_type="gra")


# Supperior Panel Letters
panel_top = ggarrange(pM, pN, pO,ncol = 3,
                       labels = c("M", "N", "O"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.98,
                       hjust = 1, vjust = 1)

# Inferior Panel Letters
panel_bot = ggarrange(pP, pQ, pR,ncol = 3,
                       labels = c("P", "Q", "R"),
                       font.label = list(size = 12, face = "bold"),
                       label.x = 0.98, label.y = 0.98,
                       hjust = 1, vjust = 1)


final_plot_exotic = ggarrange(
  panel_top,
  text_grob("Spatial distance (Δ km)", size=12),
  panel_bot,
  text_grob("Fragmentation distance (Δ LDI)", size=12),
  ncol=1,
  heights=c(1, 0.08, 1, 0.08)
)
final_plot_exotic


######          Análisis of unique species per landscape

# unify peri-urban and urban
df$Mosaico2 <- ifelse(df$Mosaic %in% c("Peri-Urban", "Urban"), "Urbanized", df$Mosaic)

#extract species
species_cols <- names(df)[which(names(df) == "Frag_LDI") + 1 : (ncol(df) - which(names(df) == "Frag_LDI"))]

# Function
get_species_list <- function(data, group_name) {
  data %>%
    filter(Mosaico2 == group_name) %>%
    select(all_of(species_cols)) %>%
    summarise(across(everything(), ~ sum(. > 0, na.rm = TRUE))) %>%
    select(where(~ . > 0)) %>%
    names()
}

sp_foothill   <- get_species_list(df, "Foothill")
sp_urbanized  <- get_species_list(df, "Urbanized")
sp_rural      <- get_species_list(df, "Rural")
sp_wetland    <- get_species_list(df, "Wetland")

# Unique species per Landscapes
unique_foothill   <- setdiff(sp_foothill, c(sp_urbanized, sp_rural, sp_wetland))
unique_foothill
unique_urbanized  <- setdiff(sp_urbanized, c(sp_foothill, sp_rural, sp_wetland))
unique_urbanized
unique_rural      <- setdiff(sp_rural, c(sp_foothill, sp_urbanized, sp_wetland))
unique_rural
unique_wetland    <- setdiff(sp_wetland, c(sp_foothill, sp_urbanized, sp_rural))
unique_wetland

# Number of unique species per landscape
length(unique_foothill)
length(unique_urbanized)
length(unique_rural)
length(unique_wetland)


# Correlation between spatial and fragmentation distance
mantel_geo_frag = mantel(dist_geo, dist_frag)
print(mantel_geo_frag)

# Dataframe para plot
df_geo_frag = data.frame(
  Spatial = as.vector(dist_geo),
  Fragmentation = as.vector(dist_frag)
)

# correlation figure
S_Figure_2 = ggplot(df_geo_frag, aes(x = Spatial, y = Fragmentation)) +
  geom_point(alpha = 0.3, size = 1) + 
  geom_smooth(method = "lm", se = FALSE, color = "red", linewidth = 1) +
  labs(x = "Δ km", y = "Δ LDI") +
  annotate("text",
           x = max(df_geo_frag$Spatial, na.rm = TRUE),
           y = min(df_geo_frag$Fragmentation, na.rm = TRUE),
           hjust = 1, vjust = 0,
           label = paste0("r² = ", round(mantel_geo_frag$statistic, 2),
                          "\np = ", signif(mantel_geo_frag$signif, 2)),
           size = 5) +                    
  theme_classic(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.9),
    axis.text  = element_text(color = "black", size = 14), 
    axis.title = element_text(size = 18, face = "bold")    
  )
S_Figure_2

#Rank the native species in urban

df_urban <- df_native %>%
  filter(Mosaic == "Urban")
abund_urban_natives <- df_urban %>%
  summarise(across(3:ncol(.), ~ sum(.x, na.rm = TRUE))) %>%
  pivot_longer(
    cols = everything(),
    names_to = "Species",
    values_to = "Total_abundance"
  ) %>%
  filter(Total_abundance > 0) %>%        
  arrange(desc(Total_abundance))

View(abund_urban_natives)


#===============================================================================
#Packages citations
citation("betapart")
citation("vegan")
citation("geosphere")
citation("stats")
citation("mgcv")


