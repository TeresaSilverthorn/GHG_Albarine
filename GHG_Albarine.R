#### Data and statistical analysis and visualization for GHG data collected along the Albarine River network, France in 2021 ####

# Associated with the publication: River network-scale drying impacts the spatio-temporal dynamics of greenhouse gas fluxes #

# Load in necessary packages
library(stringr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(ggpubr)
library(psych)
library(lme4) #for LMMs
library(car) #for VIF
library(MuMIn) #for model selection 
library(broom.mixed)
library(lmerTest)
library(ggpmisc) 
library(sf)
library(viridis)


#Read in all the data
#setwd("insert your working directory here")

dat <- read.csv("dat_clean.csv")

#view data
str(dat) #2035 observations of 94 variables

#load in site-average data, by GHG
CO2means <- read.csv("CO2means.csv")
CH4means <- read.csv("CH4means.csv")
N2Omeans <- read.csv("N2Omeans.csv")

str(CO2means) #280 obs. of  74 variables
str(CH4means)  #279 obs. of  81 variables 
str(N2Omeans)  #136 obs. of  74 variables 

# Table S1, a summary table of variables #
summary.dat <- dat %>%
  select(masl, Freq_Flow, Air_temp, Temp_water_C, soil_temp, sed_temp, soil_VWC, sed_VWC, pool_riffle_depth_cm, Ave.wet.width_m, Mean_Velocity..m.s.,  Discharge..l.s., pH, DO_mg_L, cond_us_cm, bedrock, boulders, cobbles, gravel, sand,  silt, Embeddedness....,  Canopy_Cover_Benth, Soil_OM, Sed_OM, Stock_Benth_g.m2 , sed_pH, sed_cond_uS_cm, DOC_ugL, DN_ugL, NO3_N_mg_l,NH4_N_ugL,  SRP_ugL) %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

write.csv (summary.dat, "C:/Users/teresa.silverthorn/Dropbox/My PC (lyp5183)/Documents/Data/R/GHGdata/summary.dat.csv")


#get CN ratio separately for sed and soil
sed_dat <- subset(dat, habitat=="Dry")

summary.sed_dat <- sed_dat %>%
  select(percent_C, percent_N) %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()

#relationship with m_to_source
cor.test(sed_dat$percent_N,sed_dat$m_to_source) 

soil_dat <- subset(dat, habitat=="Riparian")

summary.soil_dat <- soil_dat %>%
  select(percent_C, percent_N,) %>%
  psych::describe(quant=c(.25,.75)) %>%
  as_tibble(rownames="rowname")  %>%
  print()



##__________________________________________________________________________##

#### Stepwise Model selection for CO2 ####

#### For dry CO2 emissions ####

dat_CO2 <- dat %>% drop_na(CO2_mg_m2_h)  #subset the CO2 data
#Need to use data with no NAs
dat_CO2_num_dry <-  subset(dat_CO2,  habitat == "Dry")  #select dry

#Select relevant columns
dat_CO2_dry1 <- data.table(dat_CO2_num_dry[, c("campaign", "Site" , "CO2_g_m2_d", "Stock_Benth_g.m2", "CN_ratio", "percent_N" , "m_to_source" , "Freq_Flow", "sed_VWC" , "sand",  "Tsince_manual", "Sed_OM" ,"sed_temp",  "sed_pH")])

dat_CO2_dry1NA <- na.omit(dat_CO2_dry1) #remove NAs

#removing the NAs, we go from 101 to 28 observations... The large culprits are Cn ratio and percent N (only for 2/5 campaigns) #can try to fill gaps in sed temp and sed vwc

min(dat_CO2_dry1NA$CO2_g_m2_d) #-0.4031147

#remove the outlier BU01 (4.35) that we see in qqplot

dat_CO2_dry1NAout<- dat_CO2_dry1NA %>%
  filter(!CO2_g_m2_d >= 4)

#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense

#try also after scaling variables to get standardized coefficients
#log transform m to source
dat_CO2_dry1NAout$log_m_to_source <- log(dat_CO2_dry1NAout$m_to_source)

dat_CO2_dry1NAout_scaled <- dat_CO2_dry1NAout  %>% mutate_at(c("Stock_Benth_g.m2", "CN_ratio", "percent_N" , "log_m_to_source" , "Freq_Flow", "sed_VWC" , "sand",  "Tsince_manual", "Sed_OM" ,"sed_temp",  "sed_pH"), ~(scale(., center=T) %>% as.vector))


CO2_dry_lmer <- lmer(log(CO2_g_m2_d+1.4031147)~  
                       Freq_Flow*log_m_to_source + 
                       Tsince_manual*log_m_to_source+  
                       Sed_OM*log_m_to_source+ 
                       Stock_Benth_g.m2*log_m_to_source+ 
                       Sed_OM*sed_VWC + 
                       Sed_OM*sed_temp +
                       sed_VWC*sed_temp + 
                       sand+
                       sed_pH + 
                       (1 + campaign | Site), data=dat_CO2_dry1NAout_scaled) 


summary(CO2_dry_lmer) #get a lot of warnings when I use all of the m_to_source interactions
anova(CO2_dry_lmer)
plot(CO2_dry_lmer)
qqnorm(resid(CO2_dry_lmer))

vif(CO2_dry_lmer) 


options(na.action = "na.fail") # Required for dredge to run

CO2_dry_dredge<- dredge(CO2_dry_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CO2_dry_dredge)  #how many models were run 3924
head(CO2_dry_dredge)

top.CO2_dry_dredge  <- get.models(CO2_dry_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

CO2dry_mavg <- model.avg(top.CO2_dry_dredge)
mA <- summary(CO2dry_mavg)


#### For perennial-flowing CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_flow_p <-  subset(dat_CO2,  habitat == "Flowing" & drying_regime=="perennial") 

#Select relevant variables
CO2_flow_p <- data.table(dat_CO2_num_flow_p[, c("campaign", "Site" , "CO2_g_m2_d", "m_to_source" , "Stock_Benth_g.m2" , "Sed_OM", "Mean_Velocity..m.s." ,"percent.IR.up" , "Temp_water_C")])

CO2_flow_pNA <- na.omit(CO2_flow_p) #remove NAs (357 to 167, largest culprits are In sed/stock, velocity, percentIR)

min(CO2_flow_pNA$CO2_g_m2_d) #-0.01754717 #Min CO2 value in order to add to log transformation

#Scale variables in order to compare standardized coefficients
CO2_flow_pNA$log_m_to_source <- log(CO2_flow_pNA$m_to_source) #log transform
CO2_flow_pNA_scaled <- CO2_flow_pNA  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2" , "Sed_OM", "Mean_Velocity..m.s." ,"percent.IR.up" , "Temp_water_C"), ~(scale(., center=T) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CO2_flow_p_lmer <- lmer(log(CO2_g_m2_d+0.02816789 )~ 
                          log_m_to_source*Stock_Benth_g.m2 + 
                          log_m_to_source*Sed_OM + 
                          Stock_Benth_g.m2*Temp_water_C +
                          Sed_OM*Temp_water_C +
                          Mean_Velocity..m.s. + 
                          percent.IR.up +  
                          (1 + campaign | Site), data=CO2_flow_pNA_scaled) 

anova(CO2_flow_p_lmer)
plot(CO2_flow_p_lmer) #residuals improved after log transformation
qqnorm(resid(CO2_flow_p_lmer))
vif(CO2_flow_p_lmer)

options(na.action = "na.fail") # Required for dredge to run

CO2_flow_p_dredge <- dredge(CO2_flow_p_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default


nrow(CO2_flow_p_dredge)  #188 models were run
head(CO2_flow_p_dredge)

sw(CO2_flow_p_dredge)

top.CO2_flow_p_dredge  <- get.models(CO2_flow_p_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(top.CO2_flow_p_dredge) #importance of variables

CO2flow_p_mavg <- model.avg(top.CO2_flow_p_dredge)
mA_CO2_flow_p <- summary(CO2flow_p_mavg)


#### For intermittent-flowing CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_flow_i <-  subset(dat_CO2,  habitat == "Flowing" & (drying_regime=="intermittent_long" | drying_regime=="intermittent_moderate") )

#Select relevant variables
CO2_flow_i <- data.table(dat_CO2_num_flow_i[, c("campaign", "Site", "CO2_g_m2_d", "m_to_source", "Tsince_manual",  "Freq_Flow", "Stock_Benth_g.m2", "Sed_OM", "Embeddedness....", "Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C")])

CO2_flow_iNA <- na.omit(CO2_flow_i) #remove NAs (214 to 71, largest culprits are In sed/stock, velocity, tree div, algae, etc.)

min(CO2_flow_iNA$CO2_g_m2_d) #-0.003192917 #Min CO2 value in order to add to log transformation

#Scale variables in order to compare standardized coefficients
CO2_flow_iNA$log_m_to_source <- log(CO2_flow_iNA$m_to_source) #log transform
CO2_flow_iNA_scaled <- CO2_flow_iNA  %>% mutate_at(c("log_m_to_source",  "Tsince_manual",  "Freq_Flow", "Stock_Benth_g.m2", "Sed_OM", "Embeddedness....", "Mean_Velocity..m.s.", "Temp_water_C"), ~(scale(., center=T) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CO2_flow_i_lmer <- lmer(log(CO2_g_m2_d+0.03265892)~  #1.003192917  #0.03265892
                          Freq_Flow*log_m_to_source+
                          Tsince_manual*log_m_to_source + 
                          Stock_Benth_g.m2 *log_m_to_source+ 
                          Sed_OM*log_m_to_source + 
                          Temp_water_C*Stock_Benth_g.m2 +
                          Temp_water_C*Sed_OM +
                          Embeddedness.... + 
                          Mean_Velocity..m.s.+
                          (1 + campaign | Site), data=CO2_flow_iNA_scaled) 
summary(CO2_flow_i_lmer)

anova(CO2_flow_i_lmer)
plot(CO2_flow_i_lmer)  #residuals improved after log transformation
qqnorm(resid(CO2_flow_i_lmer))
vif(CO2_flow_i_lmer)

options(na.action = "na.fail") # Required for dredge to run

CO2_flow_i_dredge <- dredge(CO2_flow_i_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default


nrow(CO2_flow_i_dredge)  #1432 models were run
head(CO2_flow_i_dredge)

top.CO2_flow_i_dredge  <- get.models(CO2_flow_i_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(CO2_flow_i_dredge) #importance of variables

CO2flow_i_mavg <- model.avg(top.CO2_flow_i_dredge)
mA_CO2_flow_i <- summary(CO2flow_i_mavg)


#### For riparian CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_rip <-  subset(dat_CO2,  habitat == "Riparian") 

#Select relevant variables
CO2_rip <- data.table(dat_CO2_num_rip[, c("campaign", "Site", "CO2_g_m2_d", "soil_temp", "soil_VWC", "Soil_OM", "totalbiomass_g", "flow_state",  "Freq_Flow", "m_to_source", "percent_N")])


CO2_ripNA <- na.omit(CO2_rip) #remove NAs (643 to 312 obs, soil OM, percent N, stock etc.)

min(CO2_ripNA$CO2_g_m2_d) #-9.185739 #Min CO2 value in order to add to log transformation # after removing outlier -1.295294

#exclude the outlier of GR01, C1 (9.185738684) #2021-03-17_GR01_R2_Picarro
CO2_ripNAout<- CO2_ripNA %>%
  filter(!CO2_g_m2_d <= -9)

#try also after scaling variables to get standardized coefficients
#log transform m to source
CO2_ripNAout$log_m_to_source <- log(CO2_ripNAout$m_to_source)

CO2_ripNAout_scaled1 <- CO2_ripNAout  %>% mutate_at(c("log_m_to_source", "soil_temp", "soil_VWC", "Soil_OM",  "totalbiomass_g","Freq_Flow", "m_to_source", "percent_N"), ~(scale(., center=T) %>% as.vector))


#using variables based on hypotheses as well as top correlations and interactions
CO2_rip_lmer <- lmer(log(CO2_g_m2_d+2.295294)~  #2.295294  #1.299569
                       Freq_Flow*log_m_to_source+ 
                       flow_state *log_m_to_source+  
                       totalbiomass_g*log_m_to_source +
                       Soil_OM*log_m_to_source+
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+  
                       (1 + campaign | Site), data=CO2_ripNAout_scaled1)

summary(CO2_rip_lmer)
anova(CO2_rip_lmer)
plot(CO2_rip_lmer) 
qqnorm(resid(CO2_rip_lmer))
vif(CO2_rip_lmer)

#looks to have an outlier in the data, check Cook's
#Check outliers with Cook's Distance
cooksD <- cooks.distance(CO2_rip_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #51, 53, 54

n <- nrow(CO2_ripNAout_scaled1)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- CO2_ripNAout_scaled1[as.numeric(names_of_influential),]
CO2_ripNAout_scaled1_out <- CO2_ripNAout_scaled1 %>% anti_join(outliers)

#############

options(na.action = "na.fail") # Required for dredge to run

CO2_rip_dredge <- dredge(CO2_rip_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default


nrow(CO2_rip_dredge)  #1962 models were run
head(CO2_rip_dredge)

top.CO2_rip_dredge  <- get.models(CO2_rip_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(CO2_rip_dredge) #importance of variables

CO2rip_mavg <- model.avg(top.CO2_rip_dredge)
summary(CO2rip_mavg) 



#### CO2: Do aquatic conditions influence riparian fluxes? ####

#Plot relationship between riparian and in-stream emissions

CO2_rip_aq <- ggplot(subset(CO2means, habitat=="Riparian"& flow_state!="pool"), aes(In_stream_CO2_g_m2_d, mean_CO2_g_m2_d)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(),  legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black"))  +
  xlab(expression(`In-stream`~g~CO[2]*`-C`~m^-2*~d^-1))  +  ylab(expression(`Riparian`~g~CO[2]*`-C`~m^-2*~d^-1)) +
  geom_smooth(aes(colour=flow_state), method = "lm", linetype="dashed", se=FALSE,formula = y ~ x)  + scale_colour_manual(values=c( "#FBD260", "#5198B7")) 
#stat_regline_equation(label.y = 15) + stat_cor(label.y = 12)
CO2_rip_aq


#Run linear model of relationship between in-stream and riparian emissions

min(CO2means$mean_CO2_g_m2_d) #-0.4350508
#log(mean_CO2_g_m2_d + 1.4350508)
rip_aq_connect <- lmer(mean_CO2_g_m2_d~In_stream_CO2_g_m2_d*flow_state*m_to_source + (1 | campaign), data=subset(CO2means, habitat=="Riparian"& flow_state!="pool"))
summary(rip_aq_connect) 

r.squaredGLMM(rip_aq_connect)

plot(rip_aq_connect) #could be better, but a lot worse with log transformation
qqnorm(resid(rip_aq_connect))



##__________________________________________________________________________##

#### Stepwise Model selection for CH4 ####

#### For dry CH4 emissions ####

#Need to use data with no NAs
dat_CH4 <- dat %>% drop_na(CH4_ug_m2_h)  #subset the CH4 data
dat_CH4_num_dry <-  subset(dat_CH4,  habitat == "Dry")  #select dry
dat_CH4_dry1 <- data.table(dat_CH4_num_dry[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2",  "percent_N" , "m_to_source" , "Freq_Flow","sed_VWC" , "sed_pH", "Tsince_manual",   "Sed_OM" ,"sed_temp",  "sed_pH", "Canopy_Cover_Benth" )]) #Select relevant columns
dat_CH4_dry1NA <- na.omit(dat_CH4_dry1) #remove NAs

#remove outlier BU01
dat_CH4_dry1NA_out <- subset(dat_CH4_dry1NA, CH4_mg_m2_d>=-0.505900632)

sd(dat_CH4_num_dry$CH4_mg_m2_d)

#removing the NAs, we go from 100 to 42 observations... The large culprits are Cn ratio and percent N (only for 2/5 campaigns) #can try to fill gaps in sed temp and sed vwc

min(dat_CH4_dry1NA_out$CH4_mg_m2_d) #-0.2430066

dat_CH4_dry1NA_out$log_m_to_source <- log(dat_CH4_dry1NA_out$m_to_source)

dat_CH4_dry1NA_out_scaled <- dat_CH4_dry1NA_out  %>% mutate_at(c("Stock_Benth_g.m2",  "percent_N" , "log_m_to_source" , "Freq_Flow","sed_VWC" , "sed_pH", "Tsince_manual",   "Sed_OM" ,"sed_temp",  "sed_pH", "Canopy_Cover_Benth"), ~(scale(., center=T) %>% as.vector))


#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
#log(CH4_mg_m2_d+1.5059006) residuals not improve by log transforming
CH4_dry_lmer <- lmer(log(CH4_mg_m2_d+1.5059006) ~ 
                       Freq_Flow*log_m_to_source + 
                       Tsince_manual*log_m_to_source +  
                       Sed_OM*log_m_to_source + 
                       Stock_Benth_g.m2*log_m_to_source + 
                       sed_VWC*Sed_OM +  
                       Sed_OM*sed_temp + 
                       sed_VWC*sed_temp +
                       Canopy_Cover_Benth + 
                       sed_pH +
                       (1 + campaign | Site), data=dat_CH4_dry1NA_out_scaled) 

summary(CH4_dry_lmer) #get a lot of warnings when I use all of the m_to_source interactions
anova(CH4_dry_lmer)
plot(CH4_dry_lmer, id=0.05) #residuals look fine
qqnorm(resid(CH4_dry_lmer))


options(na.action = "na.fail") # Required for dredge to run

CH4_dry_dredge<- dredge(CH4_dry_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CH4_dry_dredge)  #how many models were run: 3924
head(CH4_dry_dredge)

top.CH4_dry_dredge  <- get.models(CH4_dry_dredge , subset=delta <= 2) #Top model is the null model with no models <2 AIC


#### For flowing perennial CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_flow <-  subset(dat_CH4,  habitat == "Flowing")  #select flowing
dat_CH4_num_flow_p <-  subset(dat_CH4_num_flow,  drying_regime == "perennial") #select perennial
dat_CH4_num_flow_p1 <- data.table(dat_CH4_num_flow_p[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" , "Freq_Flow",  "Temp_water_C", "Air_temp", "DO_mg_L", "Mean_Velocity..m.s.", "percent.IR.up" )]) #Select relevant columns
dat_CH4_num_flow_p1NA <- na.omit(dat_CH4_num_flow_p1) #remove NAs #Go from 353 to 166

dat_CH4_num_flow_p1NA$log_m_to_source <- log(dat_CH4_num_flow_p1NA$m_to_source) #log transform
dat_CH4_num_flow_p1NA_out_scaled1 <- dat_CH4_num_flow_p1NA  %>% mutate_at(c("Stock_Benth_g.m2", "Sed_OM", "log_m_to_source" , "Freq_Flow",  "Temp_water_C", "Air_temp", "DO_mg_L", "Mean_Velocity..m.s.", "percent.IR.up" ), ~(scale(., center=T) %>% as.vector))


#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
min(dat_CH4_num_flow_p1NA_out_scaled$CH4_mg_m2_d) #-0.7679176
min(subset(dat_CH4_num_flow_p1NA_out_scaled$CH4_mg_m2_d, dat_CH4_num_flow_p1NA_out_scaled$CH4_mg_m2_d > 0))  # 0.0005324891
#try smallest positive value - the negative value:
0.0005324891+0.7679176  #0.001064978

CH4_flow_p_lmer <- lmer(log(CH4_mg_m2_d+0.7684501) ~  #1.7679176  #0.7684501
                          log_m_to_source*Stock_Benth_g.m2 + 
                          log_m_to_source*Sed_OM + 
                          percent.IR.up +
                          Temp_water_C*Stock_Benth_g.m2 + 
                          Temp_water_C* Sed_OM  + 
                          DO_mg_L*Stock_Benth_g.m2 +  
                          DO_mg_L*Sed_OM + 
                          Mean_Velocity..m.s. +
                          (1 + campaign | Site), data=dat_CH4_num_flow_p1NA_out_scaled1)

summary(CH4_flow_p_lmer) #get a lot of warnings when I use all of the m_to_source interactions
anova(CH4_flow_p_lmer)
plot(CH4_flow_p_lmer) #  
qqnorm(resid(CH4_flow_p_lmer))

#Check outliers with Cook's Distance
cooksD <- cooks.distance(CH4_flow_p_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #158

n <- nrow(dat_CH4_num_flow_p1NA_out_scaled2)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_CH4_num_flow_p1NA_out_scaled2[as.numeric(names_of_influential),]
dat_CH4_num_flow_p1NA_out_scaled2_out <- dat_CH4_num_flow_p1NA_out_scaled1 %>% anti_join(outliers)



options(na.action = "na.fail") # Required for dredge to run

CH4_flow_p_dredge<- dredge(CH4_flow_p_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CH4_flow_p_dredge)  #how many models were run: 748
head(CH4_flow_p_dredge)

top.CH4_flow_p_dredge  <- get.models(CH4_flow_p_dredge , subset=delta <= 2) #The top models

sw(CH4_flow_p_dredge)  # Variable importance

CH4flow_p_mavg <- model.avg(top.CH4_flow_p_dredge)
mA_CH4_flow_p <- summary(CH4flow_p_mavg)






#### For flowing intermittent CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_flow <-  subset(dat_CH4,  habitat == "Flowing")  #select flowing
dat_CH4_num_flow_i <-  subset(dat_CH4_num_flow,  drying_regime == "intermittent_long" | drying_regime == "intermittent_moderate" ) #select intermittnet
dat_CH4_num_flow_i1 <- data.table(dat_CH4_num_flow_i[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" , "Tsince_manual",  "Freq_Flow", "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL",  "DOC_ugL")]) #Select relevant columns
dat_CH4_num_flow_i1NA <- na.omit(dat_CH4_num_flow_i1) #remove NAs #Go from 214 to 93 obs


dat_CH4_num_flow_i1NA$log_m_to_source <- log(dat_CH4_num_flow_i1NA$m_to_source)

dat_CH4_num_flow_i1NA_scaled1 <- dat_CH4_num_flow_i1NA  %>% mutate_at(c("Stock_Benth_g.m2", "Sed_OM", "log_m_to_source" ,"Tsince_manual", "Freq_Flow", "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL", "DOC_ugL"), ~(scale(., center=T) %>% as.vector))

#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
min(dat_CH4_num_flow_i1NA_scaled1$CH4_mg_m2_d) #-0.2719916
min(subset(dat_CH4_num_flow_i1NA_scaled1$CH4_mg_m2_d, dat_CH4_num_flow_i1NA_scaled1$CH4_mg_m2_d > 0))  # 0.01023246
#try smallest positive value - the negative value:
0.01023246+0.2719916  #0.2822241

CH4_flow_i_lmer <- lmer(log(CH4_mg_m2_d+1.2719916) ~  #0.2822241  #1.2719916
                          log_m_to_source*Freq_Flow +
                          log_m_to_source*Tsince_manual +
                          log_m_to_source*Sed_OM + 
                          log_m_to_source*Stock_Benth_g.m2 + 
                          Stock_Benth_g.m2* Temp_water_C + 
                          Sed_OM* Temp_water_C  + 
                          DO_mg_L + 
                          Mean_Velocity..m.s. + 
                          NH4_N_ugL+
                          DOC_ugL+ 
                          (1 + campaign | Site), data=dat_CH4_num_flow_i1NA_scaled1)

summary(CH4_flow_i_lmer)  # A lot of model warnings
anova(CH4_flow_i_lmer)
plot(CH4_flow_i_lmer)  #OK
qqnorm(resid(CH4_flow_i_lmer))

#Run the dredge function on the saturated model to determine the best model
options(na.action = "na.fail") # Required for dredge to run

CH4_flow_i_dredge<- dredge(CH4_flow_i_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CH4_flow_i_dredge)  #how many models were run: 5728
head(CH4_flow_i_dredge)

top.CH4_flow_i_dredge  <- get.models(CH4_flow_i_dredge , subset=delta <= 2) #The top models

sw(top.CH4_flow_i_dredge)  # Variable importance

CH4flow_i_mavg <- model.avg(top.CH4_flow_i_dredge)
mA_CH4_flow_i <- summary(CH4flow_i_mavg)




#### For riparian CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_rip <-  subset(dat_CH4,  habitat == "Riparian") 

#Select relevant variables
CH4_rip <- data.table(dat_CH4_num_rip[, c("campaign", "Site", "CH4_mg_m2_d", "soil_temp", "soil_VWC", "Soil_OM", "flow_state" ,   "totalbiomass_g",  "Freq_Flow", "m_to_source", "percent_N")])
CH4_ripNA <- na.omit(CH4_rip) #remove NAs (646 to 316 obs, soil OM, percent N, stock etc.)

#remove  outliers identified in residuals plot
CH4_ripNA_out <- CH4_ripNA %>%
  filter(!CH4_mg_m2_d <= -5)

min(CH4_ripNA_out$CH4_mg_m2_d) #-1.842777 #Min CH4 value in order to add to log transformation
min(CH4_ripNA$CH4_mg_m2_d)     #-1.842777 #outlier must have been removed earlier

#try also after scaling variables to get standardized coefficients
#log transform m to source
CH4_ripNA_out$log_m_to_source <- log(CH4_ripNA_out$m_to_source)

CH4_ripNA_scaled1 <- CH4_ripNA_out  %>% mutate_at(c("log_m_to_source",  "soil_temp", "soil_VWC", "Soil_OM",  "totalbiomass_g",  "Freq_Flow",  "percent_N"), ~(scale(., center=T) %>% as.vector))


# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
min(CH4_ripNA_scaled$CH4_mg_m2_d) #-1.842777
min(subset(CH4_ripNA_scaled$CH4_mg_m2_d, CH4_ripNA_scaled$CH4_mg_m2_d > 0))  # 0.02682259
#try smallest positive value - the negative value:
0.02682259+1.842777  #1.8696

CH4_rip_lmer <- lmer(log(CH4_mg_m2_d+2.842777)~  #2.842777   #1.8696
                       Freq_Flow*log_m_to_source+ 
                       flow_state *log_m_to_source+  
                       totalbiomass_g*log_m_to_source +
                       Soil_OM*log_m_to_source+ 
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+
                       (1 + campaign | Site), data=CH4_ripNA_scaled1)


summary(CH4_rip_lmer)
anova(CH4_rip_lmer)
plot(CH4_rip_lmer) #add id = 0.05 to label outliers
qqnorm(resid(CH4_rip_lmer))
vif(CH4_rip_lmer)

### run dredge 

options(na.action = "na.fail") # Required for dredge to run

CH4_rip_dredge <- dredge(CH4_rip_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CH4_rip_dredge)  #1962 models were run
head(CH4_rip_dredge)

top.CH4_rip_dredge  <- get.models(CH4_rip_dredge , subset=delta <= 2) #Only 1 top model

sw(CH4_rip_dredge) #importance of variables

CH4rip_mavg <- model.avg(top.CH4_rip_dredge) #only one top model
mA_CH4_flow_i <- summary(CH4rip_mavg)

#Run model with soil temp and soil VWC (the variables from the single top model)
CH4_rip_lmer2 <- lmer(log(CH4_mg_m2_d+2.842777)~ soil_VWC+ soil_temp + (1 + campaign | Site), data=CH4_ripNA_scaled)
summary(CH4_rip_lmer2)

r.squaredGLMM(CH4_rip_lmer2)




#### CH4: Do aquatic conditions influence riparian fluxes? ####

#Plot the riparian vs in-stream GHG emissions

CH4_rip_aq <- ggplot(subset(CH4means, habitat=="Riparian" & flow_state!="pool"), aes( In_stream_CH4_mg_m2_d, mean_CH4_mg_m2_d)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right",axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black"))  + 
  xlab(expression(`In-stream`~mg~CH[4]*`-C`~m^-2~d^-1))  +  ylab(expression(`Riparian`~mg~CH[4]*`-C`~m^-2~d^-1)) +
  geom_smooth(aes(colour=flow_state), method = "lm", se=FALSE, linetype="dashed", formula = y ~ x) + scale_colour_manual(values=c( "#FBD260", "#5198B7"))   
#+ stat_regline_equation(label.y = 15) + stat_cor(label.y = 14)
CH4_rip_aq


#Run linear model of relationship between in-stream and riparian emissions
min(CH4means$mean_CH4_mg_m2_d) # -1.398022

rip_aq_connect <- lmer(mean_CH4_mg_m2_d~ In_stream_CH4_mg_m2_d*flow_state*m_to_source + (1 | campaign), data=subset(CH4means, habitat=="Riparian"& flow_state!="pool"))
summary(rip_aq_connect) 

r.squaredGLMM(rip_aq_connect)
plot(rip_aq_connect) 
qqnorm(resid(rip_aq_connect)) #looks good



#### Stepwise Model selection for N2O ####

#### For dry N2O emissions ####

#Need to use data with no NAs
dat_N2O <- dat %>% drop_na(mg_N2ON_m2_d) 
dat_N2O_dry <-  subset(dat_N2O,  habitat == "Dry")  #select dry 49 obs
#Select relevant columns
dat_N2O_dry1 <- data.table(dat_N2O_dry[, c("campaign", "Site" , "mg_N2ON_m2_d",  "m_to_source" , "Freq_Flow",  "Stock_Benth_g.m2",  "Sed_OM" ,  "sed_VWC","sed_temp",  "percent_N" , "MaxDur_Ponded")])
dat_N2O_dry1NA <- na.omit(dat_N2O_dry1) #remove NAs 49 to 34 obs, missing values in Tsince manual make it 19 obs (too many NAs for BR01, and GR01, not possible to fill in)

#remove the outliers we see in the residual plot Al02, >0.1
dat_N2O_dry1NAout<- dat_N2O_dry1NA %>%
  filter(!mg_N2ON_m2_d >= 0.1)

min(dat_N2O_dry1NAout$mg_N2ON_m2_d) #-0.07417401


dat_N2O_dry1NAout$log_m_to_source <- log(dat_N2O_dry1NAout$m_to_source)

dat_N2O_dry1NAout_scaled1 <- dat_N2O_dry1NAout  %>% mutate_at(c("log_m_to_source" , "Freq_Flow",  "Stock_Benth_g.m2",  "Sed_OM" ,  "sed_VWC","sed_temp", "percent_N" , "MaxDur_Ponded"), ~(scale(., center=T) %>% as.vector))


min(dat_N2O_dry1NAout_scaled$mg_N2ON_m2_d) #-0.07417401
min(subset(dat_N2O_dry1NAout_scaled$mg_N2ON_m2_d, dat_N2O_dry1NAout_scaled$mg_N2ON_m2_d > 0))  # 0.00260448
#try smallest positive value - the negative value:
0.00260448+0.07417401  #0.07677849


N2O_dry_lmer <- lmer(log(mg_N2ON_m2_d+1.07417401)~     #0.07677849  #1.07417401
                       Freq_Flow*log_m_to_source + 
                       Sed_OM*log_m_to_source+     #removed Tsince bc too few data points
                       Stock_Benth_g.m2*log_m_to_source+ 
                       Sed_OM*sed_VWC + 
                       Sed_OM*sed_temp +
                       sed_VWC*sed_temp + 
                       percent_N + 
                       MaxDur_Ponded +
                       (1 + campaign | Site), data=dat_N2O_dry1NAout_scaled1) 

summary(N2O_dry_lmer) 
anova(N2O_dry_lmer)
plot(N2O_dry_lmer) #Looks good after we removed the Al02 outlier
qqnorm(resid(N2O_dry_lmer))

head(dat_N2O_dry1NAout_scaled1)


#Run the dredge function on the saturated model to determine the best model
options(na.action = "na.fail") # Required for dredge to run

N2O_dry_dredge<- dredge(N2O_dry_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(N2O_dry_dredge)  #how many models were run: 1404
head(N2O_dry_dredge)

top.N2O_dry_dredge  <- get.models(N2O_dry_dredge , subset=delta <= 2) #Top model is null model 


### For perennial flowing N2O emissions ####

#Need to use data with no NAs
dat_N2O_flow <-  subset(dat_N2O,  habitat == "Flowing")  #select flowing
dat_N2O_flow_p <-  subset(dat_N2O_flow,  drying_regime == "perennial")

dat_N2O_flow_p1 <- data.table(dat_N2O_flow_p[, c("campaign", "Site" , "mg_N2ON_m2_d", "m_to_source",  "Stock_Benth_g.m2",  "Sed_OM","Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C", "NH4_N_ugL" ,  "DN_ugL" ,"NO3_N_mg_l", "percent.IR.up")]) #Select relevant columns
dat_N2O_flow_p1NA <- na.omit(dat_N2O_flow_p1) #remove NAs  # 151, from 163


min(dat_N2O_flow_p1NA$mg_N2ON_m2_d) #-0.06030342

#Scale the variables
dat_N2O_flow_p1NA$log_m_to_source <- log(dat_N2O_flow_p1NA$m_to_source)

dat_N2O_flow_p1NA_scaled1 <- dat_N2O_flow_p1NA  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2",  "Sed_OM","Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C",  "DN_ugL" ,"NO3_N_mg_l", "percent.IR.up"), ~(scale(., center=T) %>% as.vector))


min(dat_N2O_flow_p1NA_scaled$mg_N2ON_m2_d) #-0.06030342
min(subset(dat_N2O_flow_p1NA_scaled$mg_N2ON_m2_d, dat_N2O_flow_p1NA_scaled$mg_N2ON_m2_d > 0))  # 0.0005535944
#try smallest positive value - the negative value:
0.0005535944+0.06030342  #0.06085701

#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
N2O_flow_p_lmer <- lmer(log(mg_N2ON_m2_d+0.06085701)~   #1.06030342  #0.06085701
                          log_m_to_source*Stock_Benth_g.m2 + 
                          log_m_to_source*Sed_OM + 
                          percent.IR.up + 
                          Temp_water_C*Stock_Benth_g.m2 + 
                          Temp_water_C*Sed_OM  + 
                          Mean_Velocity..m.s. +
                          DO_mg_L+
                          NO3_N_mg_l +
                          DN_ugL + 
                          (1 + campaign | Site), data=dat_N2O_flow_p1NA_scaled1_out) 

summary(N2O_flow_p_lmer) 
anova(N2O_flow_p_lmer)
plot(N2O_flow_p_lmer) 
qqnorm(resid(N2O_flow_p_lmer))

#Check outliers with Cook's Distance
cooksD <- cooks.distance(N2O_flow_p_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #12

n <- nrow(dat_N2O_flow_p1NA_scaled1)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_N2O_flow_p1NA_scaled1[as.numeric(names_of_influential),]
dat_N2O_flow_p1NA_scaled1_out <- dat_N2O_flow_p1NA_scaled1 %>% anti_join(outliers)


#Run the dredge function on the saturated model to determine the best model
options(na.action = "na.fail") # Required for dredge to run

N2O_flow_p_dredge<- dredge(N2O_flow_p_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(N2O_flow_p_dredge)  #how many models were run: 1504
head(N2O_flow_p_dredge)

top.N2O_flow_p_dredge  <- get.models(N2O_flow_p_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 away from the top model 

sw(N2O_flow_p_dredge)  # Variable importance

N2Oflow_p_mavg <- model.avg(top.N2O_flow_p_dredge) 
summary(N2Oflow_p_mavg)



#### For flowing intermittent N2O emissions ####

#Need to use data with no NAs
dat_N2O_flow <-  subset(dat_N2O,  habitat == "Flowing")  #select flowing
dat_N2O_flow_i <-  subset(dat_N2O_flow,  drying_regime == "intermittent_long" | drying_regime == "intermittent_moderate" ) #select intermittent
dat_N2O_flow_i1 <- data.table(dat_N2O_flow_i[, c("campaign", "Site" , "mg_N2ON_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" ,"Tsince_manual", "Freq_Flow",  "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL", "NO3_N_mg_l",  "DN_ugL" )]) #Select relevant columns


dat_N2O_flow_i1NA <- na.omit(dat_N2O_flow_i1) #remove NAs #Go from 278 to 243, 71 with manual tsince rewetting

dat_N2O_flow_i1NAout <- dat_N2O_flow_i1NA[-72,]


min(dat_N2O_flow_i1NAout$mg_N2ON_m2_d) #-0.04485442

#Scale the variables
dat_N2O_flow_i1NAout$log_m_to_source <- log(dat_N2O_flow_i1NAout$m_to_source)

dat_N2O_flow_i1NAout_scaled1 <- dat_N2O_flow_i1NAout  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2", "Sed_OM", "m_to_source" ,"Tsince_manual", "Freq_Flow",  "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL", "NO3_N_mg_l",  "DN_ugL"), ~(scale(., center=T) %>% as.vector))

#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
#log(mg_N2ON_m2_d+1.04485442)

min(dat_N2O_flow_i1NAout_scaled1$mg_N2ON_m2_d) #-0.04485442
min(subset(dat_N2O_flow_i1NAout_scaled1$mg_N2ON_m2_d, dat_N2O_flow_i1NAout_scaled1$mg_N2ON_m2_d > 0))  # 0.001925746
#try smallest positive value - the negative value:
0.001925746+0.04485442  #0.06085701 

N2O_flow_i_lmer <- lmer(log(mg_N2ON_m2_d+0.04678017) ~   #1.04485442  #0.04678017
                          log_m_to_source*Freq_Flow +
                          log_m_to_source*Tsince_manual +
                          log_m_to_source*Sed_OM + 
                          log_m_to_source*Stock_Benth_g.m2 + 
                          Temp_water_C*Sed_OM +
                          Temp_water_C*Stock_Benth_g.m2 +
                          Mean_Velocity..m.s. + 
                          NO3_N_mg_l +
                          DN_ugL + 
                          (1 + campaign | Site), data=dat_N2O_flow_i1NAout_scaled1_out)


summary(N2O_flow_i_lmer) 
anova(N2O_flow_i_lmer)
plot(N2O_flow_i_lmer) 
qqnorm(resid(N2O_flow_i_lmer))

#Check outliers with Cook's Distance
cooksD <- cooks.distance(N2O_flow_i_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #33

n <- nrow(dat_N2O_flow_i1NAout_scaled1)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_N2O_flow_i1NAout_scaled1[as.numeric(names_of_influential),]
dat_N2O_flow_i1NAout_scaled1_out <- dat_N2O_flow_i1NAout_scaled1 %>% anti_join(outliers)


#Run the dredge function on the saturated model to determine the best model
options(na.action = "na.fail") # Required for dredge to run

N2O_flow_i_dredge<- dredge(N2O_flow_i_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(N2O_flow_i_dredge)  #how many models were run:  2864
head(N2O_flow_i_dredge)

top.N2O_flow_i_dredge  <- get.models(N2O_flow_i_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(N2O_flow_i_dredge)  # Variable importance

N2Oflow_i_mavg <- model.avg(top.N2O_flow_i_dredge) 
summary(N2Oflow_i_mavg)



#### For riparian N2O emissions ####

#Need to use data with no NAs
dat_N2O_rip <-  subset(dat_N2O,  habitat == "Riparian") 

#Select relevant variables
dat_N2O_rip1 <- data.table(dat_N2O_rip[, c("campaign", "Site", "mg_N2ON_m2_d", "soil_temp", "soil_VWC", "Soil_OM",  "flow_state" ,   "totalbiomass_g",   "Freq_Flow", "m_to_source", "percent_N")])

dat_N2O_ripNA <- na.omit(dat_N2O_rip1) #remove NAs  280 from 323

min(dat_N2O_ripNA$mg_N2ON_m2_d) #-0.1833087 #Min N2O value in order to add to log 

#remove two outliers based on residual plot (AL07 and AL03)

#remove 2 outliers identified in residuals plot
dat_N2O_ripNAout <- dat_N2O_ripNA %>%
  filter(!mg_N2ON_m2_d >= 1)

#Scale the variables
dat_N2O_ripNAout$log_m_to_source <- log(dat_N2O_ripNAout$m_to_source)

dat_N2O_ripNAout_scaled1 <- dat_N2O_ripNAout  %>% mutate_at(c("log_m_to_source", "soil_temp", "soil_VWC", "Soil_OM",  "totalbiomass_g",   "Freq_Flow",  "percent_N"), ~(scale(., center=T) %>% as.vector))



min(dat_N2O_ripNAout_scaled$mg_N2ON_m2_d) #-0.1833087
min(subset(dat_N2O_ripNAout_scaled$mg_N2ON_m2_d, dat_N2O_ripNAout_scaled$mg_N2ON_m2_d > 0))  # 0.0004272289
#try smallest positive value - the negative value:
0.0004272289+0.1833087  #0.1837359  

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
N2O_rip_lmer <- lmer(log(mg_N2ON_m2_d+0.1837359)~  #0.1837359  #1.1833087
                       Freq_Flow*log_m_to_source+ 
                       flow_state*log_m_to_source+
                       totalbiomass_g*log_m_to_source+
                       Soil_OM*log_m_to_source+ 
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+
                       (1 + campaign | Site), data=dat_N2O_ripNAout_scaled1_out)

summary(N2O_rip_lmer)
anova(N2O_rip_lmer)
plot(N2O_rip_lmer) #add id = 0.05 to label outliers
qqnorm(resid(N2O_rip_lmer))



#Check outliers with Cook's Distance
cooksD <- cooks.distance(N2O_rip_lmer) #run the model with dat_sub first
influential <- cooksD[(cooksD > (15* mean(cooksD, na.rm = TRUE)))]
influential #236

n <- nrow(dat_N2O_ripNAout_scaled1)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 28/n, lty = 2, col = "steelblue") # add cutoff line

names_of_influential <- names(influential)
outliers <- dat_N2O_ripNAout_scaled1[as.numeric(names_of_influential),]
dat_N2O_ripNAout_scaled1_out <- dat_N2O_ripNAout_scaled1 %>% anti_join(outliers)

#########


#Run the dredge function on the saturated model to determine the best model
options(na.action = "na.fail") # Required for dredge to run

N2O_rip_dredge<- dredge(N2O_rip_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(N2O_rip_dredge)  #how many models were run:  1962
head(N2O_rip_dredge)

top.N2O_rip_dredge  <- get.models(N2O_rip_dredge , subset=delta <= 2) #The top models, which include those with an AIC >2 away from the top model 

sw(N2O_rip_dredge)  # Variable importance

N2Orip_mavg <- model.avg(top.N2O_rip_dredge)
mA_N2O_rip <- summary(N2Orip_mavg)


#### N2O: Do aquatic conditions influence riparian fluxes? ####


# Plot riparian vs in-stream N2O fluxes

N2O_rip_aq <- ggplot(subset(N2Omeans, habitat=="Riparian"& flow_state!="pool"), aes(In_stream_N2O_mg_m2_d, mean_mg_N2ON_m2_d)) + geom_point(aes(colour=flow_state), size=3.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(),  legend.justification = c("left", "top"), legend.box.just = "right",axis.line = element_line(colour = "black"), text = element_text(size = 20), axis.text = element_text(size = 20, colour="black"))   + 
  xlab(expression(`In-stream`~mg~N[2]*`O-`*N~m^-2~d^-1))  +  ylab(expression(`Riparian`~mg~N[2]*`O-`*N~m^-2~d^-1)) +  
  geom_smooth(aes(colour=flow_state), method = "lm", linetype="dashed", se=FALSE,  formula = y ~ x)  + scale_colour_manual(values=c( "#FBD260", "#5198B7"))   
N2O_rip_aq


#Run linear model of relationship between in-stream and riparian emissions

min(N2Omeans_connect$mean_N2O_mg_m2_d) # -1.398022

#remove outlier based on residual plots

N2Omeans_connect1 <- subset(N2Omeans_connect, !(Site=="AL07" & campaign=="1" & habitat=="Riparian"))

rip_aq_connect <- lmer(log(mean_mg_N2ON_m2_d+2.398022)~ In_stream_N2O_mg_m2_d*flow_state*m_to_source + (1 | campaign), data=subset(N2Omeans_connect1, habitat=="Riparian"& flow_state!="pool"))
summary(rip_aq_connect)

plot(rip_aq_connect, id=0.001)
outs <- influencePlot(rip_aq_connect)
qqnorm(resid(rip_aq_connect)) #log transforming does not improve, remove outlier

r.squaredGLMM(rip_aq_connect)

##_________________________________________________________________________##
#### Make figures for publication ####

# GHG flux by habitat using site means

#### CO2 by habitat using Site means ####
CO2means$habitat <- factor(CO2means$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

CO2_habitat_means <- ggplot(CO2means, aes(x=habitat, y=mean_CO2_g_m2_d)) +  theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.4, position = position_jitterdodge(jitter.width=2.0)) +  geom_boxplot( fill=NA, outlier.shape = NA) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(),axis.text.x=element_text(angle = 45, hjust = .9), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c( "#878B8C","#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C","#5598B7", "#FCD262" , "#527028" ))  + ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))  +  
  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) +
  annotate("text", x = 1.8, y = 17, label=sprintf('\u2191'), size=3) +
  #scale_y_continuous(trans='log2') +
  annotate("text", x = 2, y = 17, label="30",  size=3)  + ylim(-0.5, 21.7) +
  geom_bracket(
    xmin = c(1, 1,2, 3), xmax = c(2, 4,3, 4),
    y.position = c(21.1, 20.2, 18, 17), label = c(" ", " ", " ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x =1.5, y = 21.2 , label = "**") +
  annotate("text", x =2.8, y = 20.3 , label = "**") +
  annotate("text", x =2.5, y = 18.1 , label = "***") +
  annotate("text", x =3.5, y = 17.1 , label = "***") 
CO2_habitat_means 


#### CH4 by habitat using Site means ####

CH4means$habitat <- factor(CH4means$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

CH4_habitat_means <- ggplot(CH4means, aes(x=habitat, y=mean_CH4_mg_m2_d))  +   theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.4, position = position_jitterdodge(jitter.width=2.0)) + geom_boxplot(fill=NA, outlier.shape = NA) +  geom_hline(yintercept=0, linetype="dotted", color = "black", size=0.5) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.text.x=element_text(angle = 45, hjust = .9), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))  +  
  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(-1.4, 8.3) +
  annotate("text", x = 2.1, y = 7, label="11.2", size=3)  + 
  annotate("text", x = 1.8, y = 7.02, label=sprintf('\u2191'), size=3) +
  geom_bracket(
    xmin = c(2, 3), xmax = c(4,4),
    y.position = c(8, 7), label = c(" ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x =3, y = 8.1 , label = "***") +
  annotate("text", x =3.5, y = 7.1 , label = "***") 
CH4_habitat_means 

#### N2O by habitat using Site means ####
N2Omeans$habitat <- factor(N2Omeans$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

N2O_habitat_means <- ggplot(N2Omeans, aes(x=habitat, y=mean_mg_N2ON_m2_d)) +  theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.4, position = position_jitterdodge(jitter.width=2.0)) +  geom_boxplot(fill=NA, outlier.shape = NA) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.text.x=element_text(angle = 45, hjust = .9), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c( "#878B8C", "#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028", "#878B8C"))  + ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) +  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(-0.03, 1.45) +
  annotate("text", x = 1.6, y = 1.105, label=sprintf('\u2191'), size=3) +
  annotate("text", x = 2.15, y = 1.1, label= "1.7, 2.6", size=3) +
  geom_bracket(
    xmin = c(2, 2), xmax = c(3,4),
    y.position = c(1.39, 1.25), label = c(" ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x = 2.5, y = 1.40, label= "*") +
  annotate("text", x = 3, y = 1.26, label= "***") 
N2O_habitat_means 

## combine plots ##
tiff("MeanGHGbyhabitat", units="in", width=3, height=8.5, res=300)

MeanGHGbyhabitat <- ggarrange(CO2_habitat_means, CH4_habitat_means, N2O_habitat_means,
                              labels = c("A", "B", "C"),
                              ncol = 1, nrow = 3, align="v")
MeanGHGbyhabitat

dev.off()


#### CO2 flowing by drying regime ####

CO2means$drying_regime <- dplyr::recode(CO2means$drying_regime, intermittent_long = 'Int-S', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten

CO2means$drying_regime <- factor(CO2means$drying_regime, levels = c("Int-S", "Int-M", "Per")) #reorder


CO2_drying_regime  <- ggplot(subset(CO2means, habitat=="Flowing"), aes(x=drying_regime, y=mean_CO2_g_m2_d)) +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.6, position = position_jitterdodge(jitter.width=1.8)) + geom_boxplot( outlier.shape = NA, fill=NA) +  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1"))  + ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))  +  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(0.082, 32) +
  geom_bracket(
    xmin = c(1), xmax = c(3),
    y.position = c(30.5), label = c(" ", " "),
    tip.length = 0.01,  label.size = 5) +
  annotate("text", x =2, y = 30.9 , label = "*", size=5) 
CO2_drying_regime 


#### CH4 flowing by drying regime ####

CH4means$drying_regime <- dplyr::recode(CH4means$drying_regime, intermittent_long = 'Int-S', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten

CH4means$drying_regime <- factor(CH4means$drying_regime, levels = c("Int-S", "Int-M", "Per")) #reorder

CH4_drying_regime  <- ggplot(subset(CH4means, habitat=="Flowing"), aes(x=drying_regime, y=mean_CH4_mg_m2_d))  +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.6, position = position_jitterdodge(jitter.width=1.8)) + geom_boxplot( outlier.shape = NA, fill=NA) +  theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) 
CH4_drying_regime



#### N2O flowing by drying regime ####

N2Omeans$drying_regime <- dplyr::recode(N2Omeans$drying_regime, intermittent_long = 'Int-S', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten

N2Omeans$drying_regime <- factor(N2Omeans$drying_regime, levels = c("Int-S", "Int-M", "Per")) #reorder factors


N2O_drying_regime  <- ggplot(subset(N2Omeans, habitat=="Flowing"), aes(x=drying_regime, y=mean_mg_N2ON_m2_d)) +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.6, position = position_jitterdodge(jitter.width=1.8)) + geom_boxplot(outlier.shape = NA, fill=NA) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))   +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1"))  + ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) + theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) 
N2O_drying_regime 

#Flowing by drying regime CH4_drying_regime
tiff("FlowingbyDryingRegime", units="in", width=3, height=8, res=300)

FlowingbyDryingRegime <- ggarrange(CO2_drying_regime, CH4_drying_regime, N2O_drying_regime,
                                   labels = c("A", "B", "C"),
                                   ncol = 1, nrow = 3)
FlowingbyDryingRegime

dev.off()



#### Plot model coefficients for dry sediments/riparian soils ####

#Read in the model coefficients for dry sediment and riparian GHG fluxes
combined_data_dry <- read.csv("Model_coefficients_GHGdry.csv")
combined_data_rip <- read.csv("Model_coefficients_GHGrip.csv")

combined_data_dry$label <- as.factor(combined_data_dry$label)

levels(combined_data_dry$label)[levels(combined_data_dry$label) == "Sed_moist"] <- "Moisture"

# Plot the combined data using ggplot
rf_dry <- ggplot(combined_data_dry, aes(x=label, y=y, ymin=y-stderror, ymax=y+stderror, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(data=subset(combined_data_dry, label =="Moisture"), color="#FCD262", size=4, alpha=0.8) +
  geom_linerange(data=subset(combined_data_dry, label =="Moisture"), color = "#FCD262", linewidth = 1, alpha=0.4) + coord_flip() +
  scale_y_continuous(limits = c(-0.4, 0.4)) +
  ylab("Standardized coefficients") +
  #scale_y_continuous(limits = symmetric_limits) + 
  theme_bw() +
  theme(legend.position = "none", axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_line(), panel.border = element_blank(),axis.line = element_line(), plot.title = element_text(hjust = -0.66, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) + ggtitle("A. Dry sediments") 
rf_dry


#Reorder factor levels
combined_data_rip$variable <- factor(combined_data_rip$variable, levels=c("CO2", "CH4", "N2O"))
combined_data_rip$label <- factor(combined_data_rip$label, levels=c( "Flow state - Flowing",  "Moisture", "Temperature"))

# Plot the combined data using ggplot
rf_rip <- ggplot(combined_data_rip, aes(x=label, y=y, ymin=y-stderror, ymax=y+stderror, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(color="#527028", size=4, alpha=0.7) +
  geom_linerange(color = "#527028", size = 1, alpha=0.4) +
  coord_flip() + ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) + theme_bw() +
  theme(axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_line(), panel.border = element_blank(),axis.line = element_line(),  plot.title = element_text(hjust = -0.65, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) + ggtitle("B. Riparian soils") +  guides(shape = guide_legend(override.aes = aes(color = "grey"))) 
rf_rip


#### Combine riparian and dry sediment forest plots ####
tiff("rip_dry_combined", units="in", width=6, height=4, res=300)

rip_dry_combined <- ggarrange(rf_dry + rremove("xlab") , rf_rip,
                              ncol = 1, nrow = 2, align="v", common.legend = F,legend="bottom", heights= c(1, 1.5))
rip_dry_combined

dev.off()


#### Plot model coefficients for flowing conditions ####


#Load in the non-perennial and perennial coefficients
coef_np_flow <- read.csv("Model_coefficients_NPflow.csv")
coef_p_flow <- read.csv("Model_coefficients_Pflow.csv")


# Plot the combined data using ggplot
#A lot is going on, and the asterics would be too much, so only include the significant terms

#Reorder factor levels
coef_np_flow$variable <- factor(coef_np_flow$variable, levels=c("CO2", "CH4", "N2O"))
coef_np_flow$label <- factor(coef_np_flow$label, levels=c('Distance to source: OM stock', 'Distance to source', 'OM stock', 'Dissolved N', 'N-Nitrate', 'Velocity', 'Dissolved oxygen', 'Embeddedness', 'Time since rewetting', 'Water temperature'))

rf_np_flow1 <- ggplot(subset(coef_np_flow, abs(y) < 1), aes(x=label, y=y, ymin=y-stderror, ymax=y+stderror, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(color="#8BCEB9", size=4, alpha=0.7) +
  geom_linerange(color = "#8BCEB9", size = 1, alpha=0.4) +
  coord_flip() +
  ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = aes(color = "grey"))) +
  theme(axis.title.y=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_line(),  panel.border = element_blank(),axis.line = element_line(), 
        plot.title = element_text(hjust = -2, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  ggtitle("A. Non-perennial flowing") 
rf_np_flow1

rf_np_flow2 <- ggplot(subset(coef_np_flow, abs(y) > 1), aes(x=label, y=y, ymin=y-stderror, ymax=y+stderror, group=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(color="#8BCEB9", size=4, alpha=0.7, shape=15) +
  geom_linerange(color = "#8BCEB9", size = 1, alpha=0.4) +
  coord_flip() +
  ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = aes(color = "grey"))) +
  theme(axis.title.y=element_blank(),  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_line(), panel.border = element_blank(),axis.line = element_line(), plot.title = element_text(hjust = -2, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) 
rf_np_flow2


#combine the two np plot on different scales
np_combined <- ggarrange(rf_np_flow1 + rremove("xlab") , rf_np_flow2,
                         ncol = 1, nrow = 2, align="v", common.legend = T,legend="bottom", heights= c(1.5, .6))
np_combined



#reorder factor levels
coef_p_flow$variable <- factor(coef_p_flow$variable, levels=c("CO2", "CH4", "N2O"))
coef_p_flow$label <- factor(coef_p_flow$label, levels=c('Dissolved N', 'N-Nitrate', 'Velocity', 'OM stock', 'Sediment OM', '% non-perennial upstream'))

rf_p_flow <- ggplot(coef_p_flow, aes(x=label, y=y, ymin=y-stderror , ymax=y+stderror, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_point(color="#1C54A1", size=3.8, alpha=0.8) +
  geom_linerange(color = "#1C54A1", size = 1, alpha=0.4) +
  coord_flip() +
  ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) +
  theme_bw() +
  scale_shape_manual(values = c(16, 15))+
  guides(shape = guide_legend(override.aes = aes(color = "grey"))) +
  theme(axis.title.y=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.grid.major.y = element_line(),  panel.border = element_blank(),axis.line = element_line(), 
        plot.title = element_text(hjust = -1.5, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  ggtitle("B. Perennial flowing") 
rf_p_flow




#### Time series of each gas by campaign ####

# CO2 Site averages, boxplot by campaign and habitat #

CO2meansbox_camp<- ggplot() + 
  geom_boxplot(data = CO2means, position =  position_dodge2(width = 1.5, preserve = "single"),
               aes(x = as.factor(campaign), y = mean_CO2_g_m2_d, group = interaction(campaign, habitat),                    fill = habitat), width=0.85, lwd=.4, fatten = 1.5) + 
  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", ,axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) + 
  scale_y_continuous(limits = c(-1,15)) + 
  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1)) + xlab("Campaign") + scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "7"), labels=c("Mar", "May", "Jun", "Jul", "Sep", "Oct", "Nov")) +
  annotate("text", x = 2.897, y = 11.3, label=sprintf('\u2191')) +
  annotate("text", x = 2.897, y = 12.5, label="20")+
  annotate("text", x = 3.79, y = 13, label=sprintf('\u2191')) +
  annotate("text", x = 3.79, y = 14, label="19, 30")+
  #annotate("text", x = 3.75, y = 30, label="37, 62") +
  scale_fill_manual(values=c( "#D3DDDC", "#5598B7", "#FCD262" , "#527028")) 
CO2meansbox_camp 


## CH4 Site averages, boxplot by campaign and habitat ##
CH4means$habitat <- factor(CH4means$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

CH4meansbox_camp<- ggplot() + 
  geom_boxplot(data = CH4means, position =  position_dodge2(width = 1.5, preserve = "single"),
               aes(x = as.factor(campaign), y = mean_CH4_mg_m2_d, group = interaction(campaign, habitat),                    fill = habitat), width=0.85, lwd=.4, fatten = 1.5) + 
  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) + 
  scale_y_continuous(limits = c(-1.5,5)) + 
  geom_hline(yintercept=0, linetype='dashed', col = 'black') +
  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))+ xlab("Campaign") + scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "7"), labels=c("Mar", "May", "Jun", "Jul", "Sep", "Oct", "Nov")) +
  annotate("text", x = 1.895, y = 3.5, label=sprintf('\u2191')) +
  annotate("text", x = 1.895, y = 4.0, label="7") +
  annotate("text", x = 2.896, y = 3.5, label=sprintf('\u2191')) +
  annotate("text", x = 2.896, y = 4.0, label="11")+
  annotate("text", x = 3.795, y = 4.4, label=sprintf('\u2191')) +
  annotate("text", x = 3.795, y = 4.9, label="6, 7, 11")+
  scale_fill_manual(values=c( "#D3DDDC", "#5598B7", "#FCD262" , "#527028")) 
CH4meansbox_camp 

#### N2O Site averages, boxplot by campaign and habitat ####

N2Omeans$habitat <- factor(N2Omeans$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))


N2Omeansbox_camp<- ggplot() + 
  geom_boxplot(data = N2Omeans, position =  position_dodge2(width = 1.5, preserve = "single"),
               aes(x = as.factor(campaign), y = mean_mg_N2ON_m2_d, group = interaction(campaign, habitat),                    fill = habitat), width=0.85, lwd=.4, fatten = 1.5) + 
  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), axis.title.x = element_blank(), legend.position = c(.95, .95), legend.justification = c("right", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black")) + 
  scale_y_continuous(limits = c(-0.05,1.7)) + 
  ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) + xlab("Campaign") + scale_x_discrete(breaks=c("1","2","3", "4", "5", "6", "7"), labels=c("Mar", "May", "Jun", "Jul", "Sep", "Oct", "Nov")) +
  scale_fill_manual(values=c( "#D3DDDC", "#5598B7", "#FCD262" , "#527028")) +
  annotate("text", x = 1.79, y = 1.52, label="2.6")+
  annotate("text", x = 1.79, y = 1.42, label=sprintf('\u2191')) 
N2Omeansbox_camp 






##____________________________________________________________________________##

#### Combine plots ####

#### Figure 2. Combine riparian and dry sediment forest plots ####
tiff("rip_dry_combined", units="in", width=6, height=6, res=300)

rip_dry_combined <- ggarrange(rf_dry + rremove("xlab") , rf_rip,
                              ncol = 1, nrow = 2, align="v", common.legend = T,legend="bottom")
rip_dry_combined


dev.off()


#### Figure 3. Flowing GHGs by drying regime####
tiff("FlowingbyDryingRegime", units="in", width=3, height=8, res=300)

FlowingbyDryingRegime <- ggarrange(CO2_drying_regime, CH4_drying_regime, N2O_drying_regime,
                                   labels = c("A", "B", "C"),
                                   ncol = 1, nrow = 3)
FlowingbyDryingRegime

dev.off()


#### Figure 4. Non-perennial and perennial flowing covariates ####
tiff("flow_combined", units="in", width=6, height=6, res=300)

flow_combined <- ggarrange(rf_np_flow1 + rremove("xlab") , rf_np_flow2 + rremove("xlab") , rf_p_flow,
                           ncol = 1, nrow = 3, align="v", common.legend = T,legend="bottom", heights= c(2.2, .75, 2))
flow_combined

dev.off()


#### Figure S2. Time series by each gas ####
tiff("GHG_timeseries", units="in", width=8, height=11, res=300)
GHG_timeseries <- ggarrange(CO2meansbox_camp, CH4meansbox_camp, N2Omeansbox_camp,
                            labels = c("A", "B", "C"),
                            ncol = 1, nrow = 3, align="v", common.legend = T,legend="bottom")
GHG_timeseries
dev.off()

#### Figure S3. GHGs by habitat using Site means ####
tiff("MeanGHGbyhabitat", units="in", width=3, height=8.5, res=300)

MeanGHGbyhabitat <- ggarrange(CO2_habitat_means, CH4_habitat_means, N2O_habitat_means,
                              labels = c("A", "B", "C"),
                              ncol = 1, nrow = 3, align="v")
MeanGHGbyhabitat

dev.off()

