#### Data and statistical analysis and visalization for GHG data collected along the Albarine River network, France in 2021 ####

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


##________________________________________________________________________##

### Sediment Organic Matter LMM ####

#Results are in Table S4

# Run dredge model selection for sediment OM testing effect of drying, time since rewetting in interaction with distance to source, separately for dry and flowing. 

#Include top correlations from CO2 correlation table: dry (sand, benthic stock) and flowing (canopy cover benth, tsince manual, m to source, percent IR up). Too many variables, so I think we might need to remove sand... also canopy cover covaries with m to source.

#subset OM data
sed <- dat %>% drop_na(Sed_OM) %>%
  #Average by site, campaign and habitat
  dplyr::group_by(Site, campaign, habitat, drying_regime) %>% 
  dplyr::summarise(
    Sed_OM = mean(Sed_OM),
    Soil_OM = mean(Soil_OM),
    Freq_Flow = mean(Freq_Flow),
    m_to_source = mean(m_to_source),
    Tsince_manual = mean(Tsince_manual),
    percent.IR.up = mean(percent.IR.up),
    Canopy_Cover_Benth = mean(Canopy_Cover_Benth),
    sand = mean(sand) )

str(sed)  #135 obs

sed <- as.data.frame(sed)

# Add season
sed <- sed %>%
  mutate(Season = 
           case_when(
             campaign==1 ~ "Spring",
             campaign==2 ~ "Spring", 
             campaign==3 ~ "Summer", 
             campaign==4 ~ "Summer", 
             campaign==5 ~ "Autumn", 
             campaign==6 ~ "Autumn",
             campaign==7 ~ "Autumn"  ))

#exclude the singular pool value and remove the riparian value
sed <- subset(sed, habitat!="Pool" & habitat !="Riparian")  #67 obs

#Tsince manual has NAs for perennial sites, replace with a very high number, like 1000
sed$Tsince_manual <- ifelse(sed$drying_regime == "perennial", 1000, sed$Tsince_manual)

#remove the outlier from Cooks D
sed<- sed %>%
  filter(!Sed_OM ==16.61)

sed$log_m_to_source <- log(sed$m_to_source) #log transform before scaling

sed_scaled <- sed  %>% dplyr::mutate_at(c("Sed_OM" , "Freq_Flow", "log_m_to_source", "Tsince_manual",  "percent.IR.up"), ~(scale(., center=FALSE) %>% as.vector))

sed_NA_scaled <- na.omit(sed_scaled) #remove NAs


sed_lmer <- lmer(Sed_OM~ Freq_Flow*log_m_to_source+ 
                   Tsince_manual*log_m_to_source+ 
                   percent.IR.up*habitat + 
                   habitat*log_m_to_source+
                   Season +
                   (1|campaign ), data=sed_NA_scaled)
summary(sed_lmer)


plot(sed_lmer) #much better after outlier removed
qqnorm(resid(sed_lmer))
vif(sed_lmer) #without interactions shows freq flow and time since might be correlated, Tsince removed (higher VIF)

#looks to have an outlier in the data, check Cook's
cooksD <- cooks.distance(sed_lmer)
n <- nrow(sed)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue") # add cutoff line

influential_obs <- as.numeric(names(cooksD)[(cooksD > (25/n))])
influential_obs 


#Run dredge
options(na.action = "na.fail") # Required for dredge to run

sedOM_dredge<- dredge(sed_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(sedOM_dredge)  #how many models were run 184
head(sedOM_dredge)


top.sedOM_dredge  <- get.models(sedOM_dredge , subset=delta <= 2) #The top models, with delta AIC <2

SedOM_mavg <- model.avg(top.sedOM_dredge)
summary(SedOM_mavg) #Nothing significant

#### Benthic OM stock LMM ####

# Run dredge model for benthic OM  stock testing effect of drying, time since rewetting in interaction with distance to source, separately for dry and flowing. 

#Include top correlations from CO2 correlation table: dry (sand, benthic stock) and flowing (canopy cover benth, tsince manual, m to source, percent IR up)

#subset OM data
OMstock <- dat %>% drop_na(Stock_Benth_g.m2) %>%
  #Average by site, campaign and habitat
  dplyr::group_by(Site, campaign, habitat, drying_regime) %>% 
  dplyr::summarise(
    Stock_Benth_g.m2 = mean(Stock_Benth_g.m2),
    Freq_Flow = mean(Freq_Flow),
    freqflow7 = mean(freqflow7),
    freq_flow_2021 = mean(freq_flow_2021),
    m_to_source = mean(m_to_source),
    Tsince_manual = mean(Tsince_manual),
    length.IR.up = mean(length.IR.up),
    percent.IR.up = mean(percent.IR.up),
    Canopy_Cover_Benth = mean(Canopy_Cover_Benth),
    Mean_Velocity..m.s.= mean(Mean_Velocity..m.s.) )

str(OMstock)  #155 obs

OMstock <- as.data.frame(OMstock)

# Add season
OMstock <- OMstock %>%
  mutate(Season = 
           case_when(
             campaign==1 ~ "Spring",
             campaign==2 ~ "Spring", 
             campaign==3 ~ "Summer", 
             campaign==4 ~ "Summer", 
             campaign==5 ~ "Autumn", 
             campaign==6 ~ "Autumn",
             campaign==7 ~ "Autumn"  ))

#remove the outlier from residual plot
OMstock<- OMstock %>%
  filter(!Stock_Benth_g.m2 >= 70)

#exclude the singular pool value and remove the riparian value
OMstock <- subset(OMstock, habitat!="Pool" & habitat !="Riparian")  #77 obs

#Tsince manual has NAs for perennial sites, replace with a very high number, like 1000
OMstock$Tsince_manual <- ifelse(OMstock$drying_regime == "perennial", 1000, OMstock$Tsince_manual)

OMstock$log_m_to_source <- log(OMstock$m_to_source) #log transform before scaling

OMstock_scaled <- OMstock  %>% dplyr::mutate_at(c( "Freq_Flow", "log_m_to_source","Tsince_manual" , 
                                                   "percent.IR.up",  "Mean_Velocity..m.s." ), ~(scale(., center=FALSE) %>% as.vector))

OMstock_NA_scaled <- na.omit(OMstock_scaled) #Using complete cases reduces sample size from 76 to 56


OMstock_lmer <- lmer(log(Stock_Benth_g.m2)~ 
                       Freq_Flow*log(m_to_source)+ 
                       Tsince_manual*log(m_to_source)+ 
                       percent.IR.up*habitat + 
                       habitat*log(m_to_source)+ 
                       Mean_Velocity..m.s. +
                       Season*log(m_to_source) +
                       (1|campaign ), data=OMstock_NA_scaled)
summary(OMstock_lmer)

plot(OMstock_lmer, id=0.001) #outlier
qqnorm(resid(OMstock_lmer)) #improve after log transformation
vif(OMstock_lmer) #high VIF for tsince_manual, length IR up

#run dredge
options(na.action = "na.fail") # Required for dredge to run

OMstock_dredge<- dredge(OMstock_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(OMstock_dredge)  #how many models were run 368
head(OMstock_dredge)

top.OMstock_dredge  <- get.models(OMstock_dredge , subset=delta <= 2) #The top models

OMstock_mavg <- model.avg(top.OMstock_dredge)
summary(OMstock_mavg)  #Season spring and summer are significant

x <- subset(OMstock_dredge , delta <= 2)
mean(x$`R^2`)


##__________________________________________________________________________##

#### Stepwise Model selection for CO2 ####

#### For dry CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_dry <-  subset(dat_CO2,  habitat == "Dry")  #select dry
#Select relevant columns
dat_CO2_dry1 <- data.table(dat_CO2_num_dry[, c("campaign", "Site" , "CO2_g_m2_d", "Stock_Benth_g.m2", "CN_ratio", "percent_N" , "m_to_source" , "Freq_Flow", "sed_VWC" , "sand",  "Tsince_manual", "Sed_OM" ,"sed_temp",  "sed_pH")])
dat_CO2_dry1NA <- na.omit(dat_CO2_dry1) #remove NAs

#removing the NAs, we go from 100 to 42 observations... The large culprits are Cn ratio and percent N (only for 2/5 campaigns) #can try to fill gaps in sed temp and sed vwc

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

dat_CO2_dry1NAout_scaled <- dat_CO2_dry1NAout  %>% mutate_at(c("Stock_Benth_g.m2", "CN_ratio", "percent_N" , "log_m_to_source" , "Freq_Flow", "sed_VWC" , "sand",  "Tsince_manual", "Sed_OM" ,"sed_temp",  "sed_pH"), ~(scale(., center=FALSE) %>% as.vector))


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

vif(CO2_dry_lmer) #When we exclude the interactions, only canopy cover is over 10 (10.65). Consider running dredge with and without


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

#Run the top model variables

CO2_dry_lmer2 <- lmer(log(CO2_g_m2_d+1.4031147)~  
                        Freq_Flow+
                        Sed_OM*sed_VWC + 
                        sed_pH +
                        (1 + campaign | Site), data=dat_CO2_dry1NAout_scaled)
summary(CO2_dry_lmer2)


#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
CO2dry_coef <- coefTable(mA, full=FALSE) 
CO2dry_coef <- as.data.frame(CO2dry_coef)  
CO2dry_coef <- tibble::rownames_to_column(CO2dry_coef, "Term") #need to keep row names as column

names(CO2dry_coef)[names(CO2dry_coef) == "Std. Error"] <- "StdError"

CO2dry_coef <- CO2dry_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError
  )

#rename variables
CO2dry_coef$Term[CO2dry_coef$Term == "sed_pH"] <- "Sediment pH"
CO2dry_coef$Term[CO2dry_coef$Term == "sed_VWC"] <- "Sediment moisture"
CO2dry_coef$Term[CO2dry_coef$Term == "sed_temp"] <- "Sediment temperature"
CO2dry_coef$Term[CO2dry_coef$Term == "sed_temp:sed_VWC"] <- "Sediment temperature: Sediment moisture"


#CO2_dry_lmer2_df$term = reorder(as.factor(CO2_dry_lmer2_df$term), -CO2_dry_lmer2_df$p.value, FUN = mean) #reorder factor levels by p
#add p values from model average output 
#sed_temp           0.07629 . 
#sed_temp:sed_VWC   0.00321 **

#plot

my.labels <- c( "Sed_moist",
                "Sed_pH",
                "Sed_temp", 
                "Sed_temp*Sed_moist")

CO2_dry_standard <- ggplot(subset(CO2dry_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  scale_x_discrete(labels= my.labels) +
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 4.2, y = 1.5058001 , label = "**", size=5)  
CO2_dry_standard


#### For perennial-flowing CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_flow_p <-  subset(dat_CO2,  habitat == "Flowing" & drying_regime=="perennial") 

#Select relevant variables
CO2_flow_p <- data.table(dat_CO2_num_flow_p[, c("campaign", "Site" , "CO2_g_m2_d", "m_to_source" , "Stock_Benth_g.m2" , "Sed_OM", "Mean_Velocity..m.s." ,"percent.IR.up" , "Temp_water_C" )])

CO2_flow_pNA <- na.omit(CO2_flow_p) #remove NAs (357 to 167, largest culprits are In sed/stock, velocity, percentIR)

min(CO2_flow_pNA$CO2_g_m2_d) #-0.01754717 #Min CO2 value in order to add to log transformation

#Scale variables in order to compare standardized coefficients
CO2_flow_pNA$log_m_to_source <- log(CO2_flow_pNA$m_to_source) #log transform
CO2_flow_pNA_scaled <- CO2_flow_pNA  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2" , "Sed_OM", "Mean_Velocity..m.s." ,"percent.IR.up" , "Temp_water_C" ), ~(scale(., center=FALSE) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CO2_flow_p_lmer <- lmer(log(CO2_g_m2_d+1.01754717)~ 
                          log_m_to_source*Stock_Benth_g.m2 + 
                          log_m_to_source*Sed_OM + 
                          Stock_Benth_g.m2*Temp_water_C +
                          Sed_OM*Temp_water_C +
                          Mean_Velocity..m.s. + 
                          percent.IR.up +  
                          (1 + campaign | Site), data=CO2_flow_pNA_scaled) 


summary(CO2_flow_p_lmer)

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


##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
CO2flowp_coef <- coefTable(mA_CO2_flow_p, full=FALSE) 
CO2flowp_coef <- as.data.frame(CO2flowp_coef)  
CO2flowp_coef <- tibble::rownames_to_column(CO2flowp_coef, "Term") #need to keep row names as column

names(CO2flowp_coef)[names(CO2flowp_coef) == "Std. Error"] <- "StdError"

CO2flowp_coef <- CO2flowp_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError)

#rename variables
CO2flowp_coef$Term[CO2flowp_coef$Term == "percent.IR.up"] <- "% IRES upstream"
CO2flowp_coef$Term[CO2flowp_coef$Term == "Sed_OM"] <- "Sediment OM"
CO2flowp_coef$Term[CO2flowp_coef$Term == "log_m_to_source"] <- "Distance to source"

#get p values form model average output 
#percent.IR.up     3.06e-05 ***
#  Sed_OM          0.0333 * 

#plot
CO2_flow_p_standard <- ggplot(subset(CO2flowp_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 3.2, y = 0.20672, label = "*", size=5) +
  annotate("text", x = 1.2, y = 0.60366, label = "***", size=5) 
CO2_flow_p_standard


#### For intermittent-flowing CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_flow_i <-  subset(dat_CO2,  habitat == "Flowing" & (drying_regime=="intermittent_long" | drying_regime=="intermittent_moderate") )

#Select relevant variables
CO2_flow_i <- data.table(dat_CO2_num_flow_i[, c("campaign", "Site", "CO2_g_m2_d", "m_to_source", "Tsince_manual",  "Freq_Flow", "Stock_Benth_g.m2", "Sed_OM", "Embeddedness....",  "CH4_mg_m2_d", "Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C")])

CO2_flow_iNA <- na.omit(CO2_flow_i) #remove NAs (214 to 93, largest culprits are In sed/stock, velocity, tree div, algae, etc.)

min(CO2_flow_iNA$CO2_g_m2_d) #-0.003192917 #Min CO2 value in order to add to log transformation

#Scale variables in order to compare standardized coefficients
CO2_flow_iNA$log_m_to_source <- log(CO2_flow_iNA$m_to_source) #log transform
CO2_flow_iNA_scaled <- CO2_flow_iNA  %>% mutate_at(c("log_m_to_source",  "Tsince_manual",  "Freq_Flow", "Stock_Benth_g.m2", "Sed_OM", "Embeddedness....", "CH4_mg_m2_d", "Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C"), ~(scale(., center=FALSE) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CO2_flow_i_lmer <- lmer(log(CO2_g_m2_d+1.003192917)~ 
                          Freq_Flow*log(m_to_source)+
                          Tsince_manual*log(m_to_source) + 
                          Stock_Benth_g.m2 *log(m_to_source) + 
                          Sed_OM*log(m_to_source) + 
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


##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
CO2flowi_coef <- coefTable(mA_CO2_flow_i, full=FALSE) 
CO2flowi_coef <- as.data.frame(CO2flowi_coef)  
CO2flowi_coef <- tibble::rownames_to_column(CO2flowi_coef, "Term") #need to keep row names as column

names(CO2flowi_coef)[names(CO2flowi_coef) == "Std. Error"] <- "StdError"

CO2flowi_coef <- CO2flowi_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError)

#rename variables
CO2flowi_coef$Term[CO2flowi_coef$Term == "Embeddedness...."] <- "Embeddedness"
CO2flowi_coef$Term[CO2flowi_coef$Term == "Freq_Flow"] <- "Flowing frequency"
CO2flowi_coef$Term[CO2flowi_coef$Term == "Temp_water_C"] <- "Water temperature"

#get p values form model average output 
#Embeddedness....     <2e-16 ***
#  Freq_Flow          0.0199 *  
#  Temp_water_C      <2e-16 ***

#plot
CO2_flow_i_standard <- ggplot(subset(CO2flowi_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 3.2, y = 0.6392356 , label = "***", size=5) +
  annotate("text", x = 2.2, y = 0.6157967 , label = "*", size=5) +
  annotate("text", x = 1.2, y = 0.4355796 , label = "***", size=5) 
CO2_flow_i_standard



#### For riparian CO2 emissions ####

#Need to use data with no NAs
dat_CO2_num_rip <-  subset(dat_CO2,  habitat == "Riparian") 

#Select relevant variables
CO2_rip <- data.table(dat_CO2_num_rip[, c("campaign", "Site", "CO2_g_m2_d", "soil_temp", "soil_VWC", "Soil_OM","totalbiomass_g",   "Freq_Flow", "m_to_source", "percent_N", "flow_state")])


CO2_ripNA <- na.omit(CO2_rip) #remove NAs (642 to 311 obs, soil OM, percent N, stock etc.)

min(CO2_ripNA$CO2_g_m2_d) #-9.185739 #Min CO2 value in order to add to log transformation # after removing outlier -1.295294

#exclude the outlier of GR01, C1 (9.185738684) #2021-03-17_GR01_R2_Picarro

CO2_ripNAout<- CO2_ripNA %>%
  filter(!CO2_g_m2_d <= -9)

#try also after scaling variables to get standardized coefficients
#log transform m to source
CO2_ripNAout$log_m_to_source <- log(CO2_ripNAout$m_to_source)

CO2_ripNAout_scaled <- CO2_ripNAout  %>% mutate_at(c("log_m_to_source", "soil_temp", "soil_VWC", "Soil_OM", "totalbiomass_g","Freq_Flow", "m_to_source", "percent_N"), ~(scale(., center=FALSE) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CO2_rip_lmer <- lmer(log(CO2_g_m2_d+10.185739)~ 
                       Freq_Flow*log(m_to_source)+ 
                       flow_state *log(m_to_source)+  
                       totalbiomass_g*log(m_to_source) +
                       Soil_OM*log(m_to_source) +
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+  
                       (1 + campaign | Site), data=CO2_ripNAout_scaled)


summary(CO2_rip_lmer)

anova(CO2_rip_lmer)
plot(CO2_rip_lmer, id=0.05) #Removed outlier based on Cooks distance
qqnorm(resid(CO2_rip_lmer))
vif(CO2_rip_lmer)

#looks to have an outlier in the data, check Cook's
cooksD <- cooks.distance(CO2_rip_lmer)

n <- nrow(CO2_ripNA)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 20/n, lty = 2, col = "steelblue") # add cutoff line

influential_obs <- as.numeric(names(cooksD)[(cooksD > (100/n))])
influential_obs #17  18 184 # which correspond to CO2 values of -9.185738684, 3.428130457 and 17.0197525
#I think it is the -9 that is the outlier we see int he residual plot, so try running the model without out it

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
summary(CO2rip_mavg) #only one model supplied

#Run lmer with variables from the top model: soil temperature and soil moisture

CO2_rip_lmer3 <- lmer(log(CO2_g_m2_d+10.185739)~ 
                        soil_temp + soil_VWC + 
                        (1 + campaign | Site), data=CO2_ripNAout_scaled)
summary(CO2_rip_lmer3)

r.squaredGLMM(CO2_rip_lmer3) #marginal R2

#plot the standardized coefficicents from the top model in a forest plot
CO2_rip_lmer3_df <- broom.mixed::tidy(x = CO2_rip_lmer3)  #make data frame of model outputs

CO2_rip_lmer3_df <- CO2_rip_lmer3_df %>%
  dplyr::mutate(
    .data = .,
    conf.low = estimate - std.error,
    conf.high = estimate + std.error )

#rename variables
CO2_rip_lmer3_df$term[CO2_rip_lmer3_df$term == "soil_temp"] <- "Soil_temp"
CO2_rip_lmer3_df$term[CO2_rip_lmer3_df$term == "soil_VWC"] <- "Soil_moist"


#CO2_rip_lmer3_df$term = reorder(as.factor(CO2_rip_lmer3_df$term), -CO2_rip_lmer3_df$p.value, FUN = mean) #reorder factor levels by p

#plot
CO2_rip_standard <- ggplot(subset(CO2_rip_lmer3_df, effect=="fixed" & term!="(Intercept)"), aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 2.2, y = 0.24791515, label = "***", size=5)  +
  annotate("text", x = 1.2, y = -0.09032863, label = "***", size=5) 
CO2_rip_standard


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
dat_CH4_num_dry <-  subset(dat_CH4,  habitat == "Dry")  #select dry
dat_CH4_dry1 <- data.table(dat_CH4_num_dry[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2",  "percent_N" , "m_to_source" , "Freq_Flow","sed_VWC" , "sed_pH", "Tsince_manual",   "Sed_OM" ,"sed_temp",  "sed_pH", "Canopy_Cover_Benth" )]) #Select relevant columns
dat_CH4_dry1NA <- na.omit(dat_CH4_dry1) #remove NAs


#removing the NAs, we go from 100 to 42 observations... The large culprits are Cn ratio and percent N (only for 2/5 campaigns) #can try to fill gaps in sed temp and sed vwc

min(dat_CH4_dry1NA$CH4_mg_m2_d) #-0.5059006

dat_CH4_dry1NA$log_m_to_source <- log(dat_CH4_dry1NA$m_to_source)

dat_CH4_dry1NA_scaled <- dat_CH4_dry1NA  %>% mutate_at(c("Stock_Benth_g.m2",  "percent_N" , "log_m_to_source" , "Freq_Flow","sed_VWC" , "sed_pH", "Tsince_manual",   "Sed_OM" ,"sed_temp",  "sed_pH", "Canopy_Cover_Benth"), ~(scale(., center=FALSE) %>% as.vector))


#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
#log(CH4_mg_m2_d+1.5059006) residuals not improve by log transforming
CH4_dry_lmer <- lmer(CH4_mg_m2_d ~ 
                       Freq_Flow*log_m_to_source + 
                       Tsince_manual*log_m_to_source +  
                       Sed_OM*log_m_to_source + 
                       Stock_Benth_g.m2*log_m_to_source + 
                       sed_VWC*Sed_OM +  
                       Sed_OM*sed_temp + 
                       sed_VWC*sed_temp +
                       Canopy_Cover_Benth + 
                       sed_pH +
                       (1 + campaign | Site), data=dat_CH4_dry1NA_scaled) 

summary(CH4_dry_lmer) #get a lot of warnings when I use all of the m_to_source interactions
anova(CH4_dry_lmer)
plot(CH4_dry_lmer, id=0.05) #residuals look fine
qqnorm(resid(CH4_dry_lmer))



#reran using more data
CH4_dry_lmer2 <- lmer(log(CH4_mg_m2_d+5) ~  m_to_source*Tsince_manual +  m_to_source*Freq_Flow +   sed_temp *sed_VWC + (1 + campaign | Site), data=dat_CH4_dry1) 
summary(CH4_dry_lmer2)

#try using Site means
CH4_dry_lmer3 <- lmer(mean_CH4_mg_m2_d ~ Freq_Flow*log(m_to_source) + 
                        Tsince_manual*log(m_to_source) +  
                        sed_VWC*sed_temp  + (1 | campaign), data=subset(CH4means, habitat=="Dry" ) )
summary(CH4_dry_lmer3)



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

top.CH4_dry_dredge  <- get.models(CH4_dry_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(CH4_dry_dredge)  # Variable importance

CH4dry_mavg <- model.avg(top.CH4_dry_dredge)
summary(CH4dry_mavg)

## Make a random forest plot of the model averaged estimates (standardized coefficients)

CH4dry_mavg <- model.avg(top.CH4_dry_dredge)
mA <- summary(CH4dry_mavg)

# Extract the coefficients and standard errors 
CH4dry_coef <- coefTable(mA, full=FALSE) 
CH4dry_coef <- as.data.frame(CH4dry_coef)  
CH4dry_coef <- tibble::rownames_to_column(CH4dry_coef, "Term") #need to keep row names as column

names(CH4dry_coef)[names(CH4dry_coef) == "Std. Error"] <- "StdError"

CH4dry_coef <- CH4dry_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError
  )

#rename variables
CH4dry_coef$Term[CH4dry_coef$Term == "sed_pH"] <- "Sed_pH"

#CO2_dry_lmer2_df$term = reorder(as.factor(CO2_dry_lmer2_df$term), -CO2_dry_lmer2_df$p.value, FUN = mean) #reorder factor levels by p
#add p values from model average output 
#sed_temp           0.07629 . 
#sed_temp:sed_VWC   0.00321 **

#plot
CH4_dry_standard <- ggplot(subset(CH4dry_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) 
#annotate("text", x = 4.2, y = 1.5058001 , label = "**", size=5)  
CH4_dry_standard


#### For flowing perennial CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_flow <-  subset(dat_CH4,  habitat == "Flowing")  #select flowing
dat_CH4_num_flow_p <-  subset(dat_CH4_num_flow,  drying_regime == "perennial") #select perennial
dat_CH4_num_flow_p1 <- data.table(dat_CH4_num_flow_p[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" , "Freq_Flow",  "Temp_water_C", "Air_temp", "DO_mg_L", "Mean_Velocity..m.s.", "percent.IR.up" )]) #Select relevant columns
dat_CH4_num_flow_p1NA <- na.omit(dat_CH4_num_flow_p1) #remove NAs #Go from 353 to 166

#remove 2 outliers identified by Cooks D

dat_CH4_num_flow_p1NA_out<-dat_CH4_num_flow_p1NA[-c(164, 165),]

min(dat_CH4_num_flow_p1NA_out$CH4_mg_m2_d) #-0.7679176

dat_CH4_num_flow_p1NA_out$log_m_to_source <- log(dat_CH4_num_flow_p1NA_out$m_to_source)

dat_CH4_num_flow_p1NA_out_scaled <- dat_CH4_num_flow_p1NA_out  %>% mutate_at(c("Stock_Benth_g.m2", "Sed_OM", "log_m_to_source" , "Freq_Flow",  "Temp_water_C", "Air_temp", "DO_mg_L", "Mean_Velocity..m.s.", "percent.IR.up" ), ~(scale(., center=FALSE) %>% as.vector))

#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
CH4_flow_p_lmer <- lmer(log(CH4_mg_m2_d+1.7679176) ~  
                          log_m_to_source*Stock_Benth_g.m2 + 
                          log_m_to_source*Sed_OM + 
                          percent.IR.up +
                          Temp_water_C*Stock_Benth_g.m2 + 
                          Temp_water_C* Sed_OM  + 
                          DO_mg_L*Stock_Benth_g.m2 +  
                          DO_mg_L*Sed_OM + 
                          Mean_Velocity..m.s. +
                          (1 + campaign | Site), data=dat_CH4_num_flow_p1NA_out_scaled)

summary(CH4_flow_p_lmer) 
anova(CH4_flow_p_lmer)
plot(CH4_flow_p_lmer) #improve after log transformation
qqnorm(resid(CH4_flow_p_lmer))

#looks to have outliers in the data, check Cook's
cooksD <- cooks.distance(CH4_flow_p_lmer)

n <- nrow(dat_CH4_num_flow_p1NA)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 50/n, lty = 2, col = "steelblue") # add cutoff line

influential_obs <- as.numeric(names(cooksD)[(cooksD > (50/n))])
influential_obs # 164 165


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


##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
CH4flowp_coef <- coefTable(mA_CH4_flow_p, full=FALSE) 
CH4flowp_coef <- as.data.frame(CH4flowp_coef)  
CH4flowp_coef <- tibble::rownames_to_column(CH4flowp_coef, "Term") #need to keep row names as column

names(CH4flowp_coef)[names(CH4flowp_coef) == "Std. Error"] <- "StdError"

CH4flowp_coef <- CH4flowp_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError)

#rename variables
CH4flowp_coef$Term[CH4flowp_coef$Term == "Sed_OM"] <- "Sediment OM"
CH4flowp_coef$Term[CH4flowp_coef$Term == "DO_mg_L:Sed_OM"] <- "Dissolved oxygen:Sediment OM"
CH4flowp_coef$Term[CH4flowp_coef$Term == "DO_mg_L"] <- "Dissolved oxygen"

#get p values form model average output 
#Sed_OM            0.000132 ***
# DO_mg_L:Sed_OM   0.000254 ***

#plot
CH4_flow_p_standard <- ggplot(subset(CH4flowp_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 3.2, y = 1.6777     , label = "***", size=5) +
  annotate("text", x = 2.2, y = -1.6523, label = "***", size=5) 
CH4_flow_p_standard


#### For flowing intermittent CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_flow <-  subset(dat_CH4,  habitat == "Flowing")  #select flowing
dat_CH4_num_flow_i <-  subset(dat_CH4_num_flow,  drying_regime == "intermittent_long" | drying_regime == "intermittent_moderate" ) #select intermittnet
dat_CH4_num_flow_i1 <- data.table(dat_CH4_num_flow_i[, c("campaign", "Site" , "CH4_mg_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" , "Tsince_manual",  "Freq_Flow", "Air_temp", "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL", "CO2_mg_m2_h", "DOC_ugL", "Canopy_Cover_Benth", "percent.IR.up")]) #Select relevant columns
dat_CH4_num_flow_i1NA <- na.omit(dat_CH4_num_flow_i1) #remove NAs #Go from 214 to 93 obs


min(dat_CH4_num_flow_i1NA$CH4_mg_m2_d) #-0.2719916

dat_CH4_num_flow_i1NA$log_m_to_source <- log(dat_CH4_num_flow_i1NA$m_to_source)

dat_CH4_num_flow_i1NA_scaled <- dat_CH4_num_flow_i1NA  %>% mutate_at(c("Stock_Benth_g.m2", "Sed_OM", "log_m_to_source" , "Freq_Flow", "Air_temp", "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL", "CO2_mg_m2_h", "DOC_ugL", "Canopy_Cover_Benth", "percent.IR.up"), ~(scale(., center=FALSE) %>% as.vector))


#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
CH4_flow_i_lmer <- lmer(log(CH4_mg_m2_d+2.086478) ~ 
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
                          (1 + campaign | Site), data=dat_CH4_num_flow_i1NA_scaled)

summary(CH4_flow_i_lmer)  # A lot of model warnings
anova(CH4_flow_i_lmer)
plot(CH4_flow_i_lmer) 
qqnorm(resid(CH4_flow_i_lmer))


options(na.action = "na.fail") # Required for dredge to run

CH4_flow_i_dredge<- dredge(CH4_flow_i_lmer, trace = 2, extra = list(
  "R^2", "*" = function(x) {
    s <- summary(x)
    c(Rsq = s$r.squared, adjRsq = s$adj.r.squared,
      F = s$fstatistic[[1]])
  }))

options(na.action = "na.omit") # set back to default

nrow(CH4_flow_i_dredge)  #how many models were run: 5710
head(CH4_flow_i_dredge)

top.CH4_flow_i_dredge  <- get.models(CH4_flow_i_dredge , subset=delta <= 2) #The top models

sw(top.CH4_flow_i_dredge)  # Variable importance

CH4flow_i_mavg <- model.avg(top.CH4_flow_i_dredge)
mA_CH4_flow_i <- summary(CH4flow_i_mavg)

##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficients from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
CH4flowi_coef <- coefTable(mA_CH4_flow_i, full=FALSE) 
CH4flowi_coef <- as.data.frame(CH4flowi_coef)  
CH4flowi_coef <- tibble::rownames_to_column(CH4flowi_coef, "Term") #need to keep row names as column

names(CH4flowi_coef)[names(CH4flowi_coef) == "Std. Error"] <- "StdError"

CH4flowi_coef <- CH4flowi_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError)

#rename variables
CH4flowi_coef$Term[CH4flowi_coef$Term == "DO_mg_L"] <- "Dissolved oxygen"
CH4flowi_coef$Term[CH4flowi_coef$Term == "Temp_water_C"] <- "Water temperature"
CH4flowi_coef$Term[CH4flowi_coef$Term == "log_m_to_source"] <- "Distance to source"
CH4flowi_coef$Term[CH4flowi_coef$Term == "Stock_Benth_g.m2"] <- "OM stock"
CH4flowi_coef$Term[CH4flowi_coef$Term == "log_m_to_source:Stock_Benth_g.m2"] <- "Distance to source:OM stock"

#get p values from ave model output
#Temp_water_C                        2.1e-06 ***
# log_m_to_source                     < 2e-16 ***
#Stock_Benth_g.m2                    1.0e-07 ***
#log_m_to_source:Stock_Benth_g.m2    < 2e-16 ***


#plot
CH4_flow_i_standard <- ggplot(subset(CH4flowi_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 5.2, y = 0.6309 , label = "***", size=5) +
  annotate("text", x = 4.2, y = 7.5065 , label = "***", size=5) +
  annotate("text", x = 3.2, y = -9.6906 , label = "***", size=5) +
  annotate("text", x = 2.2, y = 4.2409 , label = "*", size=5) 
CH4_flow_i_standard



#### For riparian CH4 emissions ####

#Need to use data with no NAs
dat_CH4_num_rip <-  subset(dat_CH4,  habitat == "Riparian") 

#Select relevant variables
CH4_rip <- data.table(dat_CH4_num_rip[, c("campaign", "Site", "CH4_mg_m2_d", "soil_temp", "soil_VWC", "Soil_OM", "Stock_Rip_g.m2", "flow_state" ,   "totalbiomass_g",  "Freq_Flow", "m_to_source", "percent_N")])
CH4_ripNA <- na.omit(CH4_rip) #remove NAs (643 to 312 obs, soil OM, percent N, stock etc.)


#remove  outliers identified in residuals plot
CH4_ripNA_out <- CH4_ripNA %>%
  filter(!CH4_mg_m2_d <= -5)

min(CH4_ripNA_out$CH4_mg_m2_d) #-1.842777 #Min CH4 value in order to add to log transformation
min(CH4_ripNA$CH4_mg_m2_d)     #-49.74601

#try also after scaling variables to get standardized coefficients
#log transform m to source
CH4_ripNA_out$log_m_to_source <- log(CH4_ripNA_out$m_to_source)

CH4_ripNA_scaled <- CH4_ripNA_out  %>% mutate_at(c("log_m_to_source",  "soil_temp", "soil_VWC", "Soil_OM", "Stock_Rip_g.m2", "totalbiomass_g",  "Freq_Flow",  "percent_N"), ~(scale(., center=FALSE) %>% as.vector))

# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
CH4_rip_lmer <- lmer(log(CH4_mg_m2_d+2.842777)~ 
                       Freq_Flow*log_m_to_source+ 
                       flow_state *log_m_to_source+  
                       totalbiomass_g*log_m_to_source +
                       Soil_OM*log_m_to_source+ 
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+
                       (1 + campaign | Site), data=CH4_ripNA_scaled)



#looks to have outliers in the data, check Cook's--not working in removing them. try to label residuals: the two BR02 values

summary(CH4_rip_lmer)
anova(CH4_rip_lmer)
plot(CH4_rip_lmer, id=0.05) #residuals look ok
qqnorm(resid(CH4_rip_lmer))
vif(CH4_rip_lmer)


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

top.CH4_rip_dredge  <- get.models(CH4_rip_dredge , subset=delta <= 2) #The top models, which include those with an AIC <2 

sw(CH4_rip_dredge) #importance of variables

CH4rip_mavg <- model.avg(top.CH4_rip_dredge) #only one top model
mA_CH4_flow_i <- summary(CH4rip_mavg)

#Run model with soil temp and soil VWC (the variables from the single top model)
CH4_rip_lmer2 <- lmer(log(CH4_mg_m2_d+2.842777)~ soil_VWC+ soil_temp + (1 + campaign | Site), data=CH4_ripNA_scaled)
summary(CH4_rip_lmer2)

r.squaredGLMM(CH4_rip_lmer2)


##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors 
CH4_rip_lmer2_df <- broom.mixed::tidy(x = CH4_rip_lmer2)  #make data frame of model outputs

CH4_rip_lmer2_df <- CH4_rip_lmer2_df %>%
  dplyr::mutate(
    .data = .,
    conf.low = estimate - std.error,
    conf.high = estimate + std.error
  )

#rename variables
CH4_rip_lmer2_df$term[CH4_rip_lmer2_df$term == "soil_VWC"] <- "Soil_moist"
CH4_rip_lmer2_df$term[CH4_rip_lmer2_df$term == "soil_temp"] <- "Soil_temp"


#CH4_rip_lmer2_df$term = reorder(as.factor(CH4_rip_lmer2_df$term), -CH4_rip_lmer2_df$p.value, FUN = mean) #reorder factor levels by p


#plot
CH4_rip_standard <- ggplot(subset(CH4_rip_lmer2_df, effect=="fixed" & term!="(Intercept)"), aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 1.2, y = 0.32946 , label = "***", size=5) +
  annotate("text", x = 2.2, y = -0.08092 , label = "***", size=5)
CH4_rip_standard


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
dat_N2O_dry <-  subset(dat_N2O,  habitat == "Dry")  #select dry 49 obs
#Select relevant columns
dat_N2O_dry1 <- data.table(dat_N2O_dry[, c("campaign", "Site" , "mg_N2ON_m2_d",  "m_to_source" , "Freq_Flow",  "Stock_Benth_g.m2",  "Sed_OM" ,  "sed_VWC","sed_temp",  "percent_N" , "MaxDur_Ponded")])
dat_N2O_dry1NA <- na.omit(dat_N2O_dry1) #remove NAs 49 to 34 obs, missing values in Tsince manual make it 19 obs (too many NAs for BR01, and GR01, not possible to fill in)

#remove the outliers we see in the residual plot Al02, >0.1
dat_N2O_dry1NAout<- dat_N2O_dry1NA %>%
  filter(!mg_N2ON_m2_d >= 0.1)

min(dat_N2O_dry1NAout$mg_N2ON_m2_d) #-0.07417401


dat_N2O_dry1NAout$log_m_to_source <- log(dat_N2O_dry1NAout$m_to_source)

dat_N2O_dry1NAout_scaled <- dat_N2O_dry1NAout  %>% mutate_at(c("log_m_to_source" , "Freq_Flow",  "Stock_Benth_g.m2",  "Sed_OM" ,  "sed_VWC","sed_temp", "percent_N" , "MaxDur_Ponded"), ~(scale(., center=FALSE) %>% as.vector))


N2O_dry_lmer <- lmer(mg_N2ON_m2_d~     #log(mg_N2ON_m2_d+1.07417401)
                       Freq_Flow*log_m_to_source + 
                       #Tsince_manual*log(m_to_source) +  
                       Sed_OM*log_m_to_source+ 
                       Stock_Benth_g.m2*log_m_to_source+ 
                       Sed_OM*sed_VWC + 
                       Sed_OM*sed_temp +
                       sed_VWC*sed_temp + 
                       asin(sqrt(percent_N/100)) + 
                       MaxDur_Ponded +
                       (1 + campaign | Site), data=dat_N2O_dry1NAout_scaled) 

summary(N2O_dry_lmer) 
anova(N2O_dry_lmer)
plot(N2O_dry_lmer) #Looks good after we removed the Al02 outlier
qqnorm(resid(N2O_dry_lmer))



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

top.N2O_dry_dredge  <- get.models(N2O_dry_dredge , subset=delta <= 2) #The top models, which include those with an AIC >2 away from the top model 

sw(N2O_dry_dredge)  # Variable importance

N2Odry_mavg <- model.avg(top.N2O_dry_dredge)
summary(N2Odry_mavg)


N2O_dry_N_lmer <- lmer(mg_N2ON_m2_d~     #log(mg_N2ON_m2_d+1.07417401)
                         asin(sqrt(percent_N/100)) + sed_VWC +
                         (1 + campaign | Site), data=dat_N2O_dry1NAout) 

summary(N2O_dry_N_lmer ) #still not significant

## Make a random forest plot of the model averaged estimates (standardized coefficients)

N2Odry_mavg <- model.avg(top.N2O_dry_dredge)
mA <- summary(N2Odry_mavg)

# Extract the coefficients and standard errors 
N2Odry_coef <- coefTable(mA, full=FALSE) 
N2Odry_coef <- as.data.frame(N2Odry_coef)  
N2Odry_coef <- tibble::rownames_to_column(N2Odry_coef, "Term") #need to keep row names as column

names(N2Odry_coef)[names(N2Odry_coef) == "Std. Error"] <- "StdError"

N2Odry_coef <- N2Odry_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError
  )

#rename variables
N2Odry_coef$Term[N2Odry_coef$Term == "asin(sqrt(percent_N/100))"] <- "N_content"

#add p values from model average output : not sig

#plot
N2O_dry_standard <- ggplot(subset(N2Odry_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) 
#annotate("text", x = 4.2, y = 1.5058001 , label = "**", size=5)  
N2O_dry_standard


### For perennial flowing N2O emissions ####

#Need to use data with no NAs
dat_N2O_flow <-  subset(dat_N2O,  habitat == "Flowing")  #select flowing
dat_N2O_flow_p <-  subset(dat_N2O_flow,  drying_regime == "perennial")

dat_N2O_flow_p1 <- data.table(dat_N2O_flow_p[, c("campaign", "Site" , "mg_N2ON_m2_d", "m_to_source",  "Stock_Benth_g.m2",  "Sed_OM","Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C", "NH4_N_ugL" ,  "DN_ugL" , "percent.IR.up")]) #Select relevant columns
dat_N2O_flow_p1NA <- na.omit(dat_N2O_flow_p1) #remove NAs  # 151, from 163


min(dat_N2O_flow_p1NA$mg_N2ON_m2_d) #-0.06030342


#Scale the variables
dat_N2O_flow_p1NA$log_m_to_source <- log(dat_N2O_flow_p1NA$m_to_source)

dat_N2O_flow_p1NA_scaled <- dat_N2O_flow_p1NA  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2",  "Sed_OM","Mean_Velocity..m.s.", "DO_mg_L", "Temp_water_C", "NH4_N_ugL" ,  "DN_ugL" , "percent.IR.up"), ~(scale(., center=FALSE) %>% as.vector))

dat_N2O_flow_p1NA$campaign
#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
N2O_flow_p_lmer <- lmer(log(mg_N2ON_m2_d+1.06030342)~ 
                          log(m_to_source)*Stock_Benth_g.m2 + 
                          log(m_to_source)*Sed_OM + 
                          percent.IR.up + 
                          Temp_water_C*Stock_Benth_g.m2 + 
                          Temp_water_C*Sed_OM  + 
                          Mean_Velocity..m.s. +
                          DO_mg_L+
                          NH4_N_ugL +
                          DN_ugL + 
                          (1 + campaign | Site), data=dat_N2O_flow_p1NA_scaled) 

summary(N2O_flow_p_lmer) 
anova(N2O_flow_p_lmer)
plot(N2O_flow_p_lmer) #improve after log transformation
qqnorm(resid(N2O_flow_p_lmer))


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

N2Oflow_p_mavg <- model.avg(top.N2O_flow_p_dredge) #only 1 top model
summary(N2Oflow_p_mavg)

#A single model within 2 delta AIC, including DO, mean velocity, NO3
N2O_flow_p_lmer2 <- lmer(log(mg_N2ON_m2_d+1.06030342)~ 
                           Mean_Velocity..m.s. +
                           DN_ugL +
                           (1 + campaign | Site), data=dat_N2O_flow_p1NA_scaled)
summary(N2O_flow_p_lmer2)

r.squaredGLMM(N2O_flow_p_lmer2)

##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors 
N2O_flow_p_lmer2_df <- broom.mixed::tidy(x = N2O_flow_p_lmer2)  #make data frame of model outputs

N2O_flow_p_lmer2_df <- N2O_flow_p_lmer2_df %>%
  dplyr::mutate(
    .data = .,
    conf.low = estimate - std.error,
    conf.high = estimate + std.error
  )

#rename variables
N2O_flow_p_lmer2_df$term[N2O_flow_p_lmer2_df$term == "Mean_Velocity..m.s."] <- "Velocity"
N2O_flow_p_lmer2_df$term[N2O_flow_p_lmer2_df$term == "DN_ugL"] <- "Dissolved nitrogen"


#plot

N2O_flow_p_standard <- ggplot(subset(N2O_flow_p_lmer2_df, effect=="fixed" & term!="(Intercept)"), aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 2.2, y = 0.33675 , label = "***", size=5) +
  annotate("text", x = 1.2, y = 0.47151 , label = "***", size=5) 
N2O_flow_p_standard




#### For flowing intermittent N2O emissions ####

#Need to use data with no NAs
dat_N2O_flow <-  subset(dat_N2O,  habitat == "Flowing")  #select flowing
dat_N2O_flow_i <-  subset(dat_N2O_flow,  drying_regime == "intermittent_long" | drying_regime == "intermittent_moderate" ) #select intermittent
dat_N2O_flow_i1 <- data.table(dat_N2O_flow_i[, c("campaign", "Site" , "mg_N2ON_m2_d", "Stock_Benth_g.m2", "Sed_OM", "m_to_source" ,"Tsince_manual", "Freq_Flow",  "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL",  "DN_ugL" )]) #Select relevant columns
dat_N2O_flow_i1NA <- na.omit(dat_N2O_flow_i1) #remove NAs #Go from 278 to 243, 71 with manual tsince rewetting

dat_N2O_flow_i1NAout <- dat_N2O_flow_i1NA[-72,]


min(dat_N2O_flow_i1NAout$mg_N2ON_m2_d) #-0.04485442

#Scale the variables
dat_N2O_flow_i1NAout$log_m_to_source <- log(dat_N2O_flow_i1NAout$m_to_source)

dat_N2O_flow_i1NAout_scaled <- dat_N2O_flow_i1NAout  %>% mutate_at(c("log_m_to_source",  "Stock_Benth_g.m2", "Sed_OM", "m_to_source" ,"Tsince_manual", "Freq_Flow",  "Temp_water_C", "DO_mg_L", "Mean_Velocity..m.s." , "NH4_N_ugL",  "DN_ugL"), ~(scale(., center=FALSE) %>% as.vector))


#Create saturated model 
#using variables based on hypotheses as well as top correlations
#Also include any interactions that make biological sense
#log(mg_N2ON_m2_d+1.04485442)
N2O_flow_i_lmer <- lmer(log(mg_N2ON_m2_d+1.04485442) ~ 
                          log_m_to_source*Freq_Flow +
                          log_m_to_source*Tsince_manual +
                          log_m_to_source*Sed_OM + 
                          log_m_to_source*Stock_Benth_g.m2 + 
                          Temp_water_C*Sed_OM +
                          Temp_water_C*Stock_Benth_g.m2 +
                          Mean_Velocity..m.s. + 
                          NH4_N_ugL+
                          DN_ugL + 
                          (1 + campaign | Site), data=dat_N2O_flow_i1NAout_scaled)


summary(N2O_flow_i_lmer) 
anova(N2O_flow_i_lmer)
plot(N2O_flow_i_lmer) #Improved after log transformation
qqnorm(resid(N2O_flow_i_lmer))

#identify outliers with Cooks distance
cooksD <- cooks.distance(N2O_flow_i_lmer)
n <- nrow(dat_N2O_flow_i1NA)
plot(cooksD, main = "Cooks Distance for Influential Obs")
abline(h = 4/n, lty = 2, col = "steelblue") # add cutoff line
influential_obs <- as.numeric(names(cooksD)[(cooksD > (25/n))])
influential_obs #72

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

N2Oflow_i_mavg <- model.avg(top.N2O_flow_i_dredge) #Only one top model
summary(N2Oflow_i_mavg)

#Run a model with variable from single top model
N2O_flow_i_lmer2 <- lmer(log(mg_N2ON_m2_d+1.04485442) ~ 
                           Temp_water_C +
                           (1 + campaign | Site), data=dat_N2O_flow_i1NAout_scaled)

summary(N2O_flow_i_lmer2) 

r.squaredGLMM(N2O_flow_i_lmer2) #0.83 conditional R2


##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors 
N2O_flow_i_lmer2_df <- broom.mixed::tidy(x = N2O_flow_i_lmer2)  #make data frame of model outputs

N2O_flow_i_lmer2_df <- N2O_flow_i_lmer2_df %>%
  dplyr::mutate(
    .data = .,
    conf.low = estimate - std.error,
    conf.high = estimate + std.error
  )

#rename variables
N2O_flow_i_lmer2_df$term[N2O_flow_i_lmer2_df$term == "Temp_water_C"] <- "Water temperature"



#plot

N2O_flow_i_standard <- ggplot(subset(N2O_flow_i_lmer2_df, effect=="fixed" & term!="(Intercept)"), aes(x=term, y=estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) 
N2O_flow_i_standard



#### For riparian N2O emissions ####

#Need to use data with no NAs
dat_N2O_rip <-  subset(dat_N2O,  habitat == "Riparian") 

#Select relevant variables
dat_N2O_rip1 <- data.table(dat_N2O_rip[, c("campaign", "Site", "mg_N2ON_m2_d", "soil_temp", "soil_VWC", "Soil_OM", "Stock_Rip_g.m2", "flow_state" ,   "totalbiomass_g",   "Freq_Flow", "m_to_source", "percent_N")])

dat_N2O_ripNA <- na.omit(dat_N2O_rip1) #remove NAs  280 from 323

min(dat_N2O_ripNA$mg_N2ON_m2_d) #-0.1833087 #Min N2O value in order to add to log 

#remove two outliers based on residual plot (AL07 and AL03)

#remove  outliers identified in residuals plot
dat_N2O_ripNAout <- dat_N2O_ripNA %>%
  filter(!mg_N2ON_m2_d >= 1)


#Scale the variables
dat_N2O_ripNAout$log_m_to_source <- log(dat_N2O_ripNAout$m_to_source)

dat_N2O_ripNAout_scaled <- dat_N2O_ripNAout  %>% mutate_at(c("log_m_to_source", "soil_temp", "soil_VWC", "Soil_OM", "Stock_Rip_g.m2",   "totalbiomass_g",   "Freq_Flow",  "percent_N"), ~(scale(., center=FALSE) %>% as.vector))


# Create saturated model 
#using variables based on hypotheses as well as top correlations and interactions
N2O_rip_lmer <- lmer(log(mg_N2ON_m2_d+1.1833087)~ 
                       Freq_Flow*log_m_to_source+ 
                       flow_state*log_m_to_source+
                       totalbiomass_g*log_m_to_source+
                       Soil_OM*log_m_to_source+ 
                       soil_temp*soil_VWC + 
                       soil_VWC*Soil_OM + 
                       Soil_OM*soil_temp + 
                       percent_N+
                       (1 + campaign | Site), data=dat_N2O_ripNAout_scaled)

summary(N2O_rip_lmer)
anova(N2O_rip_lmer)
plot(N2O_rip_lmer)#improve after log transformation
qqnorm(resid(N2O_rip_lmer))

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

##  plot the standardized coefficients from the top model in a forest plot ##
#plot the standardized coefficicents from the averaged top models in a forest plot

# Extract the coefficients and standard errors using coef() and se.coef()
N2Orip_coef <- coefTable(mA_N2O_rip, full=FALSE) 
N2Orip_coef <- as.data.frame(N2Orip_coef)  
N2Orip_coef <- tibble::rownames_to_column(N2Orip_coef, "Term") #need to keep row names as column

names(N2Orip_coef)[names(N2Orip_coef) == "Std. Error"] <- "StdError"

N2Orip_coef <- N2Orip_coef %>%
  dplyr::mutate(
    conf.low = Estimate - StdError,
    conf.high = Estimate + StdError)

#rename variables
N2Orip_coef$Term[N2Orip_coef$Term == "soil_temp:soil_VWC"] <- "Soil_temp: Soil_moist"
N2Orip_coef$Term[N2Orip_coef$Term == "soil_temp"] <- "Soil_temp"
N2Orip_coef$Term[N2Orip_coef$Term == "soil_VWC"] <- "Soil_moist"

#get p values form model average output 
#soil_temp:soil_VWC  0.00753 **

#plot

my.labels <- c( "Soil_moist",
                "Soil_temp", 
                "Soil_temp*Soil_moist")


N2O_rip_standard <- ggplot(subset(N2Orip_coef, Term!="(Intercept)"), aes(x=Term, y=Estimate, ymin=conf.low, ymax=conf.high)) +
  geom_pointrange() + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  ylab("Standardized coefficients") +
  theme_bw()  + # use a white background 
  scale_x_discrete(labels= my.labels)   +
  theme(axis.title.y=element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  annotate("text", x = 3.2, y = 0.11567 , label = "**", size=5) 
N2O_rip_standard





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

CO2_habitat_means <- ggplot(CO2means, aes(x=habitat, y=mean_CO2_g_m2_d)) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.3, position = position_jitterdodge(jitter.width=2.0)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(),axis.text.x=element_text(angle = 45, hjust = .9), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c( "#878B8C","#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C","#5598B7", "#FCD262" , "#527028" ))  + ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))  +  
  theme(panel.border = element_rect(color = "black",  fill = NA, linewidth = 1)) +
  annotate("text", x = 1.8, y = 17, label=sprintf('\u2191')) +
  #scale_y_continuous(trans='log2') +
  annotate("text", x = 2, y = 17, label="30")  + ylim(-0.5, 21.7) +
  geom_bracket(
    xmin = c(1, 1,2, 3), xmax = c(2, 4,3, 4),
    y.position = c(21.1, 20.2, 18, 17), label = c(" ", " ", " ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x =1.5, y = 21.2 , label = "**", size=5) +
  annotate("text", x =2.5, y = 20.3 , label = "**", size=5) +
  annotate("text", x =2.5, y = 18.1 , label = "***", size=5) +
  annotate("text", x =3.5, y = 17.1 , label = "***", size=5) 
CO2_habitat_means 


#### CH4 by habitat using Site means ####

CH4means$habitat <- factor(CH4means$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

CH4_habitat_means <- ggplot(CH4means, aes(x=habitat, y=mean_CH4_mg_m2_d))  + geom_hline(yintercept=0, linetype="dashed", color = "black", size=0.5) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.15, position = position_jitterdodge(jitter.width=2.0)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.text.x=element_text(angle = 45, hjust = .9), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1))  +  
  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(-1.4, 8.3) +
  annotate("text", x = 2.1, y = 7, label="11.2")  + 
  annotate("text", x = 1.7, y = 7.02, label=sprintf('\u2191')) +
  geom_bracket(
    xmin = c(2, 3), xmax = c(4,4),
    y.position = c(8, 7), label = c(" ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x =3, y = 8.1 , label = "***", size=5) +
  annotate("text", x =3.5, y = 7.1 , label = "***", size=5) 
CH4_habitat_means 

#### N2O by habitat using Site means ####
N2Omeans$habitat <- factor(N2Omeans$habitat, levels = c("Pool", "Flowing", "Dry", "Riparian"))

N2O_habitat_means <- ggplot(N2Omeans, aes(x=habitat, y=mean_mg_N2ON_m2_d)) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=habitat), size=4, shape=16, alpha=0.3, position = position_jitterdodge(jitter.width=2.0)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.text.x=element_text(angle = 45, hjust = .9), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_fill_manual(values=c( "#878B8C", "#5598B7", "#FCD262" , "#527028")) +  scale_colour_manual(values=c("#878B8C", "#5598B7", "#FCD262" , "#527028", "#878B8C"))  + ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) +  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(-0.03, 1.45) +
  annotate("text", x = 1.6, y = 1.105, label=sprintf('\u2191')) +
  annotate("text", x = 2.15, y = 1.1, label= "1.7, 2.6") +
  geom_bracket(
    xmin = c(2, 2), xmax = c(3,4),
    y.position = c(1.39, 1.25), label = c(" ", " "),
    tip.length = 0.01, label.size = 5) +
  annotate("text", x = 2.5, y = 1.40, label= "*") +
  annotate("text", x = 3, y = 1.26, label= "**") 
N2O_habitat_means 

#### CO2 flowing by drying regime ####

CO2means$drying_regime <- dplyr::recode(CO2means$drying_regime, intermittent_long = 'Int-L', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten


CO2_drying_regime  <- ggplot(subset(CO2means, habitat=="Flowing"), aes(x=drying_regime, y=mean_CO2_g_m2_d)) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.6, position = position_jitterdodge(jitter.width=.71)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1"))  + ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))  +  theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) + ylim(0.082, 32) +
  geom_bracket(
    xmin = c(1, 2), xmax = c(3, 3),
    y.position = c(30.5, 27.5 ), label = c(" ", " "),
    tip.length = 0.01,  label.size = 5) +
  annotate("text", x =2, y = 30.9 , label = "*", size=5) +
  annotate("text", x =2.5, y = 27.9 , label = "*", size=5) 

CO2_drying_regime 


#### CH4 flowing by drying regime ####

CH4means$drying_regime <- dplyr::recode(CH4means$drying_regime, intermittent_long = 'Int-L', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten

CH4_drying_regime  <- ggplot(subset(CH4means, habitat=="Flowing"), aes(x=drying_regime, y=mean_CH4_mg_m2_d)) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.5, position = position_jitterdodge(jitter.width=.71)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1")) + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) 
CH4_drying_regime 



#### N2O flowing by drying regime ####

N2Omeans$drying_regime <- dplyr::recode(N2Omeans$drying_regime, intermittent_long = 'Int-L', 
                                        intermittent_moderate = 'Int-M',
                                        perennial = 'Per') #rename to shorten

N2O_drying_regime  <- ggplot(subset(N2Omeans, habitat=="Flowing"), aes(x=drying_regime, y=mean_mg_N2ON_m2_d)) + geom_boxplot( outlier.shape = NA) +   theme_pubr() + geom_point(aes(colour=drying_regime), size=4, shape=16, alpha=0.6, position = position_jitterdodge(jitter.width=.71)) + theme(legend.position = "none", axis.text = element_text(size = 14), axis.title.x=element_blank(), axis.title.y = element_text(size = 14))  +  scale_colour_manual(values=c("#C7E2B3",  "#3CB8C3", "#1C54A1"))  + ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) + theme(panel.border = element_rect(color = "black",  fill = NA,size = 1)) 
N2O_drying_regime 


#### Plot model coefficients for dry sediments/riparian soils ####

#Read in the model coefficients for dry sediment and riparian GHG fluxes
combined_data_dry <- read.csv("Model_coefficients_GHGdry.csv")
combined_data_rip <- read.csv("Model_coefficients_GHGrip.csv")

str(combined_data_dry) #6 obs of 17 vars
str(combined_data_rip) #7 obs of 17 vars

# Plot the dry sediment model coefficients using ggplot
rf_dry <- ggplot(combined_data_dry, aes(x=label, y=y, ymin=ymin, ymax=ymax, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_pointrange(color="#FCD262", size=0.9, alpha=0.7) +
  coord_flip() +
  ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = aes(color = "grey"))) +
  theme(axis.title.y=element_blank(), plot.title = element_text(hjust = -0.73, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  ggtitle("A. Dry sediments") +
  annotate("text", x = 1.2, y = 1.5058001 , label = "**", size=5)
rf_dry

# Plot the riparian soil model coefficients using ggplot
rf_rip <- ggplot(combined_data_rip, aes(x=label, y=y, ymin=ymin, ymax=ymax, group=variable, shape=variable)) +
  geom_hline(yintercept=0, lty=2) +
  geom_pointrange(color="#527028", size=0.9, alpha=0.7) +
  coord_flip() +
  ylab("Standardized coefficients") +
  scale_y_continuous(limits = symmetric_limits) +
  theme_bw() +
  theme(axis.title.y=element_blank(), plot.title = element_text(hjust = -0.73, face = "bold", size = 14), legend.title = element_blank(), text = element_text(size = 13), axis.text = element_text(size = 13, colour="black")) +
  ggtitle("B. Riparian soils") +
  guides(shape = guide_legend(override.aes = aes(color = "grey"))) +
  annotate("text", x = 3.2, y = 0.24791515, label = "***", size=5)  +
  annotate("text", x = 2.2, y = -0.09032863, label = "***", size=5) +
  annotate("text", x = 2.2, y = 0.32946 , label = "***", size=5) +
  annotate("text", x = 3.2, y = -0.08092 , label = "***", size=5) +
  annotate("text", x = 1.2, y = 0.11567 , label = "**", size=5) 
rf_rip

#### Relationship between flowing GHGs and top predictors  ####

## Intermittent CO2 vs % embededness
CO2_embed <- ggplot(subset(CO2_flow_iNA), aes(Embeddedness...., CO2_g_m2_d)) + geom_point(colour="#8BCEB9", size=3.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(),  legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))   + xlab("Embeddedness (%)") + #+ ylim(-2, 11) 
  stat_smooth(geom="line", color="grey", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CO2_embed

#Intermittent CO2 vs flow permanence 
CO2freqflow <- ggplot(CO2_flow_iNA, aes(Freq_Flow*100, CO2_g_m2_d)) + geom_point(colour="#8BCEB9",size=3.5, alpha=0.3)  +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.position = c(.05, .95), legend.justification = c("left", "top"), legend.box.just = "right", legend.margin = margin(3, 3, 3, 3),  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), legend.spacing.y = unit(0, "mm"), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))   + xlab("Flow permanence (%)") +
  stat_smooth(geom="line", color="grey", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CO2freqflow

#Intermittent CO2 vs water temperature
CO2_water_temp <- ggplot(CO2_flow_iNA, aes(Temp_water_C, CO2_g_m2_d)) + geom_point(colour="#8BCEB9", size=3.5, alpha=0.3) + theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),axis.ticks.x=element_blank(), legend.position = c(.75, .95), legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))   + xlab("Water temperature (\u00B0C)") +
  stat_smooth(geom="line", color="grey", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CO2_water_temp

# Intermittent CH4 vs water temperature
CH4temp <- ggplot(dat_CH4_num_flow_i1NA, aes(Temp_water_C, CH4_mg_m2_d)) + geom_point(colour="#8BCEB9", shape=17, size=3.5, alpha=0.3)  +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right",  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), legend.spacing.y = unit(0, "mm"), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("Water temperature (\u00B0C)") + stat_smooth(geom="line", color="grey", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CH4temp

#Intermittent CH4 vs distance to source 
CH4mtosource <- ggplot(dat_CH4_num_flow_i1NA, aes(m_to_source/100, CH4_mg_m2_d, )) + geom_point(colour="#8BCEB9", shape=17, size=3.5, alpha=0.3)  +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right",  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), legend.spacing.y = unit(0, "mm"), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) + xlab("Distance to source (km)") +  ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + scale_x_continuous(trans='log2')  +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CH4mtosource 

# Intermittent CH4 vs OM stock
CH4OMint <- ggplot(subset(dat_CH4_num_flow_i1NA,Stock_Benth_g.m2<15) , aes(Stock_Benth_g.m2, CH4_mg_m2_d)) + geom_point(colour="#8BCEB9", shape=17, size=3.5, alpha=0.3)  +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right", legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), legend.spacing.y = unit(0, "mm"), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab(expression(`OM stock g`~m^-2)) +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) + xlim(0, 2.5)
CH4OMint

#Intermittent N2O vs water temperature
N2O_water_temp <- ggplot(dat_N2O_flow_i1NAout, aes(Temp_water_C, mg_N2ON_m2_d)) + geom_point(colour="#8BCEB9", size=3.5, shape=15, alpha=0.3) + theme_bw() + theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1))+  xlab("Water temperature (\u00B0C)") +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x)
N2O_water_temp

#Perennial CO2 vs sediment OM
CO2_flow_pNA$position <- dplyr::case_when(
  CO2_flow_pNA$m_to_source < 1000 ~ "Headwaters",
  CO2_flow_pNA$m_to_source >= 1000 ~ "Mainstem",
  TRUE ~ NA_character_
)
median(CO2_flow_pNA$m_to_source)

CO2_flowp_OM<- ggplot(CO2_flow_pNA, aes(Sed_OM, CO2_g_m2_d)) + geom_point(aes(colour=position), size=4.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  legend.position = c(.5, .95), axis.ticks.x=element_blank(), legend.title=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))  + xlab("Sediment OM content (%)") +   stat_smooth(aes(colour=position), geom="line", size=1, alpha=0.7, method = "lm", se=FALSE,  formula = y ~ x) + scale_color_manual(values=c("#1C54A1", "#a1691c")) 
CO2_flowp_OM

#Perennial CO2 vs percent non-perennial upstream
CO2_percentIR <- ggplot(CO2_flow_pNA, aes(percent.IR.up, CO2_g_m2_d)) + geom_point(colour="#1C54A1", size=4.5, alpha=0.3) +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.position = c(.05, .99), legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(g~CO[2]*`-C`~m^-2*~d^-1))   + xlab("Non-per. reaches upstream (%)") + stat_smooth(geom="line", color="grey", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x) 
CO2_percentIR


# Perennial CH4 vs sediment OM
CH4OMper <- ggplot(dat_CH4_num_flow_p1NA_out, aes(Sed_OM, CH4_mg_m2_d)) + geom_point(colour= "#1C54A1" , shape=17, size=3.5, alpha=0.3)  +  theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.box.just = "right",  legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'), legend.spacing.y = unit(0, "mm"), axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black"))  + ylab(expression(mg~CH[4]*`-C`~m^-2~d^-1)) + xlab("Sediment OM content (%)") +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x)
CH4OMper

# Perennial N2O vs velocity
N2O_veloc <- ggplot(dat_N2O_flow_p1NA, aes(Mean_Velocity..m.s., mg_N2ON_m2_d)) + geom_point(colour="#1C54A1", shape=15,size=3.5, alpha=0.3) + theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.position = c(.45, .95), legend.justification = c("left", "top"), legend.box.just = "right",axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) +labs(x = expression(paste("Mean velocity (m/ [s^{-1} ])")))  + xlab(Mean~velocity~(m~s^-1)) +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x)
N2O_veloc

# Perennial N2O vs dissolved nitrogen
N2O_DN <- ggplot(dat_N2O_flow_p1NA, aes(DN_ugL, mg_N2ON_m2_d)) + geom_point(colour="#1C54A1", size=3.5, alpha=0.3, shape=15) + theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(), legend.title=element_blank(), axis.ticks.x=element_blank(), legend.position = c(.45, .95), legend.justification = c("left", "top"), legend.box.just = "right", axis.line = element_line(colour = "black"), text = element_text(size = 12), axis.text = element_text(size = 12, colour="black")) +  ylab(expression(mg~N[2]*`O-`*N~m^-2~d^-1)) + xlab(expression("Dissolved nitrogen"~(mu*g~L^-1))) +  stat_smooth(colour="grey", geom="line", size=1, alpha=0.5, method = "lm", se=FALSE,  formula = y ~ x)
N2O_DN

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


#### Make network maps for OM stock and sediment OM ####

# load in shapefile
albarine <- st_read("Hydro_Albarine.shp")

#### Spatial plot : dry sediment OM content ####

dryOM <- subset(dat, habitat=="Dry") %>% 
  group_by(Site, drying_regime) %>% # we group by the two values
  dplyr::summarise(
    Sed_OM = mean(Sed_OM, na.rm=T),   #do I need to use na.rm=T  ?
    latitude = mean(latitude),
    longitude = mean(longitude) )


site_locations <- st_as_sf(dryOM, coords = c("longitude", "latitude"), crs = 4326)

site_locations$drying_regime <- factor(site_locations$drying_regime, levels = c("perennial", "intermittent_moderate", "intermittent_long")) #reorder levels

site_locations$drying_regime <- dplyr::recode(site_locations$drying_regime, intermittent_long = 'Int-L', 
                                              intermittent_moderate = 'Int-M',
                                              perennial = 'Per') #rename to shorten

#plot 

tiff("dry_sedOM", units="in", width=3.5, height=3.5, res=300, type="cairo")

dry_sedOM<-ggplot() + 
  geom_sf(data = albarine) + 
  geom_sf(data = site_locations, aes(color = Sed_OM, shape= drying_regime), size = 4) +
  #scale_color_viridis_c(option="inferno", na.value="grey40",oob=scales::squish) +
  scale_colour_distiller(palette="YlOrRd", direction=1) +
  scale_shape_manual(values=c(17, 15))+
  theme_classic()+
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title="Sediment OM (%)")+ 
  theme(plot.title = element_text(vjust = -10, hjust=.25), 
        axis.text = element_text(size = 8),
        legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=8)) 
dry_sedOM

dev.off()


#### Spatial plot : flowing sediment OM content ####
flowingOM <- subset(dat, habitat=="Flowing") %>% 
  group_by(Site, drying_regime) %>% # we group by the two values
  dplyr::summarise(
    Sed_OM = mean(Sed_OM, na.rm=T),   #do I need to use na.rm=T  ?
    latitude = mean(latitude),
    longitude = mean(longitude) )


site_locations <- st_as_sf(flowingOM, coords = c("longitude", "latitude"), crs = 4326)

site_locations$drying_regime <- factor(site_locations$drying_regime, levels = c("perennial", "intermittent_moderate", "intermittent_long")) #reorder levels

site_locations$drying_regime <- dplyr::recode(site_locations$drying_regime, intermittent_long = 'Int-L', 
                                              intermittent_moderate = 'Int-M',
                                              perennial = 'Per') #rename to shorten

#plot 
tiff("flowing_OM", units="in", width=3.5, height=3.5, res=300, type="cairo")

flowing_OM<-ggplot() + 
  geom_sf(data = albarine) + 
  geom_sf(data = site_locations, aes(color = Sed_OM, shape=drying_regime), size = 4) +
  #scale_color_viridis_c(option="plasma", na.value="grey40", oob=scales::squish) +
  scale_colour_distiller(palette="YlGnBu", direction=1) +
  theme_classic()+
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title="Sediment OM (%)")+ 
  theme(plot.title = element_text(vjust = -10, hjust=.25), 
        axis.text = element_text(size = 8),
        legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=8)) 
#coord_sf(xlim = c(876115 , 903855 ), ylim = c(6534240 , 6557250 ), expand = FALSE)

flowing_OM

dev.off()

#### Spatial plot : dry sediment OM stock ####
dryOMstock <- subset(dat, habitat=="Dry") %>% 
  group_by(Site, drying_regime) %>% # we group by the two values
  dplyr::summarise(
    Benth_stock = mean(Stock_Benth_g.m2, na.rm=T),   #do I need to use na.rm=T  ?
    latitude = mean(latitude),
    longitude = mean(longitude) )


site_locations <- st_as_sf(dryOMstock, coords = c("longitude", "latitude"), crs = 4326)

site_locations$drying_regime <- factor(site_locations$drying_regime, levels = c("perennial", "intermittent_moderate", "intermittent_long")) #reorder levels

site_locations$drying_regime <- dplyr::recode(site_locations$drying_regime, intermittent_long = 'Int-L', 
                                              intermittent_moderate = 'Int-M',
                                              perennial = 'Per') #rename to shorten

#plot 
tiff("dry_OM_stock", units="in", width=3.5, height=3.5, res=300, type="cairo")

dry_OM_stock<-ggplot() + 
  geom_sf(data = albarine) + 
  geom_sf(data = site_locations, aes(color = Benth_stock, shape=drying_regime), size = 4) +
  #scale_color_viridis_c(option="plasma", na.value="grey40", oob=scales::squish) +
  scale_colour_distiller(palette="YlOrRd", direction=1) +
  scale_shape_manual(values=c(17, 15))+
  theme_classic()+
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = bquote("OM stock (g m"^{-2}*")")) +
  theme(plot.title = element_text(vjust = -10, hjust=.25), 
        axis.text = element_text(size = 8),
        legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=8)) 
dry_OM_stock

dev.off()

#### Spatial plot : flowing sediment OM stock ####
flowingOMstock <- subset(dat, habitat=="Flowing") %>% 
  group_by(Site, drying_regime) %>% # we group by the two values
  dplyr::summarise(
    Benth_stock = mean(Stock_Benth_g.m2, na.rm=T),   #do I need to use na.rm=T  ?
    latitude = mean(latitude),
    longitude = mean(longitude) )


site_locations <- st_as_sf(flowingOMstock, coords = c("longitude", "latitude"), crs = 4326)

site_locations$drying_regime <- factor(site_locations$drying_regime, levels = c("perennial", "intermittent_moderate", "intermittent_long")) #reorder levels

site_locations$drying_regime <- dplyr::recode(site_locations$drying_regime, intermittent_long = 'Int-L', 
                                              intermittent_moderate = 'Int-M',
                                              perennial = 'Per') #rename to shorten


#plot 
tiff("flowing_OM_stock", units="in", width=3.5, height=3.5, res=300, type="cairo")

flowing_OM_stock<-ggplot() + 
  geom_sf(data = albarine) + 
  geom_sf(data = site_locations, aes(color = Benth_stock, shape=drying_regime), size = 4) +
  #scale_color_viridis_c(option="plasma", na.value="grey40",oob=scales::squish) +
  scale_colour_distiller(palette="YlGnBu", direction=1) +
  theme_classic()+
  xlab("Longitude") + 
  ylab("Latitude") + 
  labs(title = bquote("OM stock (g m"^{-2}*")")) + 
  theme(plot.title = element_text(vjust = -10, hjust=.25), 
        axis.text = element_text(size = 8),
        legend.position ="bottom",
        legend.title = element_blank(),
        legend.text = element_text(size=8)) 
flowing_OM_stock

dev.off()


#### Figure S5. OM stock vs sediment OM ####

tiff("OMstock_OMcontent", units="in", width=6, height=5, res=300)

OM <- ggplot(CO2means, aes(Stock_Benth_g.m2, Sed_OM)) + geom_point(aes(colour=m_to_source/1000),size=3.5, alpha=0.3) + scale_colour_viridis()+ theme_bw() +  theme(axis.title = element_text(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), panel.border = element_blank(),  axis.ticks.x=element_blank(), legend.justification = c("left", "top"), legend.position = "bottom",  legend.box.just = "right",axis.line = element_line(colour = "black"), text = element_text(size = 15), axis.text = element_text(size = 15, colour="black"))  +  
  xlab(bquote("OM stock (g m"^{-2}*")"))  + ylab("Sediment OM content (%)") +
  geom_smooth(colour="black", method = "lm", linetype="dashed", se=FALSE,  formula = y ~ x)  +
  stat_regline_equation(label.x = 40, label.y = 8 ) + stat_cor(label.x = 40, label.y = 7 ) +
  ylim(0, 8) +
  labs(color = "km from source") +
  annotate("text", x =1, y = 7.4 , label = "17", size=4, colour="purple") +
  annotate("text", x =1, y = 7.8, label=sprintf('\u2191'), colour="purple") 
OM

dev.off()


##____________________________________________________________________________##

#### Combine plots ####

#### Figure 2. GHGs by habitat using Site means ####
tiff("MeanGHGbyhabitat", units="in", width=3, height=8.5, res=300)

MeanGHGbyhabitat <- ggarrange(CO2_habitat_means, CH4_habitat_means, N2O_habitat_means,
                              labels = c("A", "B", "C"),
                              ncol = 1, nrow = 3, align="v")
MeanGHGbyhabitat

dev.off()

#### Figure 3. Flowing GHGs by drying regime####
tiff("FlowingbyDryingRegime", units="in", width=3, height=8, res=300)

FlowingbyDryingRegime <- ggarrange(CO2_drying_regime, CH4_drying_regime, N2O_drying_regime,
                                   labels = c("A", "B", "C"),
                                   ncol = 1, nrow = 3)
FlowingbyDryingRegime

dev.off()

#### Figure 4. Combine riparian and dry sediment forest plots ####
tiff("rip_dry_combined", units="in", width=6, height=6, res=300)

rip_dry_combined <- ggarrange(rf_dry + rremove("xlab") , rf_rip,
                              ncol = 1, nrow = 2, align="v", common.legend = T,legend="bottom")
rip_dry_combined


dev.off()




#### Figure 5. Relationship between flowing GHGs and top predictors ####
tiff("flow_drivers", units="in", width=9, height=11, res=300)

flow_drivers <- ggarrange(CO2_embed, CO2freqflow, CO2_water_temp, 
                          CH4temp, CH4mtosource, CH4OMint,  N2O_water_temp,
                          CO2_flowp_OM, CO2_percentIR, CH4OMper,
                          N2O_veloc, N2O_DN,  
                          labels = c("A", "B", "C","D", "E", "F", "G", "H", "I", "J", "K", "L"),
                          ncol = 3, nrow = 4, align="hv")
flow_drivers

dev.off()

#### Figure S2. Time series by each gas ####
tiff("GHG_timeseries", units="in", width=8, height=11, res=300)
GHG_timeseries <- ggarrange(CO2meansbox_camp, CH4meansbox_camp, N2Omeansbox_camp,
                            labels = c("A", "B", "C"),
                            ncol = 1, nrow = 3, align="v", common.legend = T,legend="bottom")
GHG_timeseries
dev.off()


#### Figure S4. Spatial distribution of OM stock and OM content ####

SpatialsedOM <- ggarrange(dry_sedOM, dry_OM_stock ,
                          labels = c("A", "C"),
                          ncol = 1, nrow = 2)
SpatialsedOM

tiff("SpatialsOMdry", units="in", width=4, height=8, res=300, type="cairo")
annotate_figure(SpatialsedOM, top = text_grob("Dry conditions", color = "black", face = "bold", size = 14))
dev.off()



SpatialOMstock <- ggarrange(flowing_OM, flowing_OM_stock,
                            labels = c("B", "D"),
                            ncol = 1, nrow = 2)
SpatialOMstock
tiff("SpatialsOMflow", units="in", width=4, height=8, res=300, type="cairo")
annotate_figure(SpatialOMstock, top = text_grob("Flowing conditions", color = "black", face = "bold", size = 14))
dev.off()

