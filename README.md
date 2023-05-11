# GHG_Albarine
R script associated with analysis and visualization of data from the publication "River network-scale drying impacts the spatio-temporal dynamics of greenhouse gas fluxes"

The necessary files can be downloaded from Zenodo:

dat_clean.csv
Includes all of the raw data. 

CO2means.csv; CH4means.csv and N2Omeans.csv
Include the site-averaged GHG data and a subset of environmental data, also averaged by site.

Model_coefficients_GHGdry.csv and Model_coefficients_GHGrip.csv
Include the model coefficients from the top models identifying the drivers of dry condition and riparian soil GHG fluxes, respectively 

Column title meanings: 
"ID_unique": unique identifier for each sampling point GHG measurement
"CO2_mg_m2_h": CO2 flux
"CH4_ug_m2_h": CH4 flux
"Air_temp": Data logger airtemperature at GHG sampling time in degrees celsius                    
"pressure_calc_atm": Atmospheric pressure at each site/sampling time in atm                 
"campaign": Sampling campaign number
"Site" : Unique ID for each study site                            
"drying_regime" : Flow regime (perennial, intermittent_long_, intermittent_moderate)
"soil_temp": Average soil temperature in degrees Celsius taken from 3 measurememtns at each riparian soil GHG sampling point                        
"soil_VWC": cc  
"flow_state" : flowing conditions (dry, pool, flowing)                      
"habitat_type" : for flowing conditions, GHG measurement in pool (more lentic conditions) or riffle (more lotic conditions)                      
"pool_riffle_depth_cm"  : water depth (cm) of flowing condition measurements           
"Picarro_LGR" : device used to measure GHG (Picarro: CO2, CH4 or LGR: N2O)                     
"total_volume_L" : headspace volume in litres                  
"chamber_area_m2" : GHG chamber area (m2)                  
"sed_temp" :   Average sediment temperature in degrees Celsius taken from 3 measurements at each dry riverbed sediment GHG sampling point                      
"sed_VWC" : Average sediment volumteric water content (m3/m3) taken from 3 measurements at each dry riverbed sediment GHG sampling point                           
"point" : Identifier for each GHG sampling point with each site                           
"habitat" : Habitat type: Riparian, Flowing, Dry, Pool                          
"ecosystem"  : Aquatic or riparian habitat                      
"CO2_g_m2_d"   : CO2 flux                     
"CH4_mg_m2_d"  : CH4 flux                    
"ug_N2ON_m2_h" : N2O flux                    
"mg_N2ON_m2_d"  : N2O flux                    
"Canopy_Cover_Benth" : In-stream canopy cover (%)                
"Stock_Benth_g.m2"  : In-stream organic matter stock (g/m2)               
"Max..wetted.width..m."  : Maximum stream wetted width (m)          
"Min..wetted.width..m."  : Minimum stream wetted width (m)          
"Ave.wet.width_m"  : Average stream wetted width (m)                 
"Bankfull.at.max..wetted.width..m." : Stream bankfull width at maximum wetted width (m)
"Bankfull.at.min..wetted.width..m." : Stream bankfull width at minimum wetted width (m)
"Embeddedness...." : stream substrate embeddedness (%)                
"Discharge..l.s."  : Stream discharge (L/s)                 
"Mean_Velocity..m.s."   : Mean stream velocity (m/s)           
"MeanDur_Ponded"   : Mean number of days with ponded conditions                  
"MaxDur_Ponded"  : Maximum number of days with ponded conditions                  
"MeanDur_Dry"   : Mean number of days with dry conditions                    
"MaxDur_Dry"   : Maximum number of days with dry conditions                    
"MeanDur_Flow"  : Mean number of days with flowing conditions                     
"MaxDur_Flow"  : Maximum number of days with flowing conditions                    
"masl"  : metres above sea level                            
"dist.to.per"  : Distance (m) to the nearest perennial reach                     
"length.IR.up"  : Length (m) of non-perennial reaches upstream                    
"length.PR.up"  : Length (m) of perennial reaches upstream                   
"percent.IR.up"  : Percentage of non-perennial reaches upstrean                   
"percent.PR.up"  : Percentage of perennial reaches upstream                  
"Soil_OM" : Soil organic matter content (%)                          
"Sed_OM"  : Sediment organic matter content (%)                         
"percent_N"  : Nitrogen content (%)                       
"percent_C"  : Carbon content (%)                      
"CN_ratio"  : Carbon to nitrogen ratio                        
"pH"  : Water pH                             
"DO_mg_L"  : Dissolved oxygen (mg/L)                         
"DO_percent"  : Dissolved oxygen (%)                     
"cond_us_cm"   : Conductivity (uS/cm)                     
"Temp_water_C"   : Water temperature (degrees C)                  
"date"  : Date                            
"totalbiomass_g"  : Total dry aboveground biomass from riparian collars (g)                   
"grassbiomass_g"  : Dry grass biomass from riparian collars (g)                  
"herbsbiomass_g"  : Dry herb biomass from riparian collars (g)                 
"leavesbiomass_g"  : Dry leaf biomass from riparian collars (g)                
"mossbiomass_g"   : Dry moss biomass from riparian collars (g)                 
"woodybiomass_g"  : Dry woody biomass from riparian collars (g)                  
"latitude"  : Site coordinates latitude                       
"longitude"  : Site coordinates longitude                       
"Freq_Flow"  : Flow permanence (%)                      
"Freq_Dry"   : Dry days (%)                       
"Freq_Ponded"  : Pooled days (%)                    
"m_to_source"  : Distance to source (m)                     
"freqflow7" : Flow permanence (%) within 7 days before sampling                       
"freqflow30" : Flow permanence (%) within 30 days before sampling                       
"freqdry7"  : Dry days (%) within 7 days before sampling                       
"freqdry30"   : Dry days (%) within 30 days before sampling                      
"freqpool7"    : Pooled days (%) within 7 days before sampling                    
"freqpool30"   : Pooled days (%) within 30 days before sampling                     
"freq_flow_2021"   : Flow permanence (%) in 2021                
"freq_dry_2021"  : Dry days (%) in 2021                   
"freq_pool_2021"  : Pooled days (%) in 2021                 
"Tsince_manual"  : Days since the last flowing or dry day                   
"sed_pH"  :  Sediment pH                         
"sed_cond_uS_cm" : Sediment condictivity (uS/cm)                   
"NO3_N_mg_l"  : Nitrate-nitrogen (mg/L)                 
"DN_ugL"  : Dissolved nitrogen (ug/L)                          
"DOC_ugL"  : Dissolved organic carbon (ug/L)                        
"SRP_ugL"   :   Soluble reactive phosphorus (ug/L)                    
"NH4_N_ugL"  : Ammonium-nitrogen (ug/L)                      
"bedrock"  : Bedrock stream substrate (%)                         
"boulders"  : Boulder stream substrate (%)                       
"cobbles"   : Cobble stream substrate (%)                        
"gravel"  : Gravel stream substrate (%)                         
"sand" : Sand stream substrate (%)                             
"silt"   : Silt stream substrate (%)                           
