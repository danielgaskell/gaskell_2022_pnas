library(dplyr)
library(tidyr)
library(broom)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(scales)
library(bayfoxr)
library(tidync)
library(sf)
library(velociraptr)
library(seacarb)

pH_mode = FALSE # calculate CO3 correction using pH record (slower)

# Note: d18O_sw throughout is given in VSMOW, d18O_calcite in VPDB. A correction of -0.27 is applied in some cases to convert between the two (Hut 1987), though note that this is not strictly the correct conversion between bulk d18O_VMOW and bulk d18O_VPDB, but an internal correction related to how the original data were measured.

# fit poly to bayfox, for faster evaluation over 10,000s of points (indistinguishable from the real thing to many more sig figs than we care about)
bayfox_SSTs <- tibble(d18O_bayfox=seq(-10, 10, 0.5))
sst <- predict_seatemp(bayfox_SSTs$d18O_bayfox, d18osw = 0, prior_mean = 30.0, prior_std = 20.0)
sst_quantiles <- quantile(sst, c(0.025,0.05,0.5,0.95,0.975))
bayfox_SSTs$SST_2.5  <- sst_quantiles[,1]
bayfox_SSTs$SST_5    <- sst_quantiles[,2]
bayfox_SSTs$SST_50   <- sst_quantiles[,3]
bayfox_SSTs$SST_95   <- sst_quantiles[,4]
bayfox_SSTs$SST_97.5 <- sst_quantiles[,5]
bayfox_SSTs$sd       <- (sst_quantiles[,2] - sst_quantiles[,3]) / qnorm(0.05)

bayfox_model_2.5  <- lm(SST_2.5  ~ poly(d18O_bayfox, 4, raw=T), data=bayfox_SSTs)
bayfox_model_5    <- lm(SST_5    ~ poly(d18O_bayfox, 4, raw=T), data=bayfox_SSTs)
bayfox_model_50   <- lm(SST_50   ~ poly(d18O_bayfox, 4, raw=T), data=bayfox_SSTs)
bayfox_model_95   <- lm(SST_95   ~ poly(d18O_bayfox, 4, raw=T), data=bayfox_SSTs)
bayfox_model_97.5 <- lm(SST_97.5 ~ poly(d18O_bayfox, 4, raw=T), data=bayfox_SSTs)
bayfox_sd <- mean(bayfox_SSTs$sd)

predict_poly <- function(model, data) {
    # Quick helper function to speed up predicting 4th-order polynomials (model) from a vector (data)
    # without the overhead of generating a new dataframe each time.
    model$coefficients[2]*data + model$coefficients[3]*(data**2) + model$coefficients[4]*(data**3) + model$coefficients[5]*(data**4) + model$coefficients[1]
}

# fit CO3 effect in pH space (accurate to within <0.001 per mil over the pH range 7.5-8.5)
# mean slope of O. universa, G. bulloides, T. sacculifer, and G. ruber from Spero et al. (1999) - all surface-dwelling planktics, range of no, weak, and strong symbiosis
# converted to pH space using seacarb with culture parameters from Spero et al. (1997)
# the intercept 0.505 produces 0 offset at current seawater [CO3] (200 umol/kg in the Zeebe & Tyrrell 2019 curve)

cie <- tibble(pH = seq(7.5, 8.5, 0.05)) %>%
    mutate(CO3 = carb(9, pH, 2032e-6, T=22)$CO3 * 1e6) %>%
    mutate(d18O = -0.002525*CO3 + 0.505)
ciefit <- lm(d18O ~ poly(pH, 3, raw=T), cie)

# load Miller et al. (2020) ice stack (standard not specified but equation 2 implies VSMOW - see eq. 1 in Cramer et al. 2011)
d18O_benthic <- read.csv("miller_2020_sealevel.csv", header=T, stringsAsFactors = F) %>%
    mutate(ref = "Miller et al. (2020)")

# load and append Cretaceous data
d18O_benthic_cret <- read.csv("benthic_stack_cretaceous.csv", header=T, stringsAsFactors = F)
d18O_benthic[which(d18O_benthic$age > 58.08571), "d18O_sw"] <- NA # wipe d18O_sw from before ~58 Ma, where Miller et al. (2020) Mg/Ca stack runs out - we assume constant sea level before this, rather than constant temperature (as Miller et al. do).
d18O_sw_earliest <- mean((d18O_benthic %>% select(age, d18O_sw) %>% filter(!is.na(d18O_sw)) %>% filter(age > max(age) - 5) %>% select(d18O_sw))$d18O_sw, na.rm=T)
d18O_benthic <- merge(d18O_benthic, d18O_benthic_cret %>% select(age, d18O, source, ref), all=T)
d18O_benthic[which(is.na(d18O_benthic$d18O_sw)), "d18O_sw"] <- d18O_sw_earliest # extend Miller et al. (2020) d18O_sw record into the Cretaceous using mean value of earliest 5 Ma

# smooth benthic stack to 0.25 Ma
d18O_benthic <- d18O_benthic %>%
    mutate(age = round(age * 4) / 4) %>%
    group_by(age) %>%
    summarize(d18O = mean(d18O), d18O_sw = mean(d18O_sw))
d18O_benthic[which(d18O_benthic$age == 0), "d18O_sw"] <- 0.14 # reset earliest age to 0 so we don't bias the modern points in the dataset

# load Zeebe & Tyrrell (2019) CO3 stack
CO3_zeebe <- read.csv("zeebe_tyrrell_2019.csv", header=T)

# load Foster et al. (2017) CO2 stack
pCO2_foster <- read.csv("foster_2017_loess.csv", header=T, stringsAsFactors = F)

# load clumped compilation
D47 <- read.csv("D47_compilation.csv", header=T) %>% mutate(latzone = ceiling(abs(pallat) / 30))

# get LeGrande & Schmidt (2006) mean d18O in upper 50m at 1x1 degree resolution
# this file from LeGrande and Schmidt (2006), via https://data.giss.nasa.gov/o18data/grid.html
legrande_d18Osw <- tidync("calculated_d18O_v1_1.nc") %>%
    hyper_filter(depth = depth <= 50) %>%
    hyper_tibble() %>%
    group_by(lon, lat) %>%
    summarize(d18O = mean(d18o))

# load d18O stack
d18O_stack <- read.csv("d18O_stack.csv", header=T, stringsAsFactors = F)
d18O_stack$age <- as.numeric(d18O_stack$age)
d18O_stack$d18O <- as.numeric(d18O_stack$d18O)
d18O_stack <- d18O_stack %>%
    mutate(age_01Ma = round(age, 1),
           age_1Ma = round(age),
           age_5Ma = round(age / 5) * 5) %>%
    mutate(latzone = ceiling(abs(pallat) / 30))
d18O_stack[which(d18O_stack$age <= 145), "epoch"] <- "Lower Cretaceous"
d18O_stack[which(d18O_stack$age <= 100.5), "epoch"] <- "Upper Cretaceous"
d18O_stack[which(d18O_stack$age <= 66), "epoch"] <- "Paleocene"
d18O_stack[which(d18O_stack$age <= 56), "epoch"] <- "Early Eocene"
d18O_stack[which(d18O_stack$age <= 48), "epoch"] <- "Eocene"
d18O_stack[which(d18O_stack$age <= 33.9), "epoch"] <- "Oligocene"
d18O_stack[which(d18O_stack$age <= 23.03), "epoch"] <- "Miocene"
d18O_stack[which(d18O_stack$age <= 5.333), "epoch"] <- "Pliocene"
d18O_stack[which(d18O_stack$age <= 2.588), "epoch"] <- "Quaternary"

# merge in benthic/ice/CO3 stacks
d18O_stack$d18O_benthic <- approx(d18O_benthic$age, d18O_benthic$d18O, d18O_stack$age)$y # merge in Miller et al. (2020) record with linear interpolation
d18O_stack$d18O_ice <- approx(d18O_benthic$age, d18O_benthic$d18O_sw, d18O_stack$age)$y # merge in Miller et al. (2020) record with linear interpolation
d18O_stack$CO3_zeebe <- approx(CO3_zeebe$age, CO3_zeebe$CO3, d18O_stack$age)$y # merge in Zeebe & Tyrrell (2019) CO3 record with linear interpolation
d18O_stack$pCO2_foster_0.05 <- approx(pCO2_foster$age, pCO2_foster$low_95,  d18O_stack$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_stack$pCO2_foster_0.32 <- approx(pCO2_foster$age, pCO2_foster$low_68,  d18O_stack$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_stack$pCO2_foster_0.50 <- approx(pCO2_foster$age, pCO2_foster$pCO2,    d18O_stack$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_stack$pCO2_foster_0.68 <- approx(pCO2_foster$age, pCO2_foster$high_68, d18O_stack$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_stack$pCO2_foster_0.95 <- approx(pCO2_foster$age, pCO2_foster$high_95, d18O_stack$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_benthic$CO3_zeebe <- approx(CO3_zeebe$age, CO3_zeebe$CO3, d18O_benthic$age)$y # merge in Zeebe & Tyrrell (2019) CO3 record with linear interpolation
d18O_benthic$pCO2_foster_0.05 <- approx(pCO2_foster$age, pCO2_foster$low_95,  d18O_benthic$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_benthic$pCO2_foster_0.32 <- approx(pCO2_foster$age, pCO2_foster$low_68,  d18O_benthic$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_benthic$pCO2_foster_0.50 <- approx(pCO2_foster$age, pCO2_foster$pCO2,    d18O_benthic$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_benthic$pCO2_foster_0.68 <- approx(pCO2_foster$age, pCO2_foster$high_68, d18O_benthic$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation
d18O_benthic$pCO2_foster_0.95 <- approx(pCO2_foster$age, pCO2_foster$high_95, d18O_benthic$age)$y # merge in Foster et al. (2017) pCO2 record with linear interpolation

# omit everything earlier than the benthic stack (for speed)
d18O_stack <- d18O_stack %>% filter(!is.na(d18O_ice))

# how many points do we have?
nrow(d18O_stack %>% filter((preservation == "E" & latzone == 1) | (preservation %in% c("E","VG","G","M","P") & latzone == 3 & pallat < 0)))

# raw data plot
ggplot(d18O_stack, aes(x=age, y=d18O, color=as.factor(abs(latzone)))) +
    geom_point(size=0.1, shape=3) +
    labs(x="Age (Ma)",
         y=bquote(delta^18*O*" (\211)")) +
    guides(color=F, fill=F) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,max(d18O_stack$age))) +
    scale_y_reverse(breaks = c(-6,-4,-2,0,2,4)) +
    scale_color_manual(values = c("orange2", "darkgreen", "dodgerblue4")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("supplement_rawforpres.pdf", width=6, height=4, dpi=300)

# correct for [CO3]
# mean slope of O. universa, G. bulloides, T. sacculifer, and G. ruber from Spero et al. (1999) - all surface-dwelling planktics, range of no, weak, and strong symbiosis
# the intercept 0.505 produces 0 offset at current seawater [CO3] (200 umol/kg in the Zeebe & Tyrrell 2019 curve)
d18O_stack$d18O_CO3 <- d18O_stack$d18O - (-0.002525*d18O_stack$CO3_zeebe + 0.505)
d18O_benthic$d18O_CO3 <- d18O_benthic$d18O - (-0.002525*d18O_benthic$CO3_zeebe + 0.505)

# get benthic temperature with bayfox
d18O_benthic$d18O_bayfox            <- d18O_benthic$d18O_CO3 - d18O_benthic$d18O_sw
d18O_benthic$temp_benthic_bayfox_5  <- predict(bayfox_model_5,  newdata=d18O_benthic)
d18O_benthic$temp_benthic_bayfox    <- predict(bayfox_model_50, newdata=d18O_benthic)
d18O_benthic$temp_benthic_bayfox_95 <- predict(bayfox_model_95, newdata=d18O_benthic)
d18O_stack$temp_benthic_bayfox      <- approx(d18O_benthic$age, d18O_benthic$temp_benthic_bayfox, d18O_stack$age)$y # merge benthic temperature into d18O_stack with linear interpolation

# get benthic temperature with Miller et al. (2020)'s equation
d18O_benthic$temp_benthic <- 16.1 - 4.76*(d18O_benthic$d18O_CO3 - (d18O_benthic$d18O_sw - 0.27))
d18O_stack$temp_benthic <- approx(d18O_benthic$age, d18O_benthic$temp_benthic, d18O_stack$age)$y # merge benthic temperature into d18O_stack with linear interpolation

# fit salinity-to-d18Osw relationship from Schmidt dataset
d18O_modern <- read.csv("schmidt_d18O.csv", header=T)
d18O_modern[which(d18O_modern$d18O == -99.9), "d18O"] <- NA
d18O_modern[which(d18O_modern$salinity == -99.9), "salinity"] <- NA
d18O_modern <- d18O_modern %>% filter(depth <= 50)
d18O_modern_fit <- na.omit(d18O_modern %>% select(lat, d18O, salinity) %>% filter(salinity > 30, d18O > -5, abs(lat) < max(na.omit(abs(d18O_stack$pallat))) + 5))
model_modern_sal <- lm(d18O ~ salinity, data=d18O_modern_fit)

# load salinity data
# note: salinity is calculated relative to reference salinity of 35, which means HadCM3 will
# generally NOT pick up secular variations in salinity - actually what we want, since we're
# just reconstructing the latitudinal gradient and global d18O_sw is theoretically independent
# of global salinity.
# source: these files from Valdes et al. (2021); a wget script is provided (download_valdes.bat)
#   to obtain these files in the format expected.
gcm <- tibble()
for (i in 0:(ceiling(max(d18O_stack$age))+5)) {
    if (file.exists(paste(".\\Valdes\\", i, ".nc", sep=""))) {
        print(paste("Loading ", i, ".nc...", sep=""))
        gcm_load <- tidync(paste(".\\Valdes\\", i, ".nc", sep="")) %>%
            activate("D0,D1,D5,D3") %>%
            hyper_tibble() %>%
            filter(depth_1 <= 50) %>%
            group_by(longitude, latitude) %>%
            summarize(salinity = mean(salinity_ym_dpth), temp = mean(insitu_T_ym_dpth)) %>%
            mutate(salinity = 35 + (salinity*1000)) %>% # reference salinity here is 35, which is what HadCM3 usually uses
            mutate(longitude = if(longitude <= 180) longitude else longitude - 360) %>%
            mutate(age = i)
        if (nrow(gcm) == 0) {
            gcm <- gcm_load
        } else {
            gcm <- merge(gcm, gcm_load, all=T)
        }
    }
}
gcm_temps <- gcm %>%
    mutate(latzone = ceiling(abs(latitude+1e-6) / 30)) %>%
    group_by(latzone, age) %>%
    summarize(temp = mean(temp))
gcm_ages <- unique(gcm$age)
gcm$d18O_predicted <- predict(model_modern_sal, newdata=gcm)

# merge GCM salinity + temp into d18O stack
# -----------------------------------------
# get closest GCM age
agecol = which(colnames(d18O_stack) == "age")[1]
d18O_stack$age_gcm <- apply(d18O_stack, 1, function(x) gcm_ages[which(abs(gcm_ages - as.numeric(x[agecol])) == min(abs(gcm_ages - as.numeric(x[agecol]))))][1]) # find closest GCM age for each row
agecol = which(colnames(d18O_benthic) == "age")[1]
d18O_benthic$age_gcm <- apply(d18O_benthic, 1, function(x) gcm_ages[which(abs(gcm_ages - as.numeric(x[agecol])) == min(abs(gcm_ages - as.numeric(x[agecol]))))][1]) # find closest GCM age for each row
agegcmcol = which(colnames(d18O_stack) == "age_gcm")[1]

# pull surface d18O_sw from Valdes model results + salinity relationship
d18O_sites <- na.omit(unique(d18O_stack %>% select(pallat, pallong, age_gcm)))
gcm_latlon <- na.omit(unique(gcm %>% select(latitude, longitude, age)))
for (i in 1:nrow(d18O_sites)) {
    d18O_sites[i,"d18O_sw_Valdes"] <- mean((gcm %>% filter(age == d18O_sites[i,"age_gcm"], latitude >= d18O_sites[i,"pallat"] - 5, latitude <= d18O_sites[i,"pallat"] + 5, longitude >= d18O_sites[i,"pallong"] - 5, longitude <= d18O_sites[i,"pallong"] + 5))$d18O_predicted, na.rm=T)
    d18O_sites[i,"temp_Valdes"] <- mean((gcm %>% filter(age == d18O_sites[i,"age_gcm"], latitude >= d18O_sites[i,"pallat"] - 5, latitude <= d18O_sites[i,"pallat"] + 5, longitude >= d18O_sites[i,"pallong"] - 5, longitude <= d18O_sites[i,"pallong"] + 5))$temp, na.rm=T)
    d18O_sites[i,"temp_Valdes_highlat"] <- mean((gcm %>% filter(age == d18O_sites[i,"age_gcm"], abs(latitude) > 60))$temp, na.rm=T)
}
# merge in salinity/temp data from GCM
d18O_stack <- d18O_stack %>% select(-c(d18O_sw_Valdes))
d18O_stack <- d18O_stack %>% select(-c(temp_Valdes))
d18O_stack <- d18O_stack %>% select(-c(temp_Valdes_highlat))
d18O_stack <- merge(d18O_stack, d18O_sites, all.x=T)

# pull surface d18O_sw from CESM model results / modern
prep_cesm <- function(filename) {
    netcdf <- tidync(filename) %>%
        hyper_tibble() %>%
        full_join(tidync(paste(filename, ".gridarea.nc", sep="")) %>% hyper_tibble())
    lowlat_area <- sum((netcdf %>% filter(abs(lat) <= 30 & !is.na(TEMP)))$cell_area, na.rm=T)
    highlat_area <- sum((netcdf %>% filter(abs(lat) >= 60 & !is.na(TEMP)))$cell_area, na.rm=T)
    southern_area <- sum((netcdf %>% filter(lat <= -60 & !is.na(TEMP)))$cell_area, na.rm=T)
    netcdf %>%
        mutate(contrib_low = cell_area / lowlat_area,
               contrib_high = cell_area / highlat_area,
               contrib_southern = cell_area / southern_area,
               lon = ifelse(lon < 180, lon, lon - 360),
               d18O = (R18O - 1.0) * 1000.0) # per Zhu et al. (2020)
}
# these files are CESM runs from Zhu et al. (2020) and this paper. Original NetCDFs may be obtained from
# https://zenodo.org/record/3334893 and from the corresponding author of the paper.
# note: for processing speed these are remapped to r1080x540 geometry using CDO. Gridareas are extracted
# separately using CDO to yield the .gridarea.nc files expected by prep_cesm().
cesm_eocene_1x <- prep_cesm("..\\..\\..\\CESM\\b.e12.B1850C5CN.f19_g16.iPETM01x.02.pop.h.TEMP_SALT_R18O.2701-2800.climo.remapped.nc")
cesm_eocene_3x <- prep_cesm("..\\..\\..\\CESM\\b.e12.B1850C5CN.f19_g16.iPETM03x.03.pop.h.TEMP_SALT_R18O.2101-2200.climo.remapped.nc")
cesm_eocene_6x <- prep_cesm("..\\..\\..\\CESM\\b.e12.B1850C5CN.f19_g16.iPETM06x.09.pop.h.TEMP_SALT_R18O.2101-2200.climo.remapped.nc")
cesm_eocene_9x <- prep_cesm("..\\..\\..\\CESM\\b.e12.B1850C5CN.f19_g16.iPETM09x.02.pop.h.TEMP_SALT_R18O.2101-2200.climo.remapped.nc")
cesm_miocene_280 <- prep_cesm("..\\..\\..\\CESM\\surface_props_B.MMIOx2_C5_280_WISOon.pop.h.ANN_concat.remapped.nc")
cesm_miocene_400 <- prep_cesm("..\\..\\..\\CESM\\surface_props_B.MMIOx2_C5_400_WISOon.pop.h.ANN_concat.remapped.nc")
cesm_miocene_840 <- prep_cesm("..\\..\\..\\CESM\\surface_props_B.MMIOx2_C5_840_WISOon.pop.h.ANN_concat.remapped.nc")

cesm_bwt            <- function(cesm_run)  sum((cesm_run %>% filter(abs(lat) >= 60))$TEMP * (cesm_run %>% filter(abs(lat) >= 60))$contrib_high, na.rm=T)
cesm_bandtemp_low   <- function(cesm_run)  sum((cesm_run %>% filter(abs(lat) <= 30))$TEMP * (cesm_run %>% filter(abs(lat) <= 30))$contrib_low, na.rm=T)
cesm_bandtemp_high  <- function(cesm_run)  sum((cesm_run %>% filter(lat <= -60))$TEMP * (cesm_run %>% filter(lat <= -60))$contrib_southern, na.rm=T)
benthic_bwt         <- function(age_Ma)    mean(unlist(d18O_benthic[which(d18O_benthic$age > age_Ma - 0.5 & d18O_benthic$age < age_Ma + 0.5),"temp_benthic"]), na.rm=T)

legrande_d18Osw_bwt   <- benthic_bwt(0)
cesm_miocene_280_mean <- 0.68
cesm_miocene_400_mean <- 0.68
cesm_miocene_840_mean <- 0.68
cesm_miocene_280_bwt  <- cesm_bwt(cesm_miocene_280)
cesm_miocene_400_bwt  <- cesm_bwt(cesm_miocene_400)
cesm_miocene_840_bwt  <- cesm_bwt(cesm_miocene_840)
cesm_eocene_1x_mean   <- -1 # initialization value from Zhu et al. (2020)
cesm_eocene_3x_mean   <- -1 # initialization value from Zhu et al. (2020)
cesm_eocene_6x_mean   <- -1 # initialization value from Zhu et al. (2020)
cesm_eocene_9x_mean   <- -1 # initialization value from Zhu et al. (2020)
cesm_eocene_1x_bwt    <- cesm_bwt(cesm_eocene_1x)
cesm_eocene_3x_bwt    <- cesm_bwt(cesm_eocene_3x)
cesm_eocene_6x_bwt    <- cesm_bwt(cesm_eocene_6x)
cesm_eocene_9x_bwt    <- cesm_bwt(cesm_eocene_9x)
band_temps_modern <- tibble(bwt = legrande_d18Osw_bwt,
                            low_temp = 26.19064,
                            high_temp = -0.1950831)
band_temps_miocene <- tibble(bwt = c(cesm_miocene_280_bwt,
                                     cesm_miocene_400_bwt,
                                     cesm_miocene_840_bwt),
                             low_temp = c(cesm_bandtemp_low(cesm_miocene_280),
                                          cesm_bandtemp_low(cesm_miocene_400),
                                          cesm_bandtemp_low(cesm_miocene_840)),
                             high_temp = c(cesm_bandtemp_high(cesm_miocene_280),
                                           cesm_bandtemp_high(cesm_miocene_400),
                                           cesm_bandtemp_high(cesm_miocene_840)))
band_temps_eocene <- tibble(bwt = c(cesm_eocene_1x_bwt,
                                    cesm_eocene_3x_bwt,
                                    cesm_eocene_6x_bwt,
                                    cesm_eocene_9x_bwt),
                             low_temp = c(cesm_bandtemp_low(cesm_eocene_1x),
                                          cesm_bandtemp_low(cesm_eocene_3x),
                                          cesm_bandtemp_low(cesm_eocene_6x),
                                          cesm_bandtemp_low(cesm_eocene_9x)),
                             high_temp = c(cesm_bandtemp_high(cesm_eocene_1x),
                                           cesm_bandtemp_high(cesm_eocene_3x),
                                           cesm_bandtemp_high(cesm_eocene_6x),
                                           cesm_bandtemp_high(cesm_eocene_9x)))

d18O_sites <- na.omit(unique(d18O_stack %>% filter(depth == 1) %>% select(pallat, pallong, epoch)))
d18O_sites_full <- na.omit(unique(d18O_stack %>% filter(depth == 1) %>% select(pallat, pallong, epoch, temp_benthic_bayfox)))
all_sites_d18O_sw <- tibble()
all_sites_splines <- tibble()
supp_d18Osw <- tibble()
bwt_steps <- d18O_benthic %>% select(age, temp_benthic) %>% mutate(age_5Ma = round(age / 5) * 5) %>% group_by(age_5Ma) %>% summarize(temp_benthic = mean(temp_benthic))
for (i in 1:nrow(d18O_sites)) {
    print(paste("Getting CESM results for site", i, "of", nrow(d18O_sites)))
    
    pallat <- d18O_sites[i,"pallat"]
    pallong <- d18O_sites[i,"pallong"]
    
    # interpolate between model runs using local splines and BWT
    if (d18O_sites[i,"epoch"] %in% c("Quaternary", "Pliocene")) {
        bwt_steps_site <- bwt_steps %>% filter(age_5Ma < 5.333)
        site_legrande <- legrande_d18Osw %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_d18O_sw <- tibble(bwt = legrande_d18Osw_bwt,
                               d18O_sw = mean(site_legrande$d18O, na.rm=T),# - legrande_d18Osw_mean,
                               d18O_sw_sd = sd(site_legrande$d18O, na.rm=T),
                               temp = -999) # ignore temperature for lack of data
        band_temp <- band_temps_modern
    } else if (d18O_sites[i,"epoch"] %in% c("Miocene", "Oligocene")) {
        bwt_steps_site <- bwt_steps %>% filter(age_5Ma >= 5.333, age_5Ma < 33.9)
        site_miocene_280 <- cesm_miocene_280 %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_miocene_400 <- cesm_miocene_400 %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_miocene_840 <- cesm_miocene_840 %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_d18O_sw <- tibble(bwt = c(cesm_miocene_280_bwt,
                                       cesm_miocene_400_bwt,
                                       cesm_miocene_840_bwt),
                               d18O_sw = c(mean(site_miocene_280$d18O, na.rm=T) - cesm_miocene_280_mean,
                                           mean(site_miocene_400$d18O, na.rm=T) - cesm_miocene_400_mean,
                                           mean(site_miocene_840$d18O, na.rm=T) - cesm_miocene_840_mean),
                               d18O_sw_sd = c(sd(site_miocene_280$d18O, na.rm=T),
                                              sd(site_miocene_400$d18O, na.rm=T),
                                              sd(site_miocene_840$d18O, na.rm=T)),
                               temp = c(mean(site_miocene_280$TEMP, na.rm=T),
                                        mean(site_miocene_400$TEMP, na.rm=T),
                                        mean(site_miocene_840$TEMP, na.rm=T)))
        band_temp <- band_temps_miocene
    } else if (d18O_sites[i,"epoch"] %in% c("Eocene", "Early Eocene", "Paleocene", "Upper Cretaceous")) {
        bwt_steps_site <- bwt_steps %>% filter(age_5Ma >= 33.9)
        site_eocene_1x <- cesm_eocene_1x %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_eocene_3x <- cesm_eocene_3x %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_eocene_6x <- cesm_eocene_6x %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_eocene_9x <- cesm_eocene_9x %>% filter(lat >= pallat - 5, lat <= pallat + 5, lon >= pallong - 5, lon <= pallong + 5)
        site_d18O_sw <- tibble(bwt = c(cesm_eocene_1x_bwt,
                                       cesm_eocene_3x_bwt,
                                       cesm_eocene_6x_bwt,
                                       cesm_eocene_9x_bwt),
                               d18O_sw = c(mean(site_eocene_1x$d18O, na.rm=T) - cesm_eocene_1x_mean,
                                           mean(site_eocene_3x$d18O, na.rm=T) - cesm_eocene_3x_mean,
                                           mean(site_eocene_6x$d18O, na.rm=T) - cesm_eocene_6x_mean,
                                           mean(site_eocene_9x$d18O, na.rm=T) - cesm_eocene_9x_mean),
                               d18O_sw_sd = c(sd(site_eocene_1x$d18O, na.rm=T),
                                              sd(site_eocene_3x$d18O, na.rm=T),
                                              sd(site_eocene_6x$d18O, na.rm=T),
                                              sd(site_eocene_9x$d18O, na.rm=T)),
                               temp = c(mean(site_eocene_1x$TEMP, na.rm=T),
                                        mean(site_eocene_3x$TEMP, na.rm=T),
                                        mean(site_eocene_6x$TEMP, na.rm=T),
                                        mean(site_eocene_9x$TEMP, na.rm=T)))
        band_temp <- band_temps_eocene
    }
    
    # natural spline interpolation using benthic temperature for BWT
    # merge splines into all_sites_d18O_sw
    site_d18O_sw <- site_d18O_sw %>% mutate(lat = pallat, lon = pallong, epoch = d18O_sites[i,"epoch"])
    spline_seq <- seq(0, max(site_d18O_sw$bwt + 3), 0.5)
    site_splines <- tibble(bwt = spline_seq,
                           d18O_sw = spline(site_d18O_sw$bwt, site_d18O_sw$d18O_sw, xout=spline_seq, method="natural")$y,
                           lat = pallat, lon = pallong, epoch = d18O_sites[i,"epoch"])
    if (site_splines$d18O_sw[1] > 2.5) { # exclude sites with extreme splines
        exclude_site <- TRUE
    } else {
        exclude_site <- FALSE
    }
    site_splines$exclude <- exclude_site
    if (!nrow(all_sites_d18O_sw)) {
        all_sites_d18O_sw <- site_d18O_sw
        all_sites_splines <- site_splines
    } else {
        all_sites_d18O_sw <- all_sites_d18O_sw %>% full_join(site_d18O_sw)
        all_sites_splines <- all_sites_splines %>% full_join(site_splines)
    }

    # save interpolated values to site stack (if not excluded)
    # NOTE: using temp_benthic_bayfox as comparison point for high-lat SST rather than Miller et al. (2020) equation - more apples-to-apples
    if (nrow(na.omit(site_d18O_sw)) & exclude_site == FALSE) {
        indices <- which(d18O_sites_full$pallat == pallat & d18O_sites_full$pallong == pallong & d18O_sites_full$epoch == d18O_sites[i,"epoch"])
        d18O_sites_full[indices,"d18O_sw_CESM"]        <- spline(site_d18O_sw$bwt, site_d18O_sw$d18O_sw,    xout=d18O_sites_full[indices,"temp_benthic_bayfox"], method="natural")$y
        d18O_sites_full[indices,"d18O_sw_CESM_sd"]     <- spline(site_d18O_sw$bwt, site_d18O_sw$d18O_sw_sd, xout=d18O_sites_full[indices,"temp_benthic_bayfox"], method="natural")$y
        d18O_sites_full[indices,"temp_CESM"]           <- spline(site_d18O_sw$bwt, site_d18O_sw$temp,       xout=d18O_sites_full[indices,"temp_benthic_bayfox"], method="natural")$y
        d18O_sites_full[indices,"band_temp_low_CESM"]  <- spline(band_temp$bwt,    band_temp$low_temp,      xout=d18O_sites_full[indices,"temp_benthic_bayfox"], method="natural")$y
        d18O_sites_full[indices,"band_temp_high_CESM"] <- spline(band_temp$bwt,    band_temp$high_temp,     xout=d18O_sites_full[indices,"temp_benthic_bayfox"], method="natural")$y
        supp_d18Osw_site <- tibble(pallat = pallat, pallong = pallong, age = bwt_steps_site$age_5Ma, d18O_sw = spline(site_d18O_sw$bwt, site_d18O_sw$d18O_sw, xout=bwt_steps_site$temp_benthic, method="natural")$y)
        if (nrow(supp_d18Osw)) {
            supp_d18Osw <- supp_d18Osw %>% full_join(supp_d18Osw_site)
        } else {
            supp_d18Osw <- supp_d18Osw_site
        }
    }
}
# merge in d18O data from CESM runs / modern
d18O_stack <- d18O_stack %>% select(-c(d18O_sw_CESM))
d18O_stack <- d18O_stack %>% select(-c(d18O_sw_CESM_sd))
d18O_stack <- d18O_stack %>% select(-c(temp_CESM))
d18O_stack <- d18O_stack %>% select(-c(band_temp_low_CESM))
d18O_stack <- d18O_stack %>% select(-c(band_temp_high_CESM))
d18O_stack <- merge(d18O_stack, d18O_sites_full, all.x=T, by=c("pallat", "pallong", "epoch", "temp_benthic_bayfox"))

# save supplementary d18Osw table
write.csv(na.omit(unique(supp_d18Osw)), "supplement_d18Osw.csv")

# plot supplementary spatial-bias figure
band_temp_residuals <- d18O_stack %>%
    filter(pallat <= -60, temp_CESM > -900, preservation %in% c("E","VG","G","M","P")) %>%
    group_by(pallat, pallong, age) %>%
    summarize(temp_CESM = median(temp_CESM, na.rm=T),
              band_temp_CESM = median(band_temp_high_CESM, na.rm=T)) %>%
    mutate(latzone = 3) %>%
    full_join(d18O_stack %>%
                  filter(abs(pallat) <= 30, temp_CESM > -900, preservation %in% c("E")) %>%
                  group_by(pallat, pallong, age) %>%
                  summarize(temp_CESM = median(temp_CESM, na.rm=T),
                            band_temp_CESM = median(band_temp_low_CESM, na.rm=T)) %>%
                  mutate(latzone = 1)) %>%
    mutate(residual = temp_CESM - band_temp_CESM)
ggplot(band_temp_residuals, aes(x=age, y=residual)) +
    facet_grid(rows=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`3`="High-Latitude"))) +
    geom_point() +
    labs(x = "Age (Ma)",
         y = "Residuals (°C)") +
    theme_bw()
ggsave("supplement_cesmresiduals.png", width=6, height=4, dpi=150)

# show spatial-bias statistics
cells <- which(d18O_stack$pallat <= -60 & d18O_stack$preservation %in% c("E","VG","G","M","P"))
d18O_stack[cells,"band_temp_residual"] <- d18O_stack[cells,"temp_CESM"] - d18O_stack[cells,"band_temp_high_CESM"]
mean(abs(d18O_stack[cells,"band_temp_residual"]), na.rm=T)
sd(abs(d18O_stack[cells,"band_temp_residual"]), na.rm=T)

cells <- which(abs(d18O_stack$pallat) <= 30 & d18O_stack$preservation %in% c("E"))
d18O_stack[cells,"band_temp_residual"] <- d18O_stack[cells,"temp_CESM"] - d18O_stack[cells,"band_temp_low_CESM"]
mean(abs(d18O_stack[cells,"band_temp_residual"]), na.rm=T)
sd(abs(d18O_stack[cells,"band_temp_residual"]), na.rm=T)

# what is the relative variability of d18O_sw in the Southern Ocean vs the Tropical Pacific?
sd(unlist(cesm_eocene_1x[which(cesm_eocene_1x$lat < -60), "d18O"])) / sd(unlist(cesm_eocene_1x[which(abs(cesm_eocene_1x$lat) <= 30 & (cesm_eocene_1x$lon > 90 | cesm_eocene_1x$lon < -90)), "d18O"]))
sd(unlist(cesm_eocene_3x[which(cesm_eocene_3x$lat < -60), "d18O"])) / sd(unlist(cesm_eocene_3x[which(abs(cesm_eocene_3x$lat) <= 30 & (cesm_eocene_3x$lon > 90 | cesm_eocene_3x$lon < -90)), "d18O"]))
sd(unlist(cesm_eocene_6x[which(cesm_eocene_6x$lat < -60), "d18O"])) / sd(unlist(cesm_eocene_6x[which(abs(cesm_eocene_6x$lat) <= 30 & (cesm_eocene_6x$lon > 90 | cesm_eocene_6x$lon < -90)), "d18O"]))
sd(unlist(cesm_eocene_9x[which(cesm_eocene_9x$lat < -60), "d18O"])) / sd(unlist(cesm_eocene_9x[which(abs(cesm_eocene_9x$lat) <= 30 & (cesm_eocene_9x$lon > 90 | cesm_eocene_9x$lon < -90)), "d18O"]))

# plot supplementary splines figure
all_sites_d18O_sw <- all_sites_d18O_sw %>% mutate(latzone = ceiling(abs(lat) / 30), epoch = ifelse(epoch %in% c("Quaternary","Pliocene"), "Modern", ifelse(epoch %in% c("Miocene","Oligocene"), "Miocene CESM", "Eocene CESM")))
all_sites_splines <- all_sites_splines %>% mutate(latzone = ceiling(abs(lat) / 30), epoch = ifelse(epoch %in% c("Quaternary","Pliocene"), "Modern", ifelse(epoch %in% c("Miocene","Oligocene"), "Miocene CESM", "Eocene CESM")))
ggplot() +
    facet_grid(rows=vars(epoch), cols=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`2`="Mid-Latitude",`3`="High-Latitude"))) +
    geom_line(data=all_sites_splines, aes(x=bwt, y=d18O_sw, group=lon, color=exclude)) +
    geom_point(data=all_sites_d18O_sw, aes(x=bwt, y=d18O_sw, group=lon), size=0.5) +
    scale_color_manual(values = c("gray50", "red")) +
    labs(x="Southern High-Latitude Temperature (°C)",
         y=bquote("Mean "*delta^18*O['sw']*" at site (10x10°)")) +
    guides(color=F) +
    theme_bw()
ggsave("cesm_splines.png", width=6, height=6, dpi=150)

# run general splines estimate from all data
get_spline_points <- function(dataset, d18O_mean, bwt, modelname) {
    dataset %>%
        #mutate(lat = round((lat+5)/10)*10-5, lon = round((lon+5)/10)*10-5) %>%
        #mutate(lat = round((lat+5)/10)*10-5, lon = round(lon/20)*20) %>%
        mutate(lat = round((lat+2.5)/5)*5-2.5, lon = round(lon/5)*5) %>%
        group_by(lat, lon) %>%
        summarize(d18O = mean(d18O) - d18O_mean) %>%
        mutate(bwt = bwt, model = modelname)
}
spline_points <- get_spline_points(legrande_d18Osw, 0, legrande_d18Osw_bwt, "Modern") %>%
    full_join(get_spline_points(cesm_miocene_280, cesm_miocene_280_mean, cesm_miocene_280_bwt, "CESM Miocene")) %>%
    full_join(get_spline_points(cesm_miocene_400, cesm_miocene_400_mean, cesm_miocene_400_bwt, "CESM Miocene")) %>%
    full_join(get_spline_points(cesm_miocene_840, cesm_miocene_840_mean, cesm_miocene_840_bwt, "CESM Miocene")) %>%
    full_join(get_spline_points(cesm_eocene_1x, cesm_eocene_1x_mean, cesm_eocene_1x_bwt, "CESM Eocene")) %>%
    full_join(get_spline_points(cesm_eocene_3x, cesm_eocene_3x_mean, cesm_eocene_3x_bwt, "CESM Eocene")) %>%
    full_join(get_spline_points(cesm_eocene_6x, cesm_eocene_6x_mean, cesm_eocene_6x_bwt, "CESM Eocene")) %>%
    full_join(get_spline_points(cesm_eocene_9x, cesm_eocene_9x_mean, cesm_eocene_9x_bwt, "CESM Eocene"))
spline_points_filtered <- spline_points %>%
    mutate(latzone = ceiling(abs(lat+0.001) / 30) * ifelse(abs(lat) > 30, -sign(lat), 1)) %>%
    filter(latzone > 0)
lowlat_cesm_poly  <- lm(d18O ~ poly(bwt, 3), data=spline_points_filtered %>% filter(latzone == 1))
midlat_cesm_poly  <- lm(d18O ~ poly(bwt, 3), data=spline_points_filtered %>% filter(latzone == 2))
highlat_cesm_poly <- lm(d18O ~ poly(bwt, 3), data=spline_points_filtered %>% filter(latzone == 3))
ggplot(spline_points_filtered %>% filter(latzone == 3), aes(x=bwt, y=d18O, group=latzone, color=lat)) +
    geom_smooth(method="lm", formula=y~poly(x,3)) +
    geom_point() +
    theme_bw()
general_poly <- lm(d18O ~ poly(bwt, 2, raw=T) * poly(lat, 2, raw=T), data=spline_points %>% filter(lat < 30))
summary(general_poly)
spline_points$d18O_predicted <- 0.0378*spline_points$bwt - 0.00196*spline_points$bwt^2 + 0.000751*spline_points$lat - 0.000142*spline_points$lat^2 - 0.000387*spline_points$bwt*spline_points$lat + 3.05e-05*spline_points$bwt^2*spline_points$lat - 1.57e-05*spline_points$bwt*spline_points$lat^2 + 2.02e-07*spline_points$bwt^2*spline_points$lat^2 + 0.427
spline_points$d18O_residuals <- spline_points$d18O - spline_points$d18O_predicted
ggplot(spline_points, aes(x=lat, y=d18O_residuals, color=model)) +
    geom_point(alpha=0.25,shape=1) +
    labs(title=bquote("Polynomial "*delta^18*O['sw']*" residuals"),
         x="Decimal Latitude",
         y=bquote(delta^18*O['sw']*" residuals"),
         color="Source") +
    scale_x_continuous(breaks=c(-90,-60,-30,0,30,30,60,90)) +
    scale_y_continuous(breaks=c(-12,-10,-8,-6,-4,-2,0,2,4)) +
    theme_bw() +
    theme(legend.position = "bottom")
ggsave("supplement_polyresiduals.png", width=4.5, height=4, dpi=150)

print(paste(format(general_poly$coefficients[2], digits=3), "T + ",
            format(general_poly$coefficients[3], digits=3), "T^2 + ",
            format(general_poly$coefficients[4], digits=3), "L + ",
            format(general_poly$coefficients[5], digits=3), "L^2 + ",
            format(general_poly$coefficients[6], digits=3), "TL + ",
            format(general_poly$coefficients[7], digits=3), "T^2*L + ",
            format(general_poly$coefficients[8], digits=3), "T*L^2 + ",
            format(general_poly$coefficients[9], digits=3), "T^2*L^2 + ",
            format(general_poly$coefficients[1], digits=3), sep=""))
sd(unlist(spline_points[which(spline_points$lat < 30),'d18O_residuals']), na.rm=T)
sd(unlist(spline_points[which(spline_points$lat > 30),'d18O_residuals']), na.rm=T)

# apply model fits for d18O_sw - CESM version
d18O_stack$d18O_sw <- d18O_stack$d18O_ice + d18O_stack$d18O_sw_CESM
d18O_stack[which(d18O_stack$d18O_sw < -5),"d18O_sw"] <- NA # throw out extreme values (NOTE THIS!)

# SST conversion: Bemis et al. (1998) (O. universa eqs. 1 and 2 suggested by paper for SST)
d18O_stack$SST_Bemis_LL <- 16.5 - 4.80*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_sw - 0.27))
d18O_stack$SST_Bemis_HL <- 14.9 - 4.80*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_sw - 0.27))
d18O_stack$SST_Bemis_mean <- (d18O_stack$SST_Bemis_HL + d18O_stack$SST_Bemis_LL) / 2
d18O_stack$SST_Bemis_no_spatial_LL <- 16.5 - 4.80*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_ice - 0.27))
d18O_stack$SST_Bemis_no_spatial_HL <- 14.9 - 4.80*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_ice - 0.27))
d18O_stack$SST_Bemis_no_spatial_mean <- (d18O_stack$SST_Bemis_no_spatial_HL + d18O_stack$SST_Bemis_no_spatial_LL) / 2
d18O_stack$SST_Bemis_Valdes_LL <- 16.5 - 4.80*(d18O_stack$d18O_CO3 - ((d18O_stack$d18O_ice + d18O_stack$d18O_sw_Valdes) - 0.27))
d18O_stack$SST_Bemis_Valdes_HL <- 14.9 - 4.80*(d18O_stack$d18O_CO3 - ((d18O_stack$d18O_ice + d18O_stack$d18O_sw_Valdes) - 0.27))
d18O_stack$SST_Bemis_Valdes_mean <- (d18O_stack$SST_Bemis_Valdes_HL + d18O_stack$SST_Bemis_Valdes_LL) / 2

# SST conversion: Kim & O'Neil (1997) inorganic - this formulation from Barras et al. (2010)
d18O_stack$SST_inorg <- 16.1 - 4.64*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_sw - 0.27)) + 0.09*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_sw - 0.27))^2
d18O_stack$SST_inorg_no_spatial <- 16.1 - 4.64*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_ice - 0.27)) + 0.09*(d18O_stack$d18O_CO3 - (d18O_stack$d18O_ice - 0.27))^2
d18O_stack$SST_inorg_Valdes <- 16.1 - 4.64*(d18O_stack$d18O_CO3 - ((d18O_stack$d18O_ice + d18O_stack$d18O_sw_Valdes) - 0.27)) + 0.09*(d18O_stack$d18O_CO3 - ((d18O_stack$d18O_ice + d18O_stack$d18O_sw_Valdes) - 0.27))^2

# Monte Carlo bayfox SST conversion with error margins
mc_samples <- 500
# assign a unique ID to each record
d18O_stack <- d18O_stack %>%
    group_by(pallat, pallong, ref) %>%
    mutate(record_id = cur_group_id()) %>%
    ungroup()
# create randomized offset lists for each record
sample_offsets <- d18O_stack %>%
    group_by(record_id) %>%
    summarize(d18O_sw_CESM_sd_mean = max(mean(d18O_sw_CESM_sd, na.rm=T), 0)) %>%
    select(record_id, d18O_sw_CESM_sd_mean)
# replace NA values for CESM sd with mean sd
na_replace_mean <- mean(sample_offsets$d18O_sw_CESM_sd_mean, na.rm=T)
sample_offsets[which(is.na(sample_offsets$d18O_sw_CESM_sd_mean)),"d18O_sw_CESM_sd_mean"] <- na_replace_mean
# populate randomized offset lists
sample_offsets <- sample_offsets %>%
    mutate(sample_d18O_sw  = list(rnorm(mc_samples, 0, d18O_sw_CESM_sd_mean)),
           sample_d18O_ice = list(rnorm(mc_samples, 0, 0.13)), # Miller et al. (2020) provides various error estimates - this is equivalent to their SD of +-10m sea level as shown in Fig. 1.
           sample_CO3      = list(rnorm(mc_samples, 0, 14))) # Zeebe & Tyrrell (2019) do not clearly list their errors for each timestep, but per 5.3.3 the error increases with age to a maximum SD of ~9 at 100 Ma; +5 umol/kg to add ~+/- 0.1 C error due to uncertainties in latitudinal pH gradients.
# populate pCO2 offsets based on mean & sd of pCO2 record during each record
optim_func <- function(x, vals, quants) sum(abs(vals - qbeta(quants, x[1], x[2])*x[3])^2)
sample_offsets <- left_join(sample_offsets, d18O_stack %>%
    group_by(record_id) %>%
    summarize(pCO2_foster_0.05_mean = mean(pCO2_foster_0.05),
              pCO2_foster_0.32_mean = mean(pCO2_foster_0.32),
              pCO2_foster_0.50_mean = mean(pCO2_foster_0.50),
              pCO2_foster_0.68_mean = mean(pCO2_foster_0.68),
              pCO2_foster_0.95_mean = mean(pCO2_foster_0.95)) %>%
    rowwise() %>%
    mutate(pCO2_params = list(optim(c(7, 3, 350), optim_func, gr = NULL,
                                    c(pCO2_foster_0.05_mean, pCO2_foster_0.32_mean, pCO2_foster_0.50_mean, pCO2_foster_0.68_mean, pCO2_foster_0.95_mean),
                                    c(0.05, 0.32, 0.50, 0.68, 0.95))$par)) %>%
    mutate(sample_pCO2 = list(rbeta(mc_samples, unlist(pCO2_params)[1], unlist(pCO2_params)[2])*unlist(pCO2_params)[3])) %>%
    mutate(pCO2_foster_0.05_fit = quantile(unlist(sample_pCO2), c(0.05)),
           pCO2_foster_0.32_fit = quantile(unlist(sample_pCO2), c(0.32)),
           pCO2_foster_0.50_fit = quantile(unlist(sample_pCO2), c(0.50)),
           pCO2_foster_0.68_fit = quantile(unlist(sample_pCO2), c(0.68)),
           pCO2_foster_0.95_fit = quantile(unlist(sample_pCO2), c(0.95))) %>%
    select(record_id, sample_pCO2)
)
# create species-specific CIE offsets
species_offsets <- d18O_stack %>%
    select(fullsp) %>%
    unique() %>%
    mutate(sample_CO3_slope = list(rnorm(mc_samples, -0.002525, 0.0014))) # randomize CIE slope using mean and SD of the 4 species considered
# merge offset lists into d18O_stack and combine variables
d18O_stack <- d18O_stack %>% select(-c(sample_d18O_sw))
d18O_stack <- d18O_stack %>% select(-c(sample_d18O_ice))
d18O_stack <- d18O_stack %>% select(-c(sample_CO3))
d18O_stack <- d18O_stack %>% select(-c(sample_CO3_slope))
d18O_stack <- d18O_stack %>% select(-c(sample_pCO2))
d18O_stack <- d18O_stack %>%
    left_join(sample_offsets, by=c("record_id")) %>%
    left_join(species_offsets, by=c("fullsp")) %>%
    rowwise() %>%
    mutate(sample_d18O     = list(rep(d18O, mc_samples)), # bayfox already accounts for analytical and population errors, so we just create a duplicated list of raw d18O to simplify the plotting algorithm later
           sample_d18O_sw  = list(d18O_sw_CESM + sample_d18O_sw),
           sample_d18O_ice = list(d18O_ice + sample_d18O_ice),
           sample_CO3      = list(CO3_zeebe + sample_CO3))
if (pH_mode == FALSE) {
    # calculate SSTs using bayfox and the carbonate-ion effect derived from the [CO3] record
    d18O_stack <- d18O_stack %>%
        rowwise() %>%
        mutate(sample_d18O_CO3 = list(sample_d18O - (sample_CO3_slope*sample_CO3 - sample_CO3_slope*200))) %>% # correct for CO3 as above; intercept calculation yields 0 offset at modern CO3 = 200 umol/kg
        mutate(sample_d18O_bayfox = list(sample_d18O_CO3 - sample_d18O_sw - sample_d18O_ice),
               sample_d18O_bayfox_noCO3 = list(sample_d18O - sample_d18O_sw - sample_d18O_ice),
               sample_d18O_bayfox_nospatial = list(sample_d18O - sample_d18O_ice))
    d18O_stack <- d18O_stack %>%
        rowwise() %>%
        mutate(sample_SST           = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox)), bayfox_sd)),
               sample_SST_noCO3     = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox_noCO3)), bayfox_sd)),
               sample_SST_nospatial = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox_nospatial)), bayfox_sd)))
} else {
    # calculate SSTs using bayfox and the carbonate-ion derived from the pCO2 (pH) record
    d18O_stack <- d18O_stack %>% # calculate initial guess (no pH effect)
        rowwise() %>%
        mutate(sample_d18O_bayfox = list(sample_d18O - sample_d18O_sw - sample_d18O_ice)) %>%
        mutate(sample_SST = list(pmin(pmax(predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox)), 0), 50))) # capped at 0-50 (valid range for seacarb)
    # iterate pH and temperature to converge on the unique solution for each sample
    max_error = 1
    mean_error = 1
    while (max_error > 0.5 || mean_error > 0.01) {
        # get pH solution
        d18O_stack <- d18O_stack %>%
            rowwise() %>%
            filter(!is.na(d18O_sw_CESM), !is.na(d18O)) %>% # required because seacarb dislikes NAs
            mutate(sample_pH = list(carb(23, sample_pCO2, sample_CO3*1e-6, T=sample_SST)$pH)) %>%
            mutate(sample_d18O_CO3 = list(sample_d18O - (ciefit$coefficients[1] + ciefit$coefficients[2]*sample_pH + ciefit$coefficients[3]*sample_pH^2 + ciefit$coefficients[4]*sample_pH^3))) %>% # correct for pH using transformation fit above
            mutate(sample_d18O_bayfox = list(sample_d18O_CO3 - sample_d18O_sw - sample_d18O_ice),
                   sample_d18O_bayfox_noCO3 = list(sample_d18O - sample_d18O_sw - sample_d18O_ice),
                   sample_d18O_bayfox_nospatial = list(sample_d18O - sample_d18O_ice))
        # get temperature solution (w/o bayfox error, so we can check for convergence)
        d18O_stack <- d18O_stack %>%
            rowwise() %>%
            mutate(sample_SST_out = list(pmin(pmax(predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox)), 0), 50)))
        max_error = max(abs(unlist(d18O_stack$sample_SST) - unlist(d18O_stack$sample_SST_out)), na.rm=T)
        mean_error = mean(abs(unlist(d18O_stack$sample_SST) - unlist(d18O_stack$sample_SST_out)), na.rm=T)
        d18O_stack <- d18O_stack %>% mutate(sample_SST = list(sample_SST_out))
    }
    # run final temperature solution including bayfox error
    d18O_stack <- d18O_stack %>%
        rowwise() %>%
        mutate(sample_SST           = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox)), bayfox_sd)),
               sample_SST_noCO3     = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox_noCO3)), bayfox_sd)),
               sample_SST_nospatial = list(rnorm(mc_samples, predict_poly(bayfox_model_50, unlist(sample_d18O_bayfox_nospatial)), bayfox_sd)))
}
# get quantiles for plotting
d18O_stack <- d18O_stack %>%
    mutate(SST_montecarlo_2.5  = quantile(unlist(sample_SST), c(0.025), na.rm=T)[1],
           SST_montecarlo_25   = quantile(unlist(sample_SST), c(0.25),  na.rm=T)[1],
           SST_montecarlo_75   = quantile(unlist(sample_SST), c(0.75),  na.rm=T)[1],
           SST_montecarlo_97.5 = quantile(unlist(sample_SST), c(0.975), na.rm=T)[1],
           SST_montecarlo_noCO3_2.5  = quantile(unlist(sample_SST_noCO3), c(0.025), na.rm=T)[1],
           SST_montecarlo_noCO3_25   = quantile(unlist(sample_SST_noCO3), c(0.25),  na.rm=T)[1],
           SST_montecarlo_noCO3_75   = quantile(unlist(sample_SST_noCO3), c(0.75),  na.rm=T)[1],
           SST_montecarlo_noCO3_97.5 = quantile(unlist(sample_SST_noCO3), c(0.975), na.rm=T)[1],
           SST_montecarlo_nospatial_2.5  = quantile(unlist(sample_SST_nospatial), c(0.025), na.rm=T)[1],
           SST_montecarlo_nospatial_25   = quantile(unlist(sample_SST_nospatial), c(0.25),  na.rm=T)[1],
           SST_montecarlo_nospatial_75   = quantile(unlist(sample_SST_nospatial), c(0.75),  na.rm=T)[1],
           SST_montecarlo_nospatial_97.5 = quantile(unlist(sample_SST_nospatial), c(0.975), na.rm=T)[1])

# treat D47 as a preservation type for Fig. 1
d18O_stack <- d18O_stack %>% filter(!(preservation == "D"))
d18O_stack <- d18O_stack %>% full_join(D47 %>%
                                           rowwise() %>%
                                           mutate(sample_SST = list(rnorm(mc_samples, SST, SST_error / 2)), # treat +/- error as 2 standard deviations
                                                  sample_SST_nospatial = list(rnorm(mc_samples, SST, SST_error / 2)), # treat +/- error as 2 standard deviations
                                                  sample_d18O = list(rep(NA, mc_samples)),
                                                  latzone = ceiling(abs(pallat) / 30)) %>%
                                           filter(pallat <= 30) %>%
                                           select(site, age, latzone, pallat, pallong, ref, sample_SST, sample_SST_nospatial) %>%
                                           mutate(preservation = "D",
                                                  depth = 1,
                                                  temp_benthic = approx(d18O_benthic$age, d18O_benthic$temp_benthic, age)$y)) # merge benthic temperature into new data with linear interpolation

# SST conversion: bayfox
d18O_stack$d18O_bayfox <- d18O_stack$d18O_CO3 - d18O_stack$d18O_sw
d18O_stack$SST_Bayfox_2.5  <- predict(bayfox_model_2.5,  newdata=d18O_stack)
d18O_stack$SST_Bayfox_5    <- predict(bayfox_model_5,    newdata=d18O_stack)
d18O_stack$SST_Bayfox_50   <- predict(bayfox_model_50,   newdata=d18O_stack)
d18O_stack$SST_Bayfox_95   <- predict(bayfox_model_95,   newdata=d18O_stack)
d18O_stack$SST_Bayfox_97.5 <- predict(bayfox_model_97.5, newdata=d18O_stack)
d18O_stack$d18O_bayfox <- d18O_stack$d18O - d18O_stack$d18O_sw
d18O_stack$SST_Bayfox_noCO3_5  <- predict(bayfox_model_5,  newdata=d18O_stack)
d18O_stack$SST_Bayfox_noCO3_50 <- predict(bayfox_model_50, newdata=d18O_stack)
d18O_stack$SST_Bayfox_noCO3_95 <- predict(bayfox_model_95, newdata=d18O_stack)
d18O_stack$d18O_bayfox <- d18O_stack$d18O - d18O_stack$d18O_ice
d18O_stack$SST_Bayfox_no_spatial_5  <- predict(bayfox_model_5,  newdata=d18O_stack)
d18O_stack$SST_Bayfox_no_spatial_50 <- predict(bayfox_model_50, newdata=d18O_stack)
d18O_stack$SST_Bayfox_no_spatial_95 <- predict(bayfox_model_95, newdata=d18O_stack)
d18O_stack$d18O_bayfox <- d18O_stack$d18O - (d18O_stack$d18O_ice + d18O_stack$d18O_sw_Valdes)
d18O_stack$SST_Bayfox_Valdes_5  <- predict(bayfox_model_5,  newdata=d18O_stack)
d18O_stack$SST_Bayfox_Valdes_50 <- predict(bayfox_model_50, newdata=d18O_stack)
d18O_stack$SST_Bayfox_Valdes_95 <- predict(bayfox_model_95, newdata=d18O_stack)

#===============================================================================

# get Monte Carlo probabilities for a given binned dataset
make_d18O_ribbon_mc <- function(data_in, binsize, latzone_filter, preservation_filter, merge_and_lm = F) {
    prefiltered <- data_in %>% filter(depth == 1, latzone %in% latzone_filter, preservation %in% preservation_filter, latzone < 3 | pallat < 0)
    lm_intercepts <- rep(NA, mc_samples)
    lm_slopes <- rep(NA, mc_samples)
    d18O_benthic_min <- round(min(prefiltered$d18O_benthic, na.rm=T) * (1/binsize)) / (1/binsize)
    d18O_benthic_max <- round(max(prefiltered$d18O_benthic, na.rm=T) * (1/binsize)) / (1/binsize)
    
    # bootstrap sampling and extract MC samples of relevant points
    run_list <- data_in %>% # warning: this gets memory-heavy for large mc_count!
        ungroup() %>%
        mutate(matches = ((depth == 1) & (latzone %in% latzone_filter) & (preservation %in% preservation_filter) & (latzone < 3 | pallat < 0))) %>%
        select(matches, preservation, latzone, sample_d18O, d18O_benthic) %>%
        mutate(row = row_number()) %>%
        unnest(cols=(sample_d18O)) %>%
        rename(d18O=sample_d18O) %>%
        group_by(row) %>%
        mutate(sample = row_number()) %>%
        group_by(sample) %>%
        slice_sample(prop=1, replace=T) %>% # effectively doing a basic bootstrap on all medians/lms - each run becomes a separate bootstrap sample. This allows us to include sampling error in the Monte Carlo error distribution.
        filter(matches==T) %>% # only finally filter out irrelevant data after bootstrapping sampling of whole dataset
        mutate(d18O_benthic_binned = round(pmax(d18O_benthic_min, pmin(d18O_benthic_max, (d18O_benthic + runif(1, -binsize, binsize)))) * (1/binsize)) / (1/binsize)) # slightly randomizes bin to approximate temporal error, capped at ends to reduce the "flattening" of the LMs that would otherwise happen.
    
    if (merge_and_lm) run_list$preservation <- "E"
    
    # summarize into bins
    run_list <- run_list %>%
        group_by(sample, d18O_benthic_binned, preservation, latzone) %>%
        summarize(d18O = median(d18O), .groups="keep")
        
    if (merge_and_lm) {
        # calculate LMs
        run_lms <- run_list %>%
            group_by(sample) %>%
            do(tidy(lm(d18O ~ d18O_benthic_binned, .))) %>%
            select(sample, term, estimate) %>%
            pivot_wider(id_cols=sample, names_from=term, values_from=estimate) %>%
            rename(lm_intercept="(Intercept)", lm_slope="d18O_benthic_binned")
        lm_intercepts <- run_lms$lm_intercept
        lm_slopes <- run_lms$lm_slope
        run_list <- run_list %>%
            left_join(run_lms, by="sample") %>%
            mutate(lm_predict = d18O_benthic_binned * lm_slope + lm_intercept)

        # save lm coefficients to d18O_benthic
        d18O_benthic[,paste("sample_lm_intercepts_", latzone_filter[1], sep="")] <<- list(rep(list(lm_intercepts), nrow(d18O_benthic)))
        d18O_benthic[,paste("sample_lm_slopes_", latzone_filter[1], sep="")] <<- list(rep(list(lm_slopes), nrow(d18O_benthic)))
        lm_intercept_2sd <- 2*sd(lm_intercepts, na.rm=T)
        lm_slope_2sd     <- 2*sd(lm_slopes, na.rm=T)
    } else {
        run_list$lm_predict <- NA
        lm_intercept_2sd <- NA
        lm_slope_2sd     <- NA
    }
    
    # summarize with mean/sd for runs and lm
    run_list %>%
        group_by(d18O_benthic_binned, preservation, latzone) %>%
        summarize(d18O_quantile_2.5  = quantile(d18O,       c(0.025), na.rm=T)[1],
                  d18O_quantile_25   = quantile(d18O,       c(0.25),  na.rm=T)[1],
                  d18O_quantile_75   = quantile(d18O,       c(0.75),  na.rm=T)[1],
                  d18O_quantile_97.5 = quantile(d18O,       c(0.975), na.rm=T)[1],
                  lm_quantile_2.5    = quantile(lm_predict, c(0.025), na.rm=T)[1],
                  lm_quantile_25     = quantile(lm_predict, c(0.25),  na.rm=T)[1],
                  lm_quantile_75     = quantile(lm_predict, c(0.75),  na.rm=T)[1],
                  lm_quantile_97.5   = quantile(lm_predict, c(0.975), na.rm=T)[1],
                  d18O = mean(d18O),
                  lm = mean(lm_predict)) %>%
        mutate(lm_intercept_2sd = lm_intercept_2sd,
               lm_slope_2sd = lm_slope_2sd)
}
make_reg_label_mc <- function(lm_model, intercept_2sd, slope_2sd) {
    paste("y = ", round(lm_model$coefficients[2], 2), "±", round(slope_2sd, 2),
          "x ", ifelse(lm_model$coefficients[1] >= 0,
                       paste("+", round(lm_model$coefficients[1], 2)),
                       paste("-", -round(lm_model$coefficients[1], 2))),
          "±", round(intercept_2sd, 1), "\n",
          "R² = ", round(summary(lm_model)$r.squared, 2), sep="")
}
make_reg_text_mc <- function(lm1, lm2, lm3, mc_summary) {
    tibble(latzone = c(1, 2, 3),
           preservation = "E",
           label = c(make_reg_label_mc(lm1, unlist(mc_summary[which(mc_summary$latzone == 1),"lm_intercept_2sd"])[1], unlist(mc_summary[which(mc_summary$latzone == 1),"lm_slope_2sd"])[1]),
                     make_reg_label_mc(lm2, unlist(mc_summary[which(mc_summary$latzone == 2),"lm_intercept_2sd"])[1], unlist(mc_summary[which(mc_summary$latzone == 2),"lm_slope_2sd"])[1]),
                     make_reg_label_mc(lm3, unlist(mc_summary[which(mc_summary$latzone == 3),"lm_intercept_2sd"])[1], unlist(mc_summary[which(mc_summary$latzone == 3),"lm_slope_2sd"])[1])))
}
plot_gradients_mc <- function(data_in, binsize, lab_title, lab_x, lab_y, breaks_x=c(-2,-1,0,1,2,3,4,5), breaks_y=c(-6,-4,-2,0,2,4), limits_x=c(-2,4.5), limits_y=c(-6,4)) {
    d18O_plot <- make_d18O_ribbon_mc(data_in, binsize, c(1,2,3), c("E","VG","G","M","P","T","D"), F)
    d18O_plot$plotshape <- 1
    d18O_plot[which(d18O_plot$latzone == 1 & d18O_plot$preservation == "E"),"plotshape"] <- 19
    d18O_plot[which(d18O_plot$latzone == 2 & d18O_plot$preservation %in% c("E","VG")),"plotshape"] <- 19
    d18O_plot[which(d18O_plot$latzone == 3 & d18O_plot$preservation %in% c("E","VG","G","M","P")),"plotshape"] <- 19
    d18O_lms <- make_d18O_ribbon_mc(data_in, binsize, c(1), c("E"), T) %>%
        full_join(make_d18O_ribbon_mc(data_in, binsize, c(2), c("E","VG"), T)) %>%
        full_join(make_d18O_ribbon_mc(data_in, binsize, c(3), c("E","VG","G","M","P"), T))
    lowlat_lm_b  <- lm(d18O ~ d18O_benthic_binned, data = d18O_lms %>% filter(latzone == 1))
    midlat_lm_b  <- lm(d18O ~ d18O_benthic_binned, data = d18O_lms %>% filter(latzone == 2))
    highlat_lm_b <- lm(d18O ~ d18O_benthic_binned, data = d18O_lms %>% filter(latzone == 3))
    reg_lines <- tibble(latzone = c(1, 1, 2, 2, 3, 3),
                        preservation = "E",
                        d18O_benthic_binned = c(min(d18O_plot[which(d18O_plot$latzone == 1),"d18O_benthic_binned"], na.rm=T),
                                                max(d18O_plot[which(d18O_plot$latzone == 1),"d18O_benthic_binned"], na.rm=T),
                                                min(d18O_plot[which(d18O_plot$latzone == 2),"d18O_benthic_binned"], na.rm=T),
                                                max(d18O_plot[which(d18O_plot$latzone == 2),"d18O_benthic_binned"], na.rm=T),
                                                min(d18O_plot[which(d18O_plot$latzone == 3),"d18O_benthic_binned"], na.rm=T),
                                                max(d18O_plot[which(d18O_plot$latzone == 3),"d18O_benthic_binned"], na.rm=T)),
                        d18O = c(predict(lowlat_lm_b,  newdata = tibble(d18O_benthic_binned = min(d18O_plot[which(d18O_plot$latzone == 1),"d18O_benthic_binned"], na.rm=T))),
                                 predict(lowlat_lm_b,  newdata = tibble(d18O_benthic_binned = max(d18O_plot[which(d18O_plot$latzone == 1),"d18O_benthic_binned"], na.rm=T))),
                                 predict(midlat_lm_b,  newdata = tibble(d18O_benthic_binned = min(d18O_plot[which(d18O_plot$latzone == 2),"d18O_benthic_binned"], na.rm=T))),
                                 predict(midlat_lm_b,  newdata = tibble(d18O_benthic_binned = max(d18O_plot[which(d18O_plot$latzone == 2),"d18O_benthic_binned"], na.rm=T))),
                                 predict(highlat_lm_b, newdata = tibble(d18O_benthic_binned = min(d18O_plot[which(d18O_plot$latzone == 3),"d18O_benthic_binned"], na.rm=T))),
                                 predict(highlat_lm_b, newdata = tibble(d18O_benthic_binned = max(d18O_plot[which(d18O_plot$latzone == 3),"d18O_benthic_binned"], na.rm=T)))))
    ggplot(d18O_plot, aes(x=d18O_benthic_binned, y=d18O, color=preservation)) +
        facet_grid(rows=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`2`="Mid-Latitude",`3`="High-Latitude"))) +
        geom_ribbon(data=d18O_lms, aes(ymin=lm_quantile_2.5, ymax=lm_quantile_97.5), color="gray50", fill="gray50", na.rm=T) +
        geom_line(data=reg_lines, size=1) +
        geom_linerange(data=d18O_plot, aes(ymin=d18O_quantile_2.5, ymax=d18O_quantile_97.5), alpha=1) +
        geom_point(size=1.5, aes(shape=as.character(plotshape))) +
        geom_text(data=make_reg_text_mc(lowlat_lm_b, midlat_lm_b, highlat_lm_b, d18O_lms), aes(label=label), x=(limits_x[1]+limits_x[2])/2, y=limits_y[2]-(limits_y[2]-limits_y[1])/15, size=3, color="black") +
        scale_color_manual(name = "Preservation", values=c("E"="black","VG"="red","G"="orange","M"="yellow","P"="gray50", "T"="blue", "D"="blue"), labels=c("E"="Excellent","VG"="Very Good","G"="Good","M"="Moderate","P"="Poor","T"="TEX86","D"="Clumped Isotopes")) +
        scale_shape_manual(values=c(1,19)) +
        labs(title=lab_title,
             x=lab_x,
             y=lab_y) +
        guides(color=F, shape=F) +
        scale_x_continuous(breaks=breaks_x, limits=limits_x) +
        scale_y_continuous(breaks=breaks_y, limits=limits_y) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              strip.background = element_blank())
}

# Gradients Figure 1, top panel: full compilation
p1 <- ggplot(d18O_stack %>% filter(preservation %in% c("E","VG","G","M","P"), latzone %in% c(1,2,3)) %>%
                 arrange(desc(factor(preservation, levels = c("E","VG","G","M","P")))),
             aes(x = age, y = d18O, color = preservation, shape = as.factor(depth), alpha = as.factor(depth))) +
    facet_grid(rows=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`2`="Mid-Latitude",`3`="High-Latitude"))) +
    geom_point(size=1.5) +
    geom_line(data = d18O_benthic %>% mutate(latzone = 1, preservation = "E", depth = 1), aes(x=age, y=d18O)) +
    geom_line(data = d18O_benthic %>% mutate(latzone = 2, preservation = "E", depth = 1), aes(x=age, y=d18O)) +
    geom_line(data = d18O_benthic %>% mutate(latzone = 3, preservation = "E", depth = 1), aes(x=age, y=d18O)) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,max(d18O_stack$age))) +
    scale_y_reverse(breaks = c(-4, -2, 0, 2, 4, 6)) +
    scale_color_manual(name = "Preservation", values=c("E"="black","VG"="red","G"="orange","M"="yellow","P"="gray50"), labels=c("E"="Excellent","VG"="Very Good","G"="Good","M"="Moderate","P"="Poor")) +
    scale_shape_manual(name = "Depth Habitat", values=c(16,4,4,4), labels=c("Mixed-layer","Deeper","Deeper","Deeper")) +
    scale_alpha_manual(values=c(1,0.25,0.25,0.25)) +
    guides(alpha = F) +
    labs(x = "Age (Ma)",
         y = expression(paste(delta^18, "O (\211)"))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          legend.position="bottom")
ggsave("gradients_fig1d_top.pdf", width=10.25, height=6, dpi=150)

# Gradients Figure 1: step-by-step methods
p1da <- plot_gradients_mc(d18O_stack %>% filter(preservation != "T" & preservation != "D"),
                       0.25,
                       expression(paste("Raw ", delta^18, "O")),
                       expression(paste("Benthic ", delta^18, "O (\211)")), 
                       expression(paste("Planktonic ", delta^18, "O (\211)")),
                       limits_x = c(-1.75,4),
                       limits_y = c(-6,5))
p1db <- plot_gradients_mc(d18O_stack %>% filter(preservation != "T" & preservation != "D") %>% select(-c("sample_d18O")) %>% rename(sample_d18O = sample_d18O_bayfox_nospatial) %>% mutate(d18O_benthic = d18O_benthic - d18O_ice),
                       0.25,
                       bquote("- sea-level "*delta^18*O['ice']),
                       expression(paste("Benthic ", delta^18, "O (\211)")), 
                       expression(paste("Planktonic ", delta^18, "O (\211)")),
                       limits_x = c(-1.75,4),
                       limits_y = c(-6,5))
p1dc <- plot_gradients_mc(d18O_stack %>% filter(preservation != "T" & preservation != "D") %>% select(-c("sample_d18O")) %>% rename(sample_d18O = sample_d18O_bayfox) %>% mutate(d18O_benthic = d18O_benthic - d18O_ice),
                       0.25,
                       bquote("- local "*delta^18*O['sw']),
                       expression(paste("Benthic ", delta^18, "O (\211)")), 
                       expression(paste("Planktonic ", delta^18, "O (\211)")),
                       limits_x = c(-1.75,4),
                       limits_y = c(-6,5))
p1dd <- plot_gradients_mc(d18O_stack %>% select(-c("sample_d18O", "d18O_benthic")) %>% rename(sample_d18O = sample_SST_nospatial, d18O_benthic = temp_benthic),
                       1,
                       "SST (sea level only)",
                       "Bottom-water temperature (°C)", 
                       "Sea-surface temperature (°C)",
                       breaks_x = c(0,5,10,15,20),
                       breaks_y = c(-10,0,10,20,30,40,50),
                       limits_x = c(-3,22),
                       limits_y = c(-10,50))
p1de <- plot_gradients_mc(d18O_stack %>% select(-c("sample_d18O", "d18O_benthic")) %>% rename(sample_d18O = sample_SST, d18O_benthic = temp_benthic),
                          1,
                          "SST (with local)",
                          "Bottom-water temperature (°C)", 
                          "Sea-surface temperature (°C)",
                          breaks_x = c(0,5,10,15,20),
                          breaks_y = c(-10,0,10,20,30,40,50),
                          limits_x = c(-3,22),
                          limits_y = c(-10,50))
p1d <- grid.arrange(p1da, p1db, p1dc, p1dd, p1de, ncol=5)
ggsave("gradients_fig1d.pdf", p1d, width=12.5, height=5.75, dpi=150)

# Gradients Figure 2: SST by age
d18O_benthic <- d18O_benthic %>%
    rowwise() %>%
    mutate(sample_SST_1 = list(unlist(sample_lm_slopes_1)*temp_benthic + unlist(sample_lm_intercepts_1)),
           sample_SST_2 = list(unlist(sample_lm_slopes_2)*temp_benthic + unlist(sample_lm_intercepts_2)),
           sample_SST_3 = list(unlist(sample_lm_slopes_3)*temp_benthic + unlist(sample_lm_intercepts_3)))
d18O_benthic <- d18O_benthic %>%
    mutate(SST_montecarlo_1_2.5  = quantile(unlist(sample_SST_1), c(0.025), na.rm=T)[1],
           SST_montecarlo_1_25   = quantile(unlist(sample_SST_1), c(0.25),  na.rm=T)[1],
           SST_montecarlo_1_50   = quantile(unlist(sample_SST_1), c(0.50),  na.rm=T)[1],
           SST_montecarlo_1_75   = quantile(unlist(sample_SST_1), c(0.75),  na.rm=T)[1],
           SST_montecarlo_1_97.5 = quantile(unlist(sample_SST_1), c(0.975), na.rm=T)[1],
           SST_montecarlo_2_2.5  = quantile(unlist(sample_SST_2), c(0.025), na.rm=T)[1],
           SST_montecarlo_2_25   = quantile(unlist(sample_SST_2), c(0.25),  na.rm=T)[1],
           SST_montecarlo_2_50   = quantile(unlist(sample_SST_2), c(0.50),  na.rm=T)[1],
           SST_montecarlo_2_75   = quantile(unlist(sample_SST_2), c(0.75),  na.rm=T)[1],
           SST_montecarlo_2_97.5 = quantile(unlist(sample_SST_2), c(0.975), na.rm=T)[1],
           SST_montecarlo_3_2.5  = quantile(unlist(sample_SST_3), c(0.025), na.rm=T)[1],
           SST_montecarlo_3_25   = quantile(unlist(sample_SST_3), c(0.25),  na.rm=T)[1],
           SST_montecarlo_3_50   = quantile(unlist(sample_SST_3), c(0.50),  na.rm=T)[1],
           SST_montecarlo_3_75   = quantile(unlist(sample_SST_3), c(0.75),  na.rm=T)[1],
           SST_montecarlo_3_97.5 = quantile(unlist(sample_SST_3), c(0.975), na.rm=T)[1],
           LTG_montecarlo_2.5    = quantile(unlist(sample_SST_1) - unlist(sample_SST_3), c(0.025), na.rm=T)[1],
           LTG_montecarlo_25     = quantile(unlist(sample_SST_1) - unlist(sample_SST_3), c(0.25),  na.rm=T)[1],
           LTG_montecarlo_50     = quantile(unlist(sample_SST_1) - unlist(sample_SST_3), c(0.50),  na.rm=T)[1],
           LTG_montecarlo_75     = quantile(unlist(sample_SST_1) - unlist(sample_SST_3), c(0.75),  na.rm=T)[1],
           LTG_montecarlo_97.5   = quantile(unlist(sample_SST_1) - unlist(sample_SST_3), c(0.975), na.rm=T)[1])
d18O_stack_filtered <- d18O_stack %>%
    filter((preservation == "E" & latzone == 1) | (latzone == 3 & pallat < 0), depth == 1)
p2d <- ggplot() +
    geom_linerange(data=d18O_stack_filtered, aes(x=age, ymin=SST_montecarlo_25, ymax=SST_montecarlo_75, color=as.factor(latzone)), alpha=0.5, size=0.5) +
    geom_linerange(data=d18O_stack_filtered, aes(x=age, ymin=SST_montecarlo_2.5, ymax=SST_montecarlo_97.5, color=as.factor(latzone)), alpha=0.1, size=0.5) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=SST_montecarlo_3_2.5,  ymax=SST_montecarlo_3_97.5),  fill="darkslategray3", alpha=0.33) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=SST_montecarlo_3_25,   ymax=SST_montecarlo_3_75),    fill="darkslategray3", alpha=0.5) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=SST_montecarlo_1_2.5,  ymax=SST_montecarlo_1_97.5),  fill="goldenrod2", alpha=0.33) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=SST_montecarlo_1_25,   ymax=SST_montecarlo_1_75),    fill="goldenrod2", alpha=0.5) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=LTG_montecarlo_2.5*2-65, ymax=LTG_montecarlo_97.5*2-65), fill="gray50", alpha=0.33) +
    geom_ribbon(data=d18O_benthic, aes(x=age, ymin=LTG_montecarlo_25*2-65,  ymax=LTG_montecarlo_75*2-65),   fill="gray50", alpha=0.5) +
    geom_line(data=d18O_benthic, aes(x=age, y=temp_benthic), size=1) +
    geom_line(data=d18O_benthic, aes(x=age, y=SST_montecarlo_3_50),  size=1, color="darkslategray3") +
    geom_line(data=d18O_benthic, aes(x=age, y=SST_montecarlo_1_50),  size=1, color="goldenrod2") +
    geom_line(data=d18O_benthic, aes(x=age, y=LTG_montecarlo_50*2-65), size=1, color="gray50") +
    geom_point(data=d18O_stack_filtered, aes(x=age, y=SST_Bayfox_50, color=as.factor(latzone)), alpha=1, size=0.5) +
    geom_linerange(data=D47 %>% filter(latzone == 3, pallat < 0), aes(x=age, ymin=SST-SST_error, ymax=SST+SST_error), alpha=0.33, size=0.5, color="black") +
    geom_linerange(data=D47 %>% filter(latzone == 1), aes(x=age, ymin=SST-SST_error, ymax=SST+SST_error), alpha=0.33, size=0.5, color="black") +
    geom_point(data=D47 %>% filter(latzone == 3, pallat < 0), aes(x=age, y=SST), size=1, shape=23, color="black", fill="darkslategray4", alpha=1) +
    geom_point(data=D47 %>% filter(latzone == 1), aes(x=age, y=SST), size=1, shape=23, color="black", fill="goldenrod4", alpha=1) +
    geom_point(data=tibble(x=0, y=26.19064), aes(x=x, y=y), shape=21, size=3, color="black", fill="goldenrod2") +
    geom_point(data=tibble(x=0, y=-0.1950831), aes(x=x, y=y), shape=21, size=3, color="black", fill="darkslategray3") +
    geom_hline(yintercept = -11) +
    labs(x="Age (Ma)",
         y="Temperature (°C)") +
    guides(color=F, fill=F) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,max(d18O_stack_filtered$age))) +
    scale_y_continuous(breaks = c(-10,0,10,20,30,40), sec.axis = sec_axis(~ (. + 65) / 2, name="LTG (°C)", breaks=c(15,20,25))) +
    scale_color_manual(name = "Latitude", labels = c("Tropical", "High-Latitude", ""), values = c("orange2", "dodgerblue4", "")) +
    scale_fill_manual(name = "Latitude", labels = c("Tropical", "High-Latitude", ""), values = c("orange1", "dodgerblue4", "")) +
    coord_cartesian(ylim = c(-35, 42), expand=T) +
    theme_bw() +
    theme(panel.grid.major = element_line(color="gray90", size=0.1),
          panel.grid.minor = element_line(color="gray90", size=0.1))
ggsave("gradients_fig2d.pdf", p2d, width=6, height=6, dpi=150)

# Gradients Figure 3: Gradient by SST (Monte Carlo version)
sample_lm_density <- tibble(temp_benthic=seq(-10,30,0.1)) %>%
    rowwise() %>%
    mutate(SST_1=list(temp_benthic*unlist(d18O_benthic[1,"sample_lm_slopes_1"]) + unlist(d18O_benthic[1,"sample_lm_intercepts_1"])),
           SST_3=list(temp_benthic*unlist(d18O_benthic[1,"sample_lm_slopes_3"]) + unlist(d18O_benthic[1,"sample_lm_intercepts_3"]))) %>%
    mutate(SST_mean=list(0.5*unlist(SST_1) + 0.366*(unlist(SST_1)+unlist(SST_3))/2 + 0.134*unlist(SST_3)),
           SST_gradient=list(unlist(SST_1)-unlist(SST_3))) %>%
    ungroup() %>%
    select(SST_mean, SST_gradient) %>%
    unnest(cols=c(SST_mean, SST_gradient)) %>%
    mutate(SST_mean_rounded = round(SST_mean)) %>%
    group_by(SST_mean_rounded) %>%
    summarize(SST_gradient_2.5  = quantile(unlist(SST_gradient), c(0.025), na.rm=T)[1],
              SST_gradient_25   = quantile(unlist(SST_gradient), c(0.25),  na.rm=T)[1],
              SST_gradient_75   = quantile(unlist(SST_gradient), c(0.75),  na.rm=T)[1],
              SST_gradient_97.5 = quantile(unlist(SST_gradient), c(0.975), na.rm=T)[1])
plot_lines <- read.csv("slope_plottable.csv", header=T) %>%
    full_join(tibble(Mean=c(20, 20), Gradient=c(22, 20), Proxy=c("Modern", "This study"), Category=c("Proxy", "Proxy")))
ggplot() +
    geom_line(data=plot_lines, aes(x=Mean, y=Gradient, color=Proxy, group=paste(Proxy, Category), linetype=Category), size=0.5) +
    geom_ribbon(data=sample_lm_density, aes(x=SST_mean_rounded, ymin=SST_gradient_2.5, ymax=SST_gradient_97.5), fill="red", alpha=0.1) +
    geom_ribbon(data=sample_lm_density, aes(x=SST_mean_rounded, ymin=SST_gradient_25, ymax=SST_gradient_75), fill="red", alpha=0.25) +
    geom_smooth(data=d18O_benthic, aes(x=0.5*SST_montecarlo_1_50 + 0.366*(SST_montecarlo_1_50+SST_montecarlo_3_50)/2 + 0.134*SST_montecarlo_3_50, y=SST_montecarlo_1_50-SST_montecarlo_3_50), method="lm", se=F, color="red") +
    geom_point(data=tibble(x=c(17.70948), y=c(25.45324)), aes(x=x, y=y), size=4) +
    labs(x="Mean global SST (°C)",
         y="Tropical SST - High-Latitude SST (°C)",
         color="Source") +
    scale_x_continuous(breaks=c(14, 16, 18, 20, 22, 24, 26, 28, 30, 32)) +
    scale_y_continuous(breaks=c(6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 32)) +
    scale_color_manual(values=c("CESM2"="yellow1", "CESM1.2"="#fde725", "CCSM4"="orange2", "COSMOS"="#74d055", "EC-Earth3"="#29af7f", "GFDL"="#238a8d", "GISS-E2"="#32648e", "HadCM3"="#453781", "HadGEM3"="#481568", "IPSL"="#5150b5", "NorESM"="#22caff", "Modern"="black", "This study"="red", "Cramwinckel compilation"="gray25", "Sijp compilation"="gray50", "Zhang compilation"="gray75")) +
    scale_linetype_manual(values=c("Proxy"="solid", "Eocene (DeepMIP)"="twodash", "Pliocene (PlioMIP2)"="dashed", "Cretaceous"="dotted")) +
    coord_cartesian(xlim=c(15.2,31.8), ylim=c(6.5,27)) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position=c(0.5,0.5),
          legend.text=element_text(size=8))
ggsave("gradients_fig3d.pdf", width=4, height=4, dpi=150)

# Gradients Figure 4: Residuals vs. time
data_residuals <- d18O_stack_filtered %>% mutate(SST = SST_Bayfox_50, SST_min = SST_Bayfox_2.5, SST_max = SST_Bayfox_97.5) %>%
    select(age, fullsp, latzone, SST, SST_min, SST_max, band_temp_residual) %>%
    mutate(genus = gsub(" .*$", "", fullsp))
data_residuals[which(data_residuals$fullsp == "Subbotina triangularis"), "genus"] <- "Subbotina triangularis"
data_residuals[which(data_residuals$fullsp == "Archaeoglobigerina australis"), "genus"] <- "Archaeoglobigerina australis"
data_residuals[which(data_residuals$fullsp == "Archaeoglobigerina bosquensis"), "genus"] <- "Archaeoglobigerina bosquensis"
omit_list <- (data_residuals %>% filter(latzone == 3, age < 80) %>% count(genus) %>% filter(n < 50))$genus
data_residuals[which(data_residuals$genus %in% omit_list), "genus"] <- "Various (high-latitude)"
data_residuals[which(data_residuals$latzone == 1), "genus"] <- "Various (tropical)"
data_residuals[which(data_residuals$latzone == 1), "SST_inferred"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_50, unlist(data_residuals[which(data_residuals$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals[which(data_residuals$latzone == 1), "SST_inferred_min"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_2.5, unlist(data_residuals[which(data_residuals$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals[which(data_residuals$latzone == 1), "SST_inferred_max"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_97.5, unlist(data_residuals[which(data_residuals$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals[which(data_residuals$latzone == 3), "SST_inferred"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_50, unlist(data_residuals[which(data_residuals$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals[which(data_residuals$latzone == 3), "SST_inferred_min"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_2.5, unlist(data_residuals[which(data_residuals$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals[which(data_residuals$latzone == 3), "SST_inferred_max"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_97.5, unlist(data_residuals[which(data_residuals$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals$SST_residual <- data_residuals$SST - data_residuals$SST_inferred
data_residuals$SST_min_residual <- data_residuals$SST_min - data_residuals$SST_inferred
data_residuals$SST_max_residual <- data_residuals$SST_max - data_residuals$SST_inferred
p4 <- ggplot(data_residuals %>% filter(latzone != 2)) +
    facet_grid(rows=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`3`="High-Latitude"))) +
    geom_ribbon(data=d18O_benthic %>% mutate(latzone = 3), aes(x=age, ymin=SST_montecarlo_3_2.5-SST_montecarlo_3_50, ymax=SST_montecarlo_3_97.5-SST_montecarlo_3_50), fill="darkslategray3", alpha=0.33) +
    geom_ribbon(data=d18O_benthic %>% mutate(latzone = 3), aes(x=age, ymin=SST_montecarlo_3_25-SST_montecarlo_3_50,  ymax=SST_montecarlo_3_75-SST_montecarlo_3_50),   fill="darkslategray3", alpha=0.5) +
    geom_ribbon(data=d18O_benthic %>% mutate(latzone = 1), aes(x=age, ymin=SST_montecarlo_1_2.5-SST_montecarlo_1_50, ymax=SST_montecarlo_1_97.5-SST_montecarlo_1_50), fill="goldenrod2", alpha=0.33) +
    geom_ribbon(data=d18O_benthic %>% mutate(latzone = 1), aes(x=age, ymin=SST_montecarlo_1_25-SST_montecarlo_1_50,  ymax=SST_montecarlo_1_75-SST_montecarlo_1_50),   fill="goldenrod2", alpha=0.5) +
    geom_point(aes(x=age, y=SST_residual, color=genus, fill=genus), size=0.5) +
    geom_hline(yintercept=0) +
    scale_color_manual(values=c("Various (tropical)"="orange2", "Various (high-latitude)"="gray25", "Acarinina"="#f3e51c", "Archaeoglobigerina australis"="#94d740", "Archaeoglobigerina bosquensis"="#cfe11d", "Chiloguembelina"="#37b877", "Globigerinatheka"="#20908c", "Neogloboquadrina"="#2f6d8e", "Subbotina triangularis"="#453680")) +
    labs(x="Age (Ma)",
         y="Residuals (°C)") +
    guides(fill=F) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,max(d18O_stack_filtered$age))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank())
ggsave("gradients_fig4.pdf", p4, width=9, height=6, dpi=150)

# regressing residuals against model spatial bias
summary(lm(SST_residual ~ band_temp_residual, data_residuals %>% filter(latzone == 3, band_temp_residual > -900)))
summary(lm(SST_residual ~ band_temp_residual, data_residuals %>% filter(latzone == 1, band_temp_residual > -900)))

# slopes of LTG:GMSST
all_lms <- head(d18O_benthic, 1) %>%
    select(sample_lm_slopes_1, sample_lm_slopes_3, sample_lm_intercepts_1, sample_lm_intercepts_3) %>%
    unnest(cols = c(sample_lm_slopes_1, sample_lm_slopes_3, sample_lm_intercepts_1, sample_lm_intercepts_3)) %>%
    rowwise() %>%
    mutate(SST_1 = list(c(sample_lm_slopes_1 * 1 + sample_lm_intercepts_1,
                          sample_lm_slopes_1 * 20 + sample_lm_intercepts_1)),
           SST_3 = list(c(sample_lm_slopes_3 * 1 + sample_lm_intercepts_3,
                          sample_lm_slopes_3 * 20 + sample_lm_intercepts_3)),
           sample = cur_group_id()) %>%
    unnest(cols = c(SST_1, SST_3)) %>%
    mutate(LTG = SST_1 - SST_3,
           MGT = 0.5*SST_1 + 0.366*(SST_1+SST_3)/2 + 0.134*SST_3) %>%
    group_by(sample) %>%
    do(tidy(lm(LTG ~ MGT, .))) %>%
    select(sample, term, estimate) %>%
    pivot_wider(id_cols=sample, names_from=term, values_from=estimate) %>%
    rename(lm_intercept="(Intercept)", lm_slope="MGT")
mean(all_lms$lm_slope)
mean(all_lms$lm_intercept)
sd(all_lms$lm_slope)*2
sd(all_lms$lm_intercept)*2

# slopes of LTG:BWT
summary(lm(y ~ x, d18O_benthic %>% mutate(x=temp_benthic, y=SST_montecarlo_1_50-SST_montecarlo_3_50)))
slope_differences = unlist(d18O_benthic[1,"sample_lm_slopes_1"]) - unlist(d18O_benthic[1,"sample_lm_slopes_3"])
mean(slope_differences)
sd(slope_differences)*2
intercept_differences = unlist(d18O_benthic[1,"sample_lm_intercepts_1"]) - unlist(d18O_benthic[1,"sample_lm_intercepts_3"])
mean(intercept_differences)
sd(intercept_differences)*2

# Polar Amplification Factor
summary(lm(SST_montecarlo_3_50 ~ x, d18O_benthic %>% mutate(x=0.5*SST_montecarlo_1_50 + 0.366*(SST_montecarlo_1_50+SST_montecarlo_3_50)/2 + 0.134*SST_montecarlo_3_50)))
sample_paf = unlist(d18O_benthic[1,"sample_lm_slopes_3"]) / (0.5*unlist(d18O_benthic[1,"sample_lm_slopes_1"]) + 0.366*(unlist(d18O_benthic[1,"sample_lm_slopes_1"])+unlist(d18O_benthic[1,"sample_lm_slopes_3"]))/2 + 0.134*unlist(d18O_benthic[1,"sample_lm_slopes_3"]))
median(sample_paf)
sd(sample_paf)*2

# Analyze intervals and residuals
max(d18O_benthic[which(d18O_benthic$age > 47.8 & d18O_benthic$age <= 56), "SST_montecarlo_1_97.5"])
min(d18O_benthic[which(d18O_benthic$age > 47.8 & d18O_benthic$age <= 56), "SST_montecarlo_1_2.5"])
mean((d18O_benthic %>% filter(age >= 49.1, age <= 53.4) %>% mutate(GMSST = 0.5*SST_montecarlo_1_50 + 0.366*(SST_montecarlo_1_50+SST_montecarlo_3_50)/2 + 0.134*SST_montecarlo_3_50))$GMSST)
mean((d18O_benthic %>% filter(age == 57) %>% mutate(GMSST = 0.5*SST_montecarlo_1_50 + 0.366*(SST_montecarlo_1_50+SST_montecarlo_3_50)/2 + 0.134*SST_montecarlo_3_50))$GMSST)
max(d18O_benthic$SST_montecarlo_1_50 - d18O_benthic$SST_montecarlo_3_50)
min(d18O_benthic$SST_montecarlo_1_50 - d18O_benthic$SST_montecarlo_3_50)
max(0.5*d18O_benthic$SST_montecarlo_1_50 + 0.366*(d18O_benthic$SST_montecarlo_1_50+d18O_benthic$SST_montecarlo_3_50)/2 + 0.134*d18O_benthic$SST_montecarlo_3_50)
min(0.5*d18O_benthic$SST_montecarlo_1_50 + 0.366*(d18O_benthic$SST_montecarlo_1_50+d18O_benthic$SST_montecarlo_3_50)/2 + 0.134*d18O_benthic$SST_montecarlo_3_50)
d18O_benthic[which(d18O_benthic$age == max(d18O_benthic$age)), "SST_montecarlo_1_97.5"]
d18O_benthic[which(d18O_benthic$age == max(d18O_benthic$age)), "SST_montecarlo_1_2.5"]

data_residuals2 <- D47 %>% mutate(symbol = "c", SST_min=SST-(SST_error*2), SST_max=SST+(SST_error*2))
data_residuals2[which(data_residuals2$latzone == 1), "SST_inferred"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_50, unlist(data_residuals2[which(data_residuals2$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2[which(data_residuals2$latzone == 1), "SST_inferred_min"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_2.5, unlist(data_residuals2[which(data_residuals2$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2[which(data_residuals2$latzone == 1), "SST_inferred_max"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_1_97.5, unlist(data_residuals2[which(data_residuals2$latzone == 1), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2[which(data_residuals2$latzone == 3), "SST_inferred"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_50, unlist(data_residuals2[which(data_residuals2$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2[which(data_residuals2$latzone == 3), "SST_inferred_min"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_2.5, unlist(data_residuals2[which(data_residuals2$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2[which(data_residuals2$latzone == 3), "SST_inferred_max"] <- approx(d18O_benthic$age, d18O_benthic$SST_montecarlo_3_97.5, unlist(data_residuals2[which(data_residuals2$latzone == 3), "age"]))$y # merge in benthic-derived records with linear interpolation
data_residuals2$SST_residual <- data_residuals2$SST - data_residuals2$SST_inferred
data_residuals2$SST_min_residual <- data_residuals2$SST_min - data_residuals2$SST_inferred
data_residuals2$SST_max_residual <- data_residuals2$SST_max - data_residuals2$SST_inferred

mean(data_residuals2[which(data_residuals2$latzone == 3 & data_residuals2$age >= 33.9 & data_residuals2$age <= 56), "SST_residual"], na.rm=T)
mean(data_residuals2[which(data_residuals2$latzone == 3 & data_residuals2$age >= 66), "SST_residual"], na.rm=T)
sd(unlist(data_residuals[which(data_residuals$latzone == 1 & data_residuals$age > 30), "SST_residual"]), na.rm=T)
sd(unlist(data_residuals[which(data_residuals$latzone == 3 & data_residuals$age > 30), "SST_residual"]), na.rm=T)
1 - (nrow(data_residuals %>% filter(latzone == 1 & (SST_min > SST_inferred_max | SST_max < SST_inferred_min))) / nrow(data_residuals %>% filter(latzone == 1))) # percentage of individual points that overlap inference envelope by at least 5% of total range (a = 0.05) - for this specific application, equivalent to a two-value T test (although it isn't always! be careful!)
1 - (nrow(data_residuals %>% filter(latzone == 3 & (SST_min > SST_inferred_max | SST_max < SST_inferred_min))) / nrow(data_residuals %>% filter(latzone == 3))) # percentage of individual points that overlap inference envelope by at least 5% of total range (a = 0.05) - for this specific application, equivalent to a two-value T test (although it isn't always! be careful!)
1 - (nrow(data_residuals2 %>% filter(latzone == 1 & (SST_min > SST_inferred_max | SST_max < SST_inferred_min))) / nrow(data_residuals2 %>% filter(latzone == 1))) # percentage of individual points that overlap inference envelope by at least 5% of total range (a = 0.05) - for this specific application, equivalent to a two-value T test (although it isn't always! be careful!)
1 - (nrow(data_residuals2 %>% filter(latzone == 3 & (SST_min > SST_inferred_max | SST_max < SST_inferred_min))) / nrow(data_residuals2 %>% filter(latzone == 3))) # percentage of individual points that overlap inference envelope by at least 5% of total range (a = 0.05) - for this specific application, equivalent to a two-value T test (although it isn't always! be careful!)

# raw R2 of benthic vs. planktic d18O
rawmodel <- lm(d18O ~ d18O_benthic, data=d18O_stack %>% filter(latzone == 3, depth == 1, age > 0, preservation %in% c("E","VG","G","M","P")))
summary(rawmodel)

# Mahalanobis distance from null model 1x + 0 to point (m,b)
xdist <- (1.07-1)/(0.1/2) # normalize to unit size - NOTE: not pulling 1.07 automatically from Fig. 1Eb!
ydist <- (0.53-0)/(1.8/2) # normalize to unit size - NOTE: not pulling 0.53 automatically from Fig. 1Eb!
mahalanobis_dist <- sqrt(xdist^2 + ydist^2) # reduced to Euclidean distance when normalized in this way and variables are independent
pchisq(mahalanobis_dist, df=1, lower.tail=FALSE) # p-value of rejecting the null model

# supplementary figures
meansst <- d18O_benthic %>% mutate(SST_mean_global = (0.5*SST_montecarlo_1_50 + 0.366*(SST_montecarlo_1_50+SST_montecarlo_3_50)/2 + 0.134*SST_montecarlo_3_50))
ggplot(meansst, aes(x=age, y=SST_mean_global, color=SST_mean_global)) +
    geom_line(size=2) +
    labs(title="Mean global SST through time",
         x="Age (Ma)",
         y="Mean global SST (°C)",
         color="SST") +
    guides(color=F) +
    scale_x_continuous(breaks=c(0,10,20,30,40,50,60,70,80,90,100)) +
    scale_y_continuous(breaks=c(16,18,20,22,24,26,28,30,32,34)) +
    scale_color_gradientn(colors=c("blue","cyan","yellow","orange","red","darkred")) +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
ggsave("supplement_meansst.svg", width=6, height=4, dpi=150)

mean(unlist(meansst[which(meansst$age >= 49.1 & meansst$age <= 53.4), "SST_mean_global"]), na.rm=T)
mean(unlist(meansst[which(meansst$age >= 57 & meansst$age <= 58), "SST_mean_global"]), na.rm=T)

min((0.5*d18O_benthic$SST_montecarlo_1_50 + 0.366*(d18O_benthic$SST_montecarlo_1_50+d18O_benthic$SST_montecarlo_3_50)/2 + 0.134*d18O_benthic$SST_montecarlo_3_50))
max((0.5*d18O_benthic$SST_montecarlo_1_50 + 0.366*(d18O_benthic$SST_montecarlo_1_50+d18O_benthic$SST_montecarlo_3_50)/2 + 0.134*d18O_benthic$SST_montecarlo_3_50))
min(d18O_benthic$temp_benthic)
max(d18O_benthic$temp_benthic)

# pH record
if (pH_mode == TRUE) {
    ggplot(d18O_benthic %>%
               mutate(pH_1_0.50 = carb(23, pCO2_foster_0.50, CO3_zeebe*1e-6, T=SST_montecarlo_1_50)$pH)) +
        geom_line(aes(x=age, y=pH_1_0.50)) +
        labs(x="Age (Ma)",
             y="pH") +
        theme_bw()
    ggsave("supplement_ph.svg", width=6, height=4, dpi=150)
}

# by site/study
d18O_stack_prefiltered <- d18O_stack %>% filter(preservation %in% c("E","VG","G","M","P"))
site_text <- d18O_stack_prefiltered %>%
    mutate(age_10Ma = round(age / 10) * 10) %>%
    group_by(ref, latzone, age_10Ma, preservation) %>%
    summarize(age = mean(age),
              d18O = mean(d18O)) %>%
    group_by(ref) %>%
    mutate(color=paste(min(age), ref)) %>% # to order colors approximately by age rather than alphabetically
    select(ref, latzone, age, d18O, preservation, color)
d18O_stack_prefiltered <- d18O_stack_prefiltered %>% left_join(site_text %>% select(ref, color))
site_colors = sample(grDevices::colors()[grep('gr(a|e)y', grDevices::colors(), invert = T)], length(unique(site_text$color))) # get a bunch of distinct colors (except grays)
ggplot() +
    facet_grid(rows=vars(preservation), cols=vars(latzone), labeller=labeller(latzone=c(`1`="Tropical",`2`="Mid-Latitude",`3`="High-Latitude"), preservation=c(`E`="Excellent",`EVG`="Very Good",`G`="Good",`M`="Moderate",`P`="Poor"))) + 
    geom_line(data=d18O_benthic, aes(x=age, y=d18O), size=0.25) +
    geom_point(data=d18O_stack_prefiltered %>% mutate(preservation = ifelse(preservation=="VG","EVG",preservation)), aes(x=age, y=d18O, color=color), alpha=1, size=0.25) +
    geom_label_repel(data=site_text %>% mutate(preservation = ifelse(preservation=="VG","EVG",preservation)), aes(x=age, y=d18O, fill=color, label=ref), size=0.75, label.padding=unit(0.05, "lines"), label.size=0, box.padding=unit(0.05, "lines"), segment.colour = NA, alpha=0.5, direction="y") +
    labs(title=expression(paste(delta^18, "O by site")),
         x="Age (Ma)",
         y=bquote(delta^18*O*" (\211)")) +
    scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100), limits=c(0,100)) +
    scale_y_reverse(breaks = c(-5,-4,-3,-2,-1,0,1,2,3,4,5)) +
    scale_fill_manual(values=site_colors) +
    scale_color_manual(values=site_colors) +
    guides(color=F, fill=F) +
    theme_bw()
ggsave("supplement_rawbyrecord.png", width=8, height=10.5, dpi=300)

# multiple methods
crossplot_method <- function(data_in, method_name) {
    binsize = 1
    lowlat_temp_lm  <- lm(SST ~ temp_benthic_binned, data = data_in %>% filter(preservation == "E", depth == 1, latzone == 1) %>% mutate(temp_benthic_binned = round(temp_benthic * (1/binsize)) / (1/binsize)) %>% group_by(temp_benthic_binned) %>% summarize(SST = median(SST)))
    midlat_temp_lm  <- lm(SST ~ temp_benthic_binned, data = data_in %>% filter(preservation %in% c("E","VG"), depth == 1, latzone == 2) %>% mutate(temp_benthic_binned = round(temp_benthic * (1/binsize)) / (1/binsize)) %>% group_by(temp_benthic_binned) %>% summarize(SST = median(SST)))
    highlat_temp_lm <- lm(SST ~ temp_benthic_binned, data = data_in %>% filter(preservation %in% c("E","VG","G","M","P"), depth == 1, latzone == 3, pallat < 0) %>% mutate(temp_benthic_binned = round(temp_benthic * (1/binsize)) / (1/binsize)) %>% group_by(temp_benthic_binned) %>% summarize(SST = median(SST)))
    timeseries_out <- d18O_benthic %>% rename(temp_benthic_binned = temp_benthic) %>% select(age, temp_benthic_binned)
    timeseries_out <- timeseries_out %>% add_predictions(lowlat_temp_lm, var="SST") %>% mutate(latzone = 1) %>%
        full_join(timeseries_out %>% add_predictions(midlat_temp_lm, var="SST") %>% mutate(latzone = 2)) %>%
        full_join(timeseries_out %>% add_predictions(highlat_temp_lm, var="SST") %>% mutate(latzone = 3)) %>%
        mutate(method = method_name)
    timeseries_out
}
crossplot_methods <- crossplot_method(d18O_stack %>% rename(SST = SST_Bayfox_50), "bayfox") %>% mutate(spatial = 1) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_HL), "Bemis et al. 1998 (HL)") %>% mutate(spatial = 1)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_LL), "Bemis et al. 1998 (LL)") %>% mutate(spatial = 1)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_mean), "Bemis et al. 1998 (mean)") %>% mutate(spatial = 1)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_inorg), "Kim & O'Neil 1997") %>% mutate(spatial = 1)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bayfox_no_spatial_50), "bayfox") %>% mutate(spatial = 0)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_no_spatial_HL), "Bemis et al. 1998 (HL)") %>% mutate(spatial = 0)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_no_spatial_LL), "Bemis et al. 1998 (LL)") %>% mutate(spatial = 0)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_no_spatial_mean), "Bemis et al. 1998 (mean)") %>% mutate(spatial = 0)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_inorg_no_spatial), "Kim & O'Neil 1997") %>% mutate(spatial = 0)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bayfox_Valdes_50), "bayfox") %>% mutate(spatial = 2)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_Valdes_HL), "Bemis et al. 1998 (HL)") %>% mutate(spatial = 2)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_Valdes_LL), "Bemis et al. 1998 (LL)") %>% mutate(spatial = 2)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_Bemis_Valdes_mean), "Bemis et al. 1998 (mean)") %>% mutate(spatial = 2)) %>%
    full_join(crossplot_method(d18O_stack %>% rename(SST = SST_inorg_Valdes), "Kim & O'Neil 1997") %>% mutate(spatial = 2))
ggplot(crossplot_methods, aes(x=age, y=SST, color=method)) +
    facet_grid(rows=vars(latzone), cols=vars(spatial), labeller=labeller(latzone=c(`1`="Tropical",`2`="Mid-Latitude",`3`="High-Latitude"), spatial=c(`0`="No local correction",`1`="CESM correction",`2`="Valdes correction"))) +
    geom_line() +
    labs(x="Age (Ma)",
         y="SST (°C)",
         color="Calibration") +
    theme_bw()
ggsave("supplement_calibrations.png", width=8, height=6, dpi=300)

# CESM vs. Valdes
ggplot(d18O_stack, aes(x=d18O_sw_CESM, y=d18O_sw_Valdes, color=as.factor(sign(pallat)))) +
    geom_point() +
    coord_cartesian(xlim=c(-2,2), ylim=c(-2,2)) +
    labs(x=bquote("CESM-inferred "*delta^18*O['sw']),
         y=bquote("Salinity-inferred "*delta^18*O['sw'])) +
    guides(color=F) +
    scale_color_manual(values=c("black", "gray50")) +
    theme_bw()
ggsave("supplement_cesmvaldes.png", width=4, height=4, dpi=150)

# CESM vs. HadCM3 vs. LeGrande & Schmidt
tindall_d18Osw <- read.csv("tindall_2010.csv", header=T)
d18Osw_crossplot <-
    legrande_d18Osw %>%
    mutate(lat = round(lat)) %>%
    group_by(lat) %>%
    summarize(d18O = mean(d18O, na.rm=T)) %>%
    mutate(era = "Modern") %>%
    full_join(cesm_miocene_280 %>%
                  mutate(lat = round(lat)) %>%
                  group_by(lat) %>%
                  summarize(d18O = mean(d18O, na.rm=T) - cesm_miocene_280_mean) %>%
                  mutate(era = "Miocene (CESM 280ppm)") %>%
                  select(lat, era, d18O)) %>%
    full_join(cesm_eocene_1x %>%
                  mutate(lat = round(lat)) %>%
                  group_by(lat) %>%
                  summarize(d18O = mean(d18O, na.rm=T) - cesm_eocene_1x_mean) %>%
                  mutate(era = "Eocene (CESM 285ppm)") %>%
                  select(lat, era, d18O)) %>%
    full_join(cesm_eocene_3x %>%
                  mutate(lat = round(lat)) %>%
                  group_by(lat) %>%
                  summarize(d18O = mean(d18O, na.rm=T) - cesm_eocene_3x_mean) %>%
                  mutate(era = "Eocene (CESM 855ppm)") %>%
                  select(lat, era, d18O)) %>%
    full_join(tindall_d18Osw %>%
                  select(latitude, eoc_surface_d18O_w) %>%
                  rename(lat = latitude, d18O = eoc_surface_d18O_w) %>%
                  mutate(d18O = d18O + 1, era = "Eocene (HadCM3)"))
ggplot(na.omit(d18Osw_crossplot), aes(x=lat, y=d18O, color=era)) +
    geom_line() +
    labs(x="Latitude",
         y=bquote("Relative "*delta^18*O['sw']),
         color="Dataset") +
    scale_x_continuous(limits=c(-90,90), breaks=c(-90,-60,-30,0,30,60,90)) +
    scale_color_manual(values=c("Modern"="black", "Miocene (CESM 280ppm)"="blue", "Eocene (CESM 285ppm)"="cornflowerblue", "Eocene (CESM 855ppm)"="orange", "Eocene (HadCM3)"="yellow2")) +
    theme_bw()
ggsave("supplement_latgradients.png", width=6, height=4, dpi=150)

# plot time steps on maps (uses Macrostrat paleogeography for reference)
epoch_ages = c("Lower Cretaceous"=(145+100.5)/2, "Upper Cretaceous"=(100.5+66)/2, "Paleocene"=(66+56)/2, "Early Eocene"=(56+48)/2, "Eocene"=(48+33.9)/2, "Oligocene"=(33.9+23.03)/2, "Miocene"=(23.03+5.333)/2, "Pliocene"=(5.333+2.588)/2, "Quaternary"=2.588/2)
for (map_age in na.omit(unique(d18O_stack$epoch))) {
    map <- downloadPaleogeography(Age=round(epoch_ages[map_age]))
    map_sites <- d18O_stack %>% filter(preservation %in% c("E","VG","G","M","P"), epoch == map_age) %>% mutate(pallat = round(pallat), pallong = round(pallong)) %>% select(pallat, pallong, ref) %>% distinct()
    print(ggplot() +
              geom_sf(data = map, color="gray75", fill="gray99") +
              geom_point(data = map_sites, aes(x=pallong, y=pallat)) +
              geom_text_repel(data = map_sites, aes(x=pallong, y=pallat, label=ref), size=1.5, max.overlaps=1000) +
              scale_y_continuous(breaks = c(-90, -60, -30, 30, 60, 90)) +
              labs(title=map_age,
                   x="Paleo-Longitude",
                   y="Paleo-Latitude") +
              theme_bw())
    ggsave(paste("d18O_SSTs_map_", map_age, ".svg", sep=""), width=6, height=4, dpi=150)
}

# raw CESM/HadCM3 output
plot_model_SST <- function(data_in, title_text) {
    ggplot(data_in, aes(x=lon, y=lat, fill=TEMP)) +
        geom_point(data=tibble(lat=c(0,0), lon=c(0,0), TEMP=c(-10,50)), size=0) + # hacky way to force scale to -10 to 45
        geom_raster() +
        scale_fill_gradientn(colors = c("#000000","#0000cc","#5cb2fd","#a9f2fb","#fde401","#fd5502","#FF0000","#ad0006","#570002"), values=rescale(c(-10,-2,5,15,20,28,32,40,50), to=c(0,1), from=c(-10,50)), breaks=c(-10,-2,5,15,22,28,32,40,45,50)) +
        geom_line(data=tibble(x=c(-180,180,-180,180), y=c(30,30,-30,-30), TEMP=c(1,1,2,2)), aes(x=x, y=y, group=TEMP)) +
        scale_y_continuous(breaks = c(-90, -60, -30, 30, 60, 90)) +
        scale_x_continuous(breaks = c(-180, -90, 0, 90, 180)) +
        coord_cartesian(xlim = c(-180, 180)) +
        labs(title=title_text,
             x="Longitude",
             y="Latitude",
             fill="SST (°C)") +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.key.width=unit(1.5,"cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
}
ggsave("supplement_cesm_eocene_1x_SST.png", plot=plot_model_SST(cesm_eocene_1x, "Eocene CESM SSTs (285 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_3x_SST.png", plot=plot_model_SST(cesm_eocene_3x, "Eocene CESM SSTs (855 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_6x_SST.png", plot=plot_model_SST(cesm_eocene_6x, "Eocene CESM SSTs (1710 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_9x_SST.png", plot=plot_model_SST(cesm_eocene_9x, "Eocene CESM SSTs (2565 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_280_SST.png", plot=plot_model_SST(cesm_miocene_280, "Miocene CESM SSTs (280 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_400_SST.png", plot=plot_model_SST(cesm_miocene_400, "Miocene CESM SSTs (400 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_840_SST.png", plot=plot_model_SST(cesm_miocene_840, "Miocene CESM SSTs (840 ppm CO2)"), width=6, height=5, dpi=150)

plot_model_d18O <- function(data_in, offset, title_text) {
    ggplot(data_in %>% mutate(d18O=d18O-offset), aes(x=lon, y=lat, fill=d18O)) +
        geom_point(data=tibble(lat=c(0,0), lon=c(0,0), d18O=c(-15,5)), size=0) + # hacky way to force scale to -10 to 45
        geom_raster() +
        scale_fill_gradientn(colors = c("black","blue","turquoise3","white","yellow","red","darkred"), values=rescale(c(-15,-2,-1,0,1,2,5), to=c(0,1), from=c(-15,5)), breaks=c(-15,-5,-2,-1,0,1,2,5)) +
        geom_line(data=tibble(x=c(-180,180,-180,180), y=c(30,30,-30,-30), d18O=c(1,1,2,2)), aes(x=x, y=y, group=d18O)) +
        scale_y_continuous(breaks = c(-90, -60, -30, 30, 60, 90)) +
        scale_x_continuous(breaks = c(-180, -90, 0, 90, 180)) +
        coord_cartesian(xlim = c(-180, 180)) +
        labs(title=title_text,
             x="Longitude",
             y="Latitude",
             fill=bquote("Relative "*delta^18*O['sw'])) +
        theme_bw() +
        theme(legend.position = "bottom",
              legend.key.width=unit(1.5,"cm"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
}
ggsave("supplement_cesm_eocene_1x_d18O.png", plot=plot_model_d18O(cesm_eocene_1x, cesm_eocene_1x_mean, "Eocene CESM d18O (285 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_3x_d18O.png", plot=plot_model_d18O(cesm_eocene_3x, cesm_eocene_3x_mean, "Eocene CESM d18O (855 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_6x_d18O.png", plot=plot_model_d18O(cesm_eocene_6x, cesm_eocene_6x_mean, "Eocene CESM d18O (1710 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_eocene_9x_d18O.png", plot=plot_model_d18O(cesm_eocene_9x, cesm_eocene_9x_mean, "Eocene CESM d18O (2565 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_280_d18O.png", plot=plot_model_d18O(cesm_miocene_280, cesm_miocene_280_mean, "Miocene CESM d18O (280 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_400_d18O.png", plot=plot_model_d18O(cesm_miocene_400, cesm_miocene_400_mean, "Miocene CESM d18O (400 ppm CO2)"), width=6, height=5, dpi=150)
ggsave("supplement_cesm_miocene_840_d18O.png", plot=plot_model_d18O(cesm_miocene_840, cesm_miocene_840_mean, "Miocene CESM d18O (840 ppm CO2)"), width=6, height=5, dpi=150)

# supplementary tables
write.csv(d18O_benthic %>%
              select(age, d18O, d18O_sw, temp_benthic, SST_montecarlo_1_2.5, SST_montecarlo_1_25, SST_montecarlo_1_50, SST_montecarlo_1_75, SST_montecarlo_1_97.5, SST_montecarlo_3_2.5, SST_montecarlo_3_25, SST_montecarlo_3_50, SST_montecarlo_3_75, SST_montecarlo_3_97.5) %>%
              rename(SST_trop_0.025 = SST_montecarlo_1_2.5,
                     SST_trop_0.25  = SST_montecarlo_1_25,
                     SST_trop_0.50  = SST_montecarlo_1_50,
                     SST_trop_0.75  = SST_montecarlo_1_75,
                     SST_trop_0.975 = SST_montecarlo_1_97.5,
                     SST_high_0.025 = SST_montecarlo_3_2.5,
                     SST_high_0.25  = SST_montecarlo_3_25,
                     SST_high_0.50  = SST_montecarlo_3_50,
                     SST_high_0.75  = SST_montecarlo_3_75,
                     SST_high_0.975 = SST_montecarlo_3_97.5),
          "supplement_benthic.csv")

d18O_stack$d18O_sw_source <- "CESM (Eocene paleogeography)"
d18O_stack[which(d18O_stack$epoch %in% c("Holocene", "Pliocene")), "d18O_sw_source"] <- "LeGrande & Schmidt (2006)"
d18O_stack[which(d18O_stack$epoch == "Miocene"), "d18O_sw_source"] <- "CESM (Miocene paleogeography)"
write.csv(d18O_stack %>% 
              filter(depth %in% c(1,2,3), preservation %in% c("E","VG","G","M","P")) %>%
              mutate(sst_2.5 = SST_montecarlo_2.5, sst_97.5 = SST_montecarlo_97.5, sst = SST_Bayfox_50, CO3 = CO3_zeebe, d18O_sw = d18O_sw_CESM) %>%
              select(fullsp, sp_raw, core, mbsf, mcd, rmcd, lat, long, pallat, pallong, latzone, preservation, d13C, d18O, depth, symbionts, ecology_ref, zone_p, zone_n, zone_m, age_raw, age_raw_scale, age, age_error_plus, age_error_minus, age_scale, ref, notes,
                     CO3, d18O_CO3, temp_benthic, d18O_ice, d18O_sw, d18O_sw_source, sst_2.5, sst, sst_97.5),
          "supplement_planktic.csv")

# calculate PAFs resulting from data in other papers
slopes_comparisons <- read.csv("slope_comparisons.csv", header=T)
slopes_comparisons %>%
    group_by(Proxy) %>%
    mutate(PAF = lm(High ~ Mean)$coefficients[2],
           slope = lm(Gradient ~ Mean)$coefficients[2],
           intercept = lm(Gradient ~ Mean)$coefficients[1],
           error = confint(lm(High ~ Mean))[4] - PAF) %>%
    summarize(PAF = median(PAF),
              slope = median(slope),
              intercept = median(intercept),
              error = median(error),
              min_gmsst = min(Mean),
              max_gmsst = max(Mean))