# code to load and wrangle data
library(tidyverse)
library(readxl)
source("helper_funs.R")

# import data
filename <- "CleanedOxygenThresholds.xlsx"
thedata <- read_excel(filename)

# create new column called p_crit_kpa, which has all p_crit measurements in the same units (kPa)
# drop all rows that do not list units
#source("helper_funs.R") # loads the functions needed to convert ml /l to kPa

thedata <- drop_na(thedata, units)
torr.2.kpa <- 0.13332
thedata$p_crit_kpa <- NA
thedata$p_crit_kpa[thedata$units == "kPa (PO2)"] <- thedata$p_crit[thedata$units == "kPa (PO2)"]
thedata$p_crit_kpa[thedata$units == "Torr"] <- torr.2.kpa * thedata$p_crit[thedata$units == "Torr"]
thedata$p_crit_kpa[thedata$units == "mm Hg"] <- torr.2.kpa * thedata$p_crit[thedata$units == "mm Hg"]


# Set salinity valuues
thedata$salinity[thedata$fresh_or_salt=="fresh"] <- 0
thedata$salinity[thedata$fresh_or_salt=="salt"] <- 34.9

# get rows with different concentration units
mg_l_index <- which(thedata$units == "mg/L")
ml_l_index <- which(thedata$units == "ml/L")
pct_index <- which(thedata$units == "% saturation")


# convert to pKa
thedata$p_crit_kpa[mg_l_index] <- DO.unit.convert(x = thedata$p_crit[mg_l_index],
                                                  DO.units.in = "mg/L",
                                                  DO.units.out = "PP",
                                                  bar.units.out = "kpa",
                                                  bar.press = 1,
                                                  bar.units.in = "atm",
                                                  salinity = thedata$salinity[mg_l_index],
                                                  salinity.units = "uS",
                                                  temp.C = thedata$temperature_c[mg_l_index])

thedata$p_crit_kpa[ml_l_index] <- DO.unit.convert( x = 1.429 * thedata$p_crit[ml_l_index],
                                                   DO.units.in = "mg/L",
                                                   DO.units.out = "PP",
                                                   bar.units.out = "kpa",
                                                   bar.press = 1,
                                                   bar.units.in = "atm",
                                                   salinity = thedata$salinity[ml_l_index],
                                                   salinity.units = "uS",
                                                   temp.C = thedata$temperature_c[ml_l_index])

thedata$p_crit_kpa[pct_index] <- DO.unit.convert( x = thedata$p_crit[pct_index],
                                                  DO.units.in = "pct",
                                                  DO.units.out = "PP",
                                                  bar.units.out = "kpa",
                                                  bar.press = 1,
                                                  bar.units.in = "atm",
                                                  salinity = thedata$salinity[pct_index],
                                                  salinity.units = "uS",
                                                  temp.C = thedata$temperature_c[pct_index])

# fit linear model
