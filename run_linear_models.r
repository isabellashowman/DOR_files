# Code to run linear mixed effects models and plot predictions

library(tidyverse)
library(lme4)
library(readxl)
library(scales)
# Load data ####
source("load_data.R")

# run model for salmon
# step 1. filter data for Oncorhynchus and na for pcrit and method
kb <- 8.617333262E-5
tref <- 12
kelvin <- function(x) x + 273.15
make_inverse_temp <- function(x) 1 / kb * ( 1 / kelvin(x) - 1/ kelvin(tref) )

# Fit model to Salmon data ####
## Extract daata #####
salmon_data <- thedata %>%
  filter(genus == "Oncorhynchus") %>%
  drop_na(p_crit_kpa) %>%
  drop_na(method) %>%
  rename(species_name = `Species Name`,
         body_mass = med_mass_g_wet) %>%
  filter(!is.na(body_mass),
         EstMethod_Metric %in% c("Lethal/LOE", "Oxyconform", "SMR")) |>
  mutate(inverse_temp = make_inverse_temp(temperature_c),
         species_name = as.factor(species_name),
         EstMethod_Metric = factor(EstMethod_Metric, levels = c("Oxyconform", 
                                                                "SMR", 
                                                                "Lethal/LOE"))
         )
  

## step 2. remove early life history stages and methods that we don't want ####
salmon_data <- salmon_data %>%
  filter(
    !tolower(life_stage_est) %in% c("egg", "early developmental stage")
  )


## Check for NAs by species ####
salmon_data |>
  group_by(species_name) |>
  summarise(
    n = n(),
    na_mass    = sum(is.na(body_mass)),
    na_temp    = sum(is.na(inverse_temp)),
    na_pcrit   = sum(is.na(p_crit_kpa)),
    na_method  = sum(is.na(EstMethod_Metric))
  )

# setp 4: run linear model with no fixed effects
model_fit <- lm(log(p_crit_kpa) ~ log(body_mass) + EstMethod_Metric + inverse_temp + species_name ,
                     data = salmon_data)

summary(model_fit)
predict_df <- salmon_data |>
  filter(species_name == "O. mykiss") |>
  distinct(species_name, EstMethod_Metric) |>
  cross_join(
    data.frame(
      inverse_temp = with(salmon_data, seq(range(inverse_temp)[1],
                                           range(inverse_temp)[2],
                                           length.out = 10)),
      body_mass = with(salmon_data, mean(body_mass, na.rm = TRUE))
    )
  )

predict_df$log_pcrit_kpa <- predict(object = model_fit,
                             newdata = predict_df)

## plot data with best fitting line ####
salmon_plot <- ggplot(data = salmon_data, aes(x = inverse_temp, y = log(p_crit_kpa), 
                               col = EstMethod_Metric, 
                               shape = species_name)) +
  geom_point( size = 3) +
  geom_line(data = predict_df, aes(x = inverse_temp, y = log_pcrit_kpa, 
                                   col = EstMethod_Metric), linewidth = 1) +
  scale_color_viridis_d(option = "plasma", begin = 0.0, end = 0.8)  +
  theme_minimal() +
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),  # add this
    axis.ticks  = element_line(color = "black")   # and this for tick marks
  ) +
  labs(x        = expression( frac(1, k[b]) ~ (frac(1, T) - frac(1, T[ref]) )),
       y        = expression(log(p[crit]) ~ (kPa)),
       color    = "Method",
       shape    = "Species",
       linetype = "Species") +
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")),
         shape    = guide_legend(keywidth = unit(1.5, "cm")))

print(salmon_plot)
# Repeat for Euphasiids ####
## Extract Data ####
euphausiid_data <-  thedata %>%
  filter(family == "Euphausiidae") %>%
  drop_na(p_crit_kpa) %>%
  drop_na(method) %>%
  rename(species_name = `Species Name`,
         body_mass = med_mass_g_wet) %>%
  filter(!is.na(body_mass),
         EstMethod_Metric %in% c("Lethal/LOE", "Oxyconform", "SMR")) |>
  mutate(inverse_temp = make_inverse_temp(temperature_c),
         species_name = as.factor(species_name),
         EstMethod_Metric = factor(EstMethod_Metric, levels = c("Oxyconform",
                                                                "SMR",
                                                                "Lethal/LOE"))
  )

## Plot Euphausiid Data
ggplot(data = euphausiid_data, aes(x = inverse_temp, y = log(p_crit_kpa), 
                               col = EstMethod_Metric, 
                               shape = species_name)) +
  geom_point( size = 3) +
  #geom_line(data = predict_df, aes(x = inverse_temp, y = log_pcrit_kpa, 
  #                                 col = EstMethod_Metric), linewidth = 1) +
  scale_color_viridis_d(option = "plasma", begin = 0.0, end = 0.8)  +
  theme_minimal() +
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),  # add this
    axis.ticks  = element_line(color = "black")   # and this for tick marks
  ) +
  labs(x        = expression(frac(1, T) - frac(1, T[ref])),
       y        = expression(log(p[crit]) ~ (kPa)),
       color    = "Method",
       shape    = "Genus",
       linetype = "Genus") +
  guides(linetype = guide_legend(keywidth = unit(1.5, "cm")),
         shape    = guide_legend(keywidth = unit(1.5, "cm")))
## Fit lm ####
euphausiid_model <- lm(log(p_crit_kpa) ~ log(body_mass) + EstMethod_Metric + inverse_temp + species_name,
                         data = euphausiid_data)

summary(euphausiid_model)


# Load hood canal data ####
hood_canal_raw <- read.csv("orca2_L1_profiles_cf4a_c7c1_8488.csv")

hood_canal <- hood_canal_raw[-1, ] |>
  # First convert characacters to numbers)
  mutate(sea_water_temperature = as.numeric(sea_water_temperature),
         mass_concentration_of_oxygen_in_sea_water = as.numeric(mass_concentration_of_oxygen_in_sea_water),
         depth = as.numeric(depth), 
         sea_water_practical_salinity  = as.numeric(sea_water_practical_salinity)) |>
  # remove missing values
  filter(!is.nan(mass_concentration_of_oxygen_in_sea_water)) |>
  # compute pO2
  mutate(oxygen_kpa = DO.unit.convert(x = mass_concentration_of_oxygen_in_sea_water,
                                      DO.units.in = "mg/L",
                                      DO.units.out = "PP",
                                      bar.units.out = "kpa",
                                      bar.press = 1,
                                      bar.units.in = "atm",
                                      salinity = sea_water_practical_salinity,
                                      salinity.units = "uS",
                                      temp.C = sea_water_temperature)
         ) |>
  # modify variables
  mutate(inverse_temp = make_inverse_temp(sea_water_temperature),
         log_kpa = log(oxygen_kpa),
         time = as.POSIXct(time, format = "%Y-%m-%dT%H:%M:%OSZ", tz = "UTC"),
         date = as.Date(time)) |>
  # remove portions of cast when rising to surface and when warming up
  group_by(cast_number) |>
  arrange(sample_time, .by_group = TRUE) |>
  filter(sample_time >= sample_time[which.min(depth)],
         sample_time <= sample_time[which.max(depth)]) |>
  ungroup()



## Extract salmon model coefficients
coefs = coef(model_fit)
salmon_intercept = coefs["(Intercept)"] + coefs["species_nameO. mykiss"]
salmon_slope = coefs["inverse_temp"]

# make new prediction data frame
predict_df <- data.frame(sea_water_temperature = seq(7, 18, length.out = 20),
                         log_kpa = numeric(20))

predict_df <- predict_df |>
  mutate(inverse_temp = make_inverse_temp(sea_water_temperature),
         log_kpa = salmon_intercept + salmon_slope * inverse_temp)


euphausiid_coefs = coef(euphausiid_model)
euph_intercept = euphausiid_coefs["(Intercept)"]
euph_slope = euphausiid_coefs["inverse_temp"]


# Plot hood canal by day [two casts per day]
hood_canal_plot <- hood_canal |>
  filter(depth <120) |>
  ggplot(aes(x = sea_water_temperature, y = log_kpa, color = depth)) + 
  geom_point() + 
  scale_color_viridis_c(option = "plasma", name = "Depth (m)", guide = guide_colorbar(reverse=TRUE)) +
  facet_wrap(vars(date)) + 
  labs(
    #x = expression(frac(1, T) - frac(1, T[ref])),
    x = "Temperature (C)",
    y = "log(Oxygen Partial Pressure) (kPa)"
  ) +
  theme_minimal() +
  theme(
    axis.title  = element_text(size = 14),
    axis.text   = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 14),
    strip.text  = element_text(size = 12),
    panel.grid  = element_blank(),
    axis.line   = element_line(color = "black"),  # add this
    axis.ticks  = element_line(color = "black")   # and this for tick marks
  ) + 
  geom_line(data = predict_df, aes(x = sea_water_temperature, y = log_kpa), inherit.aes = FALSE)
  
plot(hood_canal_plot)

# calculate pcrit for each data point

calc_pcrit <- function(inv_temp) exp( salmon_intercept + salmon_slope * inv_temp )
calc_pcrit_surface <- function(inv_temp) exp( salmon_intercept + salmon_slope * inv_temp )

depth_at_pcrit <- hood_canal |>
  mutate(pcrit = calc_pcrit(inverse_temp)) |>
  mutate(under_pcrit = case_when(
    exp(log_kpa) < pcrit ~ 1,
    exp(log_kpa) >= pcrit ~ 0)
  ) |>
  filter(under_pcrit == 1) |>
  group_by(cast_number) |>
  summarize(min(depth))
  
print(depth_at_pcrit)

# calculate depth at which pO2 < pcrit if pcrit is based on surface 
# [warmest] temperature

hood_canal |>
  group_by(cast_number) |>
  mutate(min_inverse_temp = min(inverse_temp)) |>
  ungroup() |>
  mutate(pcrit_surface = calc_pcrit(min_inverse_temp) ) |>
  mutate(under_pcrit = case_when(
    exp(log_kpa) < pcrit_surface ~ 1,
    exp(log_kpa) >= pcrit_surface ~ 0)
  ) |>
  filter(under_pcrit == 1) |>
  group_by(cast_number) |>
  summarize(min(depth))


# make a typical plot of depth vs. temperature
ggplot(data = hood_canal, aes(x = sea_water_temperature, y = -depth)) + 
  geom_point()
ggplot(data = hood_canal, aes(x = exp(log_kpa), y = -depth)) + 
  geom_point()

# get ratio of pcrit to po2
hood_canal |>
  filter(depth < 20) |>
  mutate(pcrit = calc_pcrit(inverse_temp),
         phi = exp(log_kpa) / pcrit) |>
  ggplot(aes(x = phi)) +
  geom_histogram()
  
