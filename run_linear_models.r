# Code to run linear mixed effects models and plot predictions

library(tidyverse)
library(lme4)
library(readxl)

source("load_data.R")

# run model for salmon
# step 1. filter data for Oncorhynchus and na for pcrit and method
salmon_data <- thedata %>%
  filter(genus == "Oncorhynchus") %>%
  drop_na(p_crit_kpa) %>%
  drop_na(method) %>%
  rename(species_name = `Species Name`)

# step 2. remove early life history stages and methods that we don't want
salmon_data <- salmon_data %>%
  filter(
    !tolower(life_stage_est) %in% c("egg", "early developmental stage") &
      !(method) %in% c("CardiacActivityControl", "BloodChemistry", "Exhaustion")
  )

# step 3: use mean of mass range when mean mass not available
salmon_data <- salmon_data %>%
  mutate(mean_mass_g_wet_adj  = case_when(
    is.na(mean_mass_g_wet) ~ mean(c(mass_range_upper, mass_range_lower) ),
    !is.na(mean_mass_g_wet) ~ mean_mass_g_wet)
  )

# setp 4: run linear model
model_fit <- lmer(log(p_crit_kpa) ~ log(mean_mass_g_wet_adj) + EstMethod_Metric + inverse_temperature + (1 |species_name) ,
                  data = salmon_data
                  ) 

summary(model_fit)
predict_df <- data.frame(mean_mass_g_wet_adj = with(salmon_data, mean(mean_mass_g_wet_adj, na.rm = T) ),
                     species_name = "O. mykiss",
                     inverse_temperature = with(salmon_data,seq(range(inverse_temperature) [1],
                                                                range(inverse_temperature) [2],
                                                                length.out = 10)
                                                ),
                     EstMethod_Metric = "
                     )
predict_df$log_pcrit_kpa <- predict(object = model_fit,
                             newdata = predict_df)

# plot data with best fitting line
ggplot(data = salmon_data, aes(x = inverse_temperature, y = log(p_crit_kpa), col = species_name)) +
  geom_point() +
  geom_line(data = predict_df, aes(x = inverse_temperature, y = log_pcrit_kpa),  color = "black")
         