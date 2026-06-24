# Allometric scaling

# This script is intended for allometric scaling of CTmax by fish length.
# I will do this by population, and not by year. I will provide the justification for that here, as well.

# Length data is unequally distributed, so I cannot use a least squares linear regression, which requires even distribution of data.
# Instead, I will do two transormations: 
## 1) Log
## 2) Weighted least squares

# Load in libraries ----------------------------------------------
lapply(c("tidyverse", "lme4", "ggthemes","car","emmeans","brglm2", "ggthemes", 
         "rethinking", "usethis", "here", "mgcv", "gratia", "extrafont","rstatix"), 
       require, character.only=T)

# Load in data --------------------------------------------------------
comp <- read.csv(here("data/IRA.thermaltolerance.compileddata.csv"))
head(comp)

## Wrangling -----------------------------------------------------------
comp_clean = comp %>% 
  rename(Watershed = Location)
unique(comp_clean$Watershed)

# rename locations to match the watershed name
comp_clean = comp_clean %>% 
  mutate(Watershed = case_when(
    Watershed == "Big Sur" ~ "Big Creek", # Big Creek in Big Sur
    Watershed == "Napa" ~ "Napa River", # Napa River
    Watershed == "E Austin" ~ "Russian River", # Russian River
    TRUE~ Watershed # keep Scott Creek and Carmel River the same
  ))
str(comp_clean)

comp_clean$Date = as.Date(comp$Date, format = "%m/%d/%Y") #convert date to date structure
# Treatment should be a factor
comp_clean$Treatment = as.factor(comp_clean$Treatment)

# check for NA values
sum(is.na(comp_clean))
# there is one NA value in the length column

# add a year column
comp_clean <- comp_clean %>% 
  mutate(Year = year(Date))
comp_clean$Year = as.factor(comp_clean$Year)
comp_clean$Mass_g <- as.numeric(comp_clean$Mass_g)

# Add a LocationYear column
comp_clean <- comp_clean %>% 
  unite(WatershedYear, Watershed, Year, sep = " ", remove= FALSE)

## Filtered dataset ----------------------------------------------------------------------------------------------
treat2 = comp_clean %>% 
  filter(Treatment == 2) %>% 
  filter(Watershed %in% c("Big Creek", "Carmel River", "Scott Creek"))

treat2ALL = comp_clean %>% 
  filter(Treatment == 2)

## Color palette ----------------------------------------------------------------------------------------------
fig1_palette= c(
  "Big Creek"="#744253",
  "Scott Creek"="#496A81",
  "Carmel River"="#9BC53D")

fig2_palette = c(
  "Big Creek 2024"="#C89F9C",
  "Big Creek 2025" = "#B36A5E",
  "Scott Creek 2025"="#083d77",
  "Scott Creek 2024"="#76949f",
  "Carmel River 2024" = "#9BC53D",
  "Carmel River 2025" = "chartreuse4"
)

fig3_palette= c(
  "Big Creek"="#744253",
  "Scott Creek"="#496A81",
  "Carmel River"="#9BC53D",
  "Napa River"="blue",
  "Russian River" = "purple")

# Plot length distributions -------------------------------
colnames(treat2)
size.dist = ggplot(treat2, aes(x = Length_mm, fill = WatershedYear))+
  geom_histogram(binwidth = 10, linewidth = .2, color = "black")+
  scale_fill_manual(values = fig2_palette)+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0)),
    breaks = c(0,2,4,6,8,10)
  )+
  labs(x = "Length (mm)",
       y = NULL)+
  theme_bw()+
 theme(text = element_text(family = "Times New Roman", size = 10),
       legend.position="none")+
  facet_grid(Watershed~Year)
size.dist
?geom_histogram
ggsave("output/sizedist.png", plot = size.dist, height = 3, width = 2.5, units = "in")

# any obvious outliers?
outliercheck = treat2 %>% 
  group_by(Watershed, Year) %>%
  identify_outliers(Length_mm)

## CTMax vs length plot --------------------------------------------------------
ctvslength = ggplot(treat2, aes(x= Length_mm, y = LOE.Temp_C, color = Watershed))+
  geom_point() +
  scale_color_manual(values = fig1_palette)+
  facet_wrap(~WatershedYear)
ctvslength

ctvslength2 = ggplot(treat2ALL, aes(x= Length_mm, y = LOE.Temp_C, color = Watershed))+
  geom_point() +
  scale_color_manual(values = fig3_palette)+
  geom_smooth(method=lm)
  #facet_wrap(~WatershedYear)
ctvslength2

?identify_outliers
## Log transform length plot ----------------------------------------------------
# does this solve heteroscedasticity problem?
# make a log variable
treat2 = treat2 %>% 
  mutate(log_length = log(Length_mm),
         log_loe = log(LOE.Temp_C))

outliercheck.log = treat2 %>% 
  group_by(Watershed, Year) %>%
  identify_outliers(log_length)

ctvslog.length = ggplot(treat2, aes(x= log_length, y = LOE.Temp_C, color = Watershed))+
  geom_point() +
  scale_color_manual(values = fig1_palette)+
  facet_wrap(~WatershedYear)
ctvslog.length

log.ctvslog.length = ggplot(treat2, aes(x= log_length, y = log_loe, color = Watershed))+
  geom_point() +
  scale_color_manual(values = fig1_palette)+
  facet_wrap(~WatershedYear)
log.ctvslog.length
# to do allometric scaling, you have to log transform both x and y according to standard practice in allometry



# after checking for outliers in both raw and log-transformed data, it is clear there are still length outliers. 
# to solve that problem, I am going to use a weighted linear regression for least squares

# Baseline regression: CTmax vs length --------------------------------------------------------------

baseline = ggplot(treat2, aes(y= log_loe, x = log_length,color = Watershed))+
  geom_point()+
  scale_color_manual(values = fig1_palette)+
  geom_smooth(method = lm)
baseline

basemod <- lm(log_loe ~ log_length*Watershed, data = treat2)
summary(basemod)

# before going through all the effort of weighting, maybe Huber Regression works well?
# Huber Regression ------------------------------------------------------------------------
library(MASS)
hubermod <- MASS::rlm(log_loe ~ log_length*Watershed, data = treat2)
summary(hubermod)
coef(hubermod)

# redo huber mod by watershed and year
hubermod_year <- MASS::rlm(log_loe ~ log_length*WatershedYear, data = treat2)
summary(hubermod_year)
coef(hubermod)

hubermod_all <- MASS::rlm(log_loe ~ log_length, data = treat2)
summary(hubermod_year)
coef(hubermod)

# create lines to plot on graph
line_data <- data.frame(
  line_id = c("Big Creek", "Carmel River", "Scott Creek"),
  intercepts = c(3.4797, 3.5025, 3.4164),
  slopes = c(-0.0101, -0.0141, 0.0046)
)

line_data_year <- data.frame(
  line_id = c("Big Creek 2024", "Big Creek 2025","Carmel River 2024","Carmel River 2025", "Scott Creek 2024","Scott Creek 2025"),
  intercepts = c(3.547, 3.4227,3.4436,3.4404,3.3752,3.4505),
  slopes = c(-0.0241,0.0018)
)

treat2$predicted <- predict(hubermod_year)
treat2$weights <- hubermod_year$w #extract weight

# outlier plot- which points were adjusted by the Huber robust regression line?
ggplot(treat2, aes(x = log_length, y = log_loe)) +
  geom_point(aes(color = weights)) + # Colors points based on how the model weighs them
 # geom_line(aes(y = predicted), color = "blue", size = 1) + # The robust line
  scale_color_gradient(low = "red", high = "black") + # Highlights outliers in red
  theme_minimal() +
  labs(title = "Huber Robust Regression", subtitle = "Red points indicate down-weighted outliers")


# Chat suggested this as a way to visualize the model output more clearly:
ggplot(treat2, aes(x=log_length, y=log_loe, color = WatershedYear)) +
  geom_point() +
  geom_line(aes(y = predicted), linewidth = 1.2) +
  scale_color_manual(values= fig2_palette)+
  theme_minimal()+
  labs(title = "Size regression across locations (Huber robust)")

# Huber regression slightly changed the model output, but barely. However, it is still more robust, so I will continue with this method

ggplot(treat2, aes(y= log_loe, x = log_length,color = Watershed))+
  geom_point(size = 2, alpha = 0.75)+
  scale_color_manual(values = fig1_palette)+
  geom_abline(data = line_data, 
              aes(intercept = intercepts, slope = slopes, color = line_id),
              linewidth = 1)+
  labs(x= "Log Length (mm)",
       y = "Log LOE Temp (°C)",
       title = "Huber Robust Regression on Log Transformed X and Y")+
  theme_bw()


# how do I pick what scale to do the huber model slopes at ? Chat suggests this:
hubermod_test <- MASS::rlm(log_loe ~ log_length + Watershed + Year, data = treat2)
summary(hubermod_test)


hubermod_int <- rlm(
  LOE.Temp_C~ Length_mm * Watershed + Length_mm* Year,
  data = treat2
) # tried without logging just to see what happened

summary(hubermod_int)

coef(hubermod)

coef = tidy(hubermod)

ggplot(treat2,
       aes(log_length,
           log_loe,
           color = Watershed)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE)

scalers <- coef %>% 
  mutate(term = case_when(
    term == "(Intercept)" ~ "WatershedBig Creek", 
    term == "log_length" ~ "log_length:WatershedBig Creek",
    TRUE~ term # keep Scott Creek and Carmel River the same
  ))

# first, filter dataset into two. That will be the easiest for combining, i think
pivot <- scalers %>% 
  filter(term %in% c("log_length:WatershedBig Creek","log_length:WatershedScott Creek","log_length:WatershedCarmel River"))
scalers <- scalers %>% 
  filter(term %in% c("WatershedBig Creek","WatershedScott Creek","WatershedCarmel River"))

scalers <- scalers %>% 
  mutate(term = case_when(
    term == "WatershedBig Creek" ~ "Big Creek", 
    term == "WatershedScott Creek" ~ "Scott Creek",
    term =="WatershedCarmel River"~ "Carmel River",
    TRUE~ term # keep Scott Creek and Carmel River the same
  ))
scalers = scalers %>% 
  rename(Watershed = term,
         intercept = estimate)

# don't forget  you still have to add the estimate to the reference (Big Creek)!


pivot <- pivot %>% 
  mutate(term = case_when(
    term == "log_length:WatershedBig Creek" ~ "Big Creek", 
    term == "log_length:WatershedScott Creek" ~ "Scott Creek",
    term =="log_length:WatershedCarmel River"~ "Carmel River",
    TRUE~ term # keep Scott Creek and Carmel River the same
  )) %>% 
  rename(Watershed = term,
         slope = estimate)
# don't forget  you still have to add the estimate to the reference (Big Creek)!

intercepts_corrected <- scalers %>% # chatGPT helped with this code
  mutate(
    ref_est = intercept[Watershed == "Big Creek"],
    adjusted_int = case_when(
      Watershed == "Big Creek" ~ ref_est,
      TRUE ~ ref_est + intercept
    )
  ) %>% 
  dplyr::select(Watershed,statistic, std.error, adjusted_int) %>% 
  rename(statistic.int = statistic,
         std.error.int = std.error)

slopes_corrected <- pivot %>% # chatGPT helped with this code
  mutate(
    ref_est = slope[Watershed == "Big Creek"],
    adjusted_slope = case_when(
      Watershed == "Big Creek" ~ ref_est,
      TRUE ~ ref_est + slope
    )
  ) %>% 
  dplyr::select(Watershed,statistic, std.error, adjusted_slope) %>% 
  rename(statistic.slope = statistic,
         std.error.slope = std.error)

## merged allometric scaling data frame -----------------------------------------------
str(intercepts_corrected)
slopes_intercepts <- merge(slopes_corrected, intercepts_corrected, by = "Watershed")

# Corrected CTmax ------------------------------------------------------------------------
# Based on Eric's paper, I need to use the following equation:
# CTMax(target) = CTMax*(Fork length(grand mean across all locations)/FL observed)^b
# where b is the slope for each location, and all data is untransformed

### a) merge data frame above with main LOE data frame ------------------------------------
ctmax_corrected = merge(treat2, slopes_intercepts, by = "Watershed") %>% 
  mutate(length_target = mean(Length_mm)) %>% # grand mean of fork lengths
  group_by(Watershed) %>% 
  mutate(LOE_Target = LOE.Temp_C*(Length_mm/length_target)^adjusted_slope) # create a new column for LOE_Target

### b) plot observed vs target ---------------------------------------------------------
targetvsobserved_LOE <- ggplot(ctmax_corrected, aes(x=LOE.Temp_C,y= LOE_Target, color= Watershed))+
  geom_point(size=2, alpha = 0.75)+
  scale_color_manual(values= fig1_palette)+
  geom_abline(intercept = 0, slope=1, linetype="dashed")+ #1:1 line for vis
  theme_bw()+
  facet_grid(Watershed~Year)

# the 1:1 line reveals some deviation in Carmel River 2025 in the transformation. Ask Eric!

### c) Plot violin plot with corrected CTmax values
p <- ctmax_corrected %>% 
  ggplot(aes(x = as.factor(Year), 
             y = LOE_Target, fill = WatershedYear, shape = Year)) +
  
  geom_violin(alpha = 0.75)+
  geom_jitter(width = .05,
              size = 2.5,
              alpha = 0.75) +
  
  scale_fill_manual(values = fig2_palette) +
  
  labs(
    x = "Watershed",
    y = "Allometric Scaled LOE Temp (°C)"
  ) +
  
  theme_bw() +
  theme(
    text = element_text(family = "Times New Roman", size = 10),
    axis.title = element_text(size = 10),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10),
    strip.text = element_text(size = 10),
    panel.spacing = unit(0.4, "lines"),
    legend.position = "none"
  )+
  guides(
    color = guide_legend(nrow = 2, byrow = TRUE,
                         keywidth = 0.4,
                         keyheight = 0.4)
  )+
  facet_wrap(~Watershed);p

ggsave("output/allometricallyscaledctmax.png", plot = p, height = 5, width = 5, units = "in")

##################################################################################
# Stream Temperature Data ---------------------------------------------
#######################################################################################

# i need to pull in stream temperature 
recentst = read.csv(here("data/stream7d_sum.ch1.csv")) 

# ok so. 7DADM is going to go on the x-axis, with error bars for standard deviation.I will also plot with 7DADA on the x-axis.
# on the y-axis will be the average with standard error (or should it be standard deviation? ask eric)

# calculate standard error for CTmax target data

ctmax_corrected <- ctmax_corrected %>% 
  group_by(WatershedYear) %>% # group by watershed year
  mutate(stdev_ctmaxtarget = sd(LOE_Target), # calculate standard deviation
         sterror_ctmaxtarget = sd(LOE_Target, na.rm = TRUE)/sqrt(sum(!is.na(LOE_Target)))) %>%  # calculate standard error
  ungroup()

# merge ctmax and stream temp data
ctmax_st <- merge(comp_clean, recentst, by = "WatershedYear")
colnames(ctmax_st)

# Target Plot ----------------------------------------------------------------------------
# help from Chat
ctmax_summary <- ctmax_st %>%
  group_by(WatershedYear) %>%
  mutate(
    n = sum(!is.na(LOE.Temp_C)),
    mean_LOE = mean(LOE.Temp_C, na.rm = TRUE),
    sd_LOE = sd(LOE.Temp_C, na.rm = TRUE),
    se_LOE = sd_LOE / sqrt(n),
    .groups = "drop"
  )

ctmax_summary
 
 fig2_palette = c(
  "Big Creek 2024"="#F49CBB",
  "Big Creek 2025" = "#A44A3F",
  "Scott Creek 2025"="#467599",
  "Scott Creek 2024"="#9ED8DB",
  "Carmel River 2024" = "#98CE00",
  "Carmel River 2025" = "#618B25"
) 
 
####################################################################
 # Target plots ----------------------------------------------------
######################################################################
 
 # Target plot 7DADM --------------------------------------------

targetplot <- ggplot(ctmax_summary,
       aes(x = DADM7, y = mean_LOE, color= WatershedYear)) +
  geom_point(size=4.5) +
  scale_color_manual(values= fig2_palette)+
  # add vertical error bars
  geom_errorbar(
    aes(
      ymin = mean_LOE - se_LOE,
      ymax = mean_LOE + se_LOE
    ),
    linewidth = .75)+
  # add horizontal error bars:
      geom_segment(
        aes(
          x = DADM7 - sd_DADM7,
          xend = DADM7 + sd_DADM7,
          y = mean_LOE,
          yend = mean_LOE
        ),
        linewidth = .75
      )  +
  labs(x= "7 Day Average Daily Maximum (°C)",
       y= "NOT Allometrically scaled LOE Temperature (°C)")+
 # coord_equal()+ # this scales the x- and y-axes. But also compresses the plot so that it's hard to read. Not sure what the best approach is.
  theme_bw();targetplot

ggsave("output/targetplot.png", plot = targetplot, height = 5, width = 5, units = "in")


## Target plot with 7DADA ------------------------------------------------------


targetplot_avg <- ggplot(ctmax_summary,
                     aes(x = DADA7, y = mean_LOE, color= WatershedYear)) +
  geom_point(size=4.5) +
  scale_color_manual(values= fig2_palette)+
  # add vertical error bars
  geom_errorbar(
    aes(
      ymin = mean_LOE - se_LOE,
      ymax = mean_LOE + se_LOE
    ),
    linewidth = .75)+
  # add horizontal error bars:
  geom_segment(
    aes(
      x = DADA7 - sd_DADA7,
      xend = DADA7 + sd_DADA7,
      y = mean_LOE,
      yend = mean_LOE
    ),
    linewidth = .75
  )  +
  labs(x= "7 Day Average Daily Average (°C)",
       y= "NOT Allometrically scaled LOE Temperature (°C)")+
  # coord_equal()+ # this scales the x- and y-axes. But also compresses the plot so that it's hard to read. Not sure what the best approach is.
  theme_bw();targetplot_avg

ggsave("output/targetplot_7dayavg.png", plot = targetplot_avg, height = 5, width = 5, units = "in")



# GLM -------------------------------------------------------------------------------------------------------

## 7 DADM model --------------------
mod_7DADM <- lm(LOE_Target~ Watershed*DADM7, data=ctmax_summary)
summary(mod_7DADM)
# the above model accounts for the effect of watershed, the effect of seven day average daily maximum, and how that effect might differ between watershed
# I did not include an effect of year here, because DADM7 varies between years.
# However, to be statistically correct, i think I might have to, so I will do both.

## 7 DADA model --------------------------------------------------
mod_7DADA <- lm(LOE_Target~ Watershed*DADA7, data=ctmax_summary)
summary(mod_7DADA)

# the above two models are unlikely to differ much because DADM7 and DADA7 are highly correlated with one another.
# however, this difference may become more clear with emmeans pairwise comparisons

## 7DADM model with random effect of year -----------------------------------------
mod_7DADM_int <- lmer(LOE_Target~ Watershed*DADM7 + (1|Year), data=ctmax_summary)
summary(mod_7DADM_int)

## 7DADA model with random effect of year -----------------------------------------
mod_7DADM_int <- lmer(LOE_Target~ Watershed*DADM7 + (1|Year), data=ctmax_summary); summary(mod_7DADM_int)



mod_nested <- lmer(LOE_Target~ DADM7 + (1|Watershed:Year), data=ctmax_summary)
summary(mod_nested)

# do watersheds differ by year?
modsimple <- lmer(LOE_Target ~ Watershed + (1|Year), data = ctmax_summary)
summary(modsimple)

str(ctmax_summary$Year) # check that Year is a factor
## MODEL SELECTION --------------------------------------------------------------------------
#this model makes the most sense because we want to see how ctmax varies between years
# add the random effect of Watershed, since we expect those sites to be more similar to one another than not
library(emmeans)
mod2 <- lm(LOE_Target ~ WatershedYear*DADM7, data = ctmax_summary)
summary(mod2)

mod3 <- lm(LOE_Target ~ DADM7, data = ctmax_summary); summary(mod3)

mod4 <- lm(LOE_Target ~ WatershedYear, data = ctmax_summary); summary(mod4)
mod_both <- lm(LOE_Target ~ Watershed + DADM7,
               data = ctmax_summary)
mod_ws <- lm(LOE_Target ~ Watershed, data = ctmax_summary); summary(mod_ws)
mod_year <- lm(LOE_Target~Year, data = ctmax_summary); summary(mod_year)
anova(mod4,mod_ws)
anova(mod_both,mod3)
anova(mod4, mod_both)


mod_final <- lm(LOE_Target ~ Watershed + DADM7, data = ctmax_summary); summary(mod_final)
# Instead of including effect of 
mod_final2 <- lm(LOE_Target ~ Watershed*DADM7, data = ctmax_summary); summary(mod_final2)

# mod_final is slightly better, based on higher R2 values and lower AIC values
AIC(mod_final)
AIC(mod_final2)

# ok so the output file 


# SIZE EFFECT ----------------------------------------------------------------------------------------------
ggplot(treat2, aes(x= Length_mm, y= LOE.Temp_C, color = Watershed))+
  geom_point()+
  geom_smooth(method = lm)
