# sablefish histological data only
# ben.williams@noaa.gov
# 2022-06

# Notes: dropped early legs from ak summer survey
# see Rodgveller 2018 (poor classification early in season)

# load ----
source(here::here("R", "helper.r"))
ggplot2::theme_set(
  ggplot2::theme_light() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      # axis.ticks.length = grid::unit(base_ / 2.2, "pt"),
      strip.background = ggplot2::element_rect(fill = NA, colour = NA),
      strip.text.x = element_text(colour = "black"),
      strip.text.y = element_text(colour = "black"),
      panel.border = element_rect(fill = NA),
      legend.key.size = grid::unit(0.9, "lines"),
      legend.key = ggplot2::element_rect(colour = NA, fill = NA),
      legend.background = ggplot2::element_rect(colour = NA, fill = NA),
      text = element_text(family = "Times New Roman")
    )
)

# data ----
# west coast
read_csv("data/WCGBT_ODFW_WDFW_2010_2018_sablefishmaturity_updated022720.csv") %>% 
  rename_all(tolower) %>% 
  transmute(mature = mature,
            date = mdy(date_collected),
            year = year(date),
            age = ifelse(age_yrs>40, 40, age_yrs),
            length = length_cm,
            weight = weight_kg, 
            lat = latitude_dd,
            long = longitude_dd, 
            depth = depth_m, 
            Age = factor(age),
            Year = factor(year),
            Mature = factor(mature),
            Region = factor(case_when(lat<=36 ~ "R1",
                                      lat>36 & lat <=50 ~ "R2",
                                      lat>50 & long>= -145 ~ "R3",
                                      long< -145 ~ "R4"))) -> wc

# alaska 

read.csv("data/all_3_histo_V3.csv") %>% 
  rename_all(tolower) %>% 
  # filter(leg>5) %>%
  transmute(mature = maturity_biological_ssmature, 
            # functional = maturity_functional_ssimmature,
            date = mdy(date),
            year = year(date),
            age = ifelse(age > 40, 40, age),
            length = totlength / 10,
            weight = totalwt / 1000,
            lat = start_latitude,
            long = start_longitude,
            depth = case_when(depth ==1 ~ 100,
                              TRUE ~ depth),
            Age = factor(age),
            Year = factor(year),
            Mature = factor(mature),
            Region = factor(case_when(lat<=36 ~ "R1",
                                      lat>36 & lat <=50 ~ "R2",
                                      lat>50 & long>= -145 ~ "R3",
                                      long< -145 ~ "R4"))) -> ak

# full data set
bind_rows(wc, ak) %>% 
  mutate(Year = factor(Year, levels = c("2010", "2011", "2015", "2016", "2018"))) -> dat

# filtered data for age and region info
dat %>% 
  filter(age>0, !is.na(Region)) -> base

# vonB data from Maia Kapur
# optimized the curves on the observed maturity age/lengths to get regional mean spawning lengths
load(here::here('data', 'OM_growthPars.rdata'))

as.data.frame((growthPars$Linf_yk)) %>% 
  select(ends_with("Fem")) %>% 
  tail(1) %>% 
  `names<-`(.,gsub("\\.Fem", "", names(.) )) %>% 
  mutate(a = "a") %>% 
  pivot_longer(-a, values_to = "linf", names_to = 'Region') -> linf

as.data.frame((growthPars$L1_yk)) %>% 
  select(ends_with("Fem")) %>% 
  tail(1) %>% 
  `names<-`(.,gsub("\\.Fem", "", names(.) )) %>% 
  mutate(a = "a") %>% 
  pivot_longer(-a, values_to = "l1", names_to = 'Region') %>% 
  left_join(linf) -> l1

as.data.frame((growthPars$kappa_yk)) %>% 
  select(ends_with("Fem")) %>% 
  tail(1) %>% 
  `names<-`(.,gsub("\\.Fem", "", names(.) )) %>% 
  mutate(a = "a") %>% 
  pivot_longer(-a, values_to = "kappa", names_to = 'Region') %>% 
  left_join(l1) -> kappa

as.data.frame((growthPars$sigmaG_yk)) %>% 
  select(ends_with("Fem")) %>% 
  tail(1) %>% 
  `names<-`(.,gsub("\\.Fem", "", names(.) )) %>% 
  mutate(a = "a") %>% 
  pivot_longer(-a, values_to = "sigma", names_to = 'Region') %>% 
  left_join(kappa) %>% 
  dplyr::select(-a) -> pars

# get kappa adjustment
dat %>% 
  left_join(pars) %>% 
  drop_na(c(kappa, age)) -> df

df %>% 
  group_by(Region) %>% 
  nest() %>% 
  mutate(fit = purrr::map(data, ~ nls(length ~ linf + (l1-linf) * exp(-kappa * x * age), data = .x, start = list(x=0.5)))) -> fits

f1 <- coef(summary(fits$fit[[1]]))[[1]]
f2 <- coef(summary(fits$fit[[2]]))[[1]]
f3 <- coef(summary(fits$fit[[3]]))[[1]]
f4 <- coef(summary(fits$fit[[4]]))[[1]]

expand.grid(age = 1:42, 
            Region = c('R1', 'R2', 'R3', 'R4')) %>% 
  left_join(pars) %>% 
  mutate(length = case_when(Region=='R1' ~ round(linf + (l1-linf) * exp(-kappa * f1 * age), 1),
                            Region=='R2' ~ round(linf + (l1-linf) * exp(-kappa * f2 * age), 1),
                            Region=='R3' ~ round(linf + (l1-linf) * exp(-kappa * f3 * age), 1),
                            Region=='R4' ~ round(linf + (l1-linf) * exp(-kappa * f4 * age), 1)),
         age = as.integer(age)) %>% 
  dplyr::select(Region, age, length) %>% 
  mutate(Age = factor(age)) -> vonb

vonb %>% 
  dplyr::select(Region, Age=age, length) %>% 
  tidyr::pivot_wider(names_from=Region, values_from = length) %>% 
  vroom::vroom_write(here::here('data', "lngs.csv"), ",")

# models ----
ma <- gam(mature ~ s(age, k=4), data=base, family='binomial', gamma=1.4)
mal <-gam(mature ~ s(age, k=4) + s(length, k=4), data=base, family='binomial', gamma=1.4)
mar <- gam(mature ~ s(age, k=4) + Region, data=base, family='binomial', gamma=1.4)
marl <-gam(mature ~ s(age, k=4) + s(length, k=4) + Region, data=base, family='binomial', gamma=1.4)
mard <- gam(mature ~ s(age, k=4) + s(depth, k=4) + Region, data=base, family='binomial', gamma=1.4)
marld <- gam(mature ~ s(age, k=4) + s(depth, k=4) + s(length, k=4) + Region, data=base, family='binomial', gamma=1.4)
marlld <- gam(mature ~ s(age, k=4) + s(depth, k=4) + s(length, k=4) + te(long, lat, by = Region) + Region, data=base, family='binomial', gamma=1.4)

AIC(ma, mal, mar, marl, mard, marld, marlld) %>% 
  rownames_to_column("model")  %>% 
  mutate(AIC = round(AIC)) -> aic
mutate(delta = min(AIC) - AIC) %>% 
  arrange(-delta) -> aic

vroom::vroom_write(aic, here::here('data', "aic.csv"), ",")

BIC(ma, mal, mar, marl, mard, marld, marlld) %>% 
  rownames_to_column("model") %>% 
  mutate(BIC = round(BIC)) -> bic
mutate(delta_b = min(BIC) - BIC) %>% 
  arrange(-delta_b) -> bic

vroom::vroom_write(bic, here::here('data', "bic.csv"), ",")

summary(ma)
summary(mal)
summary(mar)
summary(marl)
summary(mard)
summary(marld)
summary(marlld)

data.frame(
  model = c('ma', 'mal', 'mar', 'marl', 'mard', 'marld', 'marlld'),
  coefficients = c('age', 'age + length', 'age + region', 'age + region + length', 'age + region + depth', 'age + region + length + depth', 'age + region + length + location + depth'),
  R = c(0.612, 0.689, 0.658, 0.76, 0.663, 0.761, 0.773),
  dev = c(55.6, 64.2, 60.1, 71.2, 60.8, 71.4, 72.7)) %>% 
  arrange(model) %>% 
  left_join(aic) %>% 
  left_join(bic) %>% 
  vroom::vroom_write(here::here('data', 'gam_results.csv'), ',')


base %>% 
  group_by(Region) %>% 
  summarise(depth = median(depth, na.rm = TRUE),
            lat = median(lat, na.rm = TRUE),
            long = median(long, na.rm = TRUE),
            n = n()) -> locs

vroom::vroom_write(locs, here::here('data', "locs.csv"), ",")

vonb %>% 
  left_join(locs) %>% 
  # filter(!(age %in% c(30, 39:42))) %>% 
  mutate(Age = factor(age)) %>% 
  mutate(# mal2 = predict.gam(mal2, ., type = 'response'),
    # ma = predict.gam(ma, ., type = 'response'),
    # mal = predict.gam(mal, ., type = 'response'),
    # mar = predict.gam(mar, ., type = 'response'),
    marl = predict.gam(marl, ., type = 'response'),
    # mard = predict.gam(mard, ., type = 'response'),
    # marld = predict.gam(marld, ., type = 'response'),
    marlld = predict.gam(marlld, ., type = 'response')
  ) %>% 
  pivot_longer(-c(Region, age, length, Age, depth, lat, long, n)) %>% 
  ggplot(aes(age, value, color = name)) + 
  geom_line() + 
  facet_wrap(~Region) + 
  # coord_cartesian(xlim = c(0,30)) + 
  geom_hline(yintercept = 0.5, lty = 3) +
  scale_color_scico_d("Model", palette = 'roma') + 
  xlab('Age') + 
  ylab("Proportion mature")

ggsave(here::here('figs', 'mat2.png'), width = 6.5, height = 6.5, units = "in", dpi=200)
