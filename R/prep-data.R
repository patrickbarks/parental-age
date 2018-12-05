

## Must have the following packages installed
# install.packages(c("Matrix", "dplyr", "tidyr", "tibble", "mgcv", "popbio",
#                    "popdemo", "ggplot2", "gridExtra"))


## Working directory should be Parental_Age_Scripts


## Load relevant packages / source required functions
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(mgcv)
library(popbio)
library(popdemo)
library(ggplot2)
library(gridExtra)
source('R/functions.R')
source('R/plot-themes.R')



##### Lemna minor data --------------------------------------------------------

## read data
# src_url <- "http://datadryad.org/bitstream/handle/10255/dryad.71128/BarksAndLaird_RawData_Phase2.csv"
# dat_lemna <- read.csv(url(src_url), stringsAsFactors = FALSE)
# save(dat_lemna, file = 'dat/dat_lemna.RData')
load('dat/dat_lemna.RData')


## extract subset containing reproduction data only
# number of offspring produced each day by each frond
dat_lemna_repro <- as.matrix(dplyr::select(dat_lemna, matches('2012|2013')))


## determine dates of first and last reproduction
# get column name (i.e. date) associated with first (min) or last (max) reproduction
GetFirstRepro <- function (x) names(x)[min(which(x > 0))]
GetLastRepro <- function (x) names(x)[max(which(x > 0))]
dat_lemna$date.first.repro <- apply(dat_lemna_repro, 1, GetFirstRepro)
dat_lemna$date.last.repro <- apply(dat_lemna_repro, 1, GetLastRepro)


## convert date values to R's date format
dat_lemna$date.focal.birth <- as.Date(dat_lemna$date.focal.birth, format = '%b.%d.%Y')
dat_lemna$date.parent.birth <- as.Date(dat_lemna$date.parent.birth, format = '%b.%d.%Y')
dat_lemna$date.first.repro <- as.Date(dat_lemna$date.first.repro, format = '%b.%d.%Y')
dat_lemna$date.last.repro <- as.Date(dat_lemna$date.last.repro, format = '%b.%d.%Y')


## determine parents' age at focals' birth
dat_lemna$parent.age <- as.numeric(dat_lemna$date.focal.birth - dat_lemna$date.parent.birth)



#### s = 4 parental age classes
## set model parameters
o <- 30   # number of age classes (max age)
s <- 4    # number of parental age classes


## group parental ages into s categories of similar sample size (using quantile fn)
dat_lemna$j <- dat_lemna$parent.age %>% 
  cut(quantile(., (0:s)/s), labels = 1:s, include.lowest = TRUE) %>% 
  as.integer()


## transform dat_lemna to flat (i.e. tidy) form
dat_lemna_flat <- dat_lemna %>%
  dplyr::select(ident.focal, ident.parent, parent.age, j, date.focal.birth,
                date.last.repro, matches('2012|2013')) %>% 
  gather(date, fecund, matches('2012|2013')) %>%
  arrange(ident.focal) %>% 
  mutate(date = as.Date(date, format = '%b.%d.%Y')) %>% 
  filter(date >= date.focal.birth & date <= date.last.repro) %>%
  mutate(age = as.numeric(date - date.focal.birth + 1),
         died = ifelse(date == date.last.repro, 1, 0),
         surv = ifelse(date == date.last.repro, 0, 1)) %>% 
  dplyr::select(-contains('date')) %>% 
  as_tibble()


## tidy df with mapping between age classes and parental age classes
age_map_lemna <- dat_lemna_flat %>% 
  dplyr::select(i = parent.age, j = j) %>% 
  mutate(i = as.integer(i),
         j = as.integer(as.factor(j))) %>% 
  unique() %>% 
  arrange(i, j) %>% 
  rbind(c(as.integer(o), as.integer(s))) %>% 
  mutate(species = 'Lemna minor') %>% 
  dplyr::select(species, i, j)


## use gam models to smoothed survival and fecundity transition rates
# gam-predictions using integer age sequence (for modeling)
P_lemna_flat <- dat_lemna_flat %>% 
  group_by(j) %>% 
  do(GamFnSurv(data = ., plot = FALSE)) %>% 
  ungroup()

F_lemna_flat <- dat_lemna_flat %>% 
  group_by(j) %>% 
  do(GamFnFert(data = ., plot = FALSE)) %>% 
  ungroup()


## bind outputs together and write to files
Tr_lemna_flat <- full_join(P_lemna_flat,
                           F_lemna_flat,
                           by = c('i', 'j')) %>% 
  mutate(species = 'Lemna minor') %>% 
  dplyr::select(species, i, j, P_ij, F_ij)

## write to file
write.csv(age_map_lemna, "dat/age_map_lemna_s4.csv", row.names = FALSE)
write.csv(Tr_lemna_flat, "dat/transition_rates_lemna_s4.csv", row.names = FALSE)




#### s = 8 parental age classes
## set model parameters
o <- 30   # number of age classes (max age)
s <- 8    # number of parental age classes

## group parental ages into s categories of similar sample size (using quantile fn)
dat_lemna$j <- dat_lemna$parent.age %>% 
  cut(quantile(., (0:s)/s), labels = 1:s, include.lowest = TRUE) %>% 
  as.integer()

## transform dat_lemna to flat (i.e. tidy) form
dat_lemna_flat <- dat_lemna %>%
  dplyr::select(ident.focal, ident.parent, parent.age, j, date.focal.birth,
                date.last.repro, matches('2012|2013')) %>% 
  gather(date, fecund, matches('2012|2013')) %>%
  arrange(ident.focal) %>% 
  mutate(date = as.Date(date, format = '%b.%d.%Y')) %>% 
  filter(date >= date.focal.birth & date <= date.last.repro) %>%
  mutate(age = as.numeric(date - date.focal.birth + 1),
         died = ifelse(date == date.last.repro, 1, 0),
         surv = ifelse(date == date.last.repro, 0, 1)) %>% 
  dplyr::select(-contains('date')) %>% 
  as_tibble()

## tidy df with mapping between age classes and parental age classes
age_map_lemna <- dat_lemna_flat %>% 
  dplyr::select(i = parent.age, j = j) %>% 
  mutate(i = as.integer(i),
         j = as.integer(as.factor(j))) %>% 
  unique() %>% 
  arrange(i, j) %>% 
  rbind(c(as.integer(o), as.integer(s))) %>% 
  mutate(species = 'Lemna minor') %>% 
  dplyr::select(species, i, j)

## use gam models to smoothed survival and fecundity transition rates
# gam-predictions using high-res sequences of ages (for plotting)
P_lemna_flat_plot <- dat_lemna_flat %>%
  group_by(j) %>%
  do(GamFnSurv(data = ., plot = TRUE)) %>%
  ungroup() %>%
  filter(i <= o-1)

F_lemna_flat_plot <- dat_lemna_flat %>%
  group_by(j) %>%
  do(GamFnFert(data = ., plot = TRUE)) %>%
  ungroup()

# gam-predictions using integer age sequence (for modeling)
P_lemna_flat <- dat_lemna_flat %>% 
  group_by(j) %>% 
  do(GamFnSurv(data = ., plot = FALSE)) %>% 
  ungroup()

F_lemna_flat <- dat_lemna_flat %>% 
  group_by(j) %>% 
  do(GamFnFert(data = ., plot = FALSE)) %>% 
  ungroup()

## bind outputs together and write to files
Tr_lemna_flat <- full_join(P_lemna_flat,
                           F_lemna_flat,
                           by = c('i', 'j')) %>% 
  mutate(species = 'Lemna minor') %>% 
  dplyr::select(species, i, j, P_ij, F_ij)

Tr_lemna_flat_plot <- full_join(P_lemna_flat_plot,
                                F_lemna_flat_plot,
                                by = c('i', 'j'))

## write to file
write.csv(age_map_lemna, "dat/age_map_lemna_s8.csv", row.names = FALSE)
write.csv(Tr_lemna_flat, "dat/transition_rates_lemna_s8.csv", row.names = FALSE)
write.csv(Tr_lemna_flat_plot, "dat/transition_rates_plot_lemna_s8.csv", row.names = FALSE)








##### Simulated life cycles ----------------------------------------------------

### global variables for (almost) all simulated life cycles
s <- 9
o <- 20
i <- 1:o

P_range <- 0.6
F_range <- 0.0

P_i_sim1 <- logit_inv(11 - 0.7 * i)
P_i_sim2 <- logit_inv(2 - 0.1 * i)
P_i_sim3 <- logit_inv(0.6 + 0 * i)
P_i_sim4 <- logit_inv(-1.1 + 0.2 * i)

F_i_init1 <- 3 * logit_inv(11 - 0.7 * i)
F_i_init2 <- exp(2 - 0.1 * i)
F_i_init3 <- exp(1 + 0.0 * i)
F_i_init4 <- exp(-2 + 0.2 * i)




### Simulated life cycles used in main text
opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim1, F_i = F_i_init1, P_range = P_range, F_range = F_range, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim2, F_i = F_i_init2, P_range = P_range, F_range = F_range, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range, F_range = F_range, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim4, F_i = F_i_init4, P_range = P_range, F_range = F_range, lam = 1)

F_i_sim1 <- F_i_init1 * opt_sim1$par
F_i_sim2 <- F_i_init2 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init4 * opt_sim4$par

# generate parental age effect with initial fecundity values (for Appendix only)
Tr_sim_init1 <- SimulateParEffect(P_i_sim1, F_i_init1, s = s, P_range = P_range, F_range = F_range)
Tr_sim_init2 <- SimulateParEffect(P_i_sim2, F_i_init2, s = s, P_range = P_range, F_range = F_range)
Tr_sim_init3 <- SimulateParEffect(P_i_sim3, F_i_init3, s = s, P_range = P_range, F_range = F_range)
Tr_sim_init4 <- SimulateParEffect(P_i_sim4, F_i_init4, s = s, P_range = P_range, F_range = F_range)

# generate parental age effect with adjusted fecundity values
Tr_sim1 <- SimulateParEffect(P_i_sim1, F_i_sim1, s = s, P_range = P_range, F_range = F_range)
Tr_sim2 <- SimulateParEffect(P_i_sim2, F_i_sim2, s = s, P_range = P_range, F_range = F_range)
Tr_sim3 <- SimulateParEffect(P_i_sim3, F_i_sim3, s = s, P_range = P_range, F_range = F_range)
Tr_sim4 <- SimulateParEffect(P_i_sim4, F_i_sim4, s = s, P_range = P_range, F_range = F_range)

# full set of sim transition rates
Tr_sim_init <- rbind.data.frame(
  Tr_sim_init1 %>% mutate(sim = 'Sim 1'),
  Tr_sim_init2 %>% mutate(sim = 'Sim 2'),
  Tr_sim_init3 %>% mutate(sim = 'Sim 3'),
  Tr_sim_init4 %>% mutate(sim = 'Sim 4')
)

Tr_sim <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4')
)

## age map for sim species
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)


## write to file (parental age effect on survival)
write.csv(Tr_sim, 'dat/transition_rates_sim_main.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_main.csv', row.names = FALSE)



### Figure A2
tr_p_full <- rbind.data.frame(
  data.frame(i, P_i = P_i_sim1, type = "Sim 1"),
  data.frame(i, P_i = P_i_sim2, type = "Sim 2"),
  data.frame(i, P_i = P_i_sim3, type = "Sim 3"),
  data.frame(i, P_i = P_i_sim4, type = "Sim 4")
)

tr_f_full <- rbind.data.frame(
  data.frame(i, F_i = F_i_init1, type = "Sim 1"),
  data.frame(i, F_i = F_i_init2, type = "Sim 2"),
  data.frame(i, F_i = F_i_init3, type = "Sim 3"),
  data.frame(i, F_i = F_i_init4, type = "Sim 4")
)

pa2_1 <- ggplot(tr_p_full) +
  geom_line(aes(i, P_i)) +
  facet_wrap(~ type, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab(expression(italic(P[i]))) +
  tt_a234 + theme(axis.text.x = element_blank(), axis.title.x = element_blank())

pa2_2 <- ggplot(tr_f_full) +
  geom_line(aes(i, F_i)) +
  facet_wrap(~ type, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_log10(labels = LabelFn) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[i]))) +
  tt_a234 + theme(strip.text.x = element_blank())

# plots to grobs
ga2_1 <- ggplotGrob(pa2_1)
ga2_2 <- ggplotGrob(pa2_2)

# arrange panels
ga2 <- arrangeGrob(rbind(ga2_1, ga2_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 4, width = 6, dpi = 150)
# grid.arrange(ga2)

# write to file
ggsave('img/Fig_A2.png', ga2, height = 4, width = 6, units = 'in', dpi = 300)



## Figure A3
lwd_1 <- 2.1

df_lab <- data.frame(sim = "Sim 1", label = "Parental age")
df_arrow <- data.frame(sim = "Sim 1", x1 = 0.4, x2 = 10, y1 = 0.16, y2 = 0.16)

pa3_1 <- ggplot(Tr_sim_init) +
  geom_line(aes(i, P_ij, col = j, group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.4, 0.22, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
             axis.title.x = element_blank())

pa3_2 <- ggplot(filter(Tr_sim_init, j == 1)) +
  geom_line(aes(i, F_ij), lwd = lwd_1, col = col_low) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_log10(labels = LabelFn) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())


# plots to grobs
ga3_1 <- ggplotGrob(pa3_1)
ga3_2 <- ggplotGrob(pa3_2)

# arrange panels
ga3 <- arrangeGrob(rbind(ga3_1, ga3_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 4, width = 6, dpi = 150)
# grid.arrange(ga3)

# write to file
ggsave('img/Fig_A3.png', ga3, height = 4, width = 6, units = 'in', dpi = 300)



### Figure A4
pa4_1 <- ggplot(Tr_sim) +
  geom_line(aes(i, P_ij, col = j, group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.4, 0.22, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank())

pa4_2 <- ggplot(filter(Tr_sim, j == 1)) +
  geom_line(aes(i, F_ij), lwd = lwd_1, col = col_low) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_log10(labels = LabelFn) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())


# plots to grobs
ga4_1 <- ggplotGrob(pa4_1)
ga4_2 <- ggplotGrob(pa4_2)

# arrange panels
ga4 <- arrangeGrob(rbind(ga4_1, ga4_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 4, width = 6, dpi = 150)
# grid.arrange(ga4)

# write to file
ggsave('img/Fig_A4.png', ga4, height = 4, width = 6, units = 'in', dpi = 300)





########## Transition rates for sensitivity analyses (Appendix S1 Section 7) ---

### Scenario a: vary the population growth rate (lambda) -----------------------
opt_wrapper <- function(lambda) {
  # create life cycles using inits for Sim #3 from main text as baseline
  optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
        P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range, F_range = F_range,
        lam = lambda)
}

opt_sim1 <- opt_wrapper(lambda = 0.5)
opt_sim2 <- opt_wrapper(lambda = 0.8)
opt_sim3 <- opt_wrapper(lambda = 1.1)
opt_sim4 <- opt_wrapper(lambda = 1.4)

F_i_sim1 <- F_i_init3 * opt_sim1$par
F_i_sim2 <- F_i_init3 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init3 * opt_sim4$par

Tr_sim1 <- SimulateParEffect(P_i_sim3, F_i_sim1, s = s, P_range = P_range, F_range = F_range)
Tr_sim2 <- SimulateParEffect(P_i_sim3, F_i_sim2, s = s, P_range = P_range, F_range = F_range)
Tr_sim3 <- SimulateParEffect(P_i_sim3, F_i_sim3, s = s, P_range = P_range, F_range = F_range)
Tr_sim4 <- SimulateParEffect(P_i_sim3, F_i_sim4, s = s, P_range = P_range, F_range = F_range)

## full set of sim transition rates
Tr_sim_scen_a <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'lambda == 0.5'),
  Tr_sim2 %>% mutate(sim = 'lambda == 0.8'),
  Tr_sim3 %>% mutate(sim = 'lambda == 1.1'),
  Tr_sim4 %>% mutate(sim = 'lambda == 1.4')
)

## age map for sim species
map_ij_sim_scen_a <- Tr_sim_scen_a %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

# write to file
write.csv(Tr_sim_scen_a, 'dat/transition_rates_sim_scen_a.csv', row.names = F)
write.csv(map_ij_sim_scen_a, 'dat/map_ij_sim_scen_a.csv', row.names = F)




### Scenario b: s = 4 parental age classes -------------------------------------
s4 <- 4

opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s4, P_i = P_i_sim1, F_i = F_i_init1, P_range = P_range, F_range = F_range, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s4, P_i = P_i_sim2, F_i = F_i_init2, P_range = P_range, F_range = F_range, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s4, P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range, F_range = F_range, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s4, P_i = P_i_sim4, F_i = F_i_init4, P_range = P_range, F_range = F_range, lam = 1)

F_i_sim1 <- F_i_init1 * opt_sim1$par
F_i_sim2 <- F_i_init2 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init4 * opt_sim4$par

Tr_sim1 <- SimulateParEffect(P_i_sim1, F_i_sim1, s = s4, P_range = P_range, F_range = F_range)
Tr_sim2 <- SimulateParEffect(P_i_sim2, F_i_sim2, s = s4, P_range = P_range, F_range = F_range)
Tr_sim3 <- SimulateParEffect(P_i_sim3, F_i_sim3, s = s4, P_range = P_range, F_range = F_range)
Tr_sim4 <- SimulateParEffect(P_i_sim4, F_i_sim4, s = s4, P_range = P_range, F_range = F_range)

# full set of sim transition rates
Tr_sim_scen_b <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4')
)

## age map for sim species
map_ij_sim_scen_b <- Tr_sim_scen_b %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s4)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

# write to file
write.csv(Tr_sim_scen_b, 'dat/transition_rates_sim_scen_b.csv', row.names = F)
write.csv(map_ij_sim_scen_b, 'dat/map_ij_sim_scen_b.csv', row.names = F)




### Scenario c: positive parental age effect -----------------------------------
opt_sim1 <- optim(0.5, opt_fn_rev, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim1, F_i = F_i_init1, P_range = P_range, F_range = F_range, lam = 1)
opt_sim2 <- optim(0.5, opt_fn_rev, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim2, F_i = F_i_init2, P_range = P_range, F_range = F_range, lam = 1)
opt_sim3 <- optim(0.5, opt_fn_rev, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range, F_range = F_range, lam = 1)
opt_sim4 <- optim(0.5, opt_fn_rev, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim4, F_i = F_i_init4, P_range = P_range, F_range = F_range, lam = 1)

F_i_sim1 <- F_i_init1 * opt_sim1$par
F_i_sim2 <- F_i_init2 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init4 * opt_sim4$par

Tr_sim1 <- SimulateParEffectRev(P_i_sim1, F_i_sim1, s = s, P_range = P_range, F_range = F_range)
Tr_sim2 <- SimulateParEffectRev(P_i_sim2, F_i_sim2, s = s, P_range = P_range, F_range = F_range)
Tr_sim3 <- SimulateParEffectRev(P_i_sim3, F_i_sim3, s = s, P_range = P_range, F_range = F_range)
Tr_sim4 <- SimulateParEffectRev(P_i_sim4, F_i_sim4, s = s, P_range = P_range, F_range = F_range)

## full set of sim transition rates
Tr_sim_scen_c <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4')
)

## age map for sim species
map_ij_sim_scen_c <- Tr_sim_scen_c %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

## write to file
write.csv(Tr_sim_scen_c, 'dat/transition_rates_sim_scen_c.csv', row.names = F)
write.csv(map_ij_sim_scen_c, 'dat/map_ij_sim_scen_c.csv', row.names = F)




### Scenario d: parental age effect on fecundity rather than survival ----------
P_range_scen_d <- 0.0
F_range_scen_d <- 0.4

opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim1, F_i = F_i_init1, P_range = P_range_scen_d, F_range = F_range_scen_d, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim2, F_i = F_i_init2, P_range = P_range_scen_d, F_range = F_range_scen_d, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range_scen_d, F_range = F_range_scen_d, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim4, F_i = F_i_init4, P_range = P_range_scen_d, F_range = F_range_scen_d, lam = 1)

F_i_sim1 <- F_i_init1 * opt_sim1$par
F_i_sim2 <- F_i_init2 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init4 * opt_sim4$par

Tr_sim1 <- SimulateParEffect(P_i_sim1, F_i_sim1, s = s, P_range = P_range_scen_d, F_range = F_range_scen_d)
Tr_sim2 <- SimulateParEffect(P_i_sim2, F_i_sim2, s = s, P_range = P_range_scen_d, F_range = F_range_scen_d)
Tr_sim3 <- SimulateParEffect(P_i_sim3, F_i_sim3, s = s, P_range = P_range_scen_d, F_range = F_range_scen_d)
Tr_sim4 <- SimulateParEffect(P_i_sim4, F_i_sim4, s = s, P_range = P_range_scen_d, F_range = F_range_scen_d)

## full set of sim transition rates
Tr_sim_scen_d <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4')
)

## age map for sim species
map_ij_sim_scen_d <- Tr_sim_scen_d %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

## write to file (parental age effect on survival)
write.csv(Tr_sim_scen_d, 'dat/transition_rates_sim_scen_d.csv', row.names = FALSE)
write.csv(map_ij_sim_scen_d, 'dat/map_ij_sim_scen_d.csv', row.names = FALSE)




### Scenario e: stasis loop at the maximum age class ---------------------------
opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim1, F_i = F_i_init1, P_range = P_range, F_range = F_range, lam = 1, stasis = TRUE)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim2, F_i = F_i_init2, P_range = P_range, F_range = F_range, lam = 1, stasis = TRUE)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim3, F_i = F_i_init3, P_range = P_range, F_range = F_range, lam = 1, stasis = TRUE)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim4, F_i = F_i_init4, P_range = P_range, F_range = F_range, lam = 1, stasis = TRUE)

F_i_sim1 <- F_i_init1 * opt_sim1$par
F_i_sim2 <- F_i_init2 * opt_sim2$par
F_i_sim3 <- F_i_init3 * opt_sim3$par
F_i_sim4 <- F_i_init4 * opt_sim4$par

Tr_sim1 <- SimulateParEffect(P_i_sim1, F_i_sim1, s = s, P_range = P_range, F_range = F_range)
Tr_sim2 <- SimulateParEffect(P_i_sim2, F_i_sim2, s = s, P_range = P_range, F_range = F_range)
Tr_sim3 <- SimulateParEffect(P_i_sim3, F_i_sim3, s = s, P_range = P_range, F_range = F_range)
Tr_sim4 <- SimulateParEffect(P_i_sim4, F_i_sim4, s = s, P_range = P_range, F_range = F_range)

# full set of sim transition rates
Tr_sim_scen_e <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4')
)

## age map for sim species
map_ij_sim_scen_e <- Tr_sim_scen_e %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

# write to file
write.csv(Tr_sim_scen_e, 'dat/transition_rates_sim_scen_e.csv', row.names = F)
write.csv(map_ij_sim_scen_e, 'dat/map_ij_sim_scen_e.csv', row.names = F)


