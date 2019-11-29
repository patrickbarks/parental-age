

## Must have the following packages installed
# install.packages(c("Matrix", "dplyr", "tidyr", "tibble", "mgcv", "popbio",
#                    "popdemo", "ggplot2", "gridExtra"))


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



#### Lemna minor data ----------------------------------------------------------

## read data
# src_url <- "https://datadryad.org/stash/downloads/file_stream/36269"
# dat_lemna <- read.csv(src_url, stringsAsFactors = FALSE)
# save(dat_lemna, file = "dat/BarksAndLaird_RawData_Phase2.RData")
load("dat/BarksAndLaird_RawData_Phase2.RData")


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



### plot empirical transtion rates for L minor
df_means <- dat_lemna_flat %>% 
  group_by(age, j) %>% 
  summarize(n = n(), surv = mean(surv), fert = mean(fecund)) %>% 
  ungroup() %>% 
  mutate(j = factor(j, labels = c("Parental age 1",
                                  "Parental age 2",
                                  "Parental age 3",
                                  "Parental age 4",
                                  "Parental age 5",
                                  "Parental age 6",
                                  "Parental age 7",
                                  "Parental age 8")))


fig_a3_1 <- ggplot(filter(df_means, age <= 30), aes(x = age, y = surv, size = n)) +
  geom_point(shape = 19, col = "grey30", alpha = 0.4) +
  geom_point(shape = 1, alpha = 0.6) +
  facet_wrap(~ j, ncol = 2, dir = "v") +
  scale_y_continuous(limits = c(0.35, 1.05), breaks = seq(0.4, 1, 0.2)) +
  scale_size(range = c(1, 3), name = "N", guid = F) +
  labs(x = "Age (days)", y = "Survival") +
  tt2 +
  theme(axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        strip.text = element_text(size = 8))

fig_a3_2 <- ggplot(filter(df_means, age <= 30), aes(x = age, y = fert, size = n)) +
  geom_point(shape = 19, col = "grey30", alpha = 0.4) +
  geom_point(shape = 1, alpha = 0.6) +
  facet_wrap(~ j, ncol = 2, dir = "v") +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.5)) +
  scale_size(range = c(1, 3), name = "Sample\nsize", breaks = c(10, 20, 40, 80)) +
  labs(x = "Age (days)", y = "Fecundity") +
  tt2 +
  theme(axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        legend.margin = margin(0, 5, 0, 0),
        strip.text = element_text(size = 8),
        legend.title = element_text(size = 9.5))

fig_a3 <- arrangeGrob(ggplotGrob(fig_a3_1),
                      ggplotGrob(fig_a3_2),
                      nrow = 1, widths = c(1, 1.27))

quartz(height = 5, width = 6.25, dpi = 150)
grid.arrange(fig_a3)

# write to file
# ggsave('img/Fig_A3.png', fig_a3, height = 5, width = 6.25, units = 'in', dpi = 300)



### compare gam models for survival
mods1 <- gam(surv ~ s(age),
             method = "REML", data = dat_lemna_flat, family = "binomial")

mods2 <- gam(surv ~ s(age) + s(j, bs = "re"),
             method = "REML", data = dat_lemna_flat, family = "binomial")

mods3 <- gam(surv ~ s(age) + s(age, j, bs = "fs"),
             method = "REML", data = dat_lemna_flat, family = "binomial")

mods4 <- gam(surv ~ s(age, bs = "tp") + s(age, by = j, m = 1, bs = "tp") + s(j, bs = "re", k = 8),
             method = "REML", data = dat_lemna_flat, family = "binomial")

mods5 <- gam(surv ~ s(age, j, bs = "fs"),
             method = "REML", data = dat_lemna_flat, family = "binomial")

mods6 <- gam(surv ~ j + s(age, by = j),
             method = "REML", data = dat_lemna_flat, family = "binomial")


dfs_aic <- as_tibble(AIC(mods1, mods2, mods3, mods4, mods5, mods6)) %>% 
  mutate(model = paste0("p", 1:6)) %>% 
  mutate(delta_aic = paste0("ΔAIC = ", sprintf("%.1f", AIC - min(AIC)))) %>% 
  mutate(model = factor(model, labels = c("No parental age effect",
                                          "Random intercept",
                                          "Group smoother: common shape,\ncommon wiggliness",
                                          "Group smoother: common shape,\nindependent wiggliness",
                                          "Group smoother: independent shape,\ncommon wiggliness",
                                          "Group smoother: independent shape,\nindependent wiggliness")))

dfs <- expand.grid(j = 1:8, age = seq(1, 30, 0.1)) %>% 
  as_tibble() %>% 
  mutate(p1 = predict(mods1, newdata = ., type = "response"),
         p2 = predict(mods2, newdata = ., type = "response"),
         p3 = predict(mods3, newdata = ., type = "response"),
         p4 = predict(mods4, newdata = ., type = "response"),
         p5 = predict(mods5, newdata = ., type = "response"),
         p6 = predict(mods5, newdata = ., type = "response")) %>% 
  gather(model, pred, p1:p6) %>% 
  mutate(model = factor(model, labels = c("No parental age effect",
                                          "Random intercept",
                                          "Group smoother: common shape,\ncommon wiggliness",
                                          "Group smoother: common shape,\nindependent wiggliness",
                                          "Group smoother: independent shape,\ncommon wiggliness",
                                          "Group smoother: independent shape,\nindependent wiggliness")))


fig_a4 <- ggplot(dfs, aes(age, pred, col = j, group = j)) +
  geom_line(size = 2) +
  geom_text(data = dfs_aic, inherit.aes = FALSE, aes(x = 30, y = 0.98, label = delta_aic),
            hjust = 1, size = 3) +
  facet_wrap(~ model, dir = "v", ncol = 2) +
  scale_color_gradient(low = "grey5", high = "grey80", name = "Parental\nage group") +
  labs(x = "Age (days)", y = "Predicted survival rate") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text = element_text(size = 8.2),
        axis.title.y = element_text(margin = margin(0, 7, 0, 0)),
        strip.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8.2))

dev.off()
quartz(height = 6, width = 6.25, dpi = 150)
print(fig_a4)

# write to file
# ggsave('img/Fig_A4.png', fig_a4, height = 6, width = 6.25, units = 'in', dpi = 300)


### compare gam models for fecundity
mod1 <- gam(fecund ~ s(age),
            method = "REML", data = dat_lemna_flat, family = "poisson")

mod2 <- gam(fecund ~ s(age) + s(j, bs = "re"),
            method = "REML", data = dat_lemna_flat, family = "poisson")

mod3 <- gam(fecund ~ s(age) + s(age, j, bs = "fs"),
            method = "REML", data = dat_lemna_flat, family = "poisson")

mod4 <- gam(fecund ~ s(age) + s(age, by = j, m = 1, bs = "tp") + s(j, bs = "re"),
            method = "REML", data = dat_lemna_flat, family = "poisson")

mod5 <- gam(fecund ~ s(age, j, bs = "fs"),
            method = "REML", data = dat_lemna_flat, family = "poisson")

mod6 <- gam(fecund ~ j + s(age, by = j),
            method = "REML", data = dat_lemna_flat, family = "poisson")



df_aic <- as_tibble(AIC(mod1, mod2, mod3, mod4, mod5, mod6)) %>% 
  mutate(model = paste0("p", 1:6)) %>% 
  mutate(delta_aic = paste0("ΔAIC = ", sprintf("%.1f", AIC - min(AIC)))) %>% 
  mutate(model = factor(model, labels = c("No parental age effect",
                                          "Random intercept",
                                          "Group smoother: common shape,\ncommon wiggliness",
                                          "Group smoother: common shape,\nindependent wiggliness",
                                          "Group smoother: independent shape,\ncommon wiggliness",
                                          "Group smoother: independent shape,\nindependent wiggliness")))

df <- expand.grid(j = 1:8, age = seq(1, 30, 0.1)) %>%
  as_tibble() %>% 
  mutate(p1 = predict(mod1, newdata = ., type = "response"),
         p2 = predict(mod2, newdata = ., type = "response"),
         p3 = predict(mod3, newdata = ., type = "response"),
         p4 = predict(mod4, newdata = ., type = "response"),
         p5 = predict(mod5, newdata = ., type = "response"),
         p6 = predict(mod6, newdata = ., type = "response")) %>% 
  gather(model, pred, p1:p6) %>% 
  mutate(model = factor(model, labels = c("No parental age effect",
                                          "Random intercept",
                                          "Group smoother: common shape,\ncommon wiggliness",
                                          "Group smoother: common shape,\nindependent wiggliness",
                                          "Group smoother: independent shape,\ncommon wiggliness",
                                          "Group smoother: independent shape,\nindependent wiggliness")))


fig_a5 <- ggplot(df, aes(age, pred, col = j, group = j)) +
  geom_line(size = 2) +
  geom_text(data = df_aic, inherit.aes = FALSE, aes(x = 30, y = 0.98, label = delta_aic),
            hjust = 1, size = 3) +
  facet_wrap(~ model, dir = "v", ncol = 2) +
  scale_color_gradient(low = "grey5", high = "grey80", name = "Parental\nage group") +
  labs(x = "Age (days)", y = "Predicted fecundity") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text = element_text(size = 8.2),
        axis.title.y = element_text(margin = margin(0, 7, 0, 0)),
        strip.text = element_text(size = 9),
        legend.title = element_text(size = 9),
        legend.text = element_text(size = 8.2))

dev.off()
quartz(height = 6, width = 6.25, dpi = 150)
print(fig_a5)

# write to file
# ggsave('img/Fig_A5.png', fig_a5, height = 6, width = 6.25, units = 'in', dpi = 300)


### arrange all transition rate
Tr_lemna_flat <- expand.grid(age = 1:30, j = 1:8) %>% 
  as_tibble() %>% 
  mutate(P_ij = predict(mods5, newdata = ., type = "response")) %>% 
  mutate(F_ij = predict(mod3, newdata = ., type = "response")) %>% 
  rename(i = age)

Tr_lemna_flat_plot <- expand.grid(age = seq(1, 30, 0.1), j = 1:8) %>% 
  as_tibble() %>% 
  mutate(P_ij = predict(mods5, newdata = ., type = "response")) %>% 
  mutate(F_ij = predict(mod3, newdata = ., type = "response")) %>% 
  rename(i = age)



### write to file
write.csv(age_map_lemna, "dat/age_map_lemna_s8.csv", row.names = FALSE)
write.csv(Tr_lemna_flat, "dat/transition_rates_lemna_s8.csv", row.names = FALSE)
write.csv(Tr_lemna_flat_plot, "dat/transition_rates_plot_lemna_s8.csv", row.names = FALSE)





#### Simulated life cycles -----------------------------------------------------

## global variables for (almost) all simulated life cycles
s <- 10
o <- 20
i <- 1:o

## simulations in main text
P_i_sim <- rep(0.5, o)
F_i_init <- rep(1, o)

P_range <- 0

F_range1 <- 0.6
F_range2 <- 0.3
F_range3 <- 0.0
F_range4 <- -0.3
F_range5 <- -0.6


opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                  P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                  F_range = F_range1, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                  P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                  F_range = F_range2, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                  P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                  F_range = F_range3, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                  P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                  F_range = F_range4, lam = 1)
opt_sim5 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                  P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                  F_range = F_range5, lam = 1)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim4$par

# generate parental age effect with initial fecundity values (for Appendix only)
Tr_sim_init1 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range,
                                  F_range = F_range1)
Tr_sim_init2 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range,
                                  F_range = F_range2)
Tr_sim_init3 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range,
                                  F_range = F_range3)
Tr_sim_init4 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range,
                                  F_range = F_range4)
Tr_sim_init5 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range,
                                  F_range = F_range5)

# generate parental age effect with adjusted fecundity values
Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = P_range,
                             F_range = F_range1)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = P_range,
                             F_range = F_range2)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = P_range,
                             F_range = F_range3)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = P_range,
                             F_range = F_range4)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = P_range,
                             F_range = F_range5)

# full set of sim transition rates
Tr_sim_init <- rbind.data.frame(
  Tr_sim_init1 %>% mutate(sim = 'Sim 1'),
  Tr_sim_init2 %>% mutate(sim = 'Sim 2'),
  Tr_sim_init3 %>% mutate(sim = 'Sim 3'),
  Tr_sim_init4 %>% mutate(sim = 'Sim 4'),
  Tr_sim_init5 %>% mutate(sim = 'Sim 5')
)

Tr_sim <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4'),
  Tr_sim5 %>% mutate(sim = 'Sim 5')
)

## map between age classes and parental age classes
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)

## write to file
write.csv(Tr_sim, 'dat/transition_rates_sim_main.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_main.csv', row.names = FALSE)




###### Appendix figures demonstrating how we simulate parental age effects -----


## Figure A6
sim_lab <- c(
  "Sim 1" = "Strongly negative",
  "Sim 2" = "Weakly negative",
  "Sim 3" = "None",
  "Sim 4" = "Weakly positive",
  "Sim 5" = "Strongly positive"
)

tr_p_full <- as_tibble(expand.grid(i = i, P_i = P_i_sim, type = sim_lab))
tr_f_full <- as_tibble(expand.grid(i = i, F_i = F_i_init, type = sim_lab))

fig_a6_1 <- ggplot(tr_p_full) +
  geom_line(aes(i, P_i)) +
  facet_wrap(~ type, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  coord_cartesian(ylim = c(0, 1)) +
  ylab(expression(italic(P[i]))) +
  tt_a234 +
  theme(strip.text.x = element_text(size = 9),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank())

fig_a6_2 <- ggplot(tr_f_full) +
  geom_line(aes(i, F_i)) +
  facet_wrap(~ type, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0.85, 1.15)) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[i]))) +
  tt_a234 + theme(strip.text.x = element_blank())

# plots to grobs
g_a6_1 <- ggplotGrob(fig_a6_1)
g_a6_2 <- ggplotGrob(fig_a6_2)

# arrange panels
fig_a6 <- arrangeGrob(rbind(g_a6_1, g_a6_2, size = 'first'))

# print (change 'quartz' to 'window' if using Windows)
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(fig_a6)

# write figure to file
ggsave('img/Fig_A6.png', fig_a6, height = 4, width = 6, units = 'in', dpi = 300)



## Figure A7
lwd_1 <- 2.1

df_lab <- data.frame(sim = "Sim 1", label = "Parental age")
df_arrow <- data.frame(sim = "Sim 1", x1 = 0.65, x2 = 13.3, y1 = 0.14, y2 = 0.14)

fig_a7_1 <- ggplot(Tr_sim_init) +
  geom_line(aes(i, P_ij, color = rev(j), group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.6, 0.2, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal',
                                barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
             axis.title.x = element_blank(),
             strip.text.x = element_text(size = 9))

fig_a7_2 <- ggplot(Tr_sim_init) +
  geom_line(aes(i, F_ij, color = j, group = j), lwd = lwd_1) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0.3, 1.7)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())

# plots to grobs
g_a7_1 <- ggplotGrob(fig_a7_1)
g_a7_2 <- ggplotGrob(fig_a7_2)

# arrange panels
fig_a7 <- arrangeGrob(rbind(g_a7_1, g_a7_2, size = 'first'))

# print (change 'quartz' to 'window' if using Windows)
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(fig_a7)

# write figure to file
ggsave('img/Fig_A7.png', fig_a7, height = 4, width = 6, units = 'in', dpi = 300)


### Figure A8
fig_a8_1 <- ggplot(Tr_sim) +
  geom_line(aes(i, P_ij, col = rev(j), group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2),
               arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.4, 0.22, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal',
                                barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  strip.text.x = element_text(size = 9))

fig_a8_2 <- ggplot(Tr_sim) +
  geom_line(aes(i, F_ij, col = j, group = j), lwd = lwd_1) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous() +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())

# plots to grobs
g_a8_1 <- ggplotGrob(fig_a8_1)
g_a8_2 <- ggplotGrob(fig_a8_2)

# arrange panels
fig_a8 <- arrangeGrob(rbind(g_a8_1, g_a8_2, size = 'first'))

# print (change 'quartz' to 'window' if using Windows)
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(fig_a8)

# write figure to file
ggsave('img/Fig_A8.png', fig_a8, height = 4, width = 6, units = 'in', dpi = 300)




###### Transition rates for sensitivity analyses (Appendix S1 Section 7) -------

#### Scenario a: vary the population growth rate (lambda) ----------------------

## initial parameters
s <- 10
o <- 20
i <- 1:o
P_i_sim <- rep(0.5, o)
F_i_init <- rep(1, o)


## find parental-age-specific fecundity values yielding desired value of lambda
opt_wrapper_a <- function(lambda) {
  # create life cycles using inits for Sim #3 from main text as baseline
  optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
        P_i = P_i_sim, F_i = F_i_init, P_range = 0, F_range = 0.6,
        lam = lambda, stasis = TRUE)
}

opt_sim1 <- opt_wrapper_a(lambda = 0.6)
opt_sim2 <- opt_wrapper_a(lambda = 0.8)
opt_sim3 <- opt_wrapper_a(lambda = 1.0)
opt_sim4 <- opt_wrapper_a(lambda = 1.2)
opt_sim5 <- opt_wrapper_a(lambda = 1.4)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim5$par


## generate transition rates reflecting desired parental age effect
Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = 0, F_range = 0.6)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = 0, F_range = 0.6)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = 0, F_range = 0.6)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = 0, F_range = 0.6)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = 0, F_range = 0.6)


## full set of simulated transition rates
Tr_sim_scen_a <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'lambda == 0.6'),
  Tr_sim2 %>% mutate(sim = 'lambda == 0.8'),
  Tr_sim3 %>% mutate(sim = 'lambda == 1.0'),
  Tr_sim4 %>% mutate(sim = 'lambda == 1.2'),
  Tr_sim5 %>% mutate(sim = 'lambda == 1.4')
)


## map between age classes and parental age classes
map_ij_sim_scen_a <- Tr_sim_scen_a %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)


## write to file (parental age effects for varying values of lambda)
write.csv(Tr_sim_scen_a, 'dat/transition_rates_sim_scen_a.csv', row.names = F)
write.csv(map_ij_sim_scen_a, 'dat/map_ij_sim_scen_a.csv', row.names = F)



#### Scenario b: parental age effect on survival -------------------------------

## initial parameters
s <- 10
o <- 20
i <- 1:o

P_i_sim <- rep(0.5, o)
F_i_init <- rep(1, o)

F_range <- 0

P_range1 <- 0.5
P_range2 <- 0.25
P_range3 <- 0
P_range4 <- -0.25
P_range5 <- -0.5


## find parental-age-specific fecundity values yielding desired value of lambda
opt_wrapper_b <- function(P_range) {
  optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
        P_i = P_i_sim, F_i = F_i_init, P_range = P_range, F_range = F_range,
        lam = 1, stasis = TRUE)
}

opt_sim1 <- opt_wrapper_b(P_range = P_range1)
opt_sim2 <- opt_wrapper_b(P_range = P_range2)
opt_sim3 <- opt_wrapper_b(P_range = P_range3)
opt_sim4 <- opt_wrapper_b(P_range = P_range4)
opt_sim5 <- opt_wrapper_b(P_range = P_range5)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim5$par


## generate transition rates reflecting desired parental age effect
Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = P_range1,
                             F_range = F_range)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = P_range2,
                             F_range = F_range)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = P_range3,
                             F_range = F_range)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = P_range4,
                             F_range = F_range)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = P_range5,
                             F_range = F_range)


## full set of simulated transition rates
Tr_sim <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'Sim 1'),
  Tr_sim2 %>% mutate(sim = 'Sim 2'),
  Tr_sim3 %>% mutate(sim = 'Sim 3'),
  Tr_sim4 %>% mutate(sim = 'Sim 4'),
  Tr_sim5 %>% mutate(sim = 'Sim 5')
)


## map between age classes and parental age classes
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)


## write to file (parental age effect on survival)
write.csv(Tr_sim, 'dat/transition_rates_sim_scen_b.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_scen_b.csv', row.names = FALSE)




##### Can selection gradients increase with age given strongly positive --------
# parental age effect?

## initial parameters
s <- 20
o <- 20
i <- 1:o

P_range <- -0.9
F_range <- -0.9

P_i_sim <- logit_inv(-1.1 + 0.2 * i)
F_i_init <- 0.1 + 3 * i


## find parental-age-specific fecundity values yielding desired value of lambda
opt_sim <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
                 P_i = P_i_sim, F_i = F_i_init, P_range = P_range,
                 F_range = F_range, lam = 1)

F_i_sim <- F_i_init * opt_sim$par


## generate transition rates reflecting desired parental age effect
Tr_sim <- SimulateParEffect(P_i_sim, F_i_sim, s = s, P_range = P_range,
                            F_range = F_range)


## map between age classes and parental age classes (1:1 here)
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  dplyr::select(i, j)


## write to file (strongly positive parental age effect)
write.csv(Tr_sim, 'dat/transition_rates_sim_strongly_pos.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_strongly_pos.csv', row.names = FALSE)







