

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



##### Lemna minor data --------------------------------------------------------

## read data
# src_url <- "http://datadryad.org/bitstream/handle/10255/dryad.71128/BarksAndLaird_RawData_Phase2.csv"
# dat_lemna <- read.csv(url(src_url), stringsAsFactors = FALSE)
# save(dat_lemna, file = 'dat/BarksAndLaird_RawData_Phase2.RData')
load('dat/BarksAndLaird_RawData_Phase2.RData')


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


p1 <- ggplot(df, aes(age, pred, col = j, group = j)) +
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
print(p1)

# write to file
# ggsave('img/Fig_smoothfert.png', p1, height = 6, width = 6.25, units = 'in', dpi = 300)



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


p2 <- ggplot(dfs, aes(age, pred, col = j, group = j)) +
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
print(p2)

# write to file
# ggsave('img/Fig_smoothsurv.png', p2, height = 6, width = 6.25, units = 'in', dpi = 300)



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
# write.csv(age_map_lemna, "dat/age_map_lemna_s8.csv", row.names = FALSE)
# write.csv(Tr_lemna_flat, "dat/transition_rates_lemna_s8.csv", row.names = FALSE)
# write.csv(Tr_lemna_flat_plot, "dat/transition_rates_plot_lemna_s8.csv", row.names = FALSE)



# df_means <- dat_lemna_flat %>% 
#   group_by(age, j) %>% 
#   summarize(fecund = mean(fecund)) %>% 
#   ungroup() %>% 
#   filter(age <= 30) %>% 
#   mutate(j = factor(j, labels = c("Parental age group 1 (1-3)",
#                                   "Parental age group 2 (4-5)",
#                                   "Parental age group 3 (6-8)",
#                                   "Parental age group 4 (9-11)",
#                                   "Parental age group 5 (12-14)",
#                                   "Parental age group 6 (15-17)",
#                                   "Parental age group 7 (18-20)",
#                                   "Parental age group 8 (21-30)")))
# 
# df_counts <- dat_lemna_flat %>% 
#   count(age, j, fecund) %>% 
#   filter(age <= 30) %>% 
#   filter(fecund < 3) %>% 
#   mutate(j = factor(j, labels = c("Parental age group 1 (1-3)",
#                                   "Parental age group 2 (4-5)",
#                                   "Parental age group 3 (6-8)",
#                                   "Parental age group 4 (9-11)",
#                                   "Parental age group 5 (12-14)",
#                                   "Parental age group 6 (15-17)",
#                                   "Parental age group 7 (18-20)",
#                                   "Parental age group 8 (21-30)")))
# 
# p1 <- ggplot(df_counts, aes(x = age, y = fecund)) +
#   geom_point(aes(size = n), shape = 19, col = "gray70") +
#   geom_line(data = rename(df_means, h = j), aes(group = h), alpha = 0.2, size = 0.4) +
#   geom_line(data = df_means, size = 0.7) +
#   scale_y_continuous(breaks = 0:2, expand = c(0.1, 0)) +
#   scale_size(name = "Count") +
#   facet_wrap(~ j, ncol = 2, dir = "v") +
#   labs(x = "Age (days)", y = "Fecundity") +
#   tt2
# 
# dev.off()
# quartz(height = 6.5, width = 6.25, dpi = 150)
# print(p1)
# 
# # ggsave('img/Fig_lemna.png', p1, height = 6.5, width = 6.25, units = 'in', dpi = 300)




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


p1 <- ggplot(filter(df_means, age <= 30), aes(x = age, y = surv, size = n)) +
  geom_point(shape = 19, col = "grey30", alpha = 0.4) +
  geom_point(shape = 1, alpha = 0.6) +
  facet_wrap(~ j, ncol = 2, dir = "v") +
  scale_y_continuous(limits = c(0.35, 1.05), breaks = seq(0.4, 1, 0.2)) +
  scale_size(range = c(1, 3), name = "N", guid = F) +
  labs(x = "Age (days)", y = "Survival") +
  tt2 +
  theme(axis.title.y = element_text(margin = margin(0, 5, 0, 0)),
        strip.text = element_text(size = 8))

p2 <- ggplot(filter(df_means, age <= 30), aes(x = age, y = fert, size = n)) +
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


g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)

g <- arrangeGrob(g1, g2, nrow = 1, widths = c(1, 1.27))

dev.off()
quartz(height = 5, width = 6.25, dpi = 150)
grid.arrange(g)

# write to file
# ggsave('img/Fig_lemna_trans.png', g, height = 5, width = 6.25, units = 'in', dpi = 300)







##### Simulated life cycles ----------------------------------------------------

### global variables for (almost) all simulated life cycles
s <- 10
o <- 20
i <- 1:o


### New simulated life cycles for revision
P_i_sim <- rep(0.5, o)
F_i_init <- rep(1, o)

P_range1 <- 0
P_range2 <- 0
P_range3 <- 0
P_range4 <- 0
P_range5 <- 0

F_range1 <- 0.6
F_range2 <- 0.3
F_range3 <- 0.0
F_range4 <- -0.3
F_range5 <- -0.6


# seq(0.4, 1.6, length.out = 9)
# seq(0.7, 1.3, length.out = 9)


opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range1, F_range = F_range1, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range2, F_range = F_range2, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range3, F_range = F_range3, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range4, F_range = F_range4, lam = 1)
opt_sim5 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range5, F_range = F_range5, lam = 1)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim4$par

# generate parental age effect with initial fecundity values (for Appendix only)
Tr_sim_init1 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range1, F_range = F_range1)
Tr_sim_init2 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range2, F_range = F_range2)
Tr_sim_init3 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range3, F_range = F_range3)
Tr_sim_init4 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range4, F_range = F_range4)
Tr_sim_init5 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range5, F_range = F_range5)

# generate parental age effect with adjusted fecundity values
Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = P_range1, F_range = F_range1)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = P_range2, F_range = F_range2)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = P_range3, F_range = F_range3)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = P_range4, F_range = F_range4)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = P_range5, F_range = F_range5)

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

## age map for sim species
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)


## write to file (parental age effect on survival)
write.csv(Tr_sim, 'dat/transition_rates_sim_revise.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_revise.csv', row.names = FALSE)




### Figure A2
sim_lab <- c(
  "Sim 1" = "Strongly negative",
  "Sim 2" = "Weakly negative",
  "Sim 3" = "None",
  "Sim 4" = "Weakly positive",
  "Sim 5" = "Strongly positive"
)

tr_p_full <- as_tibble(expand.grid(i = i, P_i = P_i_sim, type = sim_lab))
tr_f_full <- as_tibble(expand.grid(i = i, F_i = F_i_init, type = sim_lab))

pa2_1 <- ggplot(tr_p_full) +
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

pa2_2 <- ggplot(tr_f_full) +
  geom_line(aes(i, F_i)) +
  facet_wrap(~ type, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0.85, 1.15)) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[i]))) +
  tt_a234 + theme(strip.text.x = element_blank())

# plots to grobs
ga2_1 <- ggplotGrob(pa2_1)
ga2_2 <- ggplotGrob(pa2_2)

# arrange panels
ga2 <- arrangeGrob(rbind(ga2_1, ga2_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(ga2)

# write to file
ggsave('img/Fig_A2.png', ga2, height = 4, width = 6, units = 'in', dpi = 300)





## Figure A3
lwd_1 <- 2.1

df_lab <- data.frame(sim = "Sim 1", label = "Parental age")
df_arrow <- data.frame(sim = "Sim 1", x1 = 0.65, x2 = 13.3, y1 = 0.14, y2 = 0.14)

pa3_1 <- ggplot(Tr_sim_init) +
  geom_line(aes(i, P_ij, color = rev(j), group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.6, 0.2, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
             axis.title.x = element_blank(),
             strip.text.x = element_text(size = 9))

pa3_2 <- ggplot(Tr_sim_init) +
  geom_line(aes(i, F_ij, color = j, group = j), lwd = lwd_1) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(limits = c(0.3, 1.7)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())


# plots to grobs
ga3_1 <- ggplotGrob(pa3_1)
ga3_2 <- ggplotGrob(pa3_2)

# arrange panels
ga3 <- arrangeGrob(rbind(ga3_1, ga3_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(ga3)

# write to file
ggsave('img/Fig_A3.png', ga3, height = 4, width = 6, units = 'in', dpi = 300)



### Figure A4
pa4_1 <- ggplot(Tr_sim) +
  geom_line(aes(i, P_ij, col = rev(j), group = j), lwd = lwd_1) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(0.4, 0.22, label = label), hjust = 0, size = 2.6) +
  facet_wrap(~ sim, nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(italic(P[ij]))) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2), limits = c(0, 1)) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  tt_a234 + theme(axis.text.x = element_blank(),
                  axis.title.x = element_blank(),
                  strip.text.x = element_text(size = 9))

pa4_2 <- ggplot(Tr_sim) +
  geom_line(aes(i, F_ij, col = j, group = j), lwd = lwd_1) +
  facet_wrap(~ sim, nrow = 1) +
  scale_x_continuous(limits = c(0, 20)) +
  scale_y_continuous() +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  xlab(expression(Age~(italic(i)))) +
  ylab(expression(italic(F[ij]))) +
  tt_a234 + theme(strip.text.x = element_blank())


# plots to grobs
ga4_1 <- ggplotGrob(pa4_1)
ga4_2 <- ggplotGrob(pa4_2)

# arrange panels
ga4 <- arrangeGrob(rbind(ga4_1, ga4_2, size = 'first'))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 4, width = 6, dpi = 150)
grid.arrange(ga4)

# write to file
ggsave('img/Fig_A4.png', ga4, height = 4, width = 6, units = 'in', dpi = 300)






########## Transition rates for sensitivity analyses (Appendix S1 Section 7) ---

### Scenario a: vary the population growth rate (lambda) -----------------------
opt_wrapper <- function(lambda) {
  # create life cycles using inits for Sim #3 from main text as baseline
  optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s,
        P_i = P_i_sim, F_i = F_i_init, P_range = 0, F_range = 0.6,
        lam = lambda, stasis = TRUE)
}

opt_sim1 <- opt_wrapper(lambda = 0.6)
opt_sim2 <- opt_wrapper(lambda = 0.8)
opt_sim3 <- opt_wrapper(lambda = 1.0)
opt_sim4 <- opt_wrapper(lambda = 1.2)
opt_sim5 <- opt_wrapper(lambda = 1.4)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim5$par

Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = 0, F_range = 0.6)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = 0, F_range = 0.6)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = 0, F_range = 0.6)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = 0, F_range = 0.6)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = 0, F_range = 0.6)

## full set of sim transition rates
Tr_sim_scen_a <- rbind.data.frame(
  Tr_sim1 %>% mutate(sim = 'lambda == 0.6'),
  Tr_sim2 %>% mutate(sim = 'lambda == 0.8'),
  Tr_sim3 %>% mutate(sim = 'lambda == 1.0'),
  Tr_sim4 %>% mutate(sim = 'lambda == 1.2'),
  Tr_sim5 %>% mutate(sim = 'lambda == 1.4')
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





### Scenario b: parental age effect on survival
P_i_sim <- rep(0.5, o)
F_i_init <- rep(1, o)

P_range1 <- 0.5
P_range2 <- 0.25
P_range3 <- 0
P_range4 <- -0.25
P_range5 <- -0.5

F_range1 <- 0
F_range2 <- 0
F_range3 <- 0
F_range4 <- 0
F_range5 <- 0

seq(0.5, -0.5, length.out = 10)
seq(0.25, -0.25, length.out = 10)

# seq(0.5, 1.6, length.out = 9)
# seq(0.7, 1.3, length.out = 9)

opt_sim1 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range1, F_range = F_range1, lam = 1)
opt_sim2 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range2, F_range = F_range2, lam = 1)
opt_sim3 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range3, F_range = F_range3, lam = 1)
opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range4, F_range = F_range4, lam = 1)
opt_sim5 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim, F_i = F_i_init,
                  P_range = P_range5, F_range = F_range5, lam = 1)

F_i_sim1 <- F_i_init * opt_sim1$par
F_i_sim2 <- F_i_init * opt_sim2$par
F_i_sim3 <- F_i_init * opt_sim3$par
F_i_sim4 <- F_i_init * opt_sim4$par
F_i_sim5 <- F_i_init * opt_sim5$par

# generate parental age effect with initial fecundity values (for Appendix only)
Tr_sim_init1 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range1, F_range = F_range1)
Tr_sim_init2 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range2, F_range = F_range2)
Tr_sim_init3 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range3, F_range = F_range3)
Tr_sim_init4 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range4, F_range = F_range4)
Tr_sim_init5 <- SimulateParEffect(P_i_sim, F_i_init, s = s, P_range = P_range5, F_range = F_range5)

# generate parental age effect with adjusted fecundity values
Tr_sim1 <- SimulateParEffect(P_i_sim, F_i_sim1, s = s, P_range = P_range1, F_range = F_range1)
Tr_sim2 <- SimulateParEffect(P_i_sim, F_i_sim2, s = s, P_range = P_range2, F_range = F_range2)
Tr_sim3 <- SimulateParEffect(P_i_sim, F_i_sim3, s = s, P_range = P_range3, F_range = F_range3)
Tr_sim4 <- SimulateParEffect(P_i_sim, F_i_sim4, s = s, P_range = P_range4, F_range = F_range4)
Tr_sim5 <- SimulateParEffect(P_i_sim, F_i_sim5, s = s, P_range = P_range5, F_range = F_range5)

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

## age map for sim species
map_ij_sim <- Tr_sim %>% 
  filter(j == 1) %>% 
  group_by(sim) %>% 
  mutate(j = map_ij(F_ij, s = s)) %>% 
  ungroup() %>% 
  dplyr::select(sim, i, j)


## write to file (parental age effect on survival)
write.csv(Tr_sim, 'dat/transition_rates_sim_scen_b.csv', row.names = FALSE)
write.csv(map_ij_sim, 'dat/map_ij_sim_scen_b.csv', row.names = FALSE)








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






### strongly increasing?
s <- 20
o <- 20
i <- 1:o

P_range <- -0.9
F_range <- -0.9
P_i_sim4 <- logit_inv(-1.1 + 0.2 * i)
F_i_init4 <- 0.1 + 3 * i

opt_sim4 <- optim(0.5, opt_fn, method = "Brent", lower = 0, upper = 200, s = s, P_i = P_i_sim4, F_i = F_i_init4,
                  P_range = P_range, F_range = F_range, lam = 1, stasis = FALSE)

F_i_sim4 <- F_i_init4 * opt_sim4$par
Tr_sim4 <- SimulateParEffect(P_i_sim4, F_i_sim4, s = s, P_range = P_range, F_range = F_range)

ls <- 2.5
la <- 0.9
las <- 0.9



q_i_sim4 <- map_ij(F_i_sim4, s = s)
A_par_sim <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = FALSE)
mod_out_sim <- analyze_mod(A_par_sim, stasis = FALSE)
A_par_sim$A_par[380:400,380:400]

age_df_surv <- mod_out_sim$df_age %>% 
  as_tibble() %>% 
  select(i, sens_surv_ref, sens_surv_par) %>% 
  gather(type, sens_surv, sens_surv_ref:sens_surv_par) %>% 
  mutate(type = ifelse(type == "sens_surv_ref", "ref", "par"))

age_df <- mod_out_sim$df_age %>% 
  as_tibble() %>% 
  select(i, sens_fert_ref, sens_fert_par) %>% 
  gather(type, sens_fert, sens_fert_ref:sens_fert_par) %>% 
  mutate(type = ifelse(type == "sens_fert_ref", "ref", "par")) %>% 
  left_join(age_df_surv)


df_arrow <- data.frame(x1 = 14.3, x2 = 19.9, y1 = 0.151, y2 = 0.151)



p1 <- ggplot(filter(Tr_sim4, i < 20), aes(i, P_ij, color = j, group = j)) +
  geom_line(size = ls, alpha = la) +
  geom_line(data = mod_out_sim$df_age, inherit.aes = FALSE, aes(i, surv_ref),
            size = 0.6, linetype = 2) +
  # geom_segment(data = df_arrow, inherit.aes = F, aes(x = x1, y = y1, xend = x2, yend = y2),
  #              arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  # annotate('text', label = 'Parental age', x = 19.65, y = 0.195, hjust = 1, size = 3) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  scale_color_gradient(low = col_low, high = col_upp, guide = FALSE) +
  labs(x = NULL, y = expression(paste('Survival (', italic(P), ')'))) +
  # guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal',
  #                               barwidth = unit(1.8, 'cm'),
  #                               barheight = unit(0.5, 'cm'))) +
  tt2 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0.95, 0.04),
        legend.justification = c(1, 0),
        legend.margin = margin(0, 0, 0, 0))

p2 <- ggplot(Tr_sim4, aes(i, F_ij, color = j, group = j)) +
  geom_line(size = ls, alpha = la) +
  geom_line(data = mod_out_sim$df_age, inherit.aes = FALSE, aes(i, fert_ref),
            size = 0.6, linetype = 2) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_continuous(limits = c(0, 30)) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  labs(x = NULL, y = expression(paste('Fecundity (', italic(F), ')'))) +
  tt2 +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

p3 <- ggplot(age_df, aes(i, sens_surv, linetype = type)) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.01, 0.1, 1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = FALSE) +
  labs(x = "Age", y = expression(paste('Sensitivity of λ to ', italic(P)))) +
  tt2 +
  theme(legend.position = c(0.98, 0.99),
        legend.justification = c(1, 1),
        legend.text.align = 0,
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank())

p4 <- ggplot(age_df, aes(i, sens_fert, linetype = type)) +
  geom_line(size = 0.6) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.00001, 0.001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = FALSE) +
  labs(x = "Age", y = expression(paste('Sensitivity of λ to ', italic(F)))) +
  tt2

g1 <- ggplotGrob(p1)
g2 <- ggplotGrob(p2)
g3 <- ggplotGrob(p3)
g4 <- ggplotGrob(p4)

g <- arrangeGrob(g1, g2, g3, g4, nrow = 2, heights = c(0.86, 1))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 5, width = 5.75, dpi = 150)
grid.arrange(g)


# ggsave('img/Fig_6.png', g, height = 5, width = 5.75, units = 'in', dpi = 300)







