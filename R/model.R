

## Must have the following packages installed
# install.packages(c("Matrix", "dplyr", "tidyr", "tibble", "mgcv", "popbio",
#                    "popdemo", "ggplot2", "gridExtra"))


## load relevant packages
library(Matrix)
library(dplyr)
library(tidyr)
library(tibble)
library(popbio)
library(popdemo)
library(ggplot2)
library(gridExtra)
source('R/functions.R')
source('R/plot-themes.R')



## read data
map_ij_lemna <- read.csv('dat/age_map_lemna_s8.csv', stringsAsFactors = FALSE)
Tr_lemna <- read.csv('dat/transition_rates_lemna_s8.csv', stringsAsFactors = FALSE)
Tr_lemna_plot <- read.csv('dat/transition_rates_plot_lemna_s8.csv', stringsAsFactors = FALSE)

Tr_sim <- read.csv('dat/transition_rates_sim_revise.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_revise.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')
Tr_sim5 <- filter(Tr_sim, sim == 'Sim 5')


## q_i by sim
q_i_lemna <- map_ij_lemna$j
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j
q_i_sim5 <- filter(map_ij_sim, sim == 'Sim 5')$j


## construct A_par via vec-permutation method
A_par_lemna <- build_Apar(Tr = Tr_lemna, q_i = q_i_lemna)
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = TRUE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = TRUE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = TRUE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = TRUE)
A_par_sim5 <- build_Apar(Tr = Tr_sim5, q_i = q_i_sim5, stasis = TRUE)


## analysis of A_par and A_ref
mod_out_lemna <- analyze_mod(A_par_lemna)
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = TRUE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = TRUE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = TRUE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = TRUE)
mod_out_sim5 <- analyze_mod(A_par_sim5, stasis = TRUE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4'),
  mod_out_sim5$df_age %>% mutate(sim = 'Sim 5')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4'),
  mod_out_sim5$df_par %>% mutate(sim = 'Sim 5')
)


# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  # mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))


## age-specific sensitivities, tidy
df_sens_age_lemna <- mod_out_lemna$df_age %>% 
  dplyr::select(i, ref = sens_surv_ref, par = sens_surv_par) %>% 
  gather(type, sens_surv, -i) %>% 
  mutate(sens_fert = c(mod_out_lemna$df_age$sens_fert_ref,
                       mod_out_lemna$df_age$sens_fert_par))


df_sens_age_sim <- mod_age_full %>% 
  dplyr::select(sim, i, ref = sens_surv_ref, par = sens_surv_par) %>% 
  gather(type, sens_surv, -i, -sim) %>% 
  mutate(sens_fert = c(mod_age_full$sens_fert_ref,
                       mod_age_full$sens_fert_par))


## check difference in sensitivities between par and ref, at early and late ages
sens_check <- mod_age_full %>% 
  dplyr::select(i, sim, starts_with('sens')) %>% 
  group_by(sim) %>% 
  filter(i %in% c(min(i):(min(i)+1), (max(i)-1):max(i))) %>% 
  ungroup() %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref)




##### Figure 2 -----------------------------------------------------------------
# plot settings
ls <- 2.5
la <- 0.95
las <- 0.9

df_arrow <- data.frame(x1 = 1.1, x2 = 13, y1 = 0.2, y2 = 0.2)

# panels
f2_1 <- ggplot(Tr_lemna_plot) +
  geom_line(aes(i, P_ij, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, surv_ref), linetype = 2, size = 0.7) +
  # geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2),
  #              arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  # annotate('text', label = 'Parental age', x = 1.05, y = 0.265, hjust = 0, size = 3) +
  annotate('text', label = '(a)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  # guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal',
  #                               barwidth = unit(1.7, 'cm'), barheight = unit(0.5, 'cm'))) +
  xlab(NULL) + ylab(expression(paste('Survival (', italic(P), ')'))) +
  tt2 + theme(axis.text.x = element_blank(),
              legend.position = c(0.04, 0.04),
              legend.justification = c(0, 0),
              legend.margin = margin(0, 0, 0, 0))

f2_2 <- ggplot(Tr_lemna_plot) +
  geom_line(aes(i, F_ij, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, fert_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(b)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, max(c(1, max(Tr_lemna_plot$F_ij))))) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab(' ') + ylab(expression(paste('Fecundity (', italic(F), ')'))) +
  tt2

f2_3 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, w, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, w_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(c)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(NULL) + ylab(expression(paste('Stable distribution (', italic(w), ')'))) +
  tt2 + theme(axis.text.x = element_blank())

f2_4 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, v/v[1], col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, v_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(d)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10() +
  xlab('Age (days)') + ylab(expression(paste('Reproductive value (', italic(v), ')'))) +
  tt2

f2_5 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, sens_surv, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, sens_surv_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(e)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(NULL) + ylab(expression(paste('Sensitivity of λ to ', italic(P)))) +
  tt2 + theme(axis.text.x = element_blank())

f2_6 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, sens_fert, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, sens_fert_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(f)', x = 29.5, y = Inf, hjust = 0.6, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(' ') + ylab(expression(paste('Sensitivity of λ to ', italic(F)))) +
  tt2

# plots to grobs
g2_1 <- ggplotGrob(f2_1)
g2_2 <- ggplotGrob(f2_2)
g2_3 <- ggplotGrob(f2_3)
g2_4 <- ggplotGrob(f2_4)
g2_5 <- ggplotGrob(f2_5)
g2_6 <- ggplotGrob(f2_6)

# arrange panels
g2_r1 <- cbind(g2_1, g2_3, g2_5, size = 'last')
g2_r2 <- cbind(g2_2, g2_4, g2_6, size = 'last')

g2 <- arrangeGrob(g2_r1, g2_r2, nrow = 2, heights = c(0.84, 1))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 4, width = 7.25, dpi = 150)
grid.arrange(g2)

# write to file
# ggsave('img/Fig_3.png', g2, height = 4, width = 7.25, units = 'in', dpi = 300)




##### Figure 3 -----------------------------------------------------------------
# age-specific sensitivities from A_ref and A_par
w_plot <- mod_out_lemna$df_age %>% 
  gather(type, w, w_ref:w_par) %>% 
  dplyr::select(i, type, w)

# age-specific reproductive values from A_ref and A_par
v_plot <- mod_out_lemna$df_age %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(type) %>% 
  ungroup() %>% 
  dplyr::select(i, type, v)

# create plots
f3_1 <- ggplot(w_plot) +
  geom_line(aes(i, w, linetype = type)) +
  annotate('text', label = '(a)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(breaks = c(0.1, 0.001, 0.00001), labels = LabelFn) +
  xlab(NULL) + ylab(expression(paste('Stable age distribution (', italic(w[i]), ')'))) +
  tt3 + theme(axis.text.x = element_blank())

f3_2 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type), show.legend = FALSE) +
  annotate('text', label = '(b)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(labels = LabelFn) +
  xlab('Age (days)') +
  ylab(expression(paste('Reproductive value (', italic(v[i]), ')'))) +
  tt3 + theme(plot.margin = unit(c(0.21, 0.15, 0.15, 0.15), 'cm'))

f3_3 <- ggplot(df_sens_age_lemna, aes(i, sens_surv, linetype = type)) +
  geom_line(show.legend = FALSE) +
  annotate('text', label = '(c)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(breaks = c(0.1, 0.0001, 0.0000001), labels = LabelFn) +
  coord_cartesian(ylim = c(0.00000001, 0.25)) +
  xlab(NULL) +
  ylab(expression(paste('Sensitivity of ', italic(lambda), ' to ', italic(P[i])))) +
  tt3 + theme(axis.text.x = element_blank())

f3_4 <- ggplot(df_sens_age_lemna, aes(i, sens_fert, linetype = type)) +
  geom_line(show.legend = FALSE) +
  annotate('text', label = '(d)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(breaks = c(0.1, 0.0001, 0.0000001), labels = LabelFn) +
  coord_cartesian(ylim = c(0.00000001, 0.25)) +
  xlab('Age (days)') +
  ylab(expression(paste('Sensitivity of ', italic(lambda), ' to ', italic(F[i])))) +
  tt3 + theme(plot.margin = unit(c(0.21, 0.15, 0.15, 0.15), 'cm'))

# panels to grobs
g3_1 <- ggplotGrob(f3_1)
g3_2 <- ggplotGrob(f3_2)
g3_3 <- ggplotGrob(f3_3)
g3_4 <- ggplotGrob(f3_4)

# arrange panels
g3 <- arrangeGrob(g3_1, g3_3, g3_2, g3_4, nrow = 2, heights = c(1, 1.155))

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 4.2, width = 6, dpi = 190)
grid.arrange(g3)

# write to file
# ggsave('img/Fig_4.png', g3, height = 4.2, width = 6, units = 'in', dpi = 300)




##### Figure 4 -----------------------------------------------------------------

# # plot settings
# lwd_1 <- 0.6
# lwd_2 <- 1.5
# y_expand <- c(0.05, 0.05)
# 
# df_arrow <- data.frame(sim = "Sim 1",
#                        x1 = 1.35, x2 = 11.5, y1 = 0.22, y2 = 0.22)
# 
# df_lab <- data.frame(sim = "Sim 1",
#                      label = 'Parental age',
#                      stringsAsFactors = FALSE)
# 
# # create plots (by row)
# p4_1 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
#   geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
#   geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
#   guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
#   facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
#   ylab(expression(italic('P'))) +
#   tt4 +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         legend.position = c(0, 0),
#         legend.justification = c(0, 0),
#         legend.margin = margin(0, 0, 0, 0))
# 
# p4_2 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, 0.3)) +
#   facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
#   ylab(expression(italic('F'))) +
#   tt4 +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         strip.text.x = element_blank())
# 
# p4_3 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_y_log10(expand = y_expand, breaks = c(0.00001, 0.1), labels = LabelFn) +
#   facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
#   ylab(expression(italic('w'))) +
#   tt4 + theme(axis.title.x = element_blank(),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
#               strip.text.x = element_blank())
# 
# p4_4 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_y_log10(expand = y_expand, breaks = c(0.1, 1, 10),  labels = LabelFn) +
#   facet_wrap(~ sim, nrow = 1) +
#   ylab(expression(italic('v'))) +
#   tt4 + theme(axis.title.x = element_blank(),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
#               strip.text.x = element_blank())
# 
# p4_5 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_y_log10(expand = y_expand, breaks = c(0.00001, 0.1), labels = LabelFn) +
#   facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
#   ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
#   tt4 + theme(axis.title.x = element_blank(),
#               axis.text.x = element_blank(),
#               axis.ticks.x = element_blank(),
#               strip.text.x = element_blank())
# 
# p4_6 <- ggplot() +
#   geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
#   geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
#   scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
#   scale_x_continuous(breaks = c(1, 20)) +
#   scale_y_log10(expand = y_expand, breaks = c(0.00001, 0.1), labels = LabelFn) +
#   facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
#   xlab('Age') +
#   ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
#   tt4 +
#   theme(strip.text.x = element_blank(),
#         plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))
# 
# # plots to grobs
# g4_1 <- ggplotGrob(p4_1)
# g4_2 <- ggplotGrob(p4_2)
# g4_3 <- ggplotGrob(p4_3)
# g4_4 <- ggplotGrob(p4_4)
# g4_5 <- ggplotGrob(p4_5)
# g4_6 <- ggplotGrob(p4_6)
# 
# # arrange panels
# g4 <- arrangeGrob(rbind(g4_1, g4_2, g4_3, g4_4, g4_5, g4_6, size = 'last'))
# 
# # # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(g4)
# 
# # write to file
# # ggsave('img/Fig_4.png', g4, height = 7, width = 5.6, units = 'in', dpi = 300)



##### Figure 5 -----------------------------------------------------------------
# inputs and age-specific outputs for simulated life cycles
sim_lab <- c(
  "Sim 1" = "Strongly negative",
  "Sim 2" = "Weakly negative",
  "Sim 3" = "None",
  "Sim 4" = "Weakly positive",
  "Sim 5" = "Strongly positive"
)

# plot settings
lwd_1 <- 0.6
lwd_2 <- 1.8
y_expand <- c(0.05, 0.05)


df_arrow <- data.frame(sim = "Sim 1",
                       x1 = 1.35, x2 = 11, y1 = 0.16, y2 = 0.16)

df_lab <- data.frame(sim = "Sim 1",
                     label = 'Parental age',
                     stringsAsFactors = FALSE)

p4_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(expand = c(0.08, 0)) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(paste("Survival (", italic('P'), ")"))) +
  tt5 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

p4_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(expand = c(0.08, 0), breaks = c(1, 20)) +
  scale_y_continuous(limits = c(0, 1.15), breaks = seq(0, 0.8, 0.4)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  labs(x = NULL, y = expression(paste("Fecundity (", italic("F"), ")"))) +
  tt5 +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g4_1 <- ggplotGrob(p4_1)
g4_2 <- ggplotGrob(p4_2)


# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_surv_diff) < 1e-9, "=", lab))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_fert_diff) < 1e-9, "=", lab))

v_plot <- mod_age_full %>% 
  filter(i < 20) %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()


# create plots (by row)
p5_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type), size = 0.6) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  scale_x_continuous(expand = c(0.08, 0), limits = c(1, 20)) +
  # scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Rep. value (', italic(v), ')'))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank(),
              legend.justification = c(0, 0),
              legend.position = c(1, 0),)

p5_2 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_surv, linetype = type), size = 0.6) +
  geom_text(data = sens_check_surv, aes(i, 2, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.00001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(P)))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p5_3 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_fert, linetype = type), size = 0.6) +
  geom_text(data = sens_check_fert, aes(i, 2, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), breaks = c(1, 20), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.00001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(F)))) +
  tt5 + theme(strip.text.x = element_blank())

# plots to grobs
g5_1 <- ggplotGrob(p5_1)
g5_2 <- ggplotGrob(p5_2)
g5_3 <- ggplotGrob(p5_3)

# arrange panels
g5 <- arrangeGrob(rbind(g4_1, g4_2, g5_1, g5_2, g5_3, size = 'last'))
g5$layout$clip[g5$layout$name == "panel"] <- "off"


# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 6.5, width = 6.25, dpi = 140)
grid.arrange(g5)

# write to file
# ggsave('img/Fig_5.png', g5, height = 6.5, width = 6.25, units = 'in', dpi = 300)




################################################################################
################################################################################
#### Plots for sensitivity analyses in Appendix S1 Section 7 -------------------
################################################################################
################################################################################


### global plot settings for remaining plots -----------------------------------
lwd_1 <- 0.6
lwd_2 <- 1.5
y_expand <- c(0.05, 0.05)

df_arrow <- data.frame(sim = "Sim 1",
                       x1 = 1.35, x2 = 11.5, y1 = 0.22, y2 = 0.22,
                       stringsAsFactors = FALSE)

df_lab <- data.frame(sim = "Sim 1",
                     label = 'Parental age',
                     stringsAsFactors = FALSE)
################################################################################






### Scenario a: vary the population growth rate (lambda) -----------------------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_a.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_a.csv', stringsAsFactors = FALSE)


## number of parental age classes
s <- length(unique(map_ij_sim$j))

## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'lambda == 0.6')
Tr_sim2 <- filter(Tr_sim, sim == 'lambda == 0.8')
Tr_sim3 <- filter(Tr_sim, sim == 'lambda == 1.0')
Tr_sim4 <- filter(Tr_sim, sim == 'lambda == 1.2')
Tr_sim5 <- filter(Tr_sim, sim == 'lambda == 1.4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'lambda == 0.6')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'lambda == 0.8')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'lambda == 1.0')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'lambda == 1.2')$j
q_i_sim5 <- filter(map_ij_sim, sim == 'lambda == 1.4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = TRUE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = TRUE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = TRUE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = TRUE)
A_par_sim5 <- build_Apar(Tr = Tr_sim5, q_i = q_i_sim5, stasis = TRUE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = TRUE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = TRUE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = TRUE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = TRUE)
mod_out_sim5 <- analyze_mod(A_par_sim5, stasis = TRUE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'lambda == 0.6'),
  mod_out_sim2$df_age %>% mutate(sim = 'lambda == 0.8'),
  mod_out_sim3$df_age %>% mutate(sim = 'lambda == 1.0'),
  mod_out_sim4$df_age %>% mutate(sim = 'lambda == 1.2'),
  mod_out_sim5$df_age %>% mutate(sim = 'lambda == 1.4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'lambda == 0.6'),
  mod_out_sim2$df_par %>% mutate(sim = 'lambda == 0.8'),
  mod_out_sim3$df_par %>% mutate(sim = 'lambda == 1.0'),
  mod_out_sim4$df_par %>% mutate(sim = 'lambda == 1.2'),
  mod_out_sim5$df_par %>% mutate(sim = 'lambda == 1.4')
)


## age-specific sensitivities, tidy
df_sens_age_sim <- mod_age_full %>% 
  dplyr::select(sim, i, ref = sens_surv_ref, par = sens_surv_par) %>% 
  gather(type, sens_surv, -i, -sim) %>% 
  mutate(sens_fert = c(mod_age_full$sens_fert_ref,
                       mod_age_full$sens_fert_par))


## check difference in sensitivities between par and ref, at early and late ages
sens_check <- mod_age_full %>% 
  dplyr::select(i, sim, starts_with('sens')) %>% 
  group_by(sim) %>% 
  filter(i %in% c(min(i):(min(i)+1), (max(i)-1):max(i))) %>% 
  ungroup() %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref)

mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))


## check difference in sensitivities between par and ref, at early and late ages
sens_check <- mod_age_full %>% 
  dplyr::select(i, sim, starts_with('sens')) %>% 
  group_by(sim) %>% 
  filter(i %in% c(min(i):(min(i)+1), (max(i)-1):max(i))) %>% 
  ungroup() %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref)


# inputs and age-specific outputs for simulated life cycles
lwd_1 <- 0.6
lwd_2 <- 1.8

p4_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1, labeller = label_parsed) +
  ylab(expression(paste("Survival (", italic('P'), ")"))) +
  tt5 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0),
        strip.text.x = element_text(size = 9.8))

p4_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  labs(x = NULL, y = expression(paste("Fecundity (", italic("F"), ")"))) +
  tt5 +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g4_1 <- ggplotGrob(p4_1)
g4_2 <- ggplotGrob(p4_2)


# plot settings
lwd_1 <- 0.6
lwd_2 <- 1.8
y_expand <- c(0.05, 0.05)

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_surv_diff) < 1e-9, "=", lab))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_fert_diff) < 1e-9, "=", lab))

v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# create plots (by row)
p5_1 <- ggplot(filter(v_plot, i < 20)) +
  geom_line(aes(i, v, linetype = type), size = 0.6) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  scale_x_continuous(limits = c(1, 20)) +
  scale_y_continuous(limits = c(0.2, 1.1), breaks = seq(0.4, 1, 0.3)) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Rep. value (', italic(v), ')'))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank(),
              legend.justification = c(0, 0),
              legend.position = c(1, 0),)

p5_2 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_surv, linetype = type), size = 0.6) +
  geom_text(data = sens_check_surv, aes(i, 4, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.000001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(P)))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p5_3 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_fert, linetype = type), size = 0.6) +
  geom_text(data = sens_check_fert, aes(i, 4, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), breaks = c(1, 20), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.000001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(F)))) +
  tt5 + theme(strip.text.x = element_blank())

# plots to grobs
g5_1 <- ggplotGrob(p5_1)
g5_2 <- ggplotGrob(p5_2)
g5_3 <- ggplotGrob(p5_3)

# arrange panels
g5 <- arrangeGrob(rbind(g4_1, g4_2, g5_1, g5_2, g5_3, size = 'last'))
g5$layout$clip[g5$layout$name == "panel"] <- "off"

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 6.5, width = 6.25, dpi = 140)
grid.arrange(g5)

# write to file
# ggsave('img/Fig_scen_a.png', g5, height = 6.5, width = 6.25, units = 'in', dpi = 300)




### compare slopes statistically
out <- mod_age_full %>% 
  filter(sim == "lambda == 1.4") %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref) %>% 
  as_tibble()


coef(lm(log(sens_surv_par) ~ i, out))[2] # more negative, as expected given negative parental age effect
coef(lm(log(sens_surv_ref) ~ i, out))[2]




### Scenario b: parental age effect on survival --------------------------------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_b.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_b.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')
Tr_sim5 <- filter(Tr_sim, sim == 'Sim 5')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j
q_i_sim5 <- filter(map_ij_sim, sim == 'Sim 5')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = TRUE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = TRUE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = TRUE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = TRUE)
A_par_sim5 <- build_Apar(Tr = Tr_sim5, q_i = q_i_sim5, stasis = TRUE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = TRUE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = TRUE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = TRUE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = TRUE)
mod_out_sim5 <- analyze_mod(A_par_sim5, stasis = TRUE)



## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4'),
  mod_out_sim5$df_age %>% mutate(sim = 'Sim 5')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4'),
  mod_out_sim5$df_par %>% mutate(sim = 'Sim 5')
)


## age-specific sensitivities, tidy
df_sens_age_sim <- mod_age_full %>% 
  dplyr::select(sim, i, ref = sens_surv_ref, par = sens_surv_par) %>% 
  gather(type, sens_surv, -i, -sim) %>% 
  mutate(sens_fert = c(mod_age_full$sens_fert_ref,
                       mod_age_full$sens_fert_par))

mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))


## check difference in sensitivities between par and ref, at early and late ages
sens_check <- mod_age_full %>% 
  dplyr::select(i, sim, starts_with('sens')) %>% 
  group_by(sim) %>% 
  filter(i %in% c(min(i):(min(i)+1), (max(i)-1):max(i))) %>% 
  ungroup() %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref)




### compare slopes statistically
out <- mod_age_full %>% 
  filter(sim == "Sim 5") %>% 
  mutate(sens_surv_diff = sens_surv_par - sens_surv_ref,
         sens_fert_diff = sens_fert_par - sens_fert_ref) %>% 
  as_tibble()


out$sens_surv_par %>% sum()
out$sens_surv_ref %>% sum()

coef(lm(log(sens_surv_par) ~ i, out))[2] # less negative, as expected given negative parental age effect
coef(lm(log(sens_surv_ref) ~ i, out))[2]

cor(log(out$sens_surv_par), out$i)^2
cor(log(out$sens_surv_ref), out$i)^2



# inputs and age-specific outputs for simulated life cycles
sim_lab <- c(
  "Sim 1" = "Strongly negative",
  "Sim 2" = "Weakly negative",
  "Sim 3" = "None",
  "Sim 4" = "Weakly positive",
  "Sim 5" = "Strongly positive"
)

lwd_1 <- 0.6
lwd_2 <- 1.8

p4_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1, labeller = labeller(sim = sim_lab)) +
  ylab(expression(paste("Survival (", italic('P'), ")"))) +
  tt5 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0),
        strip.text.x = element_text(size = 9.8))

p4_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_continuous(limits = c(0, 0.85), breaks = seq(0, 0.6, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  labs(x = NULL, y = expression(paste("Fecundity (", italic("F"), ")"))) +
  tt5 +
  theme(strip.text.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank())

g4_1 <- ggplotGrob(p4_1)
g4_2 <- ggplotGrob(p4_2)


# plot settings
lwd_1 <- 0.6
lwd_2 <- 1.8
y_expand <- c(0.05, 0.05)

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_surv_diff) < 1e-9, "=", lab))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+')) %>% 
  mutate(lab = ifelse(abs(sens_fert_diff) < 1e-9, "=", lab))

v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v)

# create plots (by row)
p5_1 <- ggplot(filter(v_plot, i < 20)) +
  geom_line(aes(i, v, linetype = type), size = 0.6) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  scale_x_continuous(limits = c(1, 20)) +
  scale_y_continuous(breaks = c(1, 1.5)) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Rep. value (', italic(v), ')'))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank(),
              legend.justification = c(0, 0),
              legend.position = c(1, 0),)

p5_2 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_surv, linetype = type), size = 0.6) +
  geom_text(data = sens_check_surv, aes(i, 4, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.00001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(P)))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p5_3 <- ggplot(filter(df_sens_age_sim, i < 20)) +
  geom_line(aes(i, sens_fert, linetype = type), size = 0.6) +
  geom_text(data = sens_check_fert, aes(i, 4, label = lab)) +
  scale_x_continuous(expand = c(0.08, 0), breaks = c(1, 20), limits = c(1, 20)) +
  scale_y_log10(expand = c(0.08, 0), breaks = c(0.00001, 0.1), labels = LabelFn) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(d, italic(lambda), ' / ', d, italic(F)))) +
  tt5 + theme(strip.text.x = element_blank())

# plots to grobs
g5_1 <- ggplotGrob(p5_1)
g5_2 <- ggplotGrob(p5_2)
g5_3 <- ggplotGrob(p5_3)

# arrange panels
g5 <- arrangeGrob(rbind(g4_1, g4_2, g5_1, g5_2, g5_3, size = 'last'))
g5$layout$clip[g5$layout$name == "panel"] <- "off"

# # print (change 'quartz' to 'window' if using Windows)
dev.off()
quartz(height = 6.5, width = 6.25, dpi = 140)
grid.arrange(g5)

# write to file
# ggsave('img/Fig_scen_b.png', g5, height = 6.5, width = 6.25, units = 'in', dpi = 300)




