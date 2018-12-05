

## Must have the following packages installed
# install.packages(c("Matrix", "dplyr", "tidyr", "tibble", "mgcv", "popbio",
#                    "popdemo", "ggplot2", "gridExtra"))


## Working directory should be Parental_Age_Scripts


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

Tr_sim <- read.csv('dat/transition_rates_sim_main.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_main.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')


## q_i by sim
q_i_lemna <- map_ij_lemna$j
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j


## construct A_par via vec-permutation method
A_par_lemna <- build_Apar(Tr = Tr_lemna, q_i = q_i_lemna)
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4)


## analysis of A_par and A_ref
mod_out_lemna <- analyze_mod(A_par_lemna)
mod_out_sim1 <- analyze_mod(A_par_sim1)
mod_out_sim2 <- analyze_mod(A_par_sim2)
mod_out_sim3 <- analyze_mod(A_par_sim3)
mod_out_sim4 <- analyze_mod(A_par_sim4)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4')
)


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
la <- 0.9
las <- 0.9

df_arrow <- data.frame(x1 = 1.05, x2 = 11.1, y1 = 0.2, y2 = 0.2)

# panels
f2_1 <- ggplot(Tr_lemna_plot) +
  geom_line(aes(i, P_ij, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, surv_ref), linetype = 2, size = 0.7) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  annotate('text', label = 'Parental age', x = 1, y = 0.265, hjust = 0, size = 3) +
  annotate('text', label = '(a)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, 1)) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.7, 'cm'), barheight = unit(0.5, 'cm'))) +
  xlab(NULL) + ylab(NULL) + ggtitle(expression(paste('Survival (', italic(P), ')'))) +
  tt2 + theme(axis.text.x = element_blank(), legend.position = c(0.04, 0.04), legend.justification = c(0, 0), legend.margin = margin(0, 0, 0, 0))

f2_2 <- ggplot(Tr_lemna_plot) +
  geom_line(aes(i, F_ij, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, fert_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(b)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  coord_cartesian(xlim = c(0, 30), ylim = c(0, max(c(1, max(Tr_lemna_plot$F_ij))))) +
  scale_x_continuous(breaks = seq(0, 30, 10)) +
  scale_y_continuous(breaks = seq(0, 1, 0.2)) +
  xlab(' ') + ylab(NULL) + ggtitle(expression(paste('Fecundity (', italic(F), ')'))) +
  tt2

f2_3 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, w, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, w_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(c)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(NULL) + ylab(NULL) + ggtitle(expression(paste('Stable distribution (', italic(w), ')'))) +
  tt2 + theme(axis.text.x = element_blank())

f2_4 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, v/v[1], col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, v_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(d)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10() +
  xlab('Age (days)') + ylab(NULL) + ggtitle(expression(paste('Reproductive value (', italic(v), ')'))) +
  tt2

f2_5 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, sens_surv, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, sens_surv_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(e)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(NULL) + ylab(NULL) + ggtitle(expression(paste('Sensitivity of λ to ', italic(P)))) +
  tt2 + theme(axis.text.x = element_blank())

f2_6 <- ggplot(mod_out_lemna$df_par) +
  geom_line(aes(i, sens_fert, col = j, group = j), size = ls, alpha = la) +
  geom_line(data = mod_out_lemna$df_age, aes(i, sens_fert_ref), linetype = 2, size = 0.7) +
  annotate('text', label = '(f)', x = 29.5, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_color_gradient(low = col_low, high = col_upp, guide = F) +
  scale_y_log10(breaks = c(0.1, 0.00000001), labels = LabelFn) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  xlab(' ') + ylab(NULL) + ggtitle(expression(paste('Sensitivity of λ to ', italic(F)))) +
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
# dev.off()
# quartz(height = 4.35, width = 7.25, dpi = 150)
# grid.arrange(g2)

# write to file
ggsave('img/Fig_2.png', g2, height = 4.35, width = 7.25, units = 'in', dpi = 300)




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
  xlab(NULL) + ylab(NULL) +
  ggtitle(expression(paste('Stable age distribution (', italic(w[i]), ')'))) +
  tt3 + theme(axis.text.x = element_blank())

f3_2 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type), show.legend = FALSE) +
  annotate('text', label = '(b)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(labels = LabelFn) +
  xlab('Age (days)') + ylab(NULL) +
  ggtitle(expression(paste('Reproductive value (', italic(v[i]), ')'))) +
  tt3 + theme(plot.margin = unit(c(0.21, 0.15, 0.15, 0.15), 'cm'))

f3_3 <- ggplot(df_sens_age_lemna, aes(i, sens_surv, linetype = type)) +
  geom_line(show.legend = FALSE) +
  annotate('text', label = '(c)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(breaks = c(0.1, 0.0001, 0.0000001), labels = LabelFn) +
  coord_cartesian(ylim = c(0.00000001, 0.25)) +
  xlab(NULL) + ylab(NULL) +
  ggtitle(expression(paste('Sensitivity of ', italic(lambda), ' to survival (', italic(P[i]), ')'))) +
  tt3 + theme(axis.text.x = element_blank())

f3_4 <- ggplot(df_sens_age_lemna, aes(i, sens_fert, linetype = type)) +
  geom_line(show.legend = FALSE) +
  annotate('text', label = '(d)', x = 29.7, y = Inf, hjust = 0.5, vjust = 1.9, size = 3.5) +
  scale_x_continuous(limits = c(0, 30), breaks = seq(0, 30, 10)) +
  scale_y_log10(breaks = c(0.1, 0.0001, 0.0000001), labels = LabelFn) +
  coord_cartesian(ylim = c(0.00000001, 0.25)) +
  xlab('Age (days)') + ylab(NULL) +
  ggtitle(expression(paste('Sensitivity of ', italic(lambda), ' to fecundity (', italic(F[i]), ')'))) +
  tt3 + theme(plot.margin = unit(c(0.21, 0.15, 0.15, 0.15), 'cm'))

# panels to grobs
g3_1 <- ggplotGrob(f3_1)
g3_2 <- ggplotGrob(f3_2)
g3_3 <- ggplotGrob(f3_3)
g3_4 <- ggplotGrob(f3_4)

# arrange panels
g3 <- arrangeGrob(g3_1, g3_3, g3_2, g3_4, nrow = 2, heights = c(1, 1.155))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 4.5, width = 6, dpi = 150)
# grid.arrange(g3)

# write to file
ggsave('img/Fig_3.png', g3, height = 4.5, width = 6, units = 'in', dpi = 300)




##### Figure 4 -----------------------------------------------------------------

# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# plot settings
lwd_1 <- 0.6
lwd_2 <- 1.5
y_expand <- c(0.05, 0.05)

df_arrow <- data.frame(sim = "Sim 1",
                       x1 = 1.35, x2 = 11.5, y1 = 0.22, y2 = 0.22)

df_lab <- data.frame(sim = "Sim 1",
                     label = 'Parental age',
                     stringsAsFactors = FALSE)

# create plots (by row)
p4_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('P'))) +
  tt4 +
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
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt4 +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

p4_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt4 + theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p4_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.1, 1, 10),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt4 + theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p4_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt4 + theme(axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p4_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt4 +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
g4_1 <- ggplotGrob(p4_1)
g4_2 <- ggplotGrob(p4_2)
g4_3 <- ggplotGrob(p4_3)
g4_4 <- ggplotGrob(p4_4)
g4_5 <- ggplotGrob(p4_5)
g4_6 <- ggplotGrob(p4_6)

# arrange panels
g4 <- arrangeGrob(rbind(g4_1, g4_2, g4_3, g4_4, g4_5, g4_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(g4)

# write to file
ggsave('img/Fig_4.png', g4, height = 7, width = 5.6, units = 'in', dpi = 300)




##### Figure 5 -----------------------------------------------------------------
# age-specific sensitivities from A_par and A_ref, simulated life cycles

# plot settings
lwd_1 <- 0.6
lwd_2 <- 1.8
y_expand <- c(0.05, 0.05)

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  mutate(lab = ifelse(v_par < v_ref, '-', '+'))

# create plots (by row)
p5_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  geom_text(data = v_check, aes(i, Inf, label = lab), vjust = 1.1) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank())

p5_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 1.3, label = lab)) +
  scale_y_log10(breaks = c(0.001, 0.1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt5 + theme(axis.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.ticks.x = element_blank(),
              strip.text.x = element_blank())

p5_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 1.3, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.001, 0.1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt5 + theme(strip.text.x = element_blank())

# plots to grobs
g5_1 <- ggplotGrob(p5_1)
g5_2 <- ggplotGrob(p5_2)
g5_3 <- ggplotGrob(p5_3)

# arrange panels
g5 <- arrangeGrob(rbind(g5_1, g5_2, g5_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(g5)

# write to file
ggsave('img/Fig_5.png', g5, height = 5, width = 7.25, units = 'in', dpi = 300)








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
s <- 9

## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'lambda == 0.5')
Tr_sim2 <- filter(Tr_sim, sim == 'lambda == 0.8')
Tr_sim3 <- filter(Tr_sim, sim == 'lambda == 1.1')
Tr_sim4 <- filter(Tr_sim, sim == 'lambda == 1.4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'lambda == 0.5')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'lambda == 0.8')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'lambda == 1.1')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'lambda == 1.4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = FALSE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = FALSE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = FALSE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = FALSE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = FALSE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = FALSE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = FALSE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = FALSE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'lambda == 0.5'),
  mod_out_sim2$df_age %>% mutate(sim = 'lambda == 0.8'),
  mod_out_sim3$df_age %>% mutate(sim = 'lambda == 1.1'),
  mod_out_sim4$df_age %>% mutate(sim = 'lambda == 1.4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'lambda == 0.5'),
  mod_out_sim2$df_par %>% mutate(sim = 'lambda == 0.8'),
  mod_out_sim3$df_par %>% mutate(sim = 'lambda == 1.1'),
  mod_out_sim4$df_par %>% mutate(sim = 'lambda == 1.4')
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




### Figure A5

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

# age-specific reproductive values
v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# create plots (by row)
pa5_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(, labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1, labeller = label_parsed) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())

pa5_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 9, label = lab)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa5_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 15, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  coord_cartesian(ylim = c(0.000005, 15)) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt_a_age + theme(strip.text.x = element_blank())

# plots to grobs
ga5_1 <- ggplotGrob(pa5_1)
ga5_2 <- ggplotGrob(pa5_2)
ga5_3 <- ggplotGrob(pa5_3)

# arrange panels
ga5 <- arrangeGrob(rbind(ga5_1, ga5_2, ga5_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(ga5)

# write to file
ggsave('img/Fig_A5.png', ga5, height = 5, width = 7.25, units = 'in', dpi = 300)




### Figure A6
df_arrow_scen_a <- data.frame(sim = "lambda == 0.5",
                              x1 = 1.35, x2 = 11.5, y1 = 0.22, y2 = 0.22,
                              stringsAsFactors = FALSE)

df_lab_scen_a <- data.frame(sim = "lambda == 0.5",
                            label = 'Parental age',
                            stringsAsFactors = FALSE)

# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# create plots (by row)
pa6_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow_scen_a, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab_scen_a, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, nrow = 1, labeller = label_parsed) +
  ylab(expression(italic('P'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

pa6_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, breaks = seq(0, 0.9, 0.3)) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa6_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa6_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.01, 1),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa6_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa6_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt_a_out +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
ga6_1 <- ggplotGrob(pa6_1)
ga6_2 <- ggplotGrob(pa6_2)
ga6_3 <- ggplotGrob(pa6_3)
ga6_4 <- ggplotGrob(pa6_4)
ga6_5 <- ggplotGrob(pa6_5)
ga6_6 <- ggplotGrob(pa6_6)

# arrange panels
ga6 <- arrangeGrob(rbind(ga6_1, ga6_2, ga6_3, ga6_4, ga6_5, ga6_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(ga6)

# write to file
ggsave('img/Fig_A6.png', ga6, height = 7, width = 5.6, units = 'in', dpi = 300)




### Scenario b: s = 4 parental age classes -------------------------------------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_b.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_b.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = FALSE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = FALSE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = FALSE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = FALSE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = FALSE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = FALSE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = FALSE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = FALSE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4')
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




### Figure A7

# age-specific reproductive values
v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  mutate(lab = ifelse(v_par < v_ref, '-', '+'))

# create plots (by row)
pa7_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  geom_text(data = v_check, aes(i, Inf, label = lab), vjust = 1.1) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(, labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())

pa7_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 9, label = lab)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa7_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 15, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  coord_cartesian(ylim = c(0.000005, 15)) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt_a_age + theme(strip.text.x = element_blank())

# plots to grobs
ga7_1 <- ggplotGrob(pa7_1)
ga7_2 <- ggplotGrob(pa7_2)
ga7_3 <- ggplotGrob(pa7_3)

# arrange panels
ga7 <- arrangeGrob(rbind(ga7_1, ga7_2, ga7_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(ga7)

# write to file
ggsave('img/Fig_A7.png', ga7, height = 5, width = 7.25, units = 'in', dpi = 300)




### Figure A8

# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# create plots (by row)
pa8_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('P'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

pa8_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa8_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa8_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.1, 1, 10),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa8_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa8_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt_a_out +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
ga8_1 <- ggplotGrob(pa8_1)
ga8_2 <- ggplotGrob(pa8_2)
ga8_3 <- ggplotGrob(pa8_3)
ga8_4 <- ggplotGrob(pa8_4)
ga8_5 <- ggplotGrob(pa8_5)
ga8_6 <- ggplotGrob(pa8_6)

# arrange panels
ga8 <- arrangeGrob(rbind(ga8_1, ga8_2, ga8_3, ga8_4, ga8_5, ga8_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(ga8)

# write to file
ggsave('img/Fig_A8.png', ga8, height = 7, width = 5.6, units = 'in', dpi = 300)




### Scenario c: positive parental age effect -----------------------------------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_c.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_c.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = FALSE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = FALSE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = FALSE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = FALSE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = FALSE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = FALSE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = FALSE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = FALSE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4')
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



### Fig A9

# age-specific reproductive values
v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  mutate(lab = ifelse(v_par < v_ref, '-', '+'))

# create plots (by row)
pa9_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  geom_text(data = v_check, aes(i, Inf, label = lab), vjust = 1.1) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(breaks = c(0.001, 0.1, 1, 10, 100), labels = LabelFn) +
  coord_cartesian(ylim = c(0.005, 580)) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())

pa9_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 3, label = lab)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa9_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 5, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  coord_cartesian(ylim = c(0.000001, 5)) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt_a_age + theme(strip.text.x = element_blank())

# plots to grobs
ga9_1 <- ggplotGrob(pa9_1)
ga9_2 <- ggplotGrob(pa9_2)
ga9_3 <- ggplotGrob(pa9_3)

# arrange panels
ga9 <- arrangeGrob(rbind(ga9_1, ga9_2, ga9_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(ga9)

# write to file
ggsave('img/Fig_A9.png', ga9, height = 5, width = 7.25, units = 'in', dpi = 300)



### Fig A10 

# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# create plots (by row)
pa10_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('P'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

pa10_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa10_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa10_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.1, 1, 10),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa10_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa10_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt_a_out +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
ga10_1 <- ggplotGrob(pa10_1)
ga10_2 <- ggplotGrob(pa10_2)
ga10_3 <- ggplotGrob(pa10_3)
ga10_4 <- ggplotGrob(pa10_4)
ga10_5 <- ggplotGrob(pa10_5)
ga10_6 <- ggplotGrob(pa10_6)

# arrange panels
ga10 <- arrangeGrob(rbind(ga10_1, ga10_2, ga10_3, ga10_4, ga10_5, ga10_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(ga10)

# write to file
ggsave('img/Fig_A10.png', ga10, height = 7, width = 5.6, units = 'in', dpi = 300)




### Scenario d: parental age effect on fecundity rather than survival ----------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_d.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_d.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = FALSE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = FALSE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = FALSE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = FALSE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = FALSE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = FALSE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = FALSE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = FALSE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4')
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



### Fig A11

# age-specific reproductive values
v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i)-1)) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  mutate(lab = ifelse(v_par < v_ref, '-', '+'))

# create plots (by row)
pa11_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  geom_text(data = v_check, aes(i, Inf, label = lab), vjust = 1.1) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(breaks = c(0.001, 0.1, 1, 10, 100), labels = LabelFn) +
  coord_cartesian(ylim = c(0.005, 580)) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())

pa11_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 3, label = lab)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa11_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 5, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.001, 1), labels = LabelFn) +
  coord_cartesian(ylim = c(0.000001, 5)) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt_a_age + theme(strip.text.x = element_blank())

# plots to grobs
ga11_1 <- ggplotGrob(pa11_1)
ga11_2 <- ggplotGrob(pa11_2)
ga11_3 <- ggplotGrob(pa11_3)

# arrange panels
ga11 <- arrangeGrob(rbind(ga11_1, ga11_2, ga11_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(ga11)

# write to file
ggsave('img/Fig_A11.png', ga11, height = 5, width = 7.25, units = 'in', dpi = 300)



### Fig A12

## organize data frames for plotting
# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  mutate(P_ij = ifelse(i == max(i), NA_real_, P_ij)) %>%
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# create plots (by row)
pa12_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('P'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

pa12_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa12_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa12_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.1, 1, 10),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa12_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt_a_out + theme(axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa12_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt_a_out +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
ga12_1 <- ggplotGrob(pa12_1)
ga12_2 <- ggplotGrob(pa12_2)
ga12_3 <- ggplotGrob(pa12_3)
ga12_4 <- ggplotGrob(pa12_4)
ga12_5 <- ggplotGrob(pa12_5)
ga12_6 <- ggplotGrob(pa12_6)

# arrange panels
ga12 <- arrangeGrob(rbind(ga12_1, ga12_2, ga12_3, ga12_4, ga12_5, ga12_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(ga12)

# write to file
ggsave('img/Fig_A12.png', ga12, height = 7, width = 5.6, units = 'in', dpi = 300)




### Scenario e: stasis loop at the maximum age class ---------------------------

## read data
Tr_sim <- read.csv('dat/transition_rates_sim_scen_e.csv', stringsAsFactors = FALSE)
map_ij_sim <- read.csv('dat/map_ij_sim_scen_e.csv', stringsAsFactors = FALSE)


## transition rates by sim
Tr_sim1 <- filter(Tr_sim, sim == 'Sim 1')
Tr_sim2 <- filter(Tr_sim, sim == 'Sim 2')
Tr_sim3 <- filter(Tr_sim, sim == 'Sim 3')
Tr_sim4 <- filter(Tr_sim, sim == 'Sim 4')


## q_i by sim
q_i_sim1 <- filter(map_ij_sim, sim == 'Sim 1')$j
q_i_sim2 <- filter(map_ij_sim, sim == 'Sim 2')$j
q_i_sim3 <- filter(map_ij_sim, sim == 'Sim 3')$j
q_i_sim4 <- filter(map_ij_sim, sim == 'Sim 4')$j


## construct A_par via vec-permutation method
A_par_sim1 <- build_Apar(Tr = Tr_sim1, q_i = q_i_sim1, stasis = TRUE)
A_par_sim2 <- build_Apar(Tr = Tr_sim2, q_i = q_i_sim2, stasis = TRUE)
A_par_sim3 <- build_Apar(Tr = Tr_sim3, q_i = q_i_sim3, stasis = TRUE)
A_par_sim4 <- build_Apar(Tr = Tr_sim4, q_i = q_i_sim4, stasis = TRUE)


## analysis of A_par and A_ref
mod_out_sim1 <- analyze_mod(A_par_sim1, stasis = TRUE)
mod_out_sim2 <- analyze_mod(A_par_sim2, stasis = TRUE)
mod_out_sim3 <- analyze_mod(A_par_sim3, stasis = TRUE)
mod_out_sim4 <- analyze_mod(A_par_sim4, stasis = TRUE)


## arrange model outputs in single df
# age-specific values
mod_age_full <- rbind.data.frame(
  mod_out_sim1$df_age %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_age %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_age %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_age %>% mutate(sim = 'Sim 4')
)

# age-by-parental-age values
mod_par_full <- rbind.data.frame(
  mod_out_sim1$df_par %>% mutate(sim = 'Sim 1'),
  mod_out_sim2$df_par %>% mutate(sim = 'Sim 2'),
  mod_out_sim3$df_par %>% mutate(sim = 'Sim 3'),
  mod_out_sim4$df_par %>% mutate(sim = 'Sim 4')
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




### Fig A13

# age-specific reproductive values
v_plot <- mod_age_full %>% 
  gather(type, v, v_ref:v_par) %>% 
  group_by(sim, type) %>% 
  ungroup()

# labels to indicate whether sensitivities greater for A_par or A_ref
sens_check_surv <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_surv_diff < 0, '-', '+'))

sens_check_fert <- sens_check %>% 
  group_by(sim) %>% 
  filter(i %in% c(1, max(i))) %>% 
  ungroup() %>% 
  mutate(lab = ifelse(sens_fert_diff < 0, '-', '+'))

v_check <- v_plot %>% 
  group_by(sim) %>% 
  filter(i == max(i)) %>% 
  ungroup() %>% 
  spread(type, v) %>% 
  mutate(lab = ifelse(v_par < v_ref, '-', '+'))


# create plots (by row)
pa13_1 <- ggplot(v_plot) +
  geom_line(aes(i, v, linetype = type)) +
  geom_text(data = v_check, aes(i, Inf, label = lab), vjust = 1.1) +
  scale_linetype(labels = c(expression(paste('Parental age model (', bolditalic(A)^par, ')')),
                            expression(paste('Reference model (', bolditalic(A)^ref, ')')))) +
  scale_y_log10(breaks = c(0.001, 0.1, 1, 10, 100), labels = LabelFn) +
  coord_cartesian(ylim = c(0.005, 1000)) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Reprod. value (', italic(v[i]), ')'))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank())

pa13_2 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_surv, linetype = type)) +
  geom_text(data = sens_check_surv, aes(i, 2, label = lab)) +
  scale_y_log10(breaks = c(0.0001, 0.01, 1), labels = LabelFn) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(P[i])))) +
  tt_a_age + theme(axis.text.x = element_blank(),
                   axis.title.x = element_blank(),
                   axis.ticks.x = element_blank(),
                   strip.text.x = element_blank())

pa13_3 <- ggplot(df_sens_age_sim) +
  geom_line(aes(i, sens_fert, linetype = type)) +
  geom_text(data = sens_check_fert, aes(i, 1, label = lab)) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(breaks = c(0.0001, 0.01, 1), labels = LabelFn) +
  coord_cartesian(ylim = c(0.00001, 1)) +
  scale_linetype(guide = F) +
  facet_wrap(~ sim, nrow = 1) +
  xlab('Age') +
  ylab(expression(paste('Sensitivity ', italic(lambda), ' to ', italic(F[i])))) +
  tt_a_age + theme(strip.text.x = element_blank())

# plots to grobs
ga13_1 <- ggplotGrob(pa13_1)
ga13_2 <- ggplotGrob(pa13_2)
ga13_3 <- ggplotGrob(pa13_3)

# arrange panels
ga13 <- arrangeGrob(rbind(ga13_1, ga13_2, ga13_3, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 5, width = 7.25, dpi = 150)
# grid.arrange(ga13)

# write to file
ggsave('img/Fig_A13.png', ga13, height = 5, width = 7.25, units = 'in', dpi = 300)




### Fig A14

# age-by-parental age model inputs and outputs, flat form
mod_par_full_flat <- mod_par_full %>% 
  full_join(Tr_sim, by = c('sim', 'i', 'j')) %>% 
  group_by(sim) %>% 
  ungroup() %>% 
  gather(type, val, -i, -j, -sim) %>% 
  mutate(type = factor(type, levels = c('P_ij', 'F_ij', 'w', 'v',
                                        'sens_surv', 'sens_fert')))

# create plots (by row)
pa14_1 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'P_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, surv_ref), size = lwd_1, linetype = 2) +
  geom_segment(data = df_arrow, aes(x = x1, y = y1, xend = x2, yend = y2), arrow = arrow(length = unit(0.1, 'cm'), ends = 'last', type = 'closed')) +
  geom_text(data = df_lab, aes(1.35, 0.31, label = label), hjust = 0, size = 2.6) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_continuous(expand = y_expand, limits = c(0, 1), breaks = seq(0, 0.9, 0.3)) +
  guides(color = guide_colorbar(label = F, ticks = F, direction = 'horizontal', barwidth = unit(1.5, 'cm'), barheight = unit(0.4, 'cm'))) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('P'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0))

pa14_2 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'F_ij'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('F'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa14_3 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'w'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, w_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(italic('w'))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa14_4 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'v'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, v_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.01, 1, 100),  labels = LabelFn) +
  facet_wrap(~ sim, nrow = 1) +
  ylab(expression(italic('v'))) +
  tt_a_out + 
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa14_5 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_surv'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_surv_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(P)))) +
  tt_a_out +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text.x = element_blank())

pa14_6 <- ggplot() +
  geom_line(data = filter(mod_par_full_flat, type == 'sens_fert'), aes(i, val, col = j, group = j), size = lwd_2) +
  geom_line(data = mod_age_full, aes(i, sens_fert_ref), size = lwd_1, linetype = 2) +
  scale_color_gradient(low = col_low, high = col_upp, name = NULL, guide = F) +
  scale_x_continuous(breaks = c(1, 20)) +
  scale_y_log10(expand = y_expand, breaks = c(0.000001, 0.1), labels = LabelFn) +
  facet_wrap(~ sim, scales = 'free_x', nrow = 1) +
  xlab('Age') +
  ylab(expression(paste(italic(d), italic(lambda), ' / ', italic(d), italic(F)))) +
  tt_a_out +
  theme(strip.text.x = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 5.5, 5.5), units = 'pt'))

# plots to grobs
ga14_1 <- ggplotGrob(pa14_1)
ga14_2 <- ggplotGrob(pa14_2)
ga14_3 <- ggplotGrob(pa14_3)
ga14_4 <- ggplotGrob(pa14_4)
ga14_5 <- ggplotGrob(pa14_5)
ga14_6 <- ggplotGrob(pa14_6)

# arrange panels
ga14 <- arrangeGrob(rbind(ga14_1, ga14_2, ga14_3, ga14_4, ga14_5, ga14_6, size = 'last'))

# # print (change 'quartz' to 'window' if using Windows)
# dev.off()
# quartz(height = 7, width = 5.6, dpi = 150)
# grid.arrange(ga14)

# write to file
ggsave('img/Fig_A14.png', ga14, height = 7, width = 5.6, units = 'in', dpi = 300)



