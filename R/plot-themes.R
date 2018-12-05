

library(ggplot2)
library(grid)


## global plot settings
col_low <- '#006d2c'
col_upp <- '#e5f5e0'


## ggplot themes
tt2 <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text = element_text(size = 8.2),
        axis.title = element_text(size = 11),
        plot.title = element_text(size = 10, vjust = -0.8),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, unit = 'cm')),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), 'cm'))


tt3 <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        plot.title = element_text(size = 10, vjust = -1),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.text = element_text(size = 8.5),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(margin = margin(0.1, 0, 0, 0, unit = 'cm')),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.01, 0.01),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.4),
        legend.text.align = 0,
        legend.key.width = unit(0.9, 'lines'),
        legend.key.height = unit(0.8, 'lines'),
        plot.margin = unit(c(0.15, 0.15, 0.15, 0.15), 'cm'))


tt4 <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 1.5, 5.5), units = 'pt'),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA))


tt5 <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, unit = 'cm')),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        plot.margin = unit(c(1.5, 5.5, 1.5, 5.5), units = 'pt'),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.76, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.4),
        legend.text.align = 0,
        legend.key.width = unit(0.9, 'lines'),
        legend.key.height = unit(0.8, 'lines'))


tt_a234 <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 1.5, 5.5), units = 'pt'),
        legend.position = c(0, 0),
        legend.justification = c(0, 0),
        legend.margin = margin(0, 0, 0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA))
        

tt_a_age <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        axis.title.x = element_text(margin = margin(0.2, 0, 0, 0, unit = 'cm')),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        plot.margin = unit(c(1.5, 5.5, 1.5, 5.5), units = 'pt'),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA),
        legend.position = c(0.76, 0),
        legend.title = element_blank(),
        legend.text = element_text(size = 7.4),
        legend.text.align = 0,
        legend.key.width = unit(0.9, 'lines'),
        legend.key.height = unit(0.8, 'lines'))


tt_a_out <- theme_bw() +
  theme(panel.grid = element_blank(),
        line = element_line(size = 0.25),
        rect = element_rect(size = 0.4),
        axis.ticks = element_line(color = 'grey50'),
        axis.text.y = element_text(angle = 90, hjust = 0.5),
        strip.background = element_blank(),
        strip.text.x = element_text(size = 10),
        strip.text.y = element_blank(),
        plot.margin = unit(c(1.5, 5.5, 1.5, 5.5), units = 'pt'),
        legend.justification = c(0, 0),
        legend.background = element_blank(),
        legend.key = element_rect(fill = NA, color = NA))

