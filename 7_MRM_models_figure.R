# Cargar librerías
library(ggplot2)
library(ggpattern)
library(dplyr)

df <- read.csv("regression_coefficient_bray.csv")

# Asegurar orden deseado
df$variable <- factor(df$variable, levels = c("Clim", "Geog", "Envm"))
df$panel <- factor(df$panel, levels = c(
  "Forest", "Stream", "Pond",
  "P. cinereus", "E. bislineata", "N. viridescens"
))

# Gráfico
(g <- ggplot(df, aes(x = variable, y = value, 
               fill = fill_color, 
               pattern = group)) +
  geom_col_pattern(
    color = "black",
    position = position_dodge(width = 0.6),
    width = 0.5,
    pattern_density = 0.5,
    pattern_fill = "white",
    pattern_spacing = 0.02,
    pattern_angle = 45,
    pattern_size = 0.1,                      # <- this makes stripes thinner!
    linewidth = 0.2,
    pattern_key_scale_factor=.5
  ) +
  geom_text(aes(label = significance), 
            position = position_dodge(width = 0.6),
            vjust = 0.1, size = 4) +
  facet_wrap(~panel, nrow = 2) +
  scale_fill_identity() +
  scale_pattern_manual(values = c("Bd-" = "none", "Bd+" = "stripe")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Regression coefficient") +
  theme_bw(base_size = 10) +
    theme(
    strip.text = element_text(face = "bold"),
    panel.spacing = unit(1, "lines"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "left"
  )
)

df2 <- read.csv("regression_coefficient_jaccard_legend.csv")

# Asegurar orden deseado
df2$variable <- factor(df2$variable, levels = c("Clim", "Geog", "Envm"))
df2$panel <- factor(df2$panel, levels = c(
  "Forest", "Stream", "Pond",
  "P. cinereus", "E. bislineata", "N. viridescens"
))

# Gráfico
(g2 <- ggplot(df2, aes(x = variable, y = value, 
                     fill = fill_color, 
                     pattern = group)) +
                  pattern = group)) +
  geom_col_pattern(
    color = "black",
    position = position_dodge(width = 0.6),
    width = 0.5,
    pattern_density = 0.5,
    pattern_fill = "white",
    pattern_spacing = 0.02,
    pattern_angle = 45,
    pattern_size = 0.1,                      # <- this makes stripes thinner!
    linewidth = 0.2,
    pattern_key_scale_factor=.5
  ) +
    geom_text(aes(label = significance), 
              position = position_dodge(width = 0.6),
              vjust = 0.1, size = 4) +
    facet_wrap(~panel, nrow = 2) +
    scale_fill_identity() +
    scale_pattern_manual(values = c("Bd-" = "none", "Bd+" = "stripe")) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = NULL, y = "Regression coefficient") +
    theme_bw(base_size = 10) +
    theme(
      strip.text = element_text(face = "bold"),
      panel.spacing = unit(1, "lines"),
      axis.text.x = element_text(face = "bold"),
      legend.position = "none"
    )
)


