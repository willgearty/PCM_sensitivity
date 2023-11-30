theme_will <- function(...) {
  theme(axis.ticks = element_line(color = "black", linewidth = 1),
        axis.line = element_blank(),
        plot.margin = unit(c(1,1,1,1), "lines"),
        axis.text = element_text(colour = "black"),
        panel.border = element_rect(fill = NA, colour = "black"),
        strip.background = element_rect(fill = "grey85", colour = "black"),
        ...)
}
