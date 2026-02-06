library(ggplot2)

# Styling theme
consistentTheme <- theme(
  axis.line = element_line(color = "black"),
  axis.title = element_text(size = 16),
  axis.text = element_text(size = 14),
  strip.text.x = element_text(size = 15),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15),
  axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
)

# Styling trajectory analysis
gg_theme_specs <- function(fs=20, hpos=0.5, vpos=0.5, angle=0) {
  theme_bw(base_size=fs) + theme(
  axis.line = element_line(color = "black"),
  axis.title = element_text(color="black"),
  axis.text = element_text(color="black"),
  strip.text.x = element_text(color="black"),
  legend.text = element_text(size=fs-3, color="black"),
  legend.title = element_text(color="black"),
  # panel.grid.major = element_line(color="black"),
  # panel.grid.minor = element_line(color="black"),
  axis.text.x = element_text(angle=angle, vjust=vpos, hjust = hpos, color="black"))
}

IMC_colors <- c(
  `CD4+ Naive`     = "#FDBD13",
  `CD4+ Activated Naive` = "#1857F4",
  `CD4+ CXCR5+ pTfh`   = "#C0EAFE",
  `CD4+ Cytotoxic`  = '#f5a993',
  `CD4+ VD2 γδ`      = "#00C3A6",
  `CD4+ Th17` = '#b7e8e0',
  `CD4+ Type 1\nIFN responsive` = '#e04a56',
  `CD4+ AP-1+ Naive`  = '#f7e8c1',
  `CD4+ Th2`       = '#a3a2a0',
  `CD4+ VD1 γδ`  = "#7AD0FF",
  `CD4+ T-regs`  = '#e6e6e6',
  `CD4+ Transitional Naive` = '#c0b3ff',
  `CD4+ Early Activated` = '#f21d95',
  `ILC2` = '#f5ef3d',
  `CD4+ iNKT`     = "#FF6133",
  `CD4+ Mitotic`     = '#CE155C'
)

IMC_colors_numbers <- c(
  `0`     = "#FDBD13",
  `14` = "#1857F4",
  `1`   = "#C0EAFE",
  `10`  = '#f5a993',
  `7`      = "#00C3A6",
  `3` = '#b7e8e0',
  `11` = '#e04a56',
  `8`  = '#f7e8c1',
  `2`       = '#a3a2a0',
  `13`  = "#7AD0FF",
  `6`  = '#e6e6e6',
  `5` = '#c0b3ff',
  `4` = '#f21d95',
  `15` = '#f5ef3d',
  `9`     = "#FF6133",
  `12`     = '#CE155C'
) 