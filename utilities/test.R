library(dplyr)
library(ggplot2)
library(magrittr)
library(DT)

fa.sizes <- read.table(snakemake@input[["faidx"]], header = F)[,1:2]
colnames(fa.sizes) <- c("contig", "size")
fa.sizes

sv <- read.table(snakemake@input[["statsfile"]], header = T)
sv <- merge(x=sv,y=fa.sizes, by="contig")
sv$length <- gsub("\\.", "1", sv$length) %>%  as.numeric()

sv

grpstats <- sv %>%
  group_by(contig, type) %>%
  summarise(count = n(), total_bp  = sum(length))


grpstats %>%
  group_by(type) %>%
  summarise(
    total = sum(count),
    avg_per_contig = mean(count),
    sd = sd(count)
)

#DT::datatable(grpstats, rownames = F, filter = "top", extensions = 'Buttons', options = list(dom = 'Brtip', buttons = c('csv', 'pdf'), scrollX = TRUE))

.i <- 0
for (i in unique(fa.sizes$contig)) {
  .i <- .i + 1
  sv.filt <- sv %>% filter(contig == i)
  sv_stats <- group_by(sv.filt, type) %>%  summarise(n = length(type))
  if (nrow(sv_stats) > 0) {
    .subtitle <- paste(paste(sv_stats$type, sv_stats$n), collapse = " | ")
  } else {
    next
  }
  plt <- sv.filt %>%
    ggplot() +
    geom_rect(alpha = 0.7, aes(xmin = position_start, xmax = position_end, ymin = 0, ymax = 1, fill = type, color = type)) +
    theme_bw() +
    theme(
      axis.text.y=element_blank(),
      axis.line.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    scale_color_manual(values=c("DEL"="#004D40", "DUP"="#FFC107", "INV"="#1E88E5","BND"="#D81B60")) +
    scale_fill_manual(values=c("DEL"="#004D40", "DUP"="#FFC107", "INV"="#1E88E5","BND"="#D81B60")) +
    xlim(1, sv.filt$size[1]) +
    xlab("Position (bp)") +
    labs(fill = "SV type", color = "SV type") +
    ggtitle(i, subtitle = .subtitle)
  print(plt)
  if (.i == 30) {
    break
  }
}
save.image("debug.rdata")
