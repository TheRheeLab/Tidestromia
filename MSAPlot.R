# RCA Alignment Plot: Grouped by BLOSUM62 similarity

# --- Install and load required packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in c("ggplot2", "Biostrings", "pwalign")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (pkg %in% c("Biostrings", "pwalign")) {
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
  library(pkg, character.only = TRUE)
}

library(grid)  # for unit() in legend sizing

# --- User-configurable settings ---
input_fasta <- "RCA_to_Plot.fasta"
output_pdf  <- "RCA_Alignment_blosum_groups_colored.pdf"
canvas_width  <- 12  # inches
canvas_height <- 9   # inches
threshold     <- 1   # BLOSUM62 score cutoff for grouping
ignore_gaps   <- TRUE

# Color map for percent-based similarity groups
cols <- c(
  "<60%"   = "#ffffff",
  "60–80%" = "#fffbb2",
  "80–99%" = "#f9c486",
  "100%"   = "#f18171"
)

# --- Read alignment and convert to long data frame ---
alignment <- readAAStringSet(input_fasta)
df <- do.call(rbind, lapply(seq_along(alignment), function(i) {
  aa <- strsplit(as.character(alignment[i]), "")[[1]]
  data.frame(
    seq_name = names(alignment)[i],
    position = seq_along(aa),
    residue  = aa,
    stringsAsFactors = FALSE
  )
}))

# --- Hide terminal gaps ---
bounds    <- lapply(split(df, df$seq_name), function(sub) {
  range(sub$position[sub$residue != "-"])
})
start_pos <- sapply(bounds, `[`, 1)
end_pos   <- sapply(bounds, `[`, 2)
df$label <- with(df, ifelse(
  residue == "-" & (position < start_pos[seq_name] | position > end_pos[seq_name]),
  NA, residue
))

# --- Compute percent of similar residues per column ---
data("BLOSUM62", package = "pwalign")
bl62 <- BLOSUM62
df$percent <- NA_real_
for (pos in unique(df$position)) {
  idx <- which(df$position == pos)
  sub <- df[idx, ]
  if (ignore_gaps) {
    keep <- sub$residue != "-"
    idx <- idx[keep]
    sub <- sub[keep, ]
  }
  if (nrow(sub) == 0) next
  uniqs <- unique(sub$residue)
  sim_matrix <- outer(univs <- uniqs, univs,
                      Vectorize(function(a,b) bl62[a,b] >= threshold)
  )
  dimnames(sim_matrix) <- list(univs, univs)
  clusters <- list()
  rem <- univs
  while (length(rem)) {
    seed <- rem[1]; comp <- seed; queue <- seed
    while (length(queue)) {
      x <- queue[1]; queue <- queue[-1]
      nbr <- univs[sim_matrix[x,]]
      new <- setdiff(nbr, comp)
      comp <- c(comp, new); queue <- c(queue, new)
    }
    clusters <- c(clusters, list(comp))
    rem <- setdiff(rem, comp)
  }
  total <- nrow(sub)
  for (grp in clusters) {
    hits <- idx[df$residue[idx] %in% grp]
    df$percent[hits] <- length(hits) / total * 100
  }
}

# --- Assign fill colors based on percent groups ---
df$fill <- with(df, ifelse(
  is.na(percent), NA,
  ifelse(percent == 100, cols["100%"],
         ifelse(percent >= 80, cols["80–99%"],
                ifelse(percent >= 60, cols["60–80%"], cols["<60%"])))
))

# --- Calculate text sizes ---
n_sites <- max(df$position)
n_seqs  <- length(unique(df$seq_name))
text_size  <- 5 * (50 / n_sites) * 0.9
label_size <- 8.4 * (5 / n_seqs) * 6

# --- Create and save the plot ---
p <- ggplot(df, aes(x = position, y = seq_name)) +
  geom_tile(aes(fill = fill), color = NA) +
  geom_text(aes(label = label), size = text_size, na.rm = TRUE) +
  labs(
    title = "RCA Alignment"
  ) +
  scale_fill_identity(
    name   = bquote(Group~"%"~"("~BLOSUM62>=.(threshold)~")"),
    breaks = cols,
    labels = names(cols),
    guide  = "legend",
    na.value = "transparent"
  ) +
  guides(fill = guide_legend(
    title.position = "top", title.hjust = 0.5,
    keywidth = unit(1, "lines"), keyheight = unit(1, "lines")
  )) +
  scale_y_discrete(limits = rev(unique(df$seq_name))) +
  scale_x_continuous(expand = c(0, 0)) +
  coord_equal(expand = FALSE) +
  theme_void() +
  theme(
    legend.position    = "top",
    legend.title.align = 0.5,
    plot.title         = element_text(size = 16, face = "bold", hjust = 0.5),
    axis.text.y        = element_text(size = label_size, hjust = 1),
    plot.margin        = margin(5, 5, 5, 100)
  )

ggsave(output_pdf, p,
       width  = canvas_width,
       height = canvas_height,
       units  = "in"
)
