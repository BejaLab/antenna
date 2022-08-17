
snakemake@source("plot_functions.R")
library(car)
library(lme4)
library(ggpubr)
library(emmeans)
library(dunn.test)
library(photobiology)

with(snakemake@input, {
    om_rgc_file <<- om_rgc
    gem_file    <<- gem
    jgi_files   <<- jgi
})
output_file <- unlist(snakemake@output)

read.data <- function(fnames) {
    lapply(fnames, read.table, header = T, sep = "\t", fill = T, quote = "") %>%
        lapply(select, "motif", "blasso", "window", "checksum", "family", "wl_mean", "wl_sd") %>%
        bind_rows
}
residues <- c("G", "F", "W")

data <- bind_rows(
    read.data(om_rgc_file),
    read.data(gem_file),
    read.data(jgi_files)
) %>% my_rhodopsin_class %>%
    filter(class %in% c("x", "p"), window %in% residues, !is.na(wl_mean), !grepl("X", blasso)) %>%
    distinct(checksum, .keep_all = T) %>%
    group_by(blasso) %>%
    filter(n() > 1) %>%
    group_by(blasso, window, family, .keep_all = T) %>%
    summarize(wl_mean = mean(wl_mean), n = n(), .groups = "drop") %>%
    mutate(window = factor(window, levels = residues)) %>%
    mutate(family = factor(family)) %>%
    mutate(color = w_length2rgb(wl_mean))

kruskal <- split(data, f = data$family) %>%
    lapply(function(x) kruskal.test(wl_mean ~ window, data = x)) %>%
    lapply(`[`, c("statistic","p.value")) %>%
    bind_rows(.id = "family")
dunn <- split(data, f = data$family) %>%
    lapply(function(x) with(x, dunn.test(wl_mean, window, method = "bh"))) %>%
    lapply(`[`, c("Z","P.adjusted","comparisons")) %>%
    bind_rows(.id = "family") %>%
    separate(comparisons, into = c("group1", "group2"), sep = " - ", remove = F) %>%
    mutate(Q.lab = case_when(P.adjusted < 0.001 ~ "***", P.adjusted < 0.01 ~ "**", P.adjusted < 0.05 ~ "*", T ~ "ns")) %>%
    mutate(y.position = 580 - (group1 == "F" | group2 == "F") * 7)

n_fun <- function(x) {
    data.frame(y = 450, label = sprintf("n = %d", length(x)))
}
p <- ggplot(data, aes(x = window, y = wl_mean)) +
    geom_jitter(aes(color = color, size = n)) +
    geom_boxplot(color = "red", alpha = 0.1) +
    scale_size_continuous(range = c(0.1, 5), breaks = c(8, 40, 200, 1000, 5000), limits = c(1, 5000)) +
    scale_colour_identity() +
    stat_summary(fun.data = n_fun, geom = "text", size = 3) +
    facet_grid(. ~ family) +
    ylim(c(450, 580)) +
    xlab("Fenestration residue") + ylab("Predicted absorption maximum") +
    theme_bw() +
    stat_pvalue_manual(data = dunn, label = "Q.lab", xmin = "group1", xmax = "group2", y.position = "y.position", size = 2)
ggsave(output_file, p, width = 5.5, height = 2.5)
