
snakemake@source("plot_functions.R")
library(readxl)
library(car)
library(lme4)

with(snakemake@input, {
    tsv_file       <<- tsv
    metadata_file  <<- metadata
    abundance_file <<- abundance
    crtZ_file      <<- crtZ
    crtZ_abundance_file <<- crtZ_abundance
})
with(snakemake@params, {
    sheet_name <<- sheet
})
with(snakemake@output, {
    map_file <<- map
    sum_file <<- sum
    dep1_file <<- dep1
    lat1_file <<- lat1
    test1_file <<- test1
    dep2_file <<- dep2
    lat2_file <<- lat2
    test2_file <<- test2
    win_wl_file <<- win_wl
})

metadata <- read_excel(metadata_file, sheet_name, .name_repair = "universal") %>%
    mutate(x = round(Longitude), y = round(Latitude)) %>%
    mutate(Ecosystem = factor("Marine"))
abundance <- read.table(abundance_file, header = T) %>%
    gather(Sample, Abundance, -OMRGC_ID) %>%
    left_join(metadata, by = c(Sample = "PANGAEA.sample.id"))
all.data <- read.table(tsv_file, header = T, sep = "\t", na.strings = c("NA", "")) %>%
    filter(!is.na(wl_mean)) %>%
    my_rhodopsin_class %>%
    left_join(abundance, by = c(record_id = "OMRGC_ID")) %>%
    mutate(Depth = suppressWarnings(as.numeric(Depth..nominal)), fenestrated = case_when(window == "G" ~ T, window %in% c("W","F") ~ F, T ~ NA))

motifs_counts <- all.data %>%
    distinct(checksum, .keep_all = T) %>%
    group_by(family, motif) %>%
    summarize(Abundance = n())
motifs_weighted <- all.data %>%
    group_by(family, motif) %>%
    summarize(Abundance = sum(Abundance))

p <- ggplot(motifs_counts, aes(fill = motif, x = family, y = Abundance)) +
    geom_bar(position = "stack", stat = "identity")
ggsave("tmp1-om-rgc.pdf", p)
p <- ggplot(motifs_weighted, aes(fill = motif, x = family, y = Abundance)) +
    geom_bar(position = "stack", stat = "identity")
ggsave("tmp2-om-rgc.pdf", p)

ratio.data <- filter(all.data, class %in% c("x", "p"), !is.na(fenestrated), !is.na(Depth)) %>%
    # filter(Polar == "Non polar") %>%
    group_by(Sample, x, y, Ecosystem, Depth, Station, Polar, fenestrated) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    spread(fenestrated, Abundance, fill = 0) %>%
    mutate(Abundance = `TRUE` + `FALSE`, ratio = `TRUE` / Abundance)
model <- lmer(ratio ~ log10(Depth + 1) + abs(y) + (1 | Station), data = ratio.data)
capture.output(Anova(model, type = "II"), file = test1_file)

p <- ggplot(ratio.data, aes(x = Depth, y = ratio, size = Abundance)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_continuous(trans = 'log10')
ggsave(dep1_file, p)
p <- ggplot(ratio.data, aes(x = abs(y), y = ratio, size = Abundance)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x)
ggsave(lat1_file, p)

weighted_mean <- function(x, w) {
    sum(w * x) / sum(w)
}
weighted_sd <- function(x, w) {
    x_mean <- weighted_mean(x, w)
    M <- sum(w > 0)
    sqrt( sum(w * (x - x_mean)^2) * M / (M-1) / sum(w) )
}

wl.data <- filter(all.data, class %in% c("x", "p"), !is.na(fenestrated), !is.na(Depth)) %>%
    group_by(x, y, Ecosystem, Depth, Station, Polar) %>%
    summarize(wl_avg = weighted_mean(wl_mean, Abundance), wl_sd = weighted_sd(wl_mean, Abundance), Abundance = sum(Abundance), .groups = "drop")
model <- lmer(wl_avg ~ log10(Depth + 1) + abs(y) + (1 | Station), data = wl.data)
capture.output(Anova(model, type = "II"), file = test2_file)

p <- ggplot(wl.data, aes(x = Depth, y = wl_avg, size = Abundance)) +
    geom_pointrange(aes(ymin = wl_avg - wl_sd, ymax = wl_avg + wl_sd)) +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_size_continuous(range = c(0.1,0.5)) +
    scale_x_continuous(trans = 'log10')
ggsave(dep2_file, p)
p <- ggplot(wl.data, aes(x = abs(y), y = wl_avg, size = Abundance)) +
    geom_pointrange(aes(ymin = wl_avg - wl_sd, ymax = wl_avg + wl_sd)) +
    scale_size_continuous(range = c(0.1,0.5)) +
    geom_smooth(method = 'lm', formula = y ~ x)
ggsave(lat2_file, p)

crtZ.abundance <- read.table(crtZ_abundance_file, header = T) %>%
    gather(Sample, Abundance, -OMRGC_ID) %>%
    left_join(metadata, by = c(Sample = "PANGAEA.sample.id"))
crtZ.data <- read.table(crtZ_file) %>%
    select(OMRGC_ID = "V1", evalue = "V5") %>%
    filter(evalue < 1e-15) %>%
    left_join(crtZ.abundance, by = "OMRGC_ID") %>%
    group_by(Sample, x, y, Station) %>%
    summarize(CrtZ.abundance = sum(Abundance), .groups = "drop") %>%
    left_join(select(ratio.data, Sample, Depth, Abundance, `TRUE`, `FALSE`), by = "Sample") %>%
    mutate(ratio = CrtZ.abundance / `TRUE`)

model <- lmer(ratio ~ log10(Depth + 1) + abs(y) + (1 | Station), data = crtZ.data)
capture.output(Anova(model, type = "II"), file = test2_file)

p <- ggplot(crtZ.data, aes(x = CrtZ.abundance, y = `TRUE`)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10')
ggsave("TRUE.pdf", p)

p <- ggplot(crtZ.data, aes(x = CrtZ.abundance, y = `FALSE`)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10')
ggsave("FALSE.pdf", p)

p <- ggplot(crtZ.data, aes(x = Depth, y = ratio, size = Abundance)) +
    geom_point() +
    geom_smooth(method = 'lm', formula = y ~ x) +
    scale_x_continuous(trans = 'log10') +
    scale_y_continuous(trans = 'log10')
ggsave("tmp.pdf", p)

prots <- distinct(all.data, record_id, .keep_all=T)
p <- ggplot(prots, aes(x = window, y = wl_mean)) + geom_boxplot()
ggsave(win_wl_file, p)

data <- group_by(all.data, x, y, Ecosystem, class, window) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    complete(class, window, nesting(x, y, Ecosystem), fill = list(Abundance = 0))

by_location <- group_by(data, x, y, Ecosystem) %>%
    filter(class %in% c("p", "x")) %>%
    summarize(Total = log10(sum(Abundance)), G = sum(Abundance[window %in% "G"], na.rm = T), WF = sum(Abundance[window %in% c("W", "F")], na.rm = T))
by_ecosystem <- group_by(data, class, window, Ecosystem) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    my_logo_matrix

width <- 16
height <- 8
p <- my_pies(by_location, width, height)
ggsave(map_file, p, height = height, width = width)

p <- my_ggplot_ecosystems()
p <- my_add_grob_to_ecosystem(p, by_ecosystem)
ggsave(sum_file, p, height = 3, width = 4)
