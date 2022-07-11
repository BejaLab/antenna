
snakemake@source("plot_functions.R")
library(readxl)

with(snakemake@input, {
    tsv_file       <<- tsv
    metadata_file  <<- metadata
    abundance_file <<- abundance
})
with(snakemake@params, {
    sheet_name <<- sheet
})
with(snakemake@output, {
    map_file <<- map
    sum_file <<- sum
})

metadata <- read_excel(metadata_file, sheet_name, .name_repair = "universal") %>%
    mutate(x = round(Longitude), y = round(Latitude)) %>%
    mutate(Ecosystem = factor("Marine"))
abundance <- read.table(abundance_file, header = T) %>%
    gather(Sample, Abundance, -OMRGC_ID) %>%
    left_join(metadata, by = c(Sample = "PANGAEA.sample.id"))
data <- read.table(tsv_file, header = T, sep = "\t", na.strings = c("NA", "")) %>%
    filter(!is.na(wl_mean)) %>%
    my_rhodopsin_class %>%
    left_join(abundance, by = c(record_id = "OMRGC_ID")) %>%
    group_by(x, y, Ecosystem, class, window) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    ungroup %>%
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
