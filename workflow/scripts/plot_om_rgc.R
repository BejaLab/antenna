
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
output_file <- unlist(snakemake@output)

metadata <- read_excel(metadata_file, sheet_name, .name_repair = "universal") %>%
    mutate(x = round(Longitude), y = round(Latitude)) %>%
    mutate(Ecosystem = "Marine")
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
    complete(class, window, nesting(x, y, Ecosystem), fill = list(Abundance = 0)) %>%
    group_by(x, y, Ecosystem) %>%
    group_split %>%
    lapply(my_logo_matrix)

p <- my_ggplot_world()
p <- Reduce(my_add_grob_to_map, data, p)
ggsave(output_file, p, height = 8, width = 16)
