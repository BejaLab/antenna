
snakemake@source("plot_functions.R")

with(snakemake@input, {
    tsv_files     <<- tsv
    metadata_file <<- metadata
})
with(snakemake@output, {
    map_file <<- map
    sum_file <<- sum
})

metadata <- read.table(metadata_file, header = T, sep = "\t", quote = "", comment.char = "", na.strings = c("", "n/a")) %>%
    filter(JGI.Data.Utilization.Status == "Unrestricted", GOLD.Analysis.Project.Type == "Metagenome Analysis") %>%
    my_ecosystem %>%
    mutate(x = round(Longitude), y = round(Latitude))
data <- lapply(tsv_files, read.table, header = T, sep = "\t", quote = "", comment.char = "") %>%
    lapply(select, "record_id", "positions", "family", "wl_mean", "wl_sd", "Locus.Tag", "Genome.ID", "NCBI.Biosample.Accession", "Scaffold.Read.Depth") %>%
    bind_rows %>%
    distinct(record_id, .keep_all = T) %>%
    group_by(Genome.ID) %>%
    arrange(-n()) %>%
    group_by(NCBI.Biosample.Accession) %>%
    filter(Genome.ID == first(Genome.ID)) %>%
    my_rhodopsin_class %>%
    left_join(metadata, by = c(Genome.ID = "taxon_oid")) %>%
    mutate(Abundance = Scaffold.Read.Depth) %>%
    filter(!is.na(x)) %>%
    group_by(x, y) %>%
    filter(n_distinct(Ecosystem) == 1) %>%
    filter(sum(class %in% c("x", "p")) > 9) %>%
    ungroup %>%
    complete(class, window, nesting(x, y, Ecosystem), fill = list(Abundance = 0)) %>%
    group_by(x, y, Ecosystem) %>%
    group_by(x, y, Ecosystem, class, window) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop")
#by_location <- group_by(data, x, y, Ecosystem) %>%
#    group_split %>%
#    lapply(my_logo_matrix)
by_location <- group_by(data, x, y, Ecosystem) %>%
    filter(class %in% c("p", "x")) %>%
    summarize(Total = log10(sum(Abundance)), G = sum(Abundance[window %in% "G"], na.rm = T), WF = sum(Abundance[window %in% c("W", "F")], na.rm = T))
by_ecosystem <- group_by(data, class, window, Ecosystem) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Ecosystem) %>%
    group_split %>%
    lapply(my_logo_matrix)

#p <- my_ggplot_world()
#p <- Reduce(my_add_grob_to_map, by_location, p)

width <- 16
height <- 8
p <- my_pies(by_location, width, height)
ggsave(map_file, p, height = height, width = width)


p <- my_ggplot_ecosystems()
p <- Reduce(my_add_grob_to_ecosystem, by_ecosystem, p)
ggsave(sum_file, p, height = 8, width = 16)
