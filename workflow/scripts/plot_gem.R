
snakemake@source("plot_functions.R")

with(snakemake@input, {
    metadata_file <<- metadata
    tsv_file  <<- tsv
    crt_files <<- crt
    colors_file <<- colors
})
with(snakemake@output, {
    tax_file <<- tax
    sum_file <<- sum
    car_file <<- car
})

aliases <- list(
    BCD = "BLH",
    Lycopene_cycl = "CrtY",
    crtI_fam = "CrtI"
)

metadata <- read.table(metadata_file, header = T, sep = "\t", comment.char = "") %>%
    filter(ecosystem_category == "Aquatic", habitat != "Groundwater") %>%
    mutate(Ecosystem.Type = ecosystem_type, Specific.Ecosystem = habitat, Ecosystem.Subtype = habitat, Isolation = habitat, x = longitude, y = latitude) %>%
    my_ecosystem

taxa <- read.table(colors_file, sep = "\t", comment.char = "") %>%
    arrange(V1) %>%
    with(setNames(V2,V1))

genes <- lapply(crt_files, read.table) %>%
    bind_rows %>%
    select(1, 3, 5) %>%
    setNames(c("record_id", "gene", "E_value")) %>%
    filter(E_value < 1e-15) %>%
    separate(record_id, into = c("genome_id", "Gene.ID"), sep = "/") %>%
    mutate(gene = recode(gene, !!!aliases)) %>%
    distinct(genome_id, gene)

all_data <- read.table(tsv_file, header = T, sep = "\t", quote = "", comment.char = "") %>%
    rename(Taxonomy = "ecosystem") %>% # fix bug
    filter(ecosystem_category == "Aquatic", !grepl("-", motif)) %>%
    mutate(Ecosystem.Type = ecosystem_type, Specific.Ecosystem = habitat, Ecosystem.Subtype = habitat, Isolation = habitat, x = longitude, y = latitude, Abundance = 1) %>%
    my_ecosystem %>% 
    my_rhodopsin_class %>%
    separate(record_id, into = c("genome_id", "Gene.ID"), sep = "/") %>%
    mutate(window_type = case_when(
        window %in% c("W", "F") ~ "FW",
        window %in% c("G") ~ "G",
        T ~ "other"
    )) %>%
    distinct(genome_id, class, window, .keep_all = T)

data <- complete(all_data, class, window, nesting(x, y, Ecosystem), fill = list(Abundance = 0)) %>%
    group_by(x, y, Ecosystem, class, window) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop")

by_ecosystem <- group_by(data, class, window, Ecosystem) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Ecosystem) %>%
    group_split %>%
    lapply(my_logo_matrix)

aliases <- list(Cyanobacteriota = "Cyanobacteria")
by_taxa <- distinct(all_data, otu_id, .keep_all = T) %>%
    separate_rows(Taxonomy, sep = ";") %>%
    filter(class %in% c("x", "p"), window_type %in% c("FW", "G")) %>%
    separate(Taxonomy, into = c("Rank", "Name"), sep = "__") %>%
    filter(Rank == "p") %>%
    group_by(class, window_type, Ecosystem, Name) %>%
    summarize(Abundance = sum(Abundance), .groups = "drop") %>%
    group_by(Name) %>%
    mutate(Name = ifelse(sum(Abundance) > 10, Name, NA)) %>%
    mutate(Name = recode(Name, !!!aliases))

tmp <- left_join(metadata, all_data, by = "genome_id", suffix = c("", ".all_data")) %>%
    distinct(genome_id, .keep_all = T) %>%
    group_by(Ecosystem) %>%
    summarize(num_px = sum(class %in% c("x", "p") & window_type %in% c("FW", "G")), num_all = n())

all_taxa <- unique(na.omit(by_taxa$Name))
missing_taxa <- all_taxa[! all_taxa %in% names(taxa)]
if (length(missing_taxa) > 0) {
    write("Taxa not in the color file:", stderr())
    write(missing_taxa, stderr())
    q()
}

by_gene <- all_data %>%
    filter(class %in% c("x", "p"), window_type %in% c("WF", "G")) %>%
    left_join(genes, by = "genome_id") %>%
    mutate(present = T) %>%
    complete(nesting(genome_id, class, window_type), gene, fill = list(present = F)) %>%
    filter(!is.na(gene)) %>%
    group_by(class, window_type, gene, present) %>%
    summarize(n = n())
p <- ggplot(by_gene, aes(x = gene, fill = present, y = n)) +
    geom_bar(position = "stack", stat = "identity") +
    facet_grid(class ~ window_type) +
    theme_bw()
ggsave(car_file, p)

p <- ggplot(by_taxa, aes(fill = Name, x = Ecosystem, y = Abundance)) +
    geom_bar(position = "stack", stat = "identity") +
    scale_fill_manual(values = taxa) +
    facet_grid(class ~ window_type) +
    theme_bw()
ggsave(tax_file, p, width = 7, height = 5)

p <- my_ggplot_ecosystems()
p <- Reduce(my_add_grob_to_ecosystem, by_ecosystem, p)
ggsave(sum_file, p, height = 8, width = 16)
