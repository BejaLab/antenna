
library(treeio)
library(ggtree)
library(ape)
library(phytools)
library(dplyr)
library(tidyr)
library(ggplot2)
library(phangorn)
library(ggnewscale)
library(castor)
library(seqinr)

if (interactive()) {
    Snakemake <- setClass("Snakemake", slots = list(input = "list", output = "list"))
    snakemake <- Snakemake(
        input = list(outgroup = "input/outgroups.fasta", tree = "analysis/phylogeny/rhodopsins.treefile", metadata = "analysis/parse/phylogeny.tsv"),
        output = list("tmp.pdf")
    )
}

with(snakemake@input, {
    outgroup_file <<- outgroup
    tree_file     <<- tree
    tsv_file      <<- tsv
    colors_file   <<- colors
    a2m_file      <<- a2m
})
with(snakemake@output, {
    output_file_small <<- small
    output_file_big   <<- big
    output_jtree      <<- jtree
})

taxa <- read.table("metadata/taxa.txt", sep = "\t", comment.char = "") %>%
    arrange(1) %>%
    with(setNames(V2, V1))

outgroups <- names(read.fasta(outgroup_file))
tsv <- read.table(tsv_file, header = T, sep = "\t", na.strings = "", fill = T) %>%
    select(-target)
metadata <- read.fasta(a2m_file, seqtype = "AA", as.string = T) %>%
    {data.frame(label = names(.), sequence = as.character(.))} %>%
    left_join(tsv, by = c(label = "record_id")) %>%
    mutate(is_outgroup = label %in% outgroups) %>%
    mutate(Alias = gsub(",.+", "", Alias)) %>%
    mutate(Activity = ifelse(grepl("\\]$", Activity), NA, gsub("[][]", "", Activity))) %>%
    mutate(Highlight = !is.na(Highlight)) %>%
    mutate(D85 = substr(motif, 1, 1), T89 = substr(motif, 2, 2), D96 = substr(motif, 3, 3), G156 = window)

tree <- read.tree(tree_file)

tree.unrooted <- as_tibble(tree) %>%
    left_join(metadata, by = "label") %>%
    `class<-`(c("tbl_tree", "data.frame")) %>%
    as.treedata
write.jtree(tree.unrooted, file = output_jtree)

tree.tib <- ape::root(tree, outgroups, edgelabel = T, resolve.root = T) %>%
    drop.tip(outgroups) %>%
    as_tibble %>%
    mutate(support = suppressWarnings(as.numeric(label))) %>%
    left_join(metadata, by = "label") %>%
    mutate(is_outgroup = ifelse(label %in% outgroups, T, NA)) %>%
    `class<-`(c("tbl_tree", "data.frame"))
tree.phylo <- as.treedata(tree.tib)@phylo

clustalx <- c(
    A = "BLUE",
    I = "BLUE",
    L = "BLUE",
    M = "BLUE",
    F = "BLUE",
    W = "BLUE",
    V = "BLUE",
    C = "BLUE",
    K = "RED",
    R = "RED",
    E = "MAGENTA",
    D = "MAGENTA",
    N = "GREEN",
    Q = "GREEN",
    S = "GREEN",
    T = "GREEN",
    G = "ORANGE",
    P = "YELLOW",
    H = "CYAN",
    Y = "CYAN"
)

clades <- filter(tree.tib, ! node %in% parent) %>%
    pull(Clade) %>%
    as.factor
hsp <- hsp_max_parsimony(tree.phylo, as.numeric(clades)) %>%
    `$`("likelihoods") %>%
    data.frame(check.names = F) %>%
    mutate(node = row_number()) %>%
    gather(hsp, prob, -node) %>%
    filter(prob > 0.9) %>%
    mutate(hsp = levels(clades)[as.numeric(hsp)])

tree <- left_join(tree.tib, hsp, by = "node") %>%
    mutate(hsp.parent = .[match(parent, node),"hsp"]) %>%
    mutate(hsp = ifelse(!is.na(hsp.parent) & hsp == hsp.parent, NA, hsp)) %>%
    as.treedata

p_small <- ggtree(tree, layout = "circular") +
    geom_highlight(mapping = aes(subset = !is.na(hsp), fill = hsp), alpha = 0.1) +
    geom_treescale() +
    geom_tiplab(mapping = aes(subset = !is.na(Alias), label = Alias, size = Highlight), offset = 0.1) +
    # geom_tippoint(mapping = aes(subset = !is.na(Activity), color = Activity), size = 1) +
    scale_size_manual(values = c(2.5, 5)) + new_scale("size") + # labels
    geom_point2(aes(subset = !is.na(support) & support >= 90, size = support >= 95), color = "darkgray") +
    scale_size_manual(values = c(0.5, 1)) + # support values
    new_scale("color") +
    geom_text2(aes(label = G156, color = G156, angle = angle - 90, x = 4.3), size = 3) +
    scale_color_manual(values = clustalx) # residues

p_big <- ggtree(tree, aes(color = Taxon), layout = "rectangular") +
    geom_highlight(mapping = aes(subset = !is.na(hsp), fill = hsp), alpha = 0.1) +
    scale_color_manual(values = taxa) + new_scale("color") +
    geom_treescale() +
    geom_tiplab(mapping = aes(label = sprintf("%s [%s]", ifelse(is.na(Alias), label, Alias), Organism)), size = 2, offset = 0.1) +
    geom_tippoint(mapping = aes(subset = !is.na(Activity), color = Activity), size = 1) +
    geom_point2(aes(subset = !is.na(support) & support >= 90, size = support >= 95, x = branch), shape = 15, color = "darkgray") +
    scale_size_manual(values = c(0.5, 1)) + # support values
    new_scale("color") +
    geom_text2(aes(label = D85, color = D85, x = 4.8), size = 1.5) +
    geom_text2(aes(label = T89, color = T89, x = 5.0), size = 1.5) +
    geom_text2(aes(label = D96, color = D96, x = 5.2), size = 1.5) +
    geom_text2(aes(label = G156, color = G156, x = 5.4), size = 2) +
    scale_color_manual(values = clustalx) # residues

ggsave(output_file_small, p_small, height = 7, width = 8)
ggsave(output_file_big,   p_big,   height = 6.5, width = 5)
