
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggseqlogo)
library(rvcheck)
library(scatterpie)
library(ggforce)
library(ggrepel)
library(grid)

res.cols <- c(
    A = 'gray',
    B = 'gray',
    C = 'gray',
    D = 'gray',
    E = 'gray',
    F = '#005bbb',
    G = '#ffd500',
    H = 'gray',
    I = '#005bbb',
    J = 'gray',
    K = 'gray',
    L = '#cc00cc',
    M = 'gray',
    N = 'gray',
    O = 'gray',
    P = 'gray',
    Q = 'gray',
    R = 'gray',
    S = 'gray',
    T = 'gray',
    U = 'gray',
    V = 'gray',
    W = '#005bbb',
    Y = '#005bbb',

    x = '#d33600',
    z = '#686868',
    p = '#800080'
)

all_ecosystems <- c("Marine" = "#1500aa", "Coastal" = "#00daff", "Freshwater" = "#01cf00", "Saline" = "#a020f0")

my_ecosystem <- function(.data) {
    coastal_specific  <- c("Estuary", "Salt marsh", "Mangrove swamp", "Microbial mats", "Salt marsh, Tidal flats", "Beach", "Coral reef", "Intertidal zone")
    coastal_subtype   <- c("inlet", "intertidal zone", "coastal", "river estuary")
    coastal_isolation <- "canal|estuarine|estuary|coastal|tidal|inlet|freshwater|intertidal"
    mutate(.data, Ecosystem.Subtype = tolower(Ecosystem.Subtype), Ecosystem = case_when(
            Ecosystem.Type == "Marine" & (Ecosystem.Subtype %in% coastal_subtype | Specific.Ecosystem %in% coastal_specific | grepl(coastal_isolation, Isolation, i = T)) ~ "Coastal",
            Ecosystem.Type == "Freshwater" ~ "Freshwater",
            Ecosystem.Type == "Marine" ~ "Marine",
            Ecosystem.Type == "Non-marine Saline and Alkaline" ~ "Saline",
            T ~ NA_character_)) %>%
        mutate(Ecosystem = factor(Ecosystem, levels = names(all_ecosystems))) %>%
        filter(!is.na(Ecosystem), !grepl("sediment", Specific.Ecosystem, i = T))
}

my_rhodopsin_class <- function(.data) {
    .data %>%
        mutate(class = case_when(
            family == "XR" & motif == "DTE" ~ "x",
            family == "PR" & motif == "DTE" ~ "p",
            T ~ "z")
        )
}

my_logo_matrix <- function(.data) {
    meta <- head(.data, 1) %>%
        select(one_of("x", "y"), Ecosystem) %>%
        as.list
    color <- with(meta, case_when(
        Ecosystem == "Freshwater" ~ all_ecosystems[["Freshwater"]],
        Ecosystem == "Marine"     ~ all_ecosystems[["Marine"]],
        Ecosystem == "Saline"     ~ all_ecosystems[["Saline"]],
        Ecosystem == "Coastal"    ~ all_ecosystems[["Coastal"]]
    ))
    cs <- make_col_scheme(chars = names(res.cols), cols = res.cols)
    win <- spread(.data, window, Abundance, fill = 0) %>%
        select(-one_of("x", "y"), -Ecosystem) %>%
        column_to_rownames("class") %>% t
    tot <- colSums(win)
    win.pct <- t(t(win) / tot)
    tot.pct <- diag(-tot / sum(tot), length(tot)) %>%
        `rownames<-`(names(tot))
    mat <- rbind(win.pct, tot.pct)
    mat[is.nan(mat)] <- 0
    p <- ggseqlogo(mat, method = 'custom', namespace = names(res.cols), col_scheme = cs) +
        theme_void() +
        coord_cartesian(ylim = c(-1,1)) +
        theme(plot.background = element_rect(fill = "#FFFFFF88", color = color))
    c(meta, list(mat = mat, plot = ggplotGrob(p)))
}

my_add_grob_to_map <- function(.p, .logo) {
    x_width <- 20
    y_width <- 20
    .p + annotation_custom(.logo$plot, xmin = .logo$x - x_width/2, xmax = .logo$x + x_width/2, ymin = .logo$y - y_width/2, ymax = .logo$y + y_width/2)
}

my_add_grob_to_ecosystem <- function(.p, .logo) {
    x <- with(.logo, as.numeric(Ecosystem))
    x_width <- 1
    .p + annotation_custom(.logo$plot, xmin = x - x_width/2, xmax = x + x_width/2, ymin = -1, ymax = 1)
}

my_ggplot_ecosystems <- function(e) {
    ecosystem_df <- data.frame(Ecosystem = factor(names(all_ecosystems), levels = names(all_ecosystems)), y = 1)
    ggplot(ecosystem_df, aes(x = Ecosystem, y = y)) +
        geom_blank() +
        ylim(c(-1, 1)) +
        theme_bw()
}

my_ggplot_world <- function() {
    world <- map_data("world")
    ggplot(world) +
        geom_map(data = world, map = world, aes(long, lat, map_id = region), color = "black", fill = "lightgray", size = 0.1) +
        coord_equal() +
        xlim(c(-185, 185))
}

get_repel_coords <- function(.data, map_g, width, height) {
    grid.newpage()
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_point(aes(x, y), data = .data) +
        geom_text_repel(aes(x, y, size = Total), size = 3, label = ".", data = .data, box.padding = 8, max.overlaps = Inf)
    panel_params <- ggplot_build(g)$layout$panel_params[[1]]
    xrg <- panel_params$x.range
    yrg <- panel_params$y.range

    textrepeltree <- ggplotGrob(g) %>%
        grid.force(draw = F) %>%
        getGrob("textrepeltree", grep = T)
    children <- childNames(textrepeltree) %>%
        grep("textrepelgrob", ., value = T)

    get_xy <- function(n) {
        grob <- getGrob(textrepeltree, n)
        data.frame(
            x.repel = xrg[1] + diff(xrg) * convertX(grob$x, "native", valueOnly = T),
            y.repel = yrg[1] + diff(yrg) * convertY(grob$y, "native", valueOnly = T)
        )
    }
    lapply(children, get_xy) %>%
        bind_rows %>%
        cbind(.data) %>%
        mutate(theta = atan2(y - y.repel, x - x.repel), x.segm = x.repel + Total * cos(theta), y.segm = y.repel + Total * sin(theta))
}

my_pies <- function(.data, width, height, labeller = identity, repel = T) {
    g <- my_ggplot_world()
    if (repel) {
        .data <- get_repel_coords(.data, g, width, height)
    } else {
        .data <- mutate(.data, x.repel = x, y.repel = y, x.segm = x, y.segm = y)
    }
    win.cols <- c(WF = res.cols[["W"]], G = res.cols[["G"]])
    g +
        geom_segment(aes(x = x.segm, y = y.segm, xend = x, yend = y, color = Ecosystem), data = .data, arrow = arrow(length = unit(0.01, "npc"))) +
        geom_scatterpie(aes(x = x.repel, y = y.repel, r = Total), data = .data, color = NA, cols = c("G", "WF")) +
        geom_circle(aes(x0 = x.repel, y0 = y.repel, r = Total, color = Ecosystem), data = .data) +
        geom_scatterpie_legend(.data$Total, x = -140, y = -70, labeller = labeller) +
	scale_fill_manual(values = win.cols) +
        scale_color_manual(values = all_ecosystems) +
        theme_bw() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}
