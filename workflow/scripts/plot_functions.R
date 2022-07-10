
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)
library(ggseqlogo)
library(scatterpie)
library(ggforce)

cols <- c(
    A = 'orange',
    B = 'black',
    C = 'purple',
    D = 'purple',
    E = 'purple',
    F = '#f89f56',
    G = '#cc00cc',
    H = 'purple',
    I = 'black',
    J = 'black',
    K = 'purple',
    L = 'purple',
    M = 'purple',
    N = 'purple',
    O = 'black',
    P = 'green',
    Q = 'purple',
    R = 'purple',
    S = 'orange',
    T = 'orange',
    U = 'black',
    V = 'orange',
    W = '#f89f56',
    Y = '#f89f56',

    x = 'red',
    z = 'black',
    p = 'green'
)

all_ecosystems <- c("Marine", "Coastal", "Freshwater", "Saline")

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
        mutate(Ecosystem = factor(Ecosystem, levels = all_ecosystems)) %>%
        filter(!is.na(Ecosystem), !grepl("sediment", Specific.Ecosystem, i = T))
}

my_rhodopsin_class <- function(.data) {
    .data %>%
        mutate(motif = substr(positions, 1, 3), window = substr(positions, 4, 4)) %>%
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
        Ecosystem == "Freshwater" ~ "green",
        Ecosystem == "Marine"     ~ "blue",
        Ecosystem == "Saline"     ~ "purple",
        Ecosystem == "Coastal"    ~ "cyan"
    ))
    cs <- make_col_scheme(chars = names(cols), cols = cols)
    win <- spread(.data, window, Abundance, fill = 0) %>%
        select(-one_of("x", "y"), -Ecosystem) %>%
        column_to_rownames("class") %>% t
    tot <- colSums(win)
    win.pct <- t(t(win) / tot)
    tot.pct <- diag(-tot / sum(tot), length(tot)) %>%
        `rownames<-`(names(tot))
    mat <- rbind(win.pct, tot.pct)
    mat[is.nan(mat)] <- 0
    p <- ggseqlogo(mat, method = 'custom', namespace = names(cols), col_scheme = cs) +
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
    ecosystem_df <- data.frame(Ecosystem = factor(all_ecosystems, levels = all_ecosystems), y = 1)
    ggplot(ecosystem_df, aes(x = Ecosystem, y = y)) +
        geom_blank() +
        ylim(c(-1, 1)) +
        theme_bw()
}

my_ggplot_world <- function() {
    world <- map_data("world")
    ggplot(world) +
        geom_map(data = world, map = world, aes(long, lat, map_id = region), color = "black", fill = "lightgray", size = 0.1) +
        coord_equal()
}

get_repel_coords <- function(.data, map_g, width, height) {
    grid.newpage()
    pushViewport(viewport(width = width, height = height))
    g <- map_g +
        geom_point(aes(x, y), data = .data) +
        geom_text_repel(aes(x, y, size = Total), label = "O", data = .data, box.padding = 8, max.overlaps = Inf)
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

my_pies <- function(.data, res.cols = c("G", "WF"), width, height) {
    .data <- by_location
    g <- my_ggplot_world()
    .data <- get_repel_coords(.data, g, width, height)
    g +
        geom_segment(aes(x = x.segm, y = y.segm, xend = x, yend = y, color = Ecosystem), data = .data, arrow = arrow(length = unit(0.01, "npc"))) +
        geom_scatterpie(aes(x = x.repel, y = y.repel, r = Total), data = .data, color = NA, cols = res.cols, alpha = .8) +
        geom_circle(aes(x0 = x.repel, y0 = y.repel, r = Total, color = Ecosystem), data = .data) +
        geom_scatterpie_legend(.data$Total, x = -140, y = -70) +
	scale_fill_manual(values = c(G = '#cc00cc', WF = '#f89f56')) +
        scale_color_manual(values = c(Freshwater = "green", "Marine" = "blue", "Saline" = "purple", "Coastal" = "cyan")) +
        theme_bw() +
        theme(axis.title.x = element_blank(), axis.title.y = element_blank())
}
