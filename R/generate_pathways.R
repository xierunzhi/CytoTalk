# CONSTANTS


HEX <- c(0:9, LETTERS[seq_len(6)])
URL_GENECARDS = "https://www.genecards.org/cgi-bin/carddisp.pl?gene="
URL_WIKIPI = "https://hagrid.dbmi.pitt.edu/wiki-pi/index.php/search?q="
FORM_NODE = paste0(
    "\"%s\" [label = \"%s\" href = \"",
    URL_GENECARDS,
    "%s\" width = %s height = %s fillcolor = \"%s\"];"
)
FORM_EDGE = paste0(
    "\"%s\" -> \"%s\" [href = \"",
    URL_WIKIPI,
    "%s+%s\" penwidth = %s];"
)
FORM_GV = trimws("
digraph {\n
pad = 0.5;
layout = fdp;
labeljust = l;
splines = true;
overlap = false;
outputorder = \"edgesfirst\";\n
node [fixedsize = true target=\"_blank\"];
node [fontname = \"Arial\" fontsize = 10 style = filled];
edge [arrowhead = none constraint = false target=\"_blank\"];\n
subgraph cluster0 {\n
rank = same;
margin = 30;
color = none;
style = filled;
overlap = portho;
fillcolor = \"#EEEEEE\";
node [shape = circle];\n
%s\n
}\n
subgraph cluster1 {\n
rank = same;
margin = 30;
color = none;
style = filled;
overlap = portho;
fillcolor = \"#EEEEEE\";
node [shape = circle];\n
%s\n
}\n
// cluster external horizontal order
%s
// cluster external
edge [color = limegreen arrowhead = normal];
%s\n
}")


# FUNCTIONS


#' @noRd
extract_best_network <- function(type_a, type_b, dir_out) {
    # format filepaths
    fpath_net <- file.path(dir_out, "IntegratedNetwork.cfg")
    fpath_pval <- file.path(dir_out, "PCSF_EdgeTestValues.txt")
    fpath_edge <- file.path(dir_out, "PCSF_EdgeOccurance.txt")
    fpath_pem <- file.path(dir_out, "GeneCellTypeSpecific.txt")
    fpath_out <- file.path(dir_out, "PCSF_Network.txt")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!all(file.exists(fpath_net, fpath_pval, fpath_edge, fpath_pem))) {
        stop(sprintf("cannot find input file(s)"))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return()
    }

    # process network cfg
    lst_net <- strsplit(readLines(fpath_net), "\t")
    lengths <- vapply(lst_net, length, integer(1))
    df_net_nodes <- as.data.frame(do.call(rbind, lst_net[lengths == 2]))
    df_net_edges <- as.data.frame(do.call(rbind, lst_net[lengths == 3]))
    names(df_net_edges) <- c("node1", "node2", "cost")

    # load in file
    df_pval <- suppressMessages(vroom::vroom(fpath_pval, progress = FALSE))
    df_edge <- suppressMessages(vroom::vroom(fpath_edge, progress = FALSE))
    df_pem <- suppressMessages(vroom::vroom(fpath_pem, progress = FALSE))

    # column to rowname
    df_pem <- tibble::column_to_rownames(df_pem, names(df_pem)[1])

    # find cell type in PEM matrix
    vec_pem_a <- as.numeric(df_pem[, names(df_pem) == type_a])
    vec_pem_b <- as.numeric(df_pem[, names(df_pem) == type_b])

    # remove artificial nodes
    df_edge <- df_edge[df_edge$node1 != "ARTI" & 0.5 <= df_edge$omega, ]

    # which has the best score?
    index <- which(df_pval$pval == min(df_pval$pval))
    opt_beta <- df_pval[[index, "beta"]]
    opt_omega <- df_pval[[index, "omega"]]

    # merge with costs
    df_net <- df_edge[df_edge$beta == opt_beta & df_edge$omega == opt_omega, ]
    df_net <- merge(df_net[, 3:4], df_net_edges)

    # add node prizes
    index1 <- match(df_net$node1, df_net_nodes[, 1])
    index2 <- match(df_net$node2, df_net_nodes[, 1])
    df_net$node1_prize <- as.numeric(df_net_nodes[index1, 2])
    df_net$node2_prize <- as.numeric(df_net_nodes[index2, 2])

    # extract cell types
    f <- function(x) {
        ifelse(rev(unlist(strsplit(x, "")))[1] == "A", type_a, type_b)
    }
    df_net$node1_type <- vapply(df_net$node1, f, character(1))
    df_net$node2_type <- vapply(df_net$node2, f, character(1))

    # simplify names
    f <- function(x) gsub("_TypA", "", gsub("_TypB", "_", x))
    df_net$node1 <- f(df_net$node1)
    df_net$node2 <- f(df_net$node2)

    # crosstalk edges
    df_net$is_ct_edge <- (df_net$node1_type != df_net$node2_type)

    # PEM score
    index <- match(gsub("_$", "", df_net$node1), toupper(rownames(df_pem)))
    df_net$node1_pem <- vec_pem_a[index]
    index <- match(gsub("_$", "", df_net$node2), toupper(rownames(df_pem)))
    df_net$node2_pem <- vec_pem_b[index]

    # reorder columns
    index <- c(
        "node1", "node2", "node1_type", "node2_type",
        "node1_prize", "node2_prize", "node1_pem", "node2_pem",
        "is_ct_edge", "cost"
    )
    df_net <- df_net[, index]

    # write out full table
    vroom::vroom_write(df_net, fpath_out, progress = FALSE)
    NULL
}

#' @noRd
write_network_sif <- function(dir_out) {
    # format dir path
    dir_out_cs <- file.path(dir_out, "cytoscape")

    # make sure it exists
    if (!dir.exists(dir_out_cs)) {
        dir.create(dir_out_cs)
    }

    # format filepaths
    fpath_net <- file.path(dir_out, "PCSF_Network.txt")
    fpath_node <- file.path(dir_out_cs, "CytoscapeNodes.txt")
    fpath_edge <- file.path(dir_out_cs, "CytoscapeEdges.txt")
    fpath_sif <- file.path(dir_out_cs, "CytoscapeNetwork.sif")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_net)) {
        stop(sprintf("cannot find input file: %s", fpath_net))
    } else if (all(file.exists(fpath_node, fpath_edge, fpath_sif))) {
        message("files already exist, continuing")
        return()
    }

    # load in file
    df_net <- suppressMessages(vroom::vroom(fpath_net, progress = FALSE))

    # format edges
    edge_types <- ifelse(df_net$is_ct_edge, "pr", "pp")
    edge_names <- sprintf("%s (%s) %s", df_net$node1, edge_types, df_net$node2)
    edges <- gsub("[()]", "", edge_names)

    # create node table
    index1 <- c("node1", "node1_type", "node1_prize", "node1_pem")
    index2 <- c("node2", "node2_type", "node2_prize", "node2_pem")
    df_node <- data.frame(rbind(
        as.matrix(df_net[, index1]),
        as.matrix(df_net[, index2])
    ))

    # naming and type conversion
    names(df_node) <- c("node", "type", "prize", "pem")
    df_node <- utils::type.convert(df_node, as.is = TRUE)

    # create edge table
    df_edge <- cbind(edge_names, df_net[, c("cost", "is_ct_edge")])
    names(df_edge) <- c("edge", "cost", "is_ct_edge")

    # write out sif
    writeLines(paste(edges, collapse = "\n"), fpath_sif)

    # write out node and edge tables
    vroom::vroom_write(df_node, fpath_node, progress = FALSE)
    vroom::vroom_write(df_edge, fpath_edge, progress = FALSE)
    NULL
}

#' @noRd
write_neighborhood_gv <- function(row, df_net_sub, dir_out, sub=TRUE) {
    # format dir path
    dir_out_gv <- file.path(dir_out, "graphviz")

    # make sure it exists
    if (!dir.exists(dir_out_gv)) {
        dir.create(dir_out_gv)
    }

    # format filepath
    fname <- paste0(paste(gsub("_", "", row[c(1, 2)]), collapse = "_"), ".gv")
    fpath <- file.path(dir_out_gv, fname)

    # prepare nodes
    index1 <- c("node1", "node1_type", "node1_prize", "node1_pem")
    index2 <- c("node2", "node2_type", "node2_prize", "node2_pem")
    df_node <- data.frame(rbind(
        as.matrix(df_net_sub[, index1]),
        as.matrix(df_net_sub[, index2])
    ))

    df_node <- df_node[!duplicated(df_node), ]
    df_node <- utils::type.convert(df_node, as.is = TRUE)
    names(df_node) <- c("node", "type", "prize", "pem")

    # string format nodes
    ew <- function(x) endsWith(x, "_")
    index_nodes <- ew(df_node$node)
    clean <- trimws(df_node$node, whitespace = "_")
    size <- 0.5 + 9.5 * df_node$prize^2
    color <- grDevices::hsv(
        ifelse(ew(df_node$node), 0.02, 0.55), (1 - df_node$pem), 1)
    nodes <- sprintf(FORM_NODE, df_node$node, clean, clean, size, size, color)

    # string format edges
    index_edges <- ew(df_net_sub$node1) + ew(df_net_sub$node2)
    c1 <- trimws(df_net_sub$node1, whitespace = "_")
    c2 <- trimws(df_net_sub$node2, whitespace = "_")
    size <-  0.75 + 3.25 * (1 - df_net_sub$cost)^2
    size <- ifelse(df_net_sub$is_ct_edge, 2, 1) * size
    style <- ifelse(df_net_sub$is_ct_edge, "dashed", "solid")
    edges <- sprintf(
        FORM_EDGE, df_net_sub$node1, df_net_sub$node2, c1, c2, size, style)

    # cell type names
    type_a <- df_node[!ew(df_node$node), "type"][1]
    type_b <- df_node[ew(df_node$node), "type"][1]

    # string format graph
    graph <- sprintf(FORM_GV,
        paste0(
            sprintf("label = \"%s\";\ntooltip = \"%s\";\n", type_a, type_a),
            paste0(nodes[!index_nodes], collapse = "\n"), "\n",
            paste0(edges[index_edges == 0], collapse = "\n"),
            collapse = "\n"
        ),
        paste0(
            sprintf("label = \"%s\";\ntooltip = \"%s\";\n", type_b, type_b),
            paste0(nodes[index_nodes], collapse = "\n"), "\n",
            paste0(edges[index_edges == 2], collapse = "\n"),
            collapse = "\n"
        ),
        ifelse(sub, sprintf(
            "%s [style = invis, constraint = true];\n",
            gsub("\\s+\\[.+", "", edges[index_edges == 1][1])
        ), ""),
        paste0(
            edges[index_edges == 1],
            collapse = "\n"
        )
    )

    # write out
    writeLines(graph, fpath)

    # return filepath
    fpath
}

render_gv <- function(fpath) {
    clean <- gsub("\\.gv", "", fpath)
    syscmd <- function(cmd) { system(cmd, ignore.stderr = TRUE) }
    syscmd(sprintf("dot %s -Tsvg -o %s.svg", fpath, clean))
    syscmd(sprintf("dot %s -Tpng -o %s.png", fpath, clean))
}

#' @noRd
write_pathways_gv <- function(dir_out, depth) {
    # format filepaths
    fpath_net <- file.path(dir_out, "PCSF_Network.txt")

    # stop if input file doesn't exist
    if (!file.exists(fpath_net)) {
        stop(sprintf("cannot find input file: %s", fpath_net))
    }

    # load in file
    df_net <- suppressMessages(vroom::vroom(
        fpath_net, progress = FALSE, na = c("", "NA", "-Inf")
    ))
    vec_pem <- df_net[ c("node1_pem", "node2_pem")]
    vec_pem[suppressWarnings(is.na(vec_pem))] <- min(vec_pem, na.rm = TRUE)
    df_net[ c("node1_pem", "node2_pem")] <- minmax(vec_pem)

    # subset to crosstalk edges
    df_lig <- df_net[df_net$is_ct_edge, ]

    if (nrow(df_lig) == 0) {
        # write out root network
        row <- c("ROOT", "ROOT")
        fpath <- write_neighborhood_gv(row, df_net, dir_out, FALSE)
        render_gv(fpath)

        # halt
        message("no ligands found in the network!")
        return()
    }

    # loop through each
    apply(df_lig, 1, function(row) {
        nodes_all <- unlist(df_net[, c(1, 2)])
        nodes_mat <- matrix(nodes_all, ncol = 2)

        nodes_sel <- c()
        nodes_new <- unlist(row[c(1, 2)])

        for (d in seq_len(depth)) {
            nodes_sel <- unique(c(nodes_sel, nodes_new))
            index <- rowSums(matrix(nodes_all %in% nodes_sel, ncol = 2)) != 0
            nodes_new <- as.vector(nodes_mat[index, ])
        }

        df_net_sub <- df_net[index, ]
        fpath <- write_neighborhood_gv(row, df_net_sub, dir_out)
        render_gv(fpath)
    })
    NULL
}


#' Generate Network Pathways
#'
#' Starting at each ligand-receptor pair in the final network, generate small
#' neighbor graphs around these pairs. Consider them to be subsets of the larger
#' graph.
#'
#' @examples \dontrun{
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_out <- "~/CytoTalk-output"
#' depth <- 3
#'
#' generate_pathways(type_a, type_b, dir_out, depth)
#' }
#'
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_out Output directory
#' @param depth How many steps out to form neighborhood?
#' @return None
#' @export
generate_pathways <- function(type_a, type_b, dir_out, depth) {
    extract_best_network(type_a, type_b, dir_out)
    write_network_sif(dir_out)
    write_pathways_gv(dir_out, depth)
    NULL
}

