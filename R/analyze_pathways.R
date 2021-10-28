#' @noRd
score_subnetwork <- function(df_node, df_edge, beta) {
    node_prize <- as.numeric(df_node[, 2])
    edge_cost <- as.numeric(df_edge[, 3])
    c(
        mean(node_prize), mean(edge_cost),
        beta * sum(node_prize) - sum(edge_cost)
    )
}

score_subnetwork_shuffle <- function(node_prize, edge_cost, beta) {
    c(
        mean(node_prize),
        mean(edge_cost),
        beta * sum(node_prize) - sum(edge_cost)
    )
}

#' @noRd
subsample_network_simple <- function(df_node, df_edge, n_edges) {
    index_edge <- sample(seq_len(nrow(df_edge)), n_edges)
    index_node <- match(unique(unlist(df_edge[index_edge, -3])), df_node[, 1])

    list(
        nodes = df_node[index_node, ],
        edges = df_edge[index_edge, ]
    )
}

#' @noRd
subsample_network_shuffle <- function(node_prize, edge_cost, n_nodes, n_edges) {
    list(
        node_prize = as.numeric(sample(node_prize, n_nodes)),
        edge_cost = as.numeric(sample(edge_cost, n_edges))
    )
}

#' @noRd
subsample_network <- function(df_node, df_edge, ct_edge, n_edges) {
    edges_all <- df_edge[, c(1, 2)]
    
    # being selection
    edges_sel <- new <- sample(ct_edge, 1)
    nodes_sel <- unlist(edges_all[new, ])

    for (i in 2:n_edges) {
        nodes_found <- (
            edges_all[, 1] %in% nodes_sel +
            edges_all[, 2] %in% nodes_sel
        )

        new <- sample(which(nodes_found == 1), 1)
        edges_sel <- c(edges_sel, new)
        nodes_sel <- unique(c(nodes_sel, unlist(edges_all[new, ])))
    }

    list(
        nodes = df_node[match(nodes_sel, df_node[, 1]), ],
        edges = df_edge[edges_sel, ]
    )
}

#' Analyze Network Pathways
#'
#' Starting at each ligand-receptor pair in the final network, analyze small
#' neighbor graphs around these pairs. Consider them to be subsets of the larger
#' graph.
#'
#' @examples {
#' type_a <- "BCells"
#' type_b <- "TCells"
#' dir_out <- "~/CytoTalk-output"
#' depth <- 3
#' ntrial <- 1000
#'
#' analyze_pathways(type_a, type_b, dir_out, depth, ntrial)
#' }
#'
#' @param type_a Cell type A
#' @param type_b Cell type B
#' @param dir_out Output directory
#' @param depth How many steps out to form neighborhood?
#' @param ntrial How many empirical simulations to run?
#' @return None
#' @export
analyze_pathways <- function(type_a, type_b, dir_out, depth, ntrial) {
    # analysis folder
    dir_out_ana <- file.path(dir_out, "analysis")
    dir_out_ptw <- file.path(dir_out, "pathways")

    # make sure it exists
    if (!dir.exists(dir_out_ana)) {
        dir.create(dir_out_ana)
    }

    # format filepaths
    fpath_net <- file.path(dir_out, "IntegratedNetwork.cfg")
    fpath_out <- file.path(dir_out_ana, "PathwayScores.txt")
    fpath_pval <- file.path(dir_out, "PCSF_EdgeTestValues.txt")
    fpath_net <- file.path(dir_out, "PCSF_Network.txt")
    fnames_sub <- dir(dir_out_ptw)

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!all(file.exists(fpath_net, fpath_pval))) {
        stop(sprintf("cannot find input file: %s", fpath_net))
    } else if (file.exists(fpath_out)) {
        message(sprintf("file already exists, continuing: %s", fpath_out))
        return(0)
    }

    # process network cfg
    lst_net <- strsplit(readLines(fpath_net), "\t")
    lengths <- vapply(lst_net, length, integer(1))
    df_net_nodes <- as.data.frame(do.call(rbind, lst_net[lengths == 2]))
    df_net_edges <- as.data.frame(do.call(rbind, lst_net[lengths == 3]))

    # load in file
    df_pval <- suppressMessages(vroom::vroom(fpath_pval, progress = FALSE))
    beta <- df_pval[order(df_pval$pval)[1], "beta", drop=TRUE]
    df_net <- suppressMessages(vroom::vroom(
        fpath_net, progress = FALSE, na = c("", "NA", "-Inf")
    ))

    # create node table
    df_net_nodes <- data.frame(rbind(
        as.matrix(df_net[, c("node1", "node1_prize")]),
        as.matrix(df_net[, c("node2", "node2_prize")])
    ))

    # create edge table
    df_net_edges <- as.data.frame(df_net[, c("node1", "node2", "cost")])

    # set column names
    names(df_net_nodes) <- c("node", "prize")
    names(df_net_edges) <- c("node1", "node2", "cost")

    # set columns numeric
    df_net_nodes$prize <- as.numeric(df_net_nodes$prize)
    df_net_edges$cost <- as.numeric(df_net_edges$cost)

    # start at ct edge
    ct_edge <- apply(df_net_edges[, c(1, 2)], 2, endsWith, "_")
    ct_edge <- which(rowSums(ct_edge) == 1)

    df_out <- NULL
    for (fname in fnames_sub) {
        fpath <- file.path(dir_out_ptw, fname)
        df_net_sub <- suppressMessages(vroom::vroom(fpath, progress = FALSE))
        df_net_sub <- data.frame(df_net_sub)
    
        # extract information on LR pair
        lrp <- strsplit(gsub("(.+?)_(.+?)\\.txt", "\\1 \\2", fname), " ")[[1]]
        edges_plain <- apply(df_net_sub[, c(1, 2)], 2, function(x) gsub("_$", "", x))
        index <- colSums(apply(edges_plain, 1, "==", lrp)) == 2
        lrp_info <- df_net_sub[index, ]
        
        # prepare nodes
        df_node <- data.frame(rbind(
            as.matrix(df_net_sub[, c("node1", "node1_prize")]),
            as.matrix(df_net_sub[, c("node2", "node2_prize")])
        ))
        df_node <- df_node[!duplicated(df_node[, 1]), ]
        n_nodes <- nrow(df_node)

        # prepare edges
        df_edge <- df_net_sub[, c("node1", "node2", "cost")]
        n_edges <- nrow(df_edge)

        # begin scores
        scores <- score_subnetwork(df_node, df_edge, beta)

        # simulate random subsets
        for (i in seq_len(ntrial)) {
            lst <- subsample_network_shuffle(
                df_net_nodes[, 2], df_net_edges[, 3], n_nodes, n_edges
            )
            score <- score_subnetwork_shuffle(lst$node_prize, lst$edge_cost, beta)
            scores <- rbind(scores, score)
        }

        # calculate p-value
        pscores <- apply(scores, 2, function(x) {
            (match(x[1], sort(x)) - 1) / (length(x) - 1)
        })
    
        # extract node score
        pprize <- 1 - pscores[1]
        # extract edge score
        pcost <- pscores[2]
        # extract potential score
        ppot <- 1 - pscores[3]

        # new row of data
        row <- data.frame(
            ligand = lrp[1],
            ligand_cell_type = lrp_info$node1_type,
            receptor = lrp[2],
            receptor_cell_type = lrp_info$node2_type,
            num_edges = n_edges,
            num_nodes = nrow(df_node),
            mean_prize = scores[1, 1],
            mean_cost = scores[1, 2],
            potential = scores[1, 3],
            pval_prize = pprize,
            pval_cost = pcost,
            pval_potential = ppot
        )

        # start or continue dataframe
        if (is.null(df_out)) {
            df_out <- row
        } else {
            df_out <- rbind(df_out, row)
        }
    }

    # order by edge pval
    df_out <- data.frame(df_out)
    if (0 < nrow(df_out)) {
        df_out <- df_out[order(df_out$pval_potential), ]
    }

    # add adjustments
    df_out$pval_prize_adj <- stats::p.adjust(df_out$pval_prize)
    df_out$pval_cost_adj <- stats::p.adjust(df_out$pval_cost)
    df_out$pval_potential_adj <- stats::p.adjust(df_out$pval_potential)

    # write out
    vroom::vroom_write(df_out, fpath_out, progress = FALSE)
    NULL
}
