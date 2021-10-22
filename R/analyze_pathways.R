#' @noRd
subsample_network <- function(df_node, df_edge, n_edges) {
    edges_all <- df_edge[, c(1, 2)]
    nodes_all <- unlist(edges_all)

    edges_sel <- c()
    edges_new <- sample(seq_len(nrow(edges_all)), 1)

    for (i in 2:n_edges) {
        edges_sel <- c(edges_sel, edges_new)
        nodes_sel <- unique(unlist(edges_all[edges_sel, ]))
        nodes_found <- matrix(nodes_all %in% nodes_sel, ncol = 2)
        edges_new <- sample(which(rowSums(nodes_found) == 1), 1)
    }

    # final update
    edges_sel <- c(edges_sel, edges_new)
    nodes_sel <- unique(unlist(edges_all[edges_sel, ]))

    list(
        nodes = df_node[match(nodes_sel, df_node[, 1]), ],
        edges = df_edge[edges_sel, ]
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
score_subnetwork <- function(df_node, df_edge) {
    c(
        mean(as.numeric(df_node[, 2])),
        mean(as.numeric(df_edge[, 3]))
    )
}

#' Analyze Network Pathways
#'
#' Starting at each ligand-receptor pair in the final network, analyze small
#' neighbor graphs around these pairs. Consider them to be subsets of the larger
#' graph.
#'
#' @examples {
#' dir_out <- "~/CytoTalk-output"
#' depth <- 3
#'
#' analyze_pathways(dir_out, depth)
#' }
#'
#' @param dir_out Output directory
#' @param depth How many steps out to form neighborhood?
#' @return None
#' @export
analyze_pathways <- function(type_a, type_b, dir_out, depth, ntrial) {
    # analysis folder
    dir_out_ana <- file.path(dir_out, "analysis")

    # format filepaths
    fpath_net <- file.path(dir_out, "IntegratedNetwork.cfg")
    fpath_out <- file.path(dir_out_ana, "PathwayScores.txt")
    fnames_sub <- dir(dir_out_ana)

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_net)) {
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

    # set column names
    names(df_net_nodes) <- c("node", "prize")
    names(df_net_edges) <- c("node1", "node2", "cost")

    df_out <- NULL
    for (fname in fnames_sub) {
        fpath <- file.path(dir_out_ana, fname)
        df_net_sub <- read.delim(fpath)
        n_edges <- nrow(df_net_sub)
        
        # prepare nodes
        df_node <- data.frame(rbind(
            as.matrix(df_net_sub[, c("node1", "node1_prize")]),
            as.matrix(df_net_sub[, c("node2", "node2_prize")])
        ))
        df_node <- df_node[!duplicated(df_node[, 1]), ]

        # prepare edges
        df_edge <- df_net_sub[, c("node1", "node2", "cost")]

        # begin scores
        scores <- score_subnetwork(df_node, df_edge)

        # simulate random subsets
        for (i in seq_len(ntrial)) {  
            lst <- subsample_network_simple(df_net_nodes, df_net_edges, n_edges)
            score <- score_subnetwork(lst$nodes, lst$edges)
            scores <- rbind(scores, score)
        }

        # calculate Z score
        zscore <- apply(scores, 2, scale)
        # calculate p-value
        pscore <- apply(zscore, 2, pnorm)

        # extract node score
        pprize <- 1 - pscore[1, 1]
        # extract edge score
        pcost <- pscore[1, 2]

        # new row of data
        row <- data.frame(
            pathway = gsub("\\..+", "", fname),
            num_edges = n_edges, num_nodes = nrow(df_node),
            pval_prize = pprize, pval_cost = pcost
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
    df_out <- df_out[order(df_out$pval_cost), ]

    # write out
    write.csv(df_out, fpath_out, row.names = FALSE, quote = FALSE)
    NULL
}
