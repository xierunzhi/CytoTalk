#' Run Prize-Collecting Steiner Tree Algorithm
#'
#' Much thanks to The Fraenkel Lab at MIT, for providing this wonderful
#' implementation of the PCST algorithm. You can view their repository [here](https://github.com/fraenkel-lab/pcst_fast).
#'
#' Given the integrated network, run a factorial design over beta and omega
#' values, generating pruned networks for each run. These networks will be
#' compared to each other later.
#'
#' @param dir_out Output directory
#' @param beta_max End point of the range 1-`beta_max`, default 100.
#' @param omega_min Start point of the range `omega_min`-`omega_max`, default 0.5.
#' @param omega_max Start point of the range `omega_min`-`omega_max`, default 0.5. Could increase to 1.5.
#' @return NULL
#' @export
run_pcst <- function(dir_out, beta_max, omega_min, omega_max) {
    # format filepaths
    fpath_net <- file.path(dir_out, "IntegratedNetwork.cfg")
    fpath_node <- file.path(dir_out, "PCSF_NodeOccurance.txt")
    fpath_edge <- file.path(dir_out, "PCSF_EdgeOccurance.txt")

    # stop if input file doesn't exist,
    # skip if output already generated
    if (!file.exists(fpath_net)) {
        stop(sprintf("cannot find input file: %s", fpath_net))
    } else if (all(file.exists(fpath_node, fpath_edge))) {
        return()
    }

    # read in network file
    lines <- strsplit(readLines(fpath_net), "\t")
    lengths <- sapply(lines, length)

    # form data frames
    df_nodes <- as.data.frame(do.call(rbind, lines[lengths == 2]))
    df_edges <- as.data.frame(do.call(rbind, lines[lengths == 3]))

    # type conversion
    df_nodes <- utils::type.convert(df_nodes, as.is = TRUE)
    df_edges <- utils::type.convert(df_edges, as.is = TRUE)

    # names
    names(df_nodes) <- c("node", "prize")
    names(df_edges) <- c("node1", "node2", "cost")

    # generate template network
    df_nodes_tmp <- data.frame(node = "ARTI", prize = 1)
    df_edges_tmp <- data.frame(node1 = "ARTI", node2 = df_nodes$node, cost = NA)

    # bind to old network
    df_nodes_tmp <- rbind(df_nodes_tmp, df_nodes)
    df_edges_tmp <- rbind(df_edges, df_edges_tmp)

    # find index-0 node indices
    edge_inds <- apply(df_edges_tmp[, 1:2], 2, function(col) {
        as.integer(match(col, df_nodes_tmp$node) - 1)
    })

    # reticulate import (could fail)
    pcst_fast <- reticulate::import("pcst_fast")

    # for each beta
    lst_beta <- list()
    for (beta in 1:beta_max) {
        # copy of template network
        df_nodes_spc <- df_nodes_tmp
        df_edges_spc <- df_edges_tmp
        index <- is.na(df_edges_tmp$cost)

        # beta as prize multiplier
        df_nodes_spc$prize <- (df_nodes_spc$prize * beta)

        lst_omega <- list()
        for (omega in seq(omega_min, omega_max, 0.1)) {
            # omega as ARTI edge cost
            df_edges_spc$cost[index] <- omega

            lst <- pcst_fast$pcst_fast(
                edge_inds,
                df_nodes_spc$prize,
                df_edges_spc$cost,
                as.integer(0),
                as.integer(1),
                "strong",
                as.integer(0)
            )

            lst$beta <- beta
            lst$omega <- round(omega, 1)

            lst_omega <- append(lst_omega, list(lst))
        }
        lst_beta <- append(lst_beta, lst_omega)
    }

    # collect the resulting lists
    df_nodes_all <- NULL
    df_edges_all <- NULL
    for (lst in lst_beta) {
        # grab run parameters
        beta <- lst$beta
        omega <- lst$omega

        # find node and edge names
        node_names <- df_nodes_tmp[lst[[1]] + 1, "node", drop=FALSE]
        edge_names <- df_edges_tmp[lst[[2]] + 1, c("node1", "node2")]

        # if there is data to add, then add it
        if (nrow(node_names) != 0) {
            df_nodes_new <- data.frame(beta, omega, node_names)
            df_nodes_all <- rbind(df_nodes_all, df_nodes_new)
        }
        if (nrow(edge_names) != 0) {
            df_edges_new <- data.frame(beta, omega, edge_names)
            df_edges_all <- rbind(df_edges_all, df_edges_new)
        }
    }

    # write out
    vroom::vroom_write(df_nodes_all, fpath_node, progress = FALSE)
    vroom::vroom_write(df_edges_all, fpath_edge, progress = FALSE)
    NULL
}
