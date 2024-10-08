library(speedyseq)
library(tidyverse)
library(vegan)
library(scales)
library(ggtext)


#' Function to calculate beta diversity stats
#' @param phyloseq_object phyloseq object
#' @param sample_param data parameter or column from the metadata
#' @return ggplot plot with adonis2 test p-value

plot_beta_diversity_with_adonis2 <- function(physeq, sample_param) {
    df_meta <- data.frame(sample_data(physeq)) %>%
        rownames_to_column(var = "sampleID")


    bray_dist <- phyloseq::distance(physeq, method = "bray")
    mod_sample_type <- betadisper(bray_dist, sample_data(physeq)[[sample_param]], type = "median")
    sequencer_ord <- ordinate(physeq, "PCoA", "bray")

    # Get Percent explained valued from the ordination
    percent_explained <- 100 * sequencer_ord$values$Eigenvalues / sum(sequencer_ord$values$Eigenvalues)

    # Adonis test using bray distance
    sam_df <- data.frame(sample_data(physeq))


    mod <- adonis2(as.formula(paste0("bray_dist", " ~ ", sample_param)), data = sam_df, permutations = 1000)


    # Create centroids dataframe
    df_centroids <- data.frame(mod_sample_type$centroids) %>%
        dplyr::select(PCoA1, PCoA2) %>%
        rownames_to_column(var = "sample_median")


    # Join the centroids dataframe with the metadata and the ordination dataframe
    df_points <- data.frame(sequencer_ord$vectors) %>%
        dplyr::select(Axis.1, Axis.2) %>%
        rownames_to_column(var = "sampleID") %>%
        inner_join(df_meta, by = c("sampleID")) %>%
        inner_join(df_centroids, by = setNames("sample_median", sample_param)) %>%
        tibble()

    #' sym function is needed so that ggplot2 can take the variable from function
    #' using !! in the ggplot functions to use that "sym" variable
    #'
    sample_param <- sym(sample_param)
    ggplot(df_points, aes(x = Axis.1, y = Axis.2, xend = PCoA1, yend = PCoA2, color = !!sample_param)) +
        geom_segment(alpha = 0.6, linewidth = 2) +
        geom_point(aes(color = !!sample_param), size = 2, alpha = 0.6) +
        geom_point(aes(x = PCoA1, y = PCoA2, color = !!sample_param, label = !!sample_param),
            shape = 12, size = 10, inherit.aes = FALSE, show.legend = FALSE, linewidth = 2
        ) +
        theme_bw() +
        stat_ellipse(type = "t", segments = 100, linewidth = 2) +
        labs(color = "sample_param") +
        xlab(glue::glue("PCoA1: {round(percent_explained[1], digits = 1)}%")) +
        ylab(glue::glue("PCoA2: {round(percent_explained[2], digits = 1)}%")) +
        geom_richtext(
            data = mod,
            size = 5,
            aes(
                label = glue::glue("Adonis2 test: p={round(mod$`Pr(>F)`[[1]],digits=4)}; R<sup>2</sup>={round(mod$R2[[1]],digits=4)}"),
                x = -Inf, y = Inf, vjust = 1.5, hjust = -1
            ), inherit.aes = FALSE
        )
}




# Function to rearrange x axis within facets

reorder_within <- function(x, by, within, fun = mean, sep = "___", ...) {
    new_x <- paste(x, within, sep = sep)
    stats::reorder(new_x, by, FUN = fun)
}

scale_x_reordered <- function(..., sep = "___") {
    reg <- paste0(sep, ".+$")
    ggplot2::scale_x_discrete(labels = function(x) gsub(reg, "", x), ...)
}







#' Function to get data frame for Top N taxa at specified level without any groupings
#'
#' @param physeq phyloseq object
#' @param no_of_taxa Integer to indicate the number of taxa for the grouping
#' @param taxa_level Taxa level to group the dataframe
#'
#' @return dataframe containing top taxa at specified taxa level


get_top_n_taxa_df <- function(physeq, no_of_taxa, taxa_level) {
    df_meta <- data.frame(sample_data(physeq)) %>%
        rownames_to_column(var = "SampleID")

    df_top_n_taxa <- speedyseq::psmelt(physeq) %>%
        group_by(Sample, {{ taxa_level }} := fct_lump(.data[[taxa_level]], n = no_of_taxa, w = Abundance)) %>%
        dplyr::summarize(count = sum(Abundance)) %>%
        mutate(prop = ((count) / sum(count)) * 100) %>%
        inner_join(df_meta, by = c("Sample" = "SampleID")) %>%
        ungroup()

    df_top_n_taxa
}


#' Function to get the top N taxa proportion by grouping
#' @param physeq phyloseq object
#' @param no_of_taxa Integer to indicate the number of taxa for the grouping
#' @param taxa_level Taxa level to group the dataframe
#' @param groups any one grouping variable in the phyloseq metadata
#'
#' @return dataframe containing top N taxa with proportion grouped by the specified variable

get_top_n_taxa_by_group_df <- function(physeq, no_of_taxa, taxa_level, groups) {
    # Glom the physeq to the taxa level
    physeq <- speedyseq::tax_glom(physeq, taxrank = taxa_level)

    # Get the meta data
    df_meta <- data.frame(sample_data(physeq)) %>%
        rownames_to_column(var = "SampleID")

    # Get the top overall taxa without grouping
    df_top_n_taxa <- speedyseq::psmelt(physeq) %>%
        group_by(Sample, {{ taxa_level }} := fct_lump(.data[[taxa_level]], n = no_of_taxa, w = Abundance)) %>%
        dplyr::summarize(count = sum(Abundance)) %>%
        mutate(prop = ((count) / sum(count)) * 100) %>%
        inner_join(df_meta, by = c("Sample" = "SampleID")) %>%
        ungroup()

    # Get the list of taxa to be filtered to be looked into grouping
    df_taxa_unique <- unique(df_top_n_taxa[[taxa_level]])

    # Get the final table with proportions
    df_taxa_by_group <- speedyseq::psmelt(physeq) %>%
        dplyr::select(Sample, {{ taxa_level }}, {{ groups }}, Abundance) %>%
        group_by(Sample, !!sym(groups), !!sym(taxa_level), .drop = FALSE) %>% # Order matters
        dplyr::summarize(cnt = sum(Abundance)) %>%
        mutate(prop = (cnt / sum(cnt)) * 100) %>%
        filter(!!sym(taxa_level) %in% df_taxa_unique)

    df_taxa_by_group
}




#' Function to get the Ancombc output dataframe
#' @param vars all the variables in order that were used for the AncomBC formula
#' @param ancombc_res AncomBC result object
#' Function to get the Ancombc output dataframe
#' @param vars all the variables in order that were used for the AncomBC formula
#' @param ancombc_res AncomBC result object
get_df_ancombc <- function(vars, ancombc_res) {
    df_lfc <- data.frame(ancombc_res$res$lfc)
    colnames(df_lfc) <- c(paste0(vars, "_LFC"))

    df_se <- data.frame(as.matrix(ancombc_res$res$se))
    colnames(df_se) <- c(paste0(vars, "_SE"))

    df_pval <- data.frame(as.matrix(ancombc_res$res$p_val))
    colnames(df_pval) <- c(paste0(vars, "_PVAL"))

    df_qval <- data.frame(as.matrix(ancombc_res$res$q_val))
    colnames(df_qval) <- c(paste0(vars, "_QVAL"))

    df_diff_abn <- data.frame(as.matrix(ancombc_res$res$diff_abn))
    colnames(df_diff_abn) <- c(paste0(vars, "_DIFF_ABN"))


    df_final_ancombc_formatted <- data.frame(
        taxa_name = rownames(df_lfc),
        df_lfc,
        df_se,
        df_pval,
        df_qval,
        df_diff_abn
    )


    df_final_ancombc_formatted
}





# Prevalence by sample type
#' Get prevalence of the taxa by sample parameters
#' @param physeq phyloseq object
#' @param sample_param column name or metadata column from the phyloseq sample_data
#' @param taxa_rank Taxa rank on which we want to calculate prevalance
get_phyloseq_prev_by_group <- function(physeq, sample_param, taxa_rank) {
    df_prev_sam <- sample_data(physeq) %>%
        group_by(.data[[sample_param]]) %>%
        dplyr::summarise(sam_sum = n())

    data.frame(speedyseq::psmelt(physeq)) %>%
        filter(Abundance > 0) %>%
        group_by(.data[[sample_param]], .data[[taxa_rank]]) %>%
        dplyr::summarise(counts = sum(Abundance > 0, na.rm = TRUE)) %>%
        inner_join(df_prev_sam, by = join_by(!!sample_param)) %>%
        mutate(prevalence = counts / sam_sum) %>%
        ungroup()
}
