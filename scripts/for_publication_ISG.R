# Load libraries
source("src/libraries.R")


# Load commonly used function
source("src/common_functions.R")



load(file = "data/processed_data/physeq_WGS_2022_10_20_Baseline.RData")

# Set metadata variable
df_metadata <- data.frame(sample_data(physeq_WGS_2022_10_20_Baseline))


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Alpha Diversity
# ────────────────────────────────────────────────────────────────────────────────────────────────────


df_isg <- estimate_richness(physeq_WGS_2022_10_20_Baseline) %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
    inner_join(df_metadata, by = c("sample_id")) %>%
    dplyr::select(estimate_isg, InvSimpson, Observed, Shannon) %>%
    pivot_longer(-estimate_isg) %>%
    filter(!is.na(estimate_isg))



df_stat_isg <- df_isg %>%
    mutate(name = factor(name, levels = c("Observed", "InvSimpson", "Shannon"))) %>%
    group_by(name) %>%
    wilcox_test(value ~ estimate_isg) %>%
    add_y_position(scales = "free_y")



fig_alpha_isg <- df_isg %>%
    mutate(name = factor(name, levels = c("Observed", "InvSimpson", "Shannon"))) %>%
    mutate(estimate_isg = factor(estimate_isg, levels = c("Low", "High"))) %T>%
    write_csv("output/figures/alpha_diversity/By_estimate_isg_score.csv") %>%
    ggplot(aes(x = estimate_isg, y = value, fill = estimate_isg)) +
    geom_boxplot(alpha = 1, width = 0.5, fatten = 4) +
    geom_point(size = 4) +
    geom_point(size = 5, shape = 21, color = "white", stroke = 2) +
    theme(aspect.ratio = 1) +
    facet_wrap(~name, scales = "free_y") +
    theme_pubclean() +
    scale_fill_manual(values = c("High" = "blue", "Low" = "red")) +
    add_pvalue(data = df_stat_isg, label = "Wilcoxon p = {p}", inherit.aes = F, label.size = 6) +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") +
    labs(fill = "Estimate ISG score") +
    labs(
        caption = glue::glue(
            "{df_stat_isg$group1[[1]]}:{df_stat_isg$n1[[1]]} \n {df_stat_isg$group2[[1]]}:{df_stat_isg$n2[[1]]}"
        )
    ) +
    stat_summary(
        fun.y = median,
        geom = "crossbar",
        color = "white",
        linewidth = 0.6,
        width = 1
    )
save(fig_alpha_isg, file = "output/figures/alpha_diversity/By_estimate_isg_score.RData")
ggsave("output/figures/alpha_diversity/By_Whole_group_estimate_isg.pdf", width = 12, height = 8)




# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Top 20 taxa
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Set palette
myPal <- tax_palette(physeq_WGS_2022_10_20_Baseline, rank = "Species", pal = "brewerPlus", n = 30)



fig_isg_tax <- physeq_WGS_2022_10_20_Baseline %>%
    speedyseq::filter_sample_data(!is.na(estimate_isg)) %>%
    speedyseq::mutate_sample_data(estimate_isg = factor(estimate_isg, levels = c("Low", "High"))) %>%
    comp_barplot(
        tax_level = "Species", n_taxa = 20,
        sample_order = "asis",
        palette = myPal, label = "pt_id"
    ) +
    facet_wrap(~estimate_isg, scales = "free_x") +
    # ggtitle("Overall Estimate ISG score category") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(text = element_text(size = 16)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save(fig_isg_tax, file = "output/figures/top_n_taxa/top_20_species_by_estimate_isg_score.RData")
ggsave("output/figures/top_n_taxa/top_20_species_by_ISG_score_category_mOS.pdf", height = 8, width = 13)

fig_isg_tax$data %>%
    dplyr::select(pt_id, estimate_isg, Species, Abundance) %>%
    write_csv("output/figures/top_n_taxa/top_20_species_by_isg.csv")



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Beta Diversity
# ────────────────────────────────────────────────────────────────────────────────────────────────────


fig_beta_isg <- physeq_WGS_2022_10_20_Baseline %>%
    speedyseq::filter_sample_data(!is.na(estimate_isg)) %>%
    plot_beta_diversity_with_adonis2(., sample_param = "estimate_isg") +
    scale_color_manual(values = c("High" = "blue", "Low" = "red")) +
    labs(color = "Estimate ISG") +
    geom_point(size = 4, shape = 21, stroke = 2) +
    theme(text = element_text(size = 20))
save(fig_beta_isg, file = "output/figures/beta_diversity/by_estimate_isg_score.RData")
ggsave("output/figures/beta_diversity/by_estimate_isg_score.pdf", width = 12, height = 7)



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Differential analysis
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# By estimate ISG score ----
physeq_WGS_2022_10_20_Baseline_glommed <- speedyseq::tax_glom(
    physeq_WGS_2022_10_20_Baseline,
    taxrank = "Species"
) %>%
    speedyseq::filter_sample_data(estimate_isg %in% c("Low", "High")) %>%
    speedyseq::mutate_sample_data(estimate_isg = factor(estimate_isg, levels = c("Low", "High")))


tax_table_asv <- data.frame(tax_table(physeq_WGS_2022_10_20_Baseline_glommed)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(1, Species)


df_prev <- get_phyloseq_prev_by_group(
    physeq_WGS_2022_10_20_Baseline_glommed,
    sample_param = c("estimate_isg"),
    taxa_rank = "Species"
) %>%
    ungroup()



df_metadata_N <- data.frame(sample_data(physeq_WGS_2022_10_20_Baseline_glommed)) %>%
    group_by(estimate_isg) %>%
    tally()


ancombc_bc_by_isg <- ancombc(
    phyloseq = physeq_WGS_2022_10_20_Baseline_glommed,
    formula = "estimate_isg", p_adj_method = "fdr",
    group = "estimate_isg", struc_zero = T, neg_lb = TRUE, tol = 1e-5,
    prv_cut = 0.50,
    max_iter = 100, conserve = T, alpha = 0.05, global = FALSE
)


df_ancombc_by_ISG <- get_df_ancombc(var = c("estimate_isg"), ancombc_res = ancombc_bc_by_isg) %>%
    inner_join(tax_table_asv, by = c("taxa_name" = "tax_id")) %>%
    mutate(estimate_isg_PVAL = case_when(
        estimate_isg_PVAL == 0 ~ 0.0001,
        TRUE ~ as.numeric(estimate_isg_PVAL)
    )) %>%
    mutate(
        ymin = estimate_isg_LFC - estimate_isg_SE,
        ymax = estimate_isg_LFC + estimate_isg_SE
    ) %>%
    # dplyr::filter(sign(ymin) == sign(ymax)) %>%
    # dplyr::filter(estimate_isg_QVAL < 0.05) %>%
    mutate(group = case_when(
        estimate_isg_LFC > 0 ~ "High",
        estimate_isg_LFC < 0 ~ "Low"
    )) %>%
    mutate(signif = case_when(
        estimate_isg_DIFF_ABN == FALSE ~ "no",
        estimate_isg_DIFF_ABN == TRUE ~ "yes"
    ))



# Alternative plot ----
fig_diff_isg <- df_ancombc_by_ISG %>%
    inner_join(df_metadata_N, by = c("group" = "estimate_isg")) %>%
    inner_join(df_prev, by = c("Species" = "Species", "group" = "estimate_isg")) %T>%
    write_csv("output/figures/ancombc/Alternative_plot_by_estimate_isg.csv") %>%
    filter(estimate_isg_QVAL < 0.05) %>%
    filter(estimate_isg_LFC > 2 | estimate_isg_LFC < -2) %>%
    ggplot(., aes(x = estimate_isg_LFC, y = reorder(Species, -estimate_isg_LFC), color = group)) +
    geom_point(aes(size = prevalence)) +
    theme_bw() +
    annotate(geom = "rect", xmin = 0, xmax = -Inf, ymin = 0, ymax = Inf, alpha = 0.09, fill = "red") +
    annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.09, fill = "blue") +
    geom_segment(aes(x = 0, xend = estimate_isg_LFC, y = Species, yend = Species, color = group), size = 2) +
    geom_text_repel(
        aes(label = ifelse(
            estimate_isg_QVAL < 0.05 & abs(estimate_isg_LFC) > 2 & sign(ymin) == sign(ymax),
            Species,
            ""
        )),
        max.overlaps = 30, nudge_x = 0.2, nudge_y = -0.1, size = 7, show.legend = FALSE
    ) +
    scale_color_manual(values = c("blue", "red")) +
    theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
    ) +
    theme(text = element_text(size = 20)) +
    geom_vline(xintercept = 0, size = 1.5, linetype = "dashed") +
    xlab("Log fold change") +
    ylab("Species") +
    labs(caption = "qval < 0.05")
save(fig_diff_isg, file = "output/figures/ancombc/Alternative_plot_by_estimate_isg.RData")
ggsave("output/figures/ancombc/Alternative_plot_by_estimate_isg.pdf", width = 14, height = 7, dpi = 300)



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Assembling the final figure for publication
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Design layout for the combined figures
design <- "
12
34
"

(fig_alpha_isg + theme(legend.position = "none")) + fig_isg_tax + (fig_beta_isg + theme(legend.position = "none")) + fig_diff_isg + plot_layout(design = design) +
    plot_layout(widths = c(2, 2)) +
    plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") +
    plot_annotation(title = "By ISG score") &
    theme(title = element_text(size = 24, face = "bold"))
ggsave("output/figures/figure_ISG.pdf", width = 30, height = 20, scale = 0.8)
