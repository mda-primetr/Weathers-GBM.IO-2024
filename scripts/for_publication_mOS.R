# Load libraries
source("src/libraries.R")


# Load commonly used function
source("src/common_functions.R")


# Load the phyloseq object
load(file = "data/processed_data/physeq_WGS_2022_10_20_Baseline.RData")

# Set metadata variable
df_metadata <- data.frame(sample_data(physeq_WGS_2022_10_20_Baseline))


# Remove IDH mutant samples ----
physeq_WGS_2022_10_20_Baseline_mod <- physeq_WGS_2022_10_20_Baseline %>%
    ps_filter(!grepl("m", pt_id))


# Figures by OS ----
# By Whole group OS 568 days ----

# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Alpha Diversity
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# By Median OS 568 days cut-off  ----

df_alpha_final_phyloseq <- estimate_richness(physeq_WGS_2022_10_20_Baseline_mod) %>%
    rownames_to_column(var = "sample_id") %>%
    mutate(sample_id = gsub("\\.", "-", sample_id)) %>%
    inner_join(df_metadata, by = c("sample_id")) %>%
    mutate(mOS_cat_whole = factor(mOS_cat_whole, levels = c("Below Median", "Above Median")))



df_stat_mOS <- df_alpha_final_phyloseq %>%
    dplyr::select(mOS_cat_whole, InvSimpson, Observed, Shannon) %>%
    pivot_longer(-mOS_cat_whole) %>%
    mutate(name = factor(name, levels = c("Observed", "InvSimpson", "Shannon"))) %>%
    group_by(name) %>%
    wilcox_test(value ~ mOS_cat_whole) %>%
    add_y_position(scales = "free_y")



fig_alpha_mOS <- df_alpha_final_phyloseq %>%
    dplyr::select(mOS_cat_whole, InvSimpson, Observed, Shannon) %>%
    pivot_longer(-mOS_cat_whole) %>%
    mutate(name = factor(name, levels = c("Observed", "InvSimpson", "Shannon"))) %T>%
    write_csv("output/figures/alpha_diversity/By_Whole_group_median_OS_Baseline.csv") %>%
    ggplot(aes(x = mOS_cat_whole, y = value, fill = mOS_cat_whole)) +
    geom_boxplot(alpha = 1, width = 0.5, fatten = 4) +
    geom_point(size = 4) +
    geom_point(size = 5, shape = 21, color = "white", stroke = 2) +
    scale_fill_manual(values = c("Above Median" = "blue", "Below Median" = "red")) +
    theme(aspect.ratio = 1) +
    facet_wrap(~name, scales = "free_y") +
    theme_pubclean() +
    add_pvalue(data = df_stat_mOS, label = "Wilcoxon p = {p}", inherit.aes = F, label.size = 6) +
    theme(text = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    xlab("") +
    labs(fill = "Overall Survival Median OS - 568 days") +
    labs(
        caption = glue::glue(
            "{df_stat_mOS$group1[[1]]}:{df_stat_mOS$n1[[1]]} \n {df_stat_mOS$group2[[1]]}:{df_stat_mOS$n2[[1]]}"
        )
    ) +
    stat_summary(
        fun.y = median,
        geom = "crossbar",
        color = "white",
        linewidth = 0.6,
        width = 1
    )
save(fig_alpha_mOS, file = "output/figures/alpha_diversity/By_Whole_group_median_OS_Baseline.RData")
ggsave("output/figures/alpha_diversity/By_Whole_group_median_OS_Baseline.pdf", width = 12, height = 8)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Top Taxa
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# Set palette
myPal <- tax_palette(physeq_WGS_2022_10_20_Baseline_mod, rank = "Species", pal = "brewerPlus", n = 30)



# Top taxa plot ---
fig_med_tax <- physeq_WGS_2022_10_20_Baseline_mod %>%
    speedyseq::mutate_sample_data(mOS_cat_whole = factor(mOS_cat_whole, levels = c("Below Median", "Above Median"))) %>%
    comp_barplot(
        tax_level = "Species", n_taxa = 20,
        sample_order = "asis",
        palette = myPal, label = "pt_id"
    ) +
    facet_wrap(~mOS_cat_whole, scales = "free_x") +
    guides(fill = guide_legend(ncol = 1)) +
    theme(text = element_text(size = 16)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
save(fig_med_tax, file = "output/figures/top_n_taxa/top_20_species_by_mOS.RData")

# Save data for supplementary table
fig_med_tax$data %>%
    dplyr::select(pt_id, mOS_cat_whole, Species, Abundance) %>%
    write_csv("output/figures/top_n_taxa/top_20_species_by_mOS.csv")


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Beta Diversity
# ────────────────────────────────────────────────────────────────────────────────────────────────────

# By whole group median OS ----
fig_beta_mOS <- plot_beta_diversity_with_adonis2(physeq_WGS_2022_10_20_Baseline_mod, sample_param = "mOS_cat_whole") +
    scale_color_manual(values = c("Above Median" = "blue", "Below Median" = "red")) +
    labs(color = "OS Category") +
    geom_point(size = 4, shape = 21, stroke = 2) +
    theme(text = element_text(size = 20))
save(fig_beta_mOS, file = "output/figures/beta_diversity/by_Median_OS_whole_Group.RData")
ggsave("output/figures/beta_diversity/by_Median_OS_whole_Group.pdf", width = 12, height = 7)



# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Differential Abundance using ANCOMBC 1.6.2
# ────────────────────────────────────────────────────────────────────────────────────────────────────


# Differential abundance -----
# Run ANCOMBC

physeq_WGS_2022_10_20_Baseline_mod_glommed <- speedyseq::tax_glom(physeq_WGS_2022_10_20_Baseline_mod, taxrank = "Species") %>%
    speedyseq::mutate_sample_data(mOS_cat_whole = factor(
        mOS_cat_whole,
        levels = c("Below Median", "Above Median")
    ))


# Get taxa table
tax_table_asv <- data.frame(tax_table(physeq_WGS_2022_10_20_Baseline_mod_glommed)) %>%
    rownames_to_column(var = "tax_id") %>%
    dplyr::select(1, Species)


# Get prevalence
df_prev <- get_phyloseq_prev_by_group(
    physeq = physeq_WGS_2022_10_20_Baseline_mod_glommed,
    sample_param = "mOS_cat_whole", taxa_rank = "Species"
)



df_metadata_N <- data.frame(sample_data(physeq_WGS_2022_10_20_Baseline_mod_glommed)) %>%
    group_by(estimate_isg) %>%
    tally()




ancombc_bc_by_mOS <- ancombc(
    phyloseq = physeq_WGS_2022_10_20_Baseline_mod_glommed,
    formula = "mOS_cat_whole", p_adj_method = "fdr",
    group = "mOS_cat_whole", struc_zero = T, neg_lb = TRUE, tol = 1e-5,
    prv_cut = 0.30,
    max_iter = 100, conserve = T, alpha = 0.05, global = FALSE
)


fig_diff_med_os <- get_df_ancombc(var = c("mOS_cat_whole"), ancombc_res = ancombc_bc_by_mOS) %>%
    inner_join(tax_table_asv, by = c("taxa_name" = "tax_id")) %>%
    mutate(mOS_cat_whole_PVAL = case_when(
        mOS_cat_whole_PVAL == 0 ~ 0.0001,
        TRUE ~ as.numeric(mOS_cat_whole_PVAL)
    )) %>%
    mutate(
        ymin = mOS_cat_whole_LFC - mOS_cat_whole_SE,
        ymax = mOS_cat_whole_LFC + mOS_cat_whole_SE
    ) %>%
    mutate(group = case_when(
        mOS_cat_whole_LFC > 0 ~ "Above Median",
        mOS_cat_whole_LFC < 0 ~ "Below Median"
    )) %T>%
    write_csv("output/figures/ancombc/Alternative_plot_by_mOS_cat_whole_568_days.csv") %>%
    filter(mOS_cat_whole_QVAL < 0.1) %>%
    filter(abs(mOS_cat_whole_LFC) > 2) %>%
    mutate(Species = gsub("s__", "", Species)) %>%
    inner_join(df_prev, by = c("Species" = "Species", "group" = "mOS_cat_whole")) %>%
    filter(mOS_cat_whole_QVAL < 0.05) %>%
    ggplot(., aes(x = mOS_cat_whole_LFC, y = reorder(Species, -mOS_cat_whole_LFC), color = group)) +
    geom_point(aes(size = prevalence)) +
    theme_bw() +
    annotate(geom = "rect", xmin = 0, xmax = -Inf, ymin = 0, ymax = Inf, alpha = 0.09, fill = "red") +
    annotate(geom = "rect", xmin = 0, xmax = Inf, ymin = 0, ymax = Inf, alpha = 0.09, fill = "blue") +
    geom_segment(aes(x = 0, xend = mOS_cat_whole_LFC, y = Species, yend = Species, color = group), size = 2) +
    geom_text_repel(
        aes(label = ifelse(
            mOS_cat_whole_QVAL < 0.1 & abs(mOS_cat_whole_LFC) > 2 & sign(ymin) == sign(ymax),
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
save(fig_diff_med_os, file = "output/figures/ancombc/Alternative_plot_by_mOS_cat_whole_568_days_alt.RData")
ggsave("output/figures/ancombc/Alternative_plot_by_mOS_cat_whole_568_days_alt.pdf", width = 15, height = 7)


# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Assemble figure for the final figure
# ────────────────────────────────────────────────────────────────────────────────────────────────────
# Assemble figure for Median OS ----
design <- "
12
34
"
(fig_alpha_mOS + theme(legend.position = "none")) + fig_med_tax + (fig_beta_mOS + theme(legend.position = "none")) + fig_diff_med_os + plot_layout(design = design) +
    plot_layout(widths = c(2, 2)) +
    plot_annotation(tag_levels = "A") + plot_layout(guides = "collect") +
    plot_annotation(title = "By Median OS") &
    theme(title = element_text(size = 24, face = "bold"))
ggsave("output/figures/figure_OS.pdf", width = 30, height = 20, scale = 0.8)
