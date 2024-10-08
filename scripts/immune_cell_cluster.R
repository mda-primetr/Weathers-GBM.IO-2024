# Load libraries
source("src/libraries.R")

# Load commonly used function
source("src/common_functions.R")


# Load the phyloseq object
load("data/processed_data/physeq_WGS_2022_10_20_Baseline.RData")
physeq_WGS_2022_10_20_Baseline


# Read the immune cluster data ----
df_immune_cluster <- read_csv("data/processed_data/df_immune_cluster.csv")


# Set metadata
df_metadata <- sample_data(physeq_WGS_2022_10_20_Baseline) %>%
    data.frame()




# Combine the immune cluster data with the physeq data ----
df_physeq_immune <- physeq_WGS_2022_10_20_Baseline %>%
    speedyseq::tax_glom(taxrank = "Species") %>%
    transform_sample_counts(function(x) (x / sum(x)) * 100) %>%
    psmelt() %>%
    dplyr::select(sample_id, pt_id, Species, Abundance, mOS) %>%
    inner_join(df_immune_cluster, by = c("pt_id"))




# list of species for selected for correlation ----
# These species were found to be differentially abundant for the ISG score analysis
lst_of_species <- c(
    "Bacteroides_thetaiotaomicron",
    "Eubacterium_siraeum",
    "Roseburia_faecis",
    "Escherichia_coli",
    "Clostridium_sp_CAG_58",
    "Eubacterium_eligens",
    "Dorea_formicigenerans",
    "Agathobaculum_butyriciproducens",
    "Roseburia_inulinivorans",
    "Roseburia_intestinalis"
)



for (i in lst_of_species) {
    df_cor_species <- df_physeq_immune %>%
        filter(Species == i) %>%
        dplyr::select(Abundance, mOS, "x_cell_deconvolution_t_cell_cd4_effector_memory":"lee_et_al_m_dc_ccr7") %>%
        pivot_longer(cols = -c(Abundance, mOS), names_to = "Category", values_to = "value")

    df_cor_species_stat <- df_cor_species %>%
        group_by(Category) %>%
        do(tidy(cor.test(.$Abundance, .$value, method = "spearman"))) %>%
        select(Category, estimate = estimate, p.value) %>%
        ungroup() %>%
        mutate(fdr_pval = p.adjust(p.value, method = "fdr")) %>%
        mutate(fdr_pval = format.pval(fdr_pval, digits = 2, eps = 0.001))

    df_cor_species %>%
        ggplot(., aes(x = Abundance, y = value)) +
        geom_point(aes(color = mOS), size = 4) +
        geom_smooth(method = "lm", se = FALSE, color = "black") +
        facet_wrap(~Category, scales = "free_y") +
        theme_bw(base_size = 19) +
        scale_color_manual(values = c("Below_Median_OS" = "red", "Above_Median_OS" = "blue")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_text(
            data = df_cor_species_stat,
            aes(label = paste("Spearman's ⍴ :", round(estimate, 2), "\n", "p.value:", fdr_pval)),
            x = Inf, y = Inf, hjust = 2, vjust = 1.5, size = 5, color = "black"
        ) +
        xlab("Abundance (%)") +
        ggtitle(i)
    ggsave(paste0("output/figures/immune_score_cluster/immune_score_component_correlation_clusters", i, ".pdf"),
        device = cairo_pdf,
        width = 22, height = 16
    )
}

# Check the proportion of the immune cell clusters by BM subtypes ----
for (i in lst_of_species) {
    df_cor_species <- df_physeq_immune %>%
        filter(Species == i) %>%
        dplyr::select(Abundance, mOS, gbm_subtypes_mesenchymal_proneural_classical) %>%
        pivot_longer(cols = -c(Abundance, mOS), names_to = "Category", values_to = "value") %>%
        dplyr::select(Abundance, value, Category, mOS)

    df_cor_stat <- df_cor_species %>%
        group_by(value) %>%
        summarize(Abundance = list(Abundance)) %>% # Create a list of Abundance values
        deframe() %>%
        zikw(., group = c("Classical", "Mesenchymal", "Proneural"), perm = T, perm.n = 100) %>%
        enframe() %>%
        unnest(cols = c(value)) %>%
        pivot_wider(names_from = name, values_from = value) %>%
        mutate(p.value = format.pval(p.value, digits = 2, eps = 0.001))

    ggplot(df_cor_species, aes(x = value, y = Abundance)) +
        geom_boxplot(width = 0.5) +
        geom_point(aes(color = mOS), size = 4) +
        scale_color_manual(values = c("Below_Median_OS" = "red", "Above_Median_OS" = "blue")) +
        geom_text(
            data = df_cor_stat,
            aes(label = paste("p.value:", p.value)),
            inherit.aes = FALSE,
            x = 2, y = Inf, hjust = 2, vjust = 1.5, size = 5, color = "black"
        ) +
        ggtitle(i) +
        xlab("GBM Subtypes") +
        labs(caption = "Kruskal-Wallis test for Zero-Inflated data") +
        theme_bw(base_size = 19) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
    ggsave(paste0("output/figures/immune_score_cluster/immune_score_component_classifiers_with_", i, ".pdf"),
        device = cairo_pdf,
        width = 10, height = 8
    )
}
