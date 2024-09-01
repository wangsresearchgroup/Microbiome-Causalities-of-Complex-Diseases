library(TwoSampleMR)
packages <- c('dagitty', 'dplyr', 'ggdag', 'ggrepel',
              'ggpubr', 'ggraph', 'stringr', 'tidyr')
pacman::p_load(char = packages)
figure_file_exts <- c('png', 'svg', 'pdf')

## MR between Gut Microbiota Abundance and Diseases

# Microbiome
microbiome <- c('Gut microbiota abundance', 'abundance in stool')

# Diseases
diseases <- c('Inflammatory Bowel Disease', 'Type 2 Diabetes',
              'Colorectal Cancer', 'Hypertension')

# Taxonomy conversion
taxonomy_conversion <- c("_sp_" = " subspecies ", "s_" = "species ",
                         "g_" = "genus ", "f_" = "family ", "o_" = "order ",
                         "c_" = "class ", "p_" = "phylum ", "k_" = "kingdom ")

# Custom function to save CSVs
save_csv <- function(df, path) {
  # Create the directory if it does not exist
  dir.create("Mendelian Randomization", showWarnings = FALSE)
  # Save to CSV
  sprintf('Mendelian Randomization/%s.csv', path) %>%
    write.csv(df, ., row.names = FALSE)
}

# List available GWASs
ao <- available_outcomes()

# Get IDs of diseases
get_id <- function(trait, searchFor) {
  trait %>%
    sapply(function(x) {
      searchFor %>%
        sapply(function(y) {
          grepl(y, x, ignore.case = TRUE) &&
            !grepl('without', x, ignore.case = TRUE) &&
            !grepl('proteinuria', x, ignore.case = TRUE)
        }) %>% any
      }) %>% which
}
microbiome_indices <- ao$trait %>% get_id(microbiome)
disease_indices <- ao$trait %>% get_id(diseases)

microbiome_ids <- ao[microbiome_indices, 'id'] %>% unlist
disease_ids <- ao[disease_indices, 'id'] %>% unlist

disease_cause_microbiome_res <- disease_cause_microbiome(microbiome_ids,
                                                         disease_ids)
disease_cause_microbiome_res %>% save_csv('Diseases vs Microbiome')

# Get instruments
exposure_dat <- extract_instruments(microbiome_ids)

# Get effects of instruments on outcome
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP,
                                    outcomes = disease_ids)

# Harmonize the exposure and outcome data
dat <- harmonise_data(exposure_dat, outcome_dat)

# Perform MR
method_list <- subset(mr_method_list(),
                      !(obj %in% c('mr_raps', 'mr_ivw_radial')))$obj
res <- mr(dat, method_list = method_list) %>%
  # BH p value correction
  mutate(pval = p.adjust(pval, 'BH'))

# Save to CSV
res %>% save_csv('Microbiome MR Results')

res <- read.csv('Mendelian Randomization/Microbiome MR Results.csv')
# Add column for disease
res <- res %>%
  mutate(disease = NA) %>%
  {for (disease in diseases) 
    .$disease[grepl(disease, .$outcome, ignore.case = TRUE)] <- disease
  .}

# Add column for microbe
res$microbe <- res$exposure %>%
  str_extract('(?<=Gut microbiota abundance \\().*?(?= id\\.\\d+\\)|\\))') %>%
  sapply(function(x) {
    x <- x %>% trim %>% gsub("  ", " ", .)
    if (str_detect(x, "^k_.*")) {
      x %>%
        str_replace_all(taxonomy_conversion) %>%
        str_replace_all("_", " ") %>%
        word(-1, sep = fixed("."))
    }
    else {
      x
    }
  })

# Remove columns with very close b, se, and pval
threshold <- 1e-3
taxonomic_ranks <- taxonomy_conversion %>% unname %>%
  gsub(" ", "", .) %>% c('unknown')
res <- res %>%
  separate(exposure, into = c("rank", "taxon"), sep = " ", extra = "merge") %>%
  mutate(rank = factor(rank, levels = taxonomic_ranks)) %>%
  mutate(group = paste(round(b / threshold) * threshold,
                       round(se / threshold) * threshold,
                       round(pval / threshold) * threshold,
                       disease, sep = "-")) %>%
  arrange(b, se, pval, disease, rank) %>%
  group_by(group) %>% slice_min(order(rank)) %>% ungroup() %>%
  distinct(group, .keep_all = TRUE) %>%
  unite(exposure, rank, taxon, sep = " ") %>%
  select(-group)

# Filter for significant microbes
b_threshold <- 2
pval_threshold <- 0.10
significant_microbes <- res %>% filter(abs(b) > b_threshold,
                                       pval < pval_threshold)

# Volcano plot
mr_volcano_plot <- ggplot(res %>% filter(pval > 1e-50),
                          aes(x = b, y = -log10(pval), color = disease)) +
  geom_point() +
  geom_text_repel(data = significant_microbes %>% filter(b < 0) %>%
                    mutate(microbe = gsub(" ", "\n", microbe)),
                  mapping = aes(label = microbe),
                  nudge_x = 0.5, nudge_y = 0.5, force = 80,
                  max.iter = 1e4, xlim = c(-5, -2), ylim = c(3, 10)) +
  geom_text_repel(data = significant_microbes %>% filter(b > 0) %>%
                    mutate(microbe = gsub(" ", "\n", microbe)),
                  mapping = aes(label = microbe),
                  nudge_x = 0.5, nudge_y = 0.5, force = 50,
                  max.iter = 1e4, xlim = c(2, 5), ylim = c(3, 10)) + theme_bw() +
  geom_vline(xintercept = c(-b_threshold, b_threshold), linetype = "dashed") +
  geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 18),
    axis.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.position = 'inside',
    legend.position.inside = c(0, 1),
    legend.background = element_blank(),
    legend.justification.inside = c(0, 1)) +
  labs(x = "Effect Size", y = expression("-log"[10]*"(p-value)"),
       color = "Gut microbiota abundance vs") +
  ggtitle("Volcano Plot of Mendelian Randomization Results")
ggsave('mr_volcano_plot.png', plot = mr_volcano_plot,
       width = 8, height = 8)

## Disease-disease bidirectional MR analysis
bidirectional_mr_diseases <- function(disease1, disease2) {
  # Get IDs
  disease1_ids <- ao$trait %>% get_id(disease1) %>% ao[., 'id'] %>% unlist
  disease2_ids <- ao$trait %>% get_id(disease2) %>% ao[., 'id'] %>% unlist
  # Get instruments
  disease1_exposure_dat <- extract_instruments(disease1_ids)
  disease2_exposure_dat <- extract_instruments(disease2_ids)
  # Get effects of instruments on outcome
  disease1_outcome_dat <- extract_outcome_data(snps=disease2_exposure_dat$SNP,
                                               outcomes = disease1_ids)
  disease2_outcome_dat <- extract_outcome_data(snps=disease1_exposure_dat$SNP,
                                               outcomes = disease2_ids)
  # Harmonize the exposure and outcome data
  disease12_dat <- harmonise_data(disease1_exposure_dat, disease2_outcome_dat)
  disease21_dat <- harmonise_data(disease2_exposure_dat, disease1_outcome_dat)
  # MR with BH p value correction
  if (disease1 == 'Inflammatory Bowel Disease' && disease2 == 'Colorectal Cancer') {
    method_list <- mr_method_list() %>%
      subset(use_by_default) %>% .$obj %>% setdiff('mr_ivw')
  }
  else {
    method_list <- subset(mr_method_list(), use_by_default)$obj
  }
  disease12_res <- mr(disease12_dat, method_list = method_list) %>%
    mutate(pval = p.adjust(pval, 'BH'))
  disease21_res <- mr(disease21_dat) %>%
    mutate(pval = p.adjust(pval, 'BH'))
  # Return MR results
  return(rbind(disease12_res, disease21_res))
}

ibd_crc_mr_results <- bidirectional_mr_diseases('Inflammatory Bowel Disease',
                                                'Colorectal Cancer')
t2d_hypertension_mr_results <- bidirectional_mr_diseases('Type 2 Diabetes',
                                                         'Hypertension')

ibd_crc_mr_results <- read.csv('Mendelian Randomization/IBD vs CRC MR Results.csv')
t2d_hypertension_mr_results <- read.csv('Mendelian Randomization/T2D vs Hypertension MR Results.csv')

ibd_crc_mr_results <- ibd_crc_mr_results %>%
  mutate(disease = ifelse(grepl("Colorectal Cancer", outcome,
                                         ignore.case = TRUE),
                                   "Colorectal Cancer vs Inflammatory Bowel Disease",
                                   "Inflammatory Bowel Disease vs Colorectal Cancer"))
t2d_hypertension_mr_results <- t2d_hypertension_mr_results %>%
  mutate(disease = ifelse(grepl("Type 2 Diabetes", outcome,
                                         ignore.case = TRUE),
                                   "Type 2 Diabetes vs Hypertension",
                                   "Hypertension vs Type 2 Diabetes"))

make_volcano_plot <- function(mr_results, b_threshold = 2, pval_threshold = 0.10) {
  # Volcano plot
  return(ggplot(mr_results, aes(x = b, y = -log10(pval),
                                           color = disease)) +
    geom_point() + theme_bw() +
    geom_vline(xintercept = c(-b_threshold, b_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 18),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.title = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.position = 'inside',
      legend.position.inside = c(0, 1),
      legend.background = element_blank(),
      legend.justification.inside = c(0, 1)) +
    labs(x = "Effect Size", y = expression("-log"[10]*"(p-value)"),
         color = "Disease") +
    ggtitle("Volcano Plot of Mendelian Randomization Results"))
}
ibd_crc_mr_results %>% make_volcano_plot %>%
  ggsave('ibd_crc_mr_volcano_plot.png', plot = ., width = 8, height = 8)
t2d_hypertension_mr_results %>% make_volcano_plot %>%
  ggsave('t2d_hypertension_mr_volcano_plot.png', plot = ., width = 8, height = 8)

# Diet
diet_ids <- ao %>% filter(grepl('diet', trait)) %>%
  filter(!grepl('questionnaire', trait)) %>% .$id
diet_exposure_dat <- extract_instruments(diet_ids)
diet_disease_outcome_dat <- extract_outcome_data(snps=diet_exposure_dat$SNP,
                                                 outcomes = disease_ids)
diet_disease_dat <- harmonise_data(diet_exposure_dat,
                                   diet_disease_outcome_dat)
diet_disease_res <- mr(diet_disease_dat) %>% mutate(pval = p.adjust(pval, 'BH'))
diet_disease_res %>% save_csv('Diet vs Diseases MR Results')

# Produce diet volcano plot
diet_mr_results <- 'Mendelian Randomization/Diet vs Diseases MR Results.csv' %>% read.csv
diet_mr_results$disease <- NA
for (i in 1:nrow(diet_mr_results)) {
  for (disease in diseases) {
    if (grepl(disease, diet_mr_results[i, 'outcome'], ignore.case = TRUE)) {
      diet_mr_results[i, 'disease'] <- sprintf('Diet vs %s', disease)
      break
    }
  }
}
diet_mr_results %>% make_volcano_plot %>%
  ggsave('diet_mr_volcano_plot.png', plot = ., width = 8, height = 8)

# Gut bacterial pathway abundance
pathway_indices <- ao$trait %>%
  grepl('Gut bacterial pathway', ., ignore.case = TRUE) %>% which
pathway_ids <- ao$id[pathway_indices]

pathway_exposure_dat <- extract_instruments(pathway_ids)

pathway_exposure_dat$exposure <- pathway_exposure_dat$exposure %>%
  str_extract('(?<=\\().*?(?= id\\.\\d+\\)|\\))') %>%
  substr(., regexpr('\\.\\.', .) + 2, nchar(.)) %>%
  gsub('L\\.', 'L-', .) %>%
  gsub('\\.\\.\\.', ', ', .) %>%
  gsub('5\\.\\.', "5'-", .) %>%
  gsub('N\\.', 'N-', .) %>%
  gsub('E..coli', 'E. coli', .) %>%
  ifelse(endsWith(., '.'),
         str_replace(., "\\.(?=[^.]*$)", ")"), .) %>%
  sapply(function(x) {
    if (endsWith(x, ')')) {
      idx <- x %>%
        gregexpr("\\.{2}", .) %>% .[[1]] %>% .[length(.)]
      substr(x, idx, idx + 1) <- ' ('
    }
    return(x)
  }) %>% unname %>%
  gsub('\\.\\.', ', ', .) %>%
  gsub("(\\d)\\.(\\d)", "\\1,\\2", .) %>%
  gsub("(\\d)\\.", "\\1-", .) %>%
  gsub('\\.', ' ', .) %>% gsub('E  coli', 'E. coli', .)

pathway_disease_outcome_dat <- extract_outcome_data(snps=pathway_exposure_dat$SNP,
                                                 outcomes = disease_ids)
pathway_disease_dat <- harmonise_data(pathway_exposure_dat,
                                      pathway_disease_outcome_dat)
pathway_disease_res <- mr(pathway_disease_dat) %>% mutate(pval = p.adjust(pval, 'BH'))
pathway_disease_res %>% save_csv('Gut Bacterial Pathways vs Diseases MR Results')

pathway_disease_res <- read.csv('Mendelian Randomization/Gut Bacterial Pathways vs Diseases MR Results.csv')
pathway_disease_res$disease <- NA
for (disease in diseases) {
  idxs <- pathway_disease_res$outcome %>% grepl(disease, ., ignore.case = TRUE)
  pathway_disease_res$disease[idxs] <- disease
}
(make_volcano_plot(pathway_disease_res) + labs(color = 'Gut Bacterial Pathways vs')) %>%
  ggsave('bacterial_pathways_mr_volcano_plot.png', plot = ., width = 8, height = 8)

# Abundance in stool
stool_indices <- ao$trait %>%
  grepl('abundance in stool', ., ignore.case = TRUE) %>% which
stool_ids <- ao$id[stool_indices]

stool_exposure_dat <- extract_instruments(stool_ids)
stool_disease_outcome_dat <- extract_outcome_data(snps=stool_exposure_dat$SNP,
                                                    outcomes = disease_ids)
stool_disease_dat <- harmonise_data(stool_exposure_dat,
                                    stool_disease_outcome_dat)
stool_disease_res <- mr(stool_disease_dat) %>% mutate(pval = p.adjust(pval, 'BH'))
stool_disease_res %>% save_csv('Abundance in Stool vs Diseases MR Results')

stool_disease_res <- read.csv('Mendelian Randomization/Abundance in Stool vs Diseases MR Results.csv')
stool_disease_res <- stool_disease_res %>%
  filter(!is.na(se), !is.na(pval)) %>%
  mutate(microbe = exposure %>% strsplit(' abundance')
         %>% sapply(function(x) x[1])) %>%
  mutate(disease = NA) %>%
  {for (disease in diseases) 
    .$disease[grepl(disease, .$outcome, ignore.case = TRUE)] <- disease
  .}

stool_significant_bacteria <- stool_disease_res %>%
  filter(abs(b) > 4, pval < 0.01)
stool_disease_mr_volcano_plot <- make_volcano_plot(stool_disease_res %>%
                                                     filter(abs(b) < 100)) +
  labs(color = 'Abundance in Stool vs') +
  geom_text_repel(data = stool_significant_bacteria %>%
                    filter(b < 0) %>%
                    mutate(microbe = gsub(" ", "\n", microbe)),
                  mapping = aes(label = microbe),
                  nudge_x = 2, nudge_y = 2, force = 100,
                  max.iter = 5e4, xlim = c(-20, -10)) +
  geom_text_repel(data = stool_significant_bacteria %>%
                    filter(b > 0) %>%
                    mutate(microbe = gsub(" ", "\n", microbe)),
                  mapping = aes(label = microbe),
                  nudge_x = 2, nudge_y = 2, force = 100,
                  max.iter = 5e4, xlim = c(10, 20))
ggsave('stool_abundance_mr_volcano_plot.png',
       plot = stool_disease_mr_volcano_plot,
       width = 8, height = 8)

# Combine both microbiome volcano plots
combined_volcano_plots <- ggarrange(mr_volcano_plot + labs(x = '', y = '', title = ''),
          stool_disease_mr_volcano_plot + labs(x = '', y = '', title = ''),
          ncol = 2, labels = c("A", "B")) %>%
  annotate_figure(
    bottom = text_grob("Effect Size", size = 14),
    left = text_grob(expression("-log"[10]*"(p-value)"), size = 14, rot = 90)
  )
for (file_ext in figure_file_exts) {
  paste0('combined_volcano_plots.', file_ext) %>%
    ggsave(plot = combined_volcano_plots, width = 16, height = 9)
}

# Create causality directed acyclic graph
columns_to_read <- c('id.exposure', 'exposure', 'outcome', 'b', 'pval')

compiled_mr_results <- c('Microbiome', 'Abundance in Stool vs Diseases') %>%
  sprintf('Mendelian Randomization/%s MR Results.csv', .) %>%
  lapply(function(x) read.csv(x) %>% .[, columns_to_read]) %>%
  do.call(rbind, .) %>%
  filter(abs(b) > b_threshold, pval < pval_threshold)
exposure_ids <- compiled_mr_results[, 'id.exposure'] %>% unique

exposure_dat <- extract_instruments(exposure_ids)
outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP,
                                    outcomes = exposure_ids)
dat <- harmonise_data(exposure_dat, outcome_dat)
res <- mr(dat) %>% mutate(pval = p.adjust(pval, 'BH'))
res %>% save_csv('Significant Exposures MR Results')

bacteria_taxa <- read.csv("Mendelian Randomization/Bacteria taxa.csv")
compiled_mr_results <- read.csv('Mendelian Randomization/Significant Exposures MR Results.csv') %>%
  filter(exposure != outcome, abs(b) > b_threshold, pval < pval_threshold) %>%
  select(all_of(columns_to_read)) %>%
  rbind(compiled_mr_results) %>%
  mutate(id.exposure = NULL) %>%
  mutate_at(vars(exposure, outcome), ~str_replace(., " \\|\\|.*", "")) %>%
  {for (disease in diseases) 
    .$outcome[grepl(disease, .$outcome, ignore.case = TRUE)] <- disease
  .} %>%
  mutate_at(vars(exposure, outcome),
            ~str_replace(., ' abundance in stool', '')) %>%
  mutate_at(vars(exposure, outcome), ~ifelse(str_detect(., "Gut microbiota abundance"),
                                             str_extract(., "(?<=\\().+?(?=\\))"),
                                             .)) %>%
  mutate_at(vars(exposure, outcome), ~ifelse(str_detect(., " id"),
                                             str_extract(., ".*(?= id)"),
                                             .)) %>%
  filter(!grepl('UBP9', exposure), !grepl('UBP9', outcome)) %>%
  mutate(exposure = str_replace(exposure,
                     'k_Bacteria.p_Firmicutes.c_Negativicutes.o_Selenomonadales',
                     'order Selenomonadales')) %>%
  mutate(exposure = str_replace(exposure,
                                'k_Bacteria.p_Firmicutes.c_Negativicutes',
                                'class Negativicutes')) %>%
  mutate_at(vars(exposure, outcome), ~ {
    match_idx <- match(., bacteria_taxa$bacteria)
    ifelse(!is.na(match_idx), paste(bacteria_taxa$taxa[match_idx],
                                    bacteria_taxa$bacteria[match_idx]) %>% trim, .)
  }) %>%
  mutate(exposure = str_replace(exposure, ' sp000432435', '')) %>%
  mutate_at(vars(exposure, outcome), ~ifelse(grepl(' \\w$', .),
                                             gsub(' \\w$', '', .), .))
save_csv(compiled_mr_results, 'Combined Microbiome MR Results')

exposures_and_outcomes <- compiled_mr_results[, c('exposure', 'outcome')] %>%
  unlist %>% unique

formatted_e_and_o <- ifelse(exposures_and_outcomes %in% diseases,
       paste0('"', exposures_and_outcomes, '" [outcome]'),
       paste0('"', exposures_and_outcomes, '" [exposure]')) %>%
  paste(collapse = ' ')

shorten_dag_arrows <- function(tidy_dag, node_size = 16, offset = 10, scale = 16 / 9){
  tidy_dag$data <- tidy_dag$data %>%
    mutate(angle = atan2(yend - y, xend - x),
           xstart_arrow = x + (offset + node_size) * cos(angle),
           ystart_arrow = y + (scale * offset + node_size) * sin(angle),
           xend_arrow = xend - (offset + node_size) * cos(angle),
           yend_arrow = yend - (scale * offset + node_size) * sin(angle))
  return(tidy_dag)
}

microbiome_dag_plot <- compiled_mr_results %>%
  mutate(relation = paste(exposure, outcome, sep = '\" -> \"')) %>%
  .$relation %>% unlist %>% paste(collapse = '\" \"') %>%
  paste(formatted_e_and_o, ., sep = ' "') %>%
  sprintf('dag { %s" }', .) %>%
  dagitty %>% ggdag(layout = 'gem') %>%
  shorten_dag_arrows + theme_dag_blank()

microbiome_dag_plot$data <- microbiome_dag_plot$data %>%
  mutate(outcome.exposure = ifelse(is.na(direction), 'Disease', 'Bacterial Abundance')) %>%
  left_join(compiled_mr_results, by = c('name' = 'exposure', 'to' = 'outcome')) %>%
  mutate(effect_direction = as.factor(ifelse(is.na(b), 'Decrease', ifelse(b > 0, "Increase", "Decrease")))) %>%
  mutate_at(vars(name, to), ~str_replace_all(., ' ', '\n')) %>%
  mutate_at(vars(name, to), ~str_replace_all(., 'Type\n2', 'Type 2')) %>%
  mutate_at(vars(name, to), ~str_replace_all(., 'bac', '-\nbac')) %>%
  mutate_at(vars(name, to), ~str_replace_all(., 'mon', '-\nmon')) %>%
  mutate_at(vars(name, to), ~str_replace_all(., 'Eu-\n', 'Eu')) %>%
  mutate_at(vars(name, to), ~str_replace_all(., 'sp.\n', 'sp. '))

microbiome_dag_plot$layers[[3]]$mapping <- aes(color = outcome.exposure)

microbiome_dag_plot <- microbiome_dag_plot +
  geom_dag_edges(edge_color = 'white') +
  geom_segment(aes(x = xstart_arrow, y = ystart_arrow, xend = xend_arrow,
                   yend = yend_arrow, color = effect_direction,
                   linewidth = scale(log2(abs(b)) - log10(pval))),
               arrow = arrow(length = unit(0.02, "npc"), type = "closed"),
               alpha = 0.5) + geom_dag_text(color = 'black') +
  scale_color_manual(values = c('Increase' = '#0072B2', 'Decrease' = '#D55E00',
                                'Disease' = '#CC79A7',
                                'Bacterial Abundance' = '#F0E442')) +
  labs(color = 'Exposure Variable Effect') + guides(linewidth = 'none') +
  theme(legend.position = 'inside', legend.position.inside = c(1, 1),
        legend.background = element_blank(),
        legend.justification.inside = c(1, 1))

microbiome_dag_plot

for (file_ext in figure_file_exts) {
  paste0('microbiome_dag.', file_ext) %>%
    ggsave(plot = microbiome_dag_plot, width = 16, height = 9)
}
