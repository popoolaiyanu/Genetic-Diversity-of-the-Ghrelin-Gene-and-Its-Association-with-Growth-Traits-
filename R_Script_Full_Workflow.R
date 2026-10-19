## ==================================================
## Genetic Diversity of the Ghrelin Gene and Its Association
## with Growth Traits in Fulani and Yoruba Ecotypes of
## Nigerian Indigenous Chickens (Exon and Intron SNPs)
##
## Full reproducible pipeline
## ==================================================

## ------------------------------------------------------------
## WORKFLOW OVERVIEW
## ------------------------------------------------------------
## PART A - Exon SNPs (SNP1615, SNP2355)
##   0. SETUP: packages, working directory and input files
##   1. DATA IMPORT, CLEANING AND RESHAPING
##   2. DATA SCREENING: extreme values
##   3. DESCRIPTIVE STATISTICS AND NORMALITY (full population)
##   4. GENETIC DIVERSITY: exon SNPs (SNP1615, SNP2355)
##   5. HARDY-WEINBERG EQUILIBRIUM: exon SNPs
##   6. GROWTH OF THE GENOTYPED BIRDS: genotype group means
##   7. EXON SNPs: GENOTYPE AND GROWTH RATE (primary analysis)
##   8. EXON SNPs: GENOTYPE AND TRAIT VALUES BY AGE (weekly tests)
##   9. EXON SNPs: LINEAR MIXED-EFFECTS MODELS
##   10. EXON SNPs: SENSITIVITY ANALYSIS (Theil-Sen growth rates)
## PART B - Intron SNPs (16 SNPs) and the combined analysis
##   11. INTRON SNPs: IMPORT AND DATA CHECK
##   12. INTRON SNPs: DIVERSITY INDICES
##   13. INTRON + EXON SNPs COMBINED: DIVERSITY INDICES
##   14. FULANI vs YORUBA: PAIRED TESTS ON THE DIVERSITY INDICES
##   15. GENOTYPE ANOVA ON GROWTH RATE (intron + exon SNPs)
##   16. GENOTYPE ANOVA AT INDIVIDUAL AGES (intron + exon SNPs)
##   17. HARDY-WEINBERG EQUILIBRIUM: all SNPs
## PART C - Figures and session information
##   18. FIGURES
##   19. SESSION INFORMATION
##
## Run the sections in order, in a fresh R session (Session > Restart R).
## Sections marked NEW have not yet been run in R; compare their results
## with the values reported in the two Results documents.

## ------------------------------------------------------------
## 0. SETUP: packages, working directory and input files
## ------------------------------------------------------------
# Packages, working folder, input file names and settings.
# install.packages(c("readxl", "tidyverse", "e1071", "lme4", "lmerTest", "emmeans"))
library(readxl)
library(tidyverse)
library(e1071)      # skewness and kurtosis
library(lme4)
library(lmerTest)   # p-values for mixed models (load after lme4)
library(emmeans)

options(width = 150, timeout = 600)
set.seed(2026)

work_dir <- "C:/Users/DELL/OneDrive/Documents/DR LYDIA"    # <- change to your folder
setwd(work_dir)

# input files (all must be in work_dir)
f_assoc_F <- "FA Association 1.xlsx"          # genotyped Fulani birds (growth + exon SNPs)
f_assoc_Y <- "YA Association 1.xlsx"          # genotyped Yoruba birds
f_pop_F   <- "Fulani Ecotype Population.xlsx" # full Fulani population (n = 75)
f_pop_Y   <- "Yoruba Ecotype Population.xlsx" # full Yoruba population (n = 108)
f_intron  <- "INTRONIC SNP.xlsx"              # intronic SNP genotypes (16 SNPs)

traits    <- c("Weight", "BodyLength", "ChestGirth", "WingLength", "Height",
               "LegLength", "ShankLength")
key_weeks <- c(4, 8, 12, 16, 20, 24)   # ages used for the weekly genotype tests
min_n     <- 5                         # minimum birds per group (intron summaries)

## ------------------------------------------------------------
## 1. DATA IMPORT, CLEANING AND RESHAPING
## ------------------------------------------------------------
# Reads the four files, applies the confirmed corrections and converts to long format.
# Feeds: Both Results documents
# NOTE: the Y54 lines in this section are new (not yet run in R)
fa <- read_excel(f_assoc_F, col_types = "text")
ya <- read_excel(f_assoc_Y, col_types = "text")
fp <- read_excel(f_pop_F,   col_types = "text")
yp <- read_excel(f_pop_Y,   col_types = "text")

assoc_wide <- bind_rows(fa, ya)   # 20 genotyped birds
pop_wide   <- bind_rows(fp, yp)   # full population

# correction confirmed from the field records: F3 week-4 weight is 168.3 g
wk4 <- grep("^\\s*Wei\\(g\\)\\s*WK4", names(assoc_wide), ignore.case = TRUE, value = TRUE)
assoc_wide[assoc_wide$`Bird ID` == "F3", wk4] <- "168.3"

# "32,2" becomes 32.2 (comma decimal); the stray "S" (bird F46, week 14) becomes NA
clean_num <- function(x) suppressWarnings(as.numeric(gsub(",", ".", trimws(x))))

make_long <- function(df) {
  df %>%
    rename(bird_id = `Bird ID`, ecotype = Ecotype) %>%
    pivot_longer(-any_of(c("bird_id", "ecotype", "SNP1615", "SNP2355")),
                 names_to = "var", values_to = "raw") %>%
    mutate(
      code  = str_extract(str_trim(var), "^[A-Za-z]+"),
      week  = as.numeric(str_extract(var, "\\d+$")),
      trait = recode(code, WH = "Weight", Wei = "Weight", BL = "BodyLength",
                     CG = "ChestGirth", WL = "WingLength", H = "Height",
                     LL = "LegLength", SL = "ShankLength"),
      value = clean_num(raw)
    ) %>%
    select(-var, -code, -raw)
}

# genotypes: "*" means no genotype recorded; such birds are excluded for that SNP
assoc_wide <- assoc_wide %>%
  mutate(across(c(SNP1615, SNP2355),
                ~ factor(na_if(trimws(.x), "*"), levels = c("AA", "AG", "GG"))))

pop_long_raw <- make_long(pop_wide)      # population data exactly as recorded
assoc_long   <- make_long(assoc_wide)    # genotyped birds

# two impossible values for bird Y54 (decimal-point errors) are set to missing:
# body length at week 14 (333.8 cm) and wing length at week 16 (255.5 cm).
# Y54 was not genotyped, so no association result is affected.
pop_long <- pop_long_raw %>%
  mutate(value = if_else((bird_id == "Y54" & trait == "BodyLength" & week == 14) |
                         (bird_id == "Y54" & trait == "WingLength" & week == 16),
                         NA_real_, value))

# ---- checks ----
print(table(pop_wide$Ecotype))                      # expect Fulani 75, Yoruba 108
print(table(assoc_wide$Ecotype))                    # expect 10 and 10
print(table(assoc_wide$Ecotype, assoc_wide$SNP1615, useNA = "ifany"))
print(table(assoc_wide$Ecotype, assoc_wide$SNP2355, useNA = "ifany"))
print(pop_long %>% distinct(trait, week) %>% count(trait))   # weight 13 ages, others 12
print(assoc_long %>% filter(bird_id == "F3", trait == "Weight", week == 4))  # 168.3

saveRDS(pop_long,    "pop_long_clean.rds")
saveRDS(assoc_wide,  "assoc_wide_clean.rds")
saveRDS(assoc_long,  "assoc_long_clean.rds")

## ------------------------------------------------------------
## 2. DATA SCREENING: extreme values
## ------------------------------------------------------------
# Lists values far from the group median; this is how the Y54 errors were found.
# Feeds: Data screening notes
# NOTE: NEW - this section has not yet been run in R
screen <- pop_long_raw %>%
  filter(!is.na(value)) %>%
  group_by(ecotype, trait, week) %>%
  mutate(med  = median(value),
         mad_ = mad(value),
         z    = (value - med) / if_else(mad_ == 0, sd(value), mad_)) %>%
  ungroup() %>%
  filter(abs(z) > 4.5) %>%
  arrange(desc(abs(z)))
print(screen %>% select(ecotype, trait, week, bird_id, value, med, z), n = 40)

## ------------------------------------------------------------
## 3. DESCRIPTIVE STATISTICS AND NORMALITY (full population)
## ------------------------------------------------------------
# Mean, median, mode, dispersion, skewness, kurtosis and Shapiro-Wilk test by ecotype, trait and age.
# Feeds: Exon document: Tables 1, 2 and A1-A7
get_mode <- function(x) {            # most repeated value (1 decimal); NA if none repeats
  x <- round(x[!is.na(x)], 1)
  if (length(x) == 0) return(NA_real_)
  tb <- table(x)
  if (max(tb) == 1) return(NA_real_)
  as.numeric(names(tb)[which.max(tb)])
}

desc <- pop_long %>%
  group_by(ecotype, trait, week) %>%
  summarise(
    n        = sum(!is.na(value)),
    mean     = mean(value, na.rm = TRUE),
    median   = median(value, na.rm = TRUE),
    mode     = get_mode(value),
    sd       = sd(value, na.rm = TRUE),
    variance = var(value, na.rm = TRUE),
    se       = sd / sqrt(n),
    cv_pct   = sd / mean * 100,
    min      = min(value, na.rm = TRUE),
    max      = max(value, na.rm = TRUE),
    range    = max - min,
    q1       = quantile(value, 0.25, na.rm = TRUE),
    q3       = quantile(value, 0.75, na.rm = TRUE),
    iqr      = q3 - q1,
    skewness = skewness(value, na.rm = TRUE, type = 2),
    kurtosis = kurtosis(value, na.rm = TRUE, type = 2),
    shapiro_p = if (n >= 3 && n <= 5000 && sd > 0)
                  shapiro.test(value[!is.na(value)])$p.value else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(desc, "descriptive_statistics_ALL.csv", row.names = FALSE)
nrow(desc)                                          # expect 170

# normality summary across ages (Table 2)
norm_summary <- desc %>%
  filter(!is.na(shapiro_p)) %>%
  group_by(ecotype, trait) %>%
  summarise(groups_tested = dplyr::n(),
            normal_p_gt_0.05 = sum(shapiro_p > 0.05),
            non_normal_p_le_0.05 = sum(shapiro_p <= 0.05),
            mean_abs_skew = round(mean(abs(skewness)), 2), .groups = "drop")
print(norm_summary, n = 20)
write.csv(norm_summary, "normality_summary.csv", row.names = FALSE)

## ------------------------------------------------------------
## 4. GENETIC DIVERSITY: exon SNPs (SNP1615, SNP2355)
## ------------------------------------------------------------
# Genotype and allele frequencies, Ho, He, Ne and PIC by ecotype and pooled.
# Feeds: Exon document: Tables 3 and 4
assoc_geno <- assoc_wide %>% rename(bird_id = `Bird ID`, ecotype = Ecotype)

snp_long <- assoc_geno %>%
  select(bird_id, ecotype, SNP1615, SNP2355) %>%
  pivot_longer(c(SNP1615, SNP2355), names_to = "snp", values_to = "genotype")
snp_pooled <- bind_rows(snp_long, mutate(snp_long, ecotype = "Pooled"))

diversity <- function(g) {
  n_total <- length(g)
  g   <- g[!is.na(g)]
  cnt <- table(factor(g, levels = c("AA", "AG", "GG")))
  n   <- sum(cnt)
  p   <- (2 * cnt[["AA"]] + cnt[["AG"]]) / (2 * n)    # frequency of allele A
  q   <- 1 - p                                        # frequency of allele G
  tibble(
    n_total = n_total, n_excluded = n_total - n, n_genotyped = n,
    AA = cnt[["AA"]], AG = cnt[["AG"]], GG = cnt[["GG"]],
    freq_AA = AA / n_genotyped, freq_AG = AG / n_genotyped, freq_GG = GG / n_genotyped,
    freq_A = p, freq_G = q,
    Ho  = AG / n_genotyped,                           # observed heterozygosity
    He  = 1 - (p^2 + q^2),                            # expected heterozygosity
    Ne  = 1 / (p^2 + q^2),                            # effective number of alleles
    PIC = 1 - (p^2 + q^2) - 2 * p^2 * q^2
  )
}

diversity_results <- snp_pooled %>%
  mutate(ecotype = factor(ecotype, levels = c("Fulani", "Yoruba", "Pooled"))) %>%
  group_by(snp, ecotype) %>%
  group_modify(~ diversity(.x$genotype)) %>%
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(diversity_results, width = Inf)
write.csv(diversity_results, "genetic_diversity_results.csv", row.names = FALSE)

## ------------------------------------------------------------
## 5. HARDY-WEINBERG EQUILIBRIUM: exon SNPs
## ------------------------------------------------------------
# Chi-square, likelihood-ratio and exact tests. The exact test function is reused in section 7.
# Feeds: Exon document: Table 5
hwe_exact <- function(nAA, nAB, nBB) {            # exact test (two alleles)
  n  <- nAA + nAB + nBB
  nA <- 2 * nAA + nAB
  nB <- 2 * nBB + nAB
  hs <- seq(nA %% 2, min(nA, nB), by = 2)
  logp <- sapply(hs, function(h) {
    aa <- (nA - h) / 2; bb <- (nB - h) / 2
    lgamma(n + 1) - lgamma(aa + 1) - lgamma(h + 1) - lgamma(bb + 1) +
      h * log(2) + lgamma(nA + 1) + lgamma(nB + 1) - lgamma(2 * n + 1)
  })
  p <- exp(logp)
  pobs <- p[hs == nAB]
  sum(p[p <= pobs * (1 + 1e-8)])
}

hwe_test <- function(nAA, nAG, nGG) {
  n   <- nAA + nAG + nGG
  p   <- (2 * nAA + nAG) / (2 * n);  q <- 1 - p
  obs <- c(nAA, nAG, nGG)
  exp_counts <- n * c(p^2, 2 * p * q, q^2)
  chi <- sum((obs - exp_counts)^2 / exp_counts)
  G   <- 2 * sum(ifelse(obs > 0, obs * log(obs / exp_counts), 0))
  tibble(
    n = n, obs_AA = nAA, obs_AG = nAG, obs_GG = nGG,
    exp_AA = exp_counts[1], exp_AG = exp_counts[2], exp_GG = exp_counts[3],
    chi2 = chi, df = 1, p_chi2 = pchisq(chi, 1, lower.tail = FALSE),
    G2 = G, p_LRT = pchisq(G, 1, lower.tail = FALSE),
    p_exact = hwe_exact(nAA, nAG, nGG),
    min_expected = min(exp_counts),
    F_is = 1 - (nAG / n) / (1 - (p^2 + q^2))
  )
}

hwe_results <- diversity_results %>%
  rowwise() %>%
  reframe(snp, ecotype, hwe_test(AA, AG, GG)) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)),
         HWE_status = ifelse(p_exact < 0.05, "Deviates", "In HWE"))
print(hwe_results, width = Inf)
write.csv(hwe_results, "hwe_results.csv", row.names = FALSE)

## ------------------------------------------------------------
## 6. GROWTH OF THE GENOTYPED BIRDS: genotype group means
## ------------------------------------------------------------
# Mean and SD of body weight of the genotyped birds, overall and by genotype (Tables 6 to 8 of the exon document).
# Feeds: Exon document: Tables 6 to 8
# NOTE: NEW - this section has not yet been run in R
assoc_snp <- assoc_long %>%
  pivot_longer(c(SNP1615, SNP2355), names_to = "snp", values_to = "genotype") %>%
  filter(!is.na(genotype), !is.na(value)) %>%
  mutate(genotype = factor(genotype, levels = c("AA", "AG", "GG")))

# Table 6: all genotyped birds, by ecotype and age
weight_genotyped <- assoc_long %>%
  filter(trait == "Weight", !is.na(value)) %>%
  group_by(ecotype, week) %>%
  summarise(n = dplyr::n(), M = mean(value), SD = sd(value), .groups = "drop")
print(weight_genotyped, n = 30)

# Tables 7 and 8 and the genotype means: every trait, age, SNP, ecotype and genotype
geno_means <- assoc_snp %>%
  group_by(snp, ecotype, genotype, trait, week) %>%
  summarise(n = dplyr::n(), mean = mean(value), sd = sd(value),
            se = sd / sqrt(n), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 2)))
write.csv(geno_means, "genotype_group_means.csv", row.names = FALSE)

## ------------------------------------------------------------
## 7. EXON SNPs: GENOTYPE AND GROWTH RATE (primary analysis)
## ------------------------------------------------------------
# Each bird's growth rate (slope on age, weeks 2 to 24) compared between genotypes.
# Feeds: Exon document: Tables 9 and 10
perm_p <- function(x, g, nperm = 20000) {         # exact permutation test, two genotypes
  g <- droplevels(factor(g)); lv <- levels(g)
  n <- length(x); n2 <- sum(g == lv[2])
  obs <- mean(x[g == lv[2]]) - mean(x[g == lv[1]])
  if (choose(n, n2) <= nperm) {
    idx <- combn(n, n2)
    d <- apply(idx, 2, function(i) mean(x[i]) - mean(x[-i]))
  } else {
    set.seed(2026)
    d <- replicate(nperm, { i <- sample(n, n2); mean(x[i]) - mean(x[-i]) })
  }
  mean(abs(d) >= abs(obs) - 1e-9)
}

slope_two <- function(d) {                        # Welch t test + permutation test
  d <- droplevels(d); lv <- levels(d$genotype)
  a <- d$slope[d$genotype == lv[1]]; b <- d$slope[d$genotype == lv[2]]
  tt <- t.test(b, a)
  tibble(ref = lv[1], comp = lv[2], n_ref = length(a), n_comp = length(b),
         mean_ref = mean(a), mean_comp = mean(b), diff = mean(b) - mean(a),
         ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
         p_welch = tt$p.value, p_perm = perm_p(d$slope, d$genotype))
}

# birds measured through week 24 (F13 has no records after week 12 and is left out)
complete_birds <- assoc_long %>%
  filter(trait == "Weight", week == 24, !is.na(value)) %>%
  distinct(bird_id) %>% pull(bird_id)

slopes_full <- assoc_snp %>%
  filter(week >= 2, bird_id %in% complete_birds) %>%
  group_by(snp, ecotype, genotype, trait, bird_id) %>%
  filter(dplyr::n() >= 10) %>%
  summarise(slope = unname(coef(lm(value ~ week))[2]), .groups = "drop")

bird_level_tests <- function(sl) {
  # SNP1615: two genotypes compared within each ecotype (14 comparisons)
  s1615 <- sl %>%
    filter(snp == "SNP1615") %>%
    group_by(ecotype, trait) %>%
    group_modify(~ slope_two(.x)) %>%
    ungroup() %>%
    mutate(p_adj_14tests = p.adjust(p_welch, "BH"))
  # SNP2355: additive effect of the G allele (AA = 0, AG = 1, GG = 2), ecotype as covariate
  s2355 <- sl %>%
    filter(snp == "SNP2355") %>%
    mutate(dose = as.numeric(genotype) - 1) %>%
    group_by(trait) %>%
    group_modify(~ {
      m  <- lm(slope ~ dose + ecotype, data = .x)
      ci <- confint(m)["dose", ]
      tibble(n = nrow(.x), effect_per_G = unname(coef(m)["dose"]),
             ci_low = ci[[1]], ci_high = ci[[2]],
             p = summary(m)$coefficients["dose", 4])
    }) %>%
    ungroup() %>%
    mutate(p_adj = p.adjust(p, "BH"))
  list(s1615 = s1615, s2355 = s2355)
}

bl <- bird_level_tests(slopes_full)
sl1615_full <- bl$s1615 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
sl2355_full <- bl$s2355 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(sl1615_full), row.names = FALSE)
print(as.data.frame(sl2355_full), row.names = FALSE)
write.csv(sl1615_full, "birdlevel_full_SNP1615.csv", row.names = FALSE)
write.csv(sl2355_full, "birdlevel_full_SNP2355.csv", row.names = FALSE)

## ------------------------------------------------------------
## 8. EXON SNPs: GENOTYPE AND TRAIT VALUES BY AGE (weekly tests)
## ------------------------------------------------------------
# Weekly comparisons: Welch t test within ecotype for SNP1615; ecotype-adjusted F test for SNP2355.
# Feeds: Exon document: Tables 11, A8 and A9
# NOTE: the sd, t, df, F and eta-squared columns are new (not yet run in R)
get_d <- function(s, e, tr) assoc_snp %>% filter(snp == s, ecotype == e, trait == tr)

two_group_weeks <- function(d, weeks) {
  map_dfr(weeks, function(w) {
    x  <- d %>% filter(week == w, !is.na(value)) %>% droplevels()
    lv <- levels(x$genotype)
    if (length(lv) < 2 || min(table(x$genotype)) < 2) return(NULL)
    a  <- x$value[x$genotype == lv[1]]; b <- x$value[x$genotype == lv[2]]
    tt <- t.test(b, a)
    tibble(week = w, ref = lv[1], comp = lv[2], n_ref = length(a), n_comp = length(b),
           mean_ref = mean(a), mean_comp = mean(b), sd_ref = sd(a), sd_comp = sd(b),
           diff = unname(tt$estimate[1] - tt$estimate[2]),
           ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
           t = unname(tt$statistic), df = unname(tt$parameter),
           p_welch = tt$p.value, p_perm = perm_p(x$value, x$genotype))
  })
}

p1615 <- expand_grid(ecotype = c("Fulani", "Yoruba"), trait = traits) %>%
  mutate(res = map2(ecotype, trait,
                    ~ two_group_weeks(get_d("SNP1615", .x, .y), key_weeks))) %>%
  unnest(res) %>%
  group_by(ecotype) %>%
  mutate(p_adj = p.adjust(p_welch, "BH")) %>%      # 42 tests per ecotype
  ungroup() %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(p1615, "table_SNP1615_weekly_tests.csv", row.names = FALSE)

p2355 <- assoc_snp %>%
  filter(snp == "SNP2355", week %in% key_weeks) %>%
  group_by(trait, week) %>%
  group_modify(~ {
    a <- anova(lm(value ~ ecotype + genotype, data = .x))
    tibble(n = nrow(.x), F = a["genotype", "F value"],
           df1 = a["genotype", "Df"], df2 = a["Residuals", "Df"],
           p_ecotype = a["ecotype", "Pr(>F)"], p_genotype = a["genotype", "Pr(>F)"],
           eta2_partial = a["genotype", "Sum Sq"] /
                          (a["genotype", "Sum Sq"] + a["Residuals", "Sum Sq"]))
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p_genotype, "BH")) %>%   # 42 tests
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(p2355, "table_SNP2355_weekly_tests.csv", row.names = FALSE)

cat("SNP1615 tests:", nrow(p1615), "| raw p<.05:", sum(p1615$p_welch < 0.05),
    "| adjusted p<.05:", sum(p1615$p_adj < 0.05), "\n")
cat("SNP2355 tests:", nrow(p2355), "| raw p<.05:", sum(p2355$p_genotype < 0.05),
    "| adjusted p<.05:", sum(p2355$p_adj < 0.05), "\n")

## ------------------------------------------------------------
## 9. EXON SNPs: LINEAR MIXED-EFFECTS MODELS
## ------------------------------------------------------------
# Supporting analysis. The interaction p-values are optimistic with so few birds; see the exon document.
# Feeds: Exon document: Tables 13, 14 and A10
# ---- body weight (log scale): models for each SNP ----
w0 <- assoc_long %>% filter(trait == "Weight", week == 0) %>% select(bird_id, w0 = value)

dat_w <- assoc_long %>%
  filter(trait == "Weight", week >= 2, !is.na(value)) %>%
  left_join(w0, by = "bird_id") %>%
  mutate(logwt = log(value),
         wk    = (week - 12) / 10,                 # age centred on week 12
         w0c   = w0 - mean(w0, na.rm = TRUE),
         dose  = as.numeric(SNP2355) - 1,          # AA = 0, AG = 1, GG = 2
         het   = as.numeric(SNP2355 == "AG"))

d_fs <- dat_w %>% filter(ecotype == "Fulani", !is.na(SNP1615)) %>% droplevels()
d_ys <- dat_w %>% filter(ecotype == "Yoruba", !is.na(SNP1615)) %>% droplevels()
d_p  <- dat_w %>% filter(!is.na(SNP2355)) %>% droplevels()

m_fs  <- lmer(logwt ~ SNP1615 * (wk + I(wk^2)) + (1 | bird_id), data = d_fs)
m_ys  <- lmer(logwt ~ SNP1615 * (wk + I(wk^2)) + (1 | bird_id), data = d_ys)
m_add <- lmer(logwt ~ dose * (wk + I(wk^2)) + ecotype + (1 | bird_id), data = d_p)
m_dom <- update(m_add, . ~ . + het * (wk + I(wk^2)))   # adds dominance
models <- list("SNP1615 Fulani" = m_fs, "SNP1615 Yoruba" = m_ys, "SNP2355 pooled" = m_add)

get_anova <- function(m, label) {
  a <- as.data.frame(anova(m, type = 3))
  a$term <- rownames(a); a$model <- label
  a %>% select(model, term, NumDF, DenDF, `F value`, `Pr(>F)`)
}
aov_tab <- bind_rows(get_anova(m_fs,  "SNP1615 Fulani"),
                     get_anova(m_ys,  "SNP1615 Yoruba"),
                     get_anova(m_add, "SNP2355 pooled")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(aov_tab)
write.csv(aov_tab, "lmm_anova_weight.csv", row.names = FALSE)
print(sapply(models, isSingular))                    # TRUE = bird variance estimated as 0
print(anova(m_add, m_dom))                           # dominance test (SNP2355)

# genotype contrasts at weeks 12 and 24 (ratio of trait values)
contrast_two <- function(m, var) {
  map_dfr(c(12, 24), function(w) {
    emm <- emmeans(m, as.formula(paste("~", var)), at = list(wk = (w - 12) / 10))
    as.data.frame(summary(pairs(emm), infer = TRUE)) %>% mutate(week = w)
  }) %>% mutate(ratio = exp(estimate), ratio_lo = exp(lower.CL), ratio_hi = exp(upper.CL))
}
dose_two <- function(m) {                             # effect of one extra G allele
  map_dfr(c(12, 24), function(w) {
    emm <- emmeans(m, ~ dose, at = list(dose = c(0, 1), wk = (w - 12) / 10))
    as.data.frame(summary(pairs(emm), infer = TRUE)) %>% mutate(week = w)
  }) %>% transmute(week, contrast = "per G allele", ratio = exp(-estimate),
                   ratio_lo = exp(-upper.CL), ratio_hi = exp(-lower.CL), p.value)
}

# residual check for the weight models
par(mfrow = c(2, 2))
for (nm in names(models)) {
  plot(fitted(models[[nm]]), resid(models[[nm]]), main = nm,
       xlab = "Fitted", ylab = "Residual")
  abline(h = 0, lty = 2)
}
dev.copy(png, "fig_lmm_residuals_weight.png", width = 1600, height = 1600, res = 200)
dev.off()
par(mfrow = c(1, 1))

# ---- the same models for all seven traits ----
run_trait <- function(tr) {
  d <- assoc_long %>%
    filter(trait == tr, week >= 2, !is.na(value)) %>%
    mutate(logv = log(value), wk = (week - 12) / 10, dose = as.numeric(SNP2355) - 1)
  d_fs <- d %>% filter(ecotype == "Fulani", !is.na(SNP1615)) %>% droplevels()
  d_ys <- d %>% filter(ecotype == "Yoruba", !is.na(SNP1615)) %>% droplevels()
  d_p  <- d %>% filter(!is.na(SNP2355)) %>% droplevels()
  fits <- list(
    S1615_Fulani = lmer(logv ~ SNP1615 * (wk + I(wk^2)) + (1 | bird_id), data = d_fs),
    S1615_Yoruba = lmer(logv ~ SNP1615 * (wk + I(wk^2)) + (1 | bird_id), data = d_ys),
    S2355_pooled = lmer(logv ~ dose * (wk + I(wk^2)) + ecotype + (1 | bird_id), data = d_p))
  an <- map_dfr(names(fits), function(nm) {
    m <- fits[[nm]]
    a <- as.data.frame(anova(m, type = 3)); a$term <- rownames(a)
    a %>% filter(str_detect(term, "SNP1615|dose")) %>%
      transmute(trait = tr, model = nm, term, F_value = `F value`,
                p = `Pr(>F)`, singular = isSingular(m))
  })
  ef <- bind_rows(
    contrast_two(fits[["S1615_Fulani"]], "SNP1615") %>%
      transmute(trait = tr, model = "S1615_Fulani", week, contrast,
                ratio, ratio_lo, ratio_hi, p.value),
    contrast_two(fits[["S1615_Yoruba"]], "SNP1615") %>%
      transmute(trait = tr, model = "S1615_Yoruba", week, contrast,
                ratio, ratio_lo, ratio_hi, p.value),
    dose_two(fits[["S2355_pooled"]]) %>%
      mutate(trait = tr, model = "S2355_pooled", .before = 1))
  list(an = an, ef = ef)
}
safe_run <- function(tr) tryCatch(run_trait(tr),
  error = function(e) { message(tr, " FAILED: ", e$message); NULL })

lmm_out <- map(setNames(traits, traits), safe_run)

an_all <- map_dfr(lmm_out, "an") %>%
  mutate(term_type = case_when(str_detect(term, "I\\(wk") ~ "GxAge2",
                               str_detect(term, ":wk")    ~ "GxAge",
                               TRUE                       ~ "Genotype")) %>%
  group_by(model, term_type) %>%
  mutate(p_adj = p.adjust(p, "BH")) %>%               # across the 7 traits
  ungroup()
ef_all <- map_dfr(lmm_out, "ef")
write.csv(an_all, "lmm_all_traits_tests.csv", row.names = FALSE)
write.csv(ef_all, "lmm_all_traits_effects.csv", row.names = FALSE)

p_tab_compact <- an_all %>%
  mutate(m = recode(model, S1615_Fulani = "F1615", S1615_Yoruba = "Y1615",
                    S2355_pooled = "P2355"),
         t = recode(term_type, Genotype = "G", GxAge = "GxA", GxAge2 = "GxA2")) %>%
  select(trait, m, t, p) %>%
  pivot_wider(names_from = c(m, t), values_from = p) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(p_tab_compact), row.names = FALSE)

## ------------------------------------------------------------
## 10. EXON SNPs: SENSITIVITY ANALYSIS (Theil-Sen growth rates)
## ------------------------------------------------------------
# Repeats the primary analysis with growth rates that resist isolated erroneous readings (exon document, Table 12 equivalent; intron-document note).
# Feeds: Exon document: Table 12
# NOTE: NEW - this section has not yet been run in R
theil_sen <- function(x, y) {            # median of all pairwise slopes
  ij <- combn(length(x), 2)
  dx <- x[ij[2, ]] - x[ij[1, ]]
  dy <- y[ij[2, ]] - y[ij[1, ]]
  median(dy[dx != 0] / dx[dx != 0])
}

slopes_ts <- assoc_snp %>%
  filter(week >= 2, bird_id %in% complete_birds) %>%
  group_by(snp, ecotype, genotype, trait, bird_id) %>%
  filter(dplyr::n() >= 10) %>%
  summarise(slope = theil_sen(week, value), .groups = "drop")

bl_ts <- bird_level_tests(slopes_ts)
ts1615 <- bl_ts$s1615 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
ts2355 <- bl_ts$s2355 %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(ts1615), row.names = FALSE)
print(as.data.frame(ts2355), row.names = FALSE)
write.csv(ts1615, "sensitivity_theilsen_SNP1615.csv", row.names = FALSE)
write.csv(ts2355, "sensitivity_theilsen_SNP2355.csv", row.names = FALSE)

## ------------------------------------------------------------
## 11. INTRON SNPs: IMPORT AND DATA CHECK
## ------------------------------------------------------------
# Reads the intronic genotypes, matches bird IDs to the exon files (FA1 = F1, YA4 = Y4) and counts genotyped birds per SNP.
# Feeds: Intron document: Table 1
f_int <- if (file.exists(f_intron)) f_intron else
  list.files(pattern = "intron.*\\.xlsx$", ignore.case = TRUE)
if (length(f_int) == 0) f_int <- file.choose()

intr <- read_excel(f_int[1]) %>%
  rename(bird_raw = 1) %>%
  mutate(bird_id = str_replace(str_replace(bird_raw, "^FA", "F"), "^YA", "Y"),
         ecotype = if_else(str_starts(bird_raw, "FA"), "Fulani", "Yoruba"))
snps <- names(intr) %>% str_subset("^SNP")          # 16 intronic SNPs

# blank cells are treated as missing genotypes
qual <- intr %>%
  pivot_longer(all_of(snps), names_to = "snp", values_to = "genotype") %>%
  group_by(snp, ecotype) %>%
  summarise(n = sum(!is.na(genotype)),
            genotypes = paste(sort(unique(na.omit(genotype))), collapse = "/"),
            .groups = "drop") %>%
  pivot_wider(names_from = ecotype, values_from = c(n, genotypes)) %>%
  arrange(as.numeric(str_extract(snp, "\\d+")))
print(qual, n = 20, width = Inf)

## ------------------------------------------------------------
## 12. INTRON SNPs: DIVERSITY INDICES
## ------------------------------------------------------------
# Allele and genotype frequencies, MAF, Ho, He, Ne and PIC for any allele pair; summary across SNPs.
# Feeds: Intron document: Tables 1, 2 and A1
il <- intr %>%
  pivot_longer(all_of(snps), names_to = "snp", values_to = "genotype") %>%
  filter(!is.na(genotype)) %>%
  select(bird_id, ecotype, snp, genotype)
il_all <- bind_rows(il, mutate(il, ecotype = "Pooled"))

div_snp <- function(g) {
  n  <- length(g)
  a  <- unlist(strsplit(g, ""))
  tb <- table(g)
  fr <- sort(table(a) / length(a), decreasing = TRUE)
  p  <- as.numeric(fr)
  het <- sum(substr(g, 1, 1) != substr(g, 2, 2))
  He  <- 1 - sum(p^2)
  PIC <- He - if (length(p) > 1)
    sum(combn(length(p), 2, function(i) 2 * p[i[1]]^2 * p[i[2]]^2)) else 0
  tibble(n = n,
         genotypes = paste0(names(tb), "=", as.integer(tb), collapse = ", "),
         allele1 = names(fr)[1], f1 = p[1],
         allele2 = if (length(p) > 1) names(fr)[2] else NA_character_,
         f2 = if (length(p) > 1) p[2] else 0,
         MAF = if (length(p) > 1) min(p) else 0,
         Ho = het / n, He = He, Ne = 1 / sum(p^2), PIC = PIC)
}

intron_div <- il_all %>%
  group_by(snp, ecotype) %>%
  group_modify(~ div_snp(.x$genotype)) %>%
  ungroup() %>%
  mutate(pos = as.numeric(str_extract(snp, "\\d+")),
         ecotype = factor(ecotype, levels = c("Fulani", "Yoruba", "Pooled"))) %>%
  arrange(pos, ecotype) %>% select(-pos) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(intron_div %>% select(snp, ecotype, n, genotypes, allele1, f1, allele2, f2,
                            Ho, He, Ne, PIC), n = 60, width = Inf)
write.csv(intron_div, "intron_diversity_by_snp.csv", row.names = FALSE)

intron_summary <- intron_div %>%
  filter(n >= min_n) %>%
  group_by(ecotype) %>%
  summarise(n_snps = dplyr::n(), n_polymorphic = sum(MAF > 0),
            across(c(MAF, Ho, He, Ne, PIC), list(mean = mean, sd = sd),
                   .names = "{.col}_{.fn}"), .groups = "drop") %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(intron_summary), row.names = FALSE)
write.csv(intron_summary, "intron_diversity_summary.csv", row.names = FALSE)

## ------------------------------------------------------------
## 13. INTRON + EXON SNPs COMBINED: DIVERSITY INDICES
## ------------------------------------------------------------
# Repeats the indices for all 18 SNPs and summarises Exon, Intron and Combined sets.
# Feeds: Intron document: Tables 2 and A1
ex <- assoc_wide %>%
  rename(bird_id = `Bird ID`, ecotype = Ecotype) %>%
  select(bird_id, ecotype, SNP1615, SNP2355) %>%
  pivot_longer(c(SNP1615, SNP2355), names_to = "snp", values_to = "genotype") %>%
  mutate(genotype = as.character(genotype)) %>%
  filter(!is.na(genotype))
print(table(ex$snp, ex$ecotype))      # SNP1615: 6 and 7; SNP2355: 8 and 9

all_long  <- bind_rows(mutate(il, region = "Intron"), mutate(ex, region = "Exon"))
all_long2 <- bind_rows(all_long, mutate(all_long, ecotype = "Pooled"))

all_div <- all_long2 %>%
  group_by(region, snp, ecotype) %>%
  group_modify(~ div_snp(.x$genotype)) %>%
  ungroup() %>%
  mutate(pos = as.numeric(str_extract(snp, "\\d+")),
         ecotype = factor(ecotype, levels = c("Fulani", "Yoruba", "Pooled"))) %>%
  arrange(pos, ecotype) %>% select(-pos) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write.csv(all_div, "combined_diversity_by_snp.csv", row.names = FALSE)

sets <- list(Exon = "Exon", Intron = "Intron", Combined = c("Exon", "Intron"))
all_summary <- purrr::map_dfr(names(sets), function(s) {
  all_div %>%
    filter(region %in% sets[[s]], n >= min_n) %>%
    group_by(ecotype) %>%
    summarise(n_snps = dplyr::n(), n_polymorphic = sum(MAF > 0),
              across(c(MAF, Ho, He, Ne, PIC), list(mean = mean, sd = sd),
                     .names = "{.col}_{.fn}"), .groups = "drop") %>%
    mutate(snp_set = s, .before = 1)
}) %>% mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(all_summary), row.names = FALSE)
write.csv(all_summary, "combined_diversity_summary.csv", row.names = FALSE)

## ------------------------------------------------------------
## 14. FULANI vs YORUBA: PAIRED TESTS ON THE DIVERSITY INDICES
## ------------------------------------------------------------
# Paired t test and Wilcoxon test across the SNPs genotyped in both ecotypes (at least min_n birds in each).
# Feeds: Intron document: Table 4
wide <- all_div %>%
  filter(ecotype %in% c("Fulani", "Yoruba"), n >= min_n) %>%
  select(region, snp, ecotype, MAF, Ho, He, Ne, PIC) %>%
  pivot_wider(names_from = ecotype, values_from = c(MAF, Ho, He, Ne, PIC)) %>%
  drop_na()
cat("SNPs in the paired test:", nrow(wide), "\n"); print(wide$snp)

paired_test <- function(d, set_label) {
  map_dfr(c("MAF", "Ho", "He", "Ne", "PIC"), function(v) {
    x  <- d[[paste0(v, "_Fulani")]]; y <- d[[paste0(v, "_Yoruba")]]
    tt <- t.test(x, y, paired = TRUE)
    wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE))
    tibble(snp_set = set_label, index = v, n_snps = length(x),
           mean_Fulani = mean(x), mean_Yoruba = mean(y), diff = mean(x - y),
           ci_low = tt$conf.int[1], ci_high = tt$conf.int[2],
           t = unname(tt$statistic), df = unname(tt$parameter), p_t = tt$p.value,
           p_wilcoxon = wt$p.value, dz = mean(x - y) / sd(x - y))
  })
}
paired_res <- bind_rows(paired_test(wide, "Intron + exon"),
                        paired_test(filter(wide, region == "Intron"), "Intron only")) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
print(as.data.frame(paired_res), row.names = FALSE)
write.csv(paired_res, "paired_tests_diversity.csv", row.names = FALSE)

## ------------------------------------------------------------
## 15. GENOTYPE ANOVA ON GROWTH RATE (intron + exon SNPs)
## ------------------------------------------------------------
# One-way ANOVA of genotype on each bird's growth rate within ecotype, with a permutation check; Benjamini-Hochberg adjustment over all tests.
# Feeds: Intron document: Tables 5 and A3
gl <- assoc_long %>% select(bird_id, ecotype, trait, week, value)

rates <- gl %>%
  filter(week >= 2, bird_id %in% complete_birds, !is.na(value)) %>%
  group_by(bird_id, ecotype, trait) %>%
  filter(dplyr::n() >= 10) %>%
  summarise(rate = unname(coef(lm(value ~ week))[2]), .groups = "drop")
cat("birds used:", n_distinct(rates$bird_id), "\n")          # expect 19

# genotypes of all 18 SNPs (16 intronic + 2 exonic), one row per bird
geno_exon <- assoc_wide %>%
  rename(bird_id = `Bird ID`) %>%
  select(bird_id, SNP1615, SNP2355) %>%
  mutate(across(c(SNP1615, SNP2355), as.character))
geno_all <- intr %>% select(bird_id, all_of(snps)) %>% left_join(geno_exon, by = "bird_id")
snp_names_all <- c(snps, "SNP1615", "SNP2355")

F_stat <- function(y, g) {
  g <- droplevels(g); k <- nlevels(g); n <- length(y)
  gm <- tapply(y, g, mean); ng <- tapply(y, g, length)
  ssb <- sum(ng * (gm - mean(y))^2)
  ssw <- sum((y - gm[as.character(g)])^2)
  (ssb / (k - 1)) / (ssw / (n - k))
}

run_anova <- function(d) {
  d <- d %>% filter(!is.na(genotype), !is.na(rate))
  tab <- table(d$genotype)
  d <- d %>% filter(genotype %in% names(tab)[tab >= 2]) %>%   # classes with >= 2 birds
    mutate(genotype = factor(genotype))
  if (nlevels(d$genotype) < 2) return(NULL)
  k <- nlevels(d$genotype); n <- nrow(d)
  Fobs <- F_stat(d$rate, d$genotype)
  gm   <- tapply(d$rate, d$genotype, mean)
  ssb  <- sum(table(d$genotype) * (gm - mean(d$rate))^2)
  ssw  <- sum((d$rate - gm[as.character(d$genotype)])^2)
  set.seed(2026)
  perm <- replicate(9999, F_stat(d$rate, sample(d$genotype)))
  tibble(n = n,
         groups = paste0(levels(d$genotype), "=", as.integer(table(d$genotype)),
                         collapse = ", "),
         means  = paste0(round(gm, 3), collapse = " / "),
         F = Fobs, df1 = k - 1, df2 = n - k,
         p = pf(Fobs, k - 1, n - k, lower.tail = FALSE),
         eta2 = ssb / (ssb + ssw),
         p_perm = (1 + sum(perm >= Fobs - 1e-9)) / 10000)
}

anova_rate <- expand_grid(snp = snp_names_all, ecotype = c("Fulani", "Yoruba"),
                          trait = traits) %>%
  mutate(res = pmap(list(snp, ecotype, trait), function(s, e, t) {
    d <- rates %>% filter(ecotype == e, trait == t) %>%
      left_join(geno_all %>% select(bird_id, genotype = all_of(s)), by = "bird_id")
    run_anova(d)
  })) %>%
  unnest(res) %>%
  mutate(p_bh = p.adjust(p, "BH"))                   # across all tests

cat("tests run:", nrow(anova_rate), "| SNP x ecotype combinations:",
    nrow(distinct(anova_rate, snp, ecotype)), "\n")  # expect 140 and 20
cat("raw p < .05:", sum(anova_rate$p < .05), "| BH p < .05:", sum(anova_rate$p_bh < .05),
    "| smallest BH p:", round(min(anova_rate$p_bh), 3), "\n")
anova_rate_out <- anova_rate %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(anova_rate_out %>% filter(p < .05) %>% arrange(p) %>%
        select(snp, ecotype, trait, groups, F, df1, df2, p, eta2, p_perm, p_bh)),
      row.names = FALSE)
write.csv(anova_rate_out, "intron_exon_genotype_anova_growthrate.csv", row.names = FALSE)

## ------------------------------------------------------------
## 16. GENOTYPE ANOVA AT INDIVIDUAL AGES (intron + exon SNPs)
## ------------------------------------------------------------
# The same ANOVA on the trait value at ages 0 (weight only), 4, 8, 12, 16, 20 and 24; at least 3 residual degrees of freedom required.
# Feeds: Intron document: Tables 6 and A4
key_ages <- c(0, 4, 8, 12, 16, 20, 24)
aged <- gl %>%
  filter(!is.na(value), week %in% key_ages, !(week == 0 & trait != "Weight"))

anova_age <- function(d) {
  d <- d %>% filter(!is.na(genotype), !is.na(value))
  tab <- table(d$genotype)
  d <- d %>% filter(genotype %in% names(tab)[tab >= 2]) %>%
    mutate(genotype = factor(genotype))
  if (nlevels(d$genotype) < 2) return(NULL)
  k <- nlevels(d$genotype); n <- nrow(d)
  if (n - k < 3) return(NULL)                         # at least 3 residual df
  gm  <- tapply(d$value, d$genotype, mean)
  ssb <- sum(table(d$genotype) * (gm - mean(d$value))^2)
  ssw <- sum((d$value - gm[as.character(d$genotype)])^2)
  if (ssw == 0) return(NULL)
  Fobs <- (ssb / (k - 1)) / (ssw / (n - k))
  tibble(n = n,
         groups = paste0(levels(d$genotype), "=", as.integer(table(d$genotype)),
                         collapse = ", "),
         F = Fobs, df1 = k - 1, df2 = n - k,
         p = pf(Fobs, k - 1, n - k, lower.tail = FALSE),
         eta2 = ssb / (ssb + ssw))
}

age_tests <- expand_grid(snp = snp_names_all, ecotype = c("Fulani", "Yoruba"),
                         trait = traits, week = key_ages) %>%
  filter(!(week == 0 & trait != "Weight")) %>%
  mutate(res = pmap(list(snp, ecotype, trait, week), function(s, e, t, w) {
    d <- aged %>% filter(ecotype == e, trait == t, week == w) %>%
      left_join(geno_all %>% select(bird_id, genotype = all_of(s)), by = "bird_id")
    anova_age(d)
  })) %>%
  unnest(res) %>%
  mutate(p_bh = p.adjust(p, "BH"))

cat("tests run:", nrow(age_tests), "| combinations:",
    nrow(distinct(age_tests, snp, ecotype)), "\n")    # expect 862 and 22
cat("raw p < .05:", sum(age_tests$p < .05), "| BH p < .05:", sum(age_tests$p_bh < .05),
    "| smallest BH p:", round(min(age_tests$p_bh), 3), "\n")
print(table(age_tests$week))
age_tests_out <- age_tests %>% mutate(across(where(is.numeric), ~ round(.x, 4)))
print(as.data.frame(age_tests_out %>% arrange(p) %>% head(10) %>%
        select(snp, ecotype, trait, week, groups, F, df1, df2, p, eta2, p_bh)),
      row.names = FALSE)
write.csv(age_tests_out, "intron_exon_genotype_anova_by_age.csv", row.names = FALSE)

## ------------------------------------------------------------
## 17. HARDY-WEINBERG EQUILIBRIUM: all SNPs
## ------------------------------------------------------------
# Chi-square, likelihood-ratio and exact tests for every SNP and group with at least min_n birds that is polymorphic.
# Feeds: Intron document: Tables 3 and A2
hwe_one <- function(g) {
  a <- unlist(strsplit(g, ""))
  alle <- names(sort(table(a), decreasing = TRUE))
  if (length(alle) < 2) return(NULL)                  # monomorphic: nothing to test
  A <- alle[1]; B <- alle[2]
  nAA <- sum(g == paste0(A, A)); nBB <- sum(g == paste0(B, B))
  n   <- length(g); nAB <- n - nAA - nBB
  p   <- (2 * nAA + nAB) / (2 * n); q <- 1 - p
  obs <- c(nAA, nAB, nBB)
  expected <- n * c(p^2, 2 * p * q, q^2)
  min_exp  <- min(expected)               # computed before expected becomes text
  chi <- sum((obs - expected)^2 / expected)
  G   <- 2 * sum(ifelse(obs > 0, obs * log(obs / expected), 0))
  tibble(n = n,
         obs = paste(obs, collapse = "/"),
         expected = paste(round(expected, 2), collapse = "/"),
         chi2 = chi,  p_chi2 = pchisq(chi, 1, lower.tail = FALSE),
         G2 = G,      p_LRT  = pchisq(G,   1, lower.tail = FALSE),
         p_exact = hwe_exact(nAA, nAB, nBB),
         min_expected = min_exp,
         F_IS = 1 - (nAB / n) / (2 * p * q))
}

hwe_res <- all_long2 %>%
  group_by(region, snp, ecotype) %>%
  group_modify(~ {
    g <- .x$genotype
    if (length(g) < min_n) return(tibble())
    out1 <- hwe_one(g)
    if (is.null(out1)) tibble() else out1
  }) %>%
  ungroup() %>%
  mutate(fam = paste(region, if_else(ecotype == "Pooled", "pooled", "by ecotype"))) %>%
  group_by(fam) %>%
  mutate(p_exact_BH = p.adjust(p_exact, "BH")) %>%    # within each family of tests
  ungroup() %>%
  mutate(pos = as.numeric(str_extract(snp, "\\d+")),
         ecotype = factor(ecotype, levels = c("Fulani", "Yoruba", "Pooled"))) %>%
  arrange(region, pos, ecotype) %>% select(-pos) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

print(as.data.frame(hwe_res %>% filter(region == "Intron", ecotype != "Pooled") %>%
        select(snp, ecotype, n, obs, expected, chi2, p_chi2, G2, p_LRT, p_exact,
               min_expected, F_IS, p_exact_BH)), row.names = FALSE)
print(as.data.frame(hwe_res %>% group_by(fam) %>%
        summarise(tests = dplyr::n(), exact_p_lt_05 = sum(p_exact < .05),
                  chi2_p_lt_05 = sum(p_chi2 < .05), LRT_p_lt_05 = sum(p_LRT < .05),
                  min_exact_p = min(p_exact), .groups = "drop")), row.names = FALSE)
write.csv(hwe_res, "hwe_intron_exon.csv", row.names = FALSE)

## ------------------------------------------------------------
## 18. FIGURES
## ------------------------------------------------------------
# R versions of the figures in the two Word documents (the documents' figures were drawn with matplotlib).
# Feeds: Figures of both Results documents
# NOTE: NEW - this section has not yet been run in R
th <- theme_classic(base_size = 11)

# Figure 1 (exon document): mean body weight by ecotype, full population
f1 <- desc %>% filter(trait == "Weight") %>%
  ggplot(aes(week, mean, linetype = ecotype, shape = ecotype)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  labs(x = "Age (weeks)", y = "Body weight (g)",
       linetype = "Ecotype", shape = "Ecotype") + th
ggsave("fig1_weight_population.png", f1, width = 6.2, height = 4.2, dpi = 300)

# Figure 2: the six linear measurements
f2 <- desc %>% filter(trait != "Weight") %>%
  ggplot(aes(week, mean, linetype = ecotype, shape = ecotype)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  facet_wrap(~ trait, scales = "free_y") +
  labs(x = "Age (weeks)", y = "Mean (cm)", linetype = "Ecotype", shape = "Ecotype") + th
ggsave("fig2_linear_traits.png", f2, width = 8, height = 5, dpi = 300)

# Figure 3: genotype counts of the exon SNPs
f3 <- snp_long %>% filter(!is.na(genotype)) %>% count(snp, ecotype, genotype) %>%
  ggplot(aes(genotype, n, fill = ecotype)) +
  geom_col(position = position_dodge(preserve = "single")) +
  facet_wrap(~ snp) +
  scale_fill_grey(start = 0.2, end = 0.7) +
  labs(x = "Genotype", y = "Number of birds", fill = "Ecotype") + th
ggsave("fig3_genotype_counts.png", f3, width = 7, height = 3.4, dpi = 300)

# Figure 4: body weight by genotype (mean +- SE)
f4 <- geno_means %>% filter(trait == "Weight") %>%
  ggplot(aes(week, mean, linetype = genotype, shape = genotype)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  facet_grid(snp ~ ecotype) +
  labs(x = "Age (weeks)", y = "Body weight (g)",
       linetype = "Genotype", shape = "Genotype") + th
ggsave("fig4_weight_by_genotype.png", f4, width = 7.4, height = 6, dpi = 300)

# Figure 5: shank length by SNP1615 genotype in Yoruba
f5 <- geno_means %>%
  filter(trait == "ShankLength", snp == "SNP1615", ecotype == "Yoruba") %>%
  ggplot(aes(week, mean, linetype = genotype, shape = genotype)) +
  geom_line() + geom_point() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.4) +
  labs(x = "Age (weeks)", y = "Shank length (cm)",
       linetype = "SNP1615", shape = "SNP1615") + th
ggsave("fig5_shank_snp1615_yoruba.png", f5, width = 6, height = 4, dpi = 300)

# Figure 1 (intron document): expected heterozygosity of every SNP
snp_order <- all_div %>% distinct(snp) %>%
  mutate(pos = as.numeric(str_extract(snp, "\\d+"))) %>% arrange(pos) %>% pull(snp)
fi1 <- all_div %>%
  filter(ecotype %in% c("Fulani", "Yoruba"), n >= 3) %>%
  mutate(snp = factor(snp, levels = snp_order), enough = n >= 5) %>%
  ggplot(aes(snp, He, shape = ecotype, alpha = enough)) +
  geom_point(size = 2.5, position = position_dodge(width = 0.4)) +
  scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.4), guide = "none") +
  labs(x = "SNP", y = "Expected heterozygosity (He)", shape = "Ecotype") + th +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
ggsave("figI1_he_by_snp.png", fi1, width = 7.2, height = 3.8, dpi = 300)

## ------------------------------------------------------------
## 19. SESSION INFORMATION
## ------------------------------------------------------------
# Records the R and package versions used, for reproducibility.
sessionInfo()

# ==================================================
## END OF SCRIPT
## ==================================================
