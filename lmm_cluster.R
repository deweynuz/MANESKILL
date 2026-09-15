# MANESKILL - cluster-aware linear mixed models
# Called by maneskill_analysis.ipynb. Do not run directly unless analysis.csv exists.
#
# Usage : Rscript lmm_cluster.R <analysis.csv> <output_dir> [prefix]
# Input : one row per participant with columns arm, team, date and the outcome variables
# Output: <output_dir>/<prefix>lmm_icc.csv       ICC (null model and adjusted for arm)
#         <output_dir>/<prefix>lmm_anova.csv     omnibus F test for arm (Kenward-Roger df)
#         <output_dir>/<prefix>lmm_emm.csv       adjusted means (EMM) per arm
#         <output_dir>/<prefix>lmm_contrasts.csv pairwise contrasts, Tukey-adjusted and raw
# Requires: lme4, lmerTest, pbkrtest, emmeans

suppressMessages({library(lme4); library(lmerTest); library(emmeans); library(pbkrtest)})
args <- commandArgs(trailingOnly = TRUE)
IN  <- if (length(args) >= 1) args[1] else "analysis.csv"
DIR <- if (length(args) >= 2) args[2] else "."
PFX <- if (length(args) >= 3) args[3] else ""
dir.create(DIR, showWarnings = FALSE, recursive = TRUE)
out <- function(name) file.path(DIR, paste0(PFX, name))

d <- read.csv(IN)
d$arm  <- factor(d$arm, levels = c("Control", "Passive", "Active"))
d$team <- factor(d$team)
d$date <- factor(d$date)
emm_options(lmer.df = "kenward-roger")

# ---- Intraclass correlation: null model (used for the design effect) and adjusted for arm ----
icc_fun <- function(y, adjusted = FALSE) {
  rhs <- if (adjusted) "arm + (1|team)" else "1 + (1|team)"
  m <- lmer(as.formula(paste(y, "~", rhs)), d, REML = TRUE)
  v <- as.data.frame(VarCorr(m)); v$vcov[1] / sum(v$vcov)
}
icc <- data.frame(outcome = c("cal_ants", "game_ants", "delta_ants", "game_stress", "pdi"))
icc$icc          <- sapply(icc$outcome, icc_fun)
icc$icc_adjusted <- sapply(icc$outcome, icc_fun, adjusted = TRUE)
write.csv(icc, out("lmm_icc.csv"), row.names = FALSE)

# ---- Mixed models: arm as fixed effect, random intercept for session (team) ----
anova_rows <- list(); emm_rows <- list(); con_rows <- list()
contrast_table <- function(e, label) {
  pt <- as.data.frame(summary(pairs(e, adjust = "tukey"), infer = TRUE))
  pr <- as.data.frame(summary(pairs(e, adjust = "none"), infer = TRUE))
  # emmeans reports "first - second"; flip the sign to report "second minus first" (e.g. Active vs Control)
  cs <- data.frame(outcome = label, contrast = as.character(pt$contrast), estimate = -pt$estimate,
                   lower = -pt$upper.CL, upper = -pt$lower.CL, df = pt$df, t = -pt$t.ratio, p_tukey = pt$p.value,
                   lower_raw = -pr$upper.CL, upper_raw = -pr$lower.CL, p_raw = pr$p.value)
  cs$contrast <- sapply(strsplit(cs$contrast, " - "), function(z) paste(z[2], "vs", z[1]))
  cs
}
fit <- function(y, cov = NULL, label = y, extra_random = NULL) {
  rhs <- paste(c("arm", cov, "(1|team)", extra_random), collapse = " + ")
  m <- lmer(as.formula(paste(y, "~", rhs)), d, REML = TRUE)
  a <- anova(m, ddf = "Kenward-Roger")
  anova_rows[[length(anova_rows) + 1]] <<- data.frame(outcome = label, F = a["arm", "F value"], df1 = a["arm", "NumDF"],
                                                      df2 = a["arm", "DenDF"], p = a["arm", "Pr(>F)"])
  e <- emmeans(m, "arm"); es <- as.data.frame(e); es$outcome <- label
  emm_rows[[length(emm_rows) + 1]] <<- es
  con_rows[[length(con_rows) + 1]] <<- contrast_table(e, label)
  invisible(m)
}
# Primary outcome: ANCOVA (phase 2 adjusted for phase 1) and change score
fit("game_ants",  "cal_ants", "ANTS total (phase 2, adjusted for phase 1)")
fit("delta_ants", NULL,       "ANTS total (change score)")
# ANTS domains: ANCOVA (Table 2) and change score (Figure 1)
fit("game_sa", "cal_sa", "Situation awareness (phase 2, adjusted)")
fit("game_dm", "cal_dm", "Decision making (phase 2, adjusted)")
fit("game_tw", "cal_tw", "Team working (phase 2, adjusted)")
fit("game_tm", "cal_tm", "Task management (phase 2, adjusted)")
fit("delta_sa", NULL, "Situation awareness (change score)")
fit("delta_dm", NULL, "Decision making (change score)")
fit("delta_tw", NULL, "Team working (change score)")
fit("delta_tm", NULL, "Task management (change score)")
# Secondary individual outcomes
fit("game_stress",  "cal_stress", "Stress NRS (phase 2, adjusted)")
fit("delta_stress", NULL,         "Stress NRS (change score)")
fit("pdi",     NULL, "Peritraumatic Distress Inventory")
fit("self_sa", NULL, "Self-rated situation awareness")
fit("self_dm", NULL, "Self-rated decision making")
fit("self_tw", NULL, "Self-rated team working")
fit("self_tm", NULL, "Self-rated task management")
# Sensitivity: additional random intercept for the session date (proxy for site)
fit("game_ants", "cal_ants", "ANTS total (sensitivity: + date random effect)", extra_random = "(1|date)")

write.csv(do.call(rbind, anova_rows), out("lmm_anova.csv"), row.names = FALSE)
write.csv(do.call(rbind, emm_rows),   out("lmm_emm.csv"),   row.names = FALSE)
write.csv(do.call(rbind, con_rows),   out("lmm_contrasts.csv"), row.names = FALSE)
cat("Mixed models done:", nrow(d), "participants,", nlevels(d$team), "sessions ->", DIR, "\n")
