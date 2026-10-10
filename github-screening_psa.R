# Six-cancer screening: deaths potentially averted and mortality-risk-only PSA
# PSA varies stage-specific 5-year index-cancer death risks; other inputs are fixed.
# Input: stage.data.xlsx in the working directory (or edit input_file).

library(readxl)
library(dplyr)
library(tidyr)
library(openxlsx)

input_file <- "8-stage.data-暴雨.xlsx"
n_sim <- 1000L
set.seed(123)
stages <- c("I", "II", "III", "IV")
scenarios <- c("China", "Multi", "China_Adherence", "Multi_Adherence")

# Stage distributions from routine practice, Chinese screening and multi-country screening.
raw <- read_excel(input_file)
stopifnot(all(c("Cancer", "Group", "Stage", "Percentage") %in% names(raw)))
dist <- raw %>%
  filter(Group %in% c("Non-screening", "Chinese post-screening", "Multi-country post-screening")) %>%
  mutate(Cancer = recode(Cancer, "Colorectal" = "Colorectum", "Esophageal" = "Oesophageal"),
         Group = recode(Group, "Non-screening" = "Patient",
                        "Chinese post-screening" = "China", "Multi-country post-screening" = "Multi"),
         Stage = as.character(Stage), Percentage = as.numeric(Percentage)) %>%
  select(Cancer, Group, Stage, Percentage)

# Parameters in cancer order: cases, stage-specific risks and SEs, observed adherence.
params <- data.frame(
  Cancer = c("Oesophageal", "Stomach", "Colorectum", "Liver", "Lung", "Female breast"),
  Num_Can = c(142213, 241318, 410087, 212249, 412115, 319230),
  DR_I = c(.170124481, .099152542, .057630736, .313653137, .076762402, .027518172),
  DR_II = c(.292622951, .223456790, .107459572, .374631268, .261718750, .059295331),
  DR_III = c(.498831776, .394371780, .257269591, .544753086, .564563463, .185448092),
  DR_IV = c(.751269036, .692798541, .677595628, .753694581, .694253859, .545454545),
  DR_I_SE = c(.013985191, .008704307, .007617452, .028252682, .006085071, .003728532),
  DR_II_SE = c(.013032458, .014647284, .007075346, .026342963, .010387366, .003997853),
  DR_III_SE = c(.012089388, .009732152, .009707284, .019590413, .009498809, .011583194),
  DR_IV_SE = c(.021863717, .013945197, .014119727, .017501765, .006746784, .032126753),
  Adherence_Rate = c(.358, .358, .1902, .7406, .330, .5136)
)

risk_cols <- paste0("DR_", stages)
se_cols <- paste0(risk_cols, "_SE")
stopifnot(all(as.matrix(params[risk_cols]) > 0 & as.matrix(params[risk_cols]) < 1),
          all(is.finite(as.matrix(params[se_cols]))), all(as.matrix(params[se_cols]) >= 0),
          all(params$Adherence_Rate >= 0 & params$Adherence_Rate <= 1))

# Require a complete four-stage distribution for every cancer and group.
expected <- expand.grid(Cancer = params$Cancer, Group = c("Patient", "China", "Multi"),
                        Stage = stages, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(dist[c("Cancer", "Group", "Stage")]),
          nrow(anti_join(expected, dist, by = c("Cancer", "Group", "Stage"))) == 0,
          all(is.finite(dist$Percentage)), all(dist$Percentage >= 0))

get_p <- function(cancer, group, normalise = TRUE) {
  x <- dist %>% filter(Cancer == cancer, Group == group) %>%
    arrange(match(Stage, stages)) %>% pull(Percentage)
  stopifnot(length(x) == 4L, sum(x) > 0)
  if (normalise) x / sum(x) else x  # Original PSA normalised; original point estimates did not.
}

# Logit-normal sampling keeps simulated risks strictly within (0, 1).
draw_risks <- function(p, se) {
  ifelse(se == 0, p, plogis(rnorm(length(p), qlogis(p), se / (p * (1 - p)))))
}

# Return the four averted-death outcomes for a given risk vector.
averted <- function(N, p0, pc, pm, adherence, risks, round_point = FALSE) {
  deaths <- function(p) N * sum(p * risks)
  if (!round_point) {
    baseline <- deaths(p0)
    return(c(China = baseline - deaths(pc), Multi = baseline - deaths(pm),
             China_Adherence = adherence * (baseline - deaths(pc)),
             Multi_Adherence = adherence * (baseline - deaths(pm))))
  }
  # Original main-analysis rounding: round each death component before subtraction.
  baseline <- round(deaths(p0))
  china <- round(deaths(pc)); multi <- round(deaths(pm))
  china_adh <- round(adherence * deaths(pc)) + round((1 - adherence) * deaths(p0))
  multi_adh <- round(adherence * deaths(pm)) + round((1 - adherence) * deaths(p0))
  c(China = baseline - china, Multi = baseline - multi,
    China_Adherence = baseline - china_adh, Multi_Adherence = baseline - multi_adh)
}

# For each simulation, sum the six cancer outcomes BEFORE estimating Overall UI.
cancers <- params$Cancer
sim <- array(NA_real_, c(n_sim, length(cancers), length(scenarios)),
             dimnames = list(NULL, cancers, scenarios))
points <- matrix(NA_real_, length(cancers), length(scenarios),
                 dimnames = list(cancers, scenarios))
baseline_deaths <- numeric(length(cancers))

for (i in seq_along(cancers)) {
  d <- params[i, ]
  p0 <- get_p(d$Cancer, "Patient")
  pc <- get_p(d$Cancer, "China")
  pm <- get_p(d$Cancer, "Multi")
  risks <- unlist(d[1, risk_cols], use.names = FALSE)
  ses <- unlist(d[1, se_cols], use.names = FALSE)
  baseline_deaths[i] <- round(d$Num_Can * sum(p0 * risks))
  points[i, ] <- averted(d$Num_Can, get_p(d$Cancer, "Patient", FALSE),
                         get_p(d$Cancer, "China", FALSE), get_p(d$Cancer, "Multi", FALSE),
                         d$Adherence_Rate, risks, TRUE)
  for (b in seq_len(n_sim)) {
    sim[b, i, ] <- averted(d$Num_Can, p0, pc, pm, d$Adherence_Rate,
                          draw_risks(risks, ses))
  }
}

summarise_sim <- function(x, cancer, point, baseline) {
  data.frame(Cancer = cancer, Scenario = scenarios,
             Cases = if (cancer == "Overall") sum(params$Num_Can) else params$Num_Can[match(cancer, cancers)],
             Baseline_Deaths = baseline, Averted_Deaths = as.numeric(point),
             Death_Reduction_Pct = 100 * as.numeric(point) / baseline,
             PSA_Mean = colMeans(x),
             UI_Lower = apply(x, 2, quantile, probs = .025, names = FALSE),
             UI_Upper = apply(x, 2, quantile, probs = .975, names = FALSE))
}

results <- bind_rows(lapply(seq_along(cancers), function(i) {
  x <- matrix(sim[, i, ], ncol = length(scenarios), dimnames = list(NULL, scenarios))
  summarise_sim(x, cancers[i], points[i, ], baseline_deaths[i])
}))
overall_sim <- apply(sim, c(1, 3), sum)
overall <- summarise_sim(overall_sim, "Overall", colSums(points), sum(baseline_deaths))
results <- bind_rows(results, overall)

# Export editable results; keep intervals numeric for downstream figures/tables.
write.xlsx(list(Results = results, Parameters = params),
           file = "screening_psa_results.xlsx", overwrite = TRUE)
print(results, row.names = FALSE)
