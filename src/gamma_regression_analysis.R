# =============================================================================
# Bayesian Gamma Regression of Glycosylated Hemoglobin (GHB)
# Computational Statistics - MSc Data Science and Engineering, IST (2025/26)
# Authors: Joao Pires, Gustavo Brazil
# =============================================================================

#----------------------Alignment and Data Treatment-----------------------------

library(tidyverse)

# Run from the repository root, e.g. Rscript src/gamma_regression_analysis.R
data_path <- "data/diabetes.txt"
dir.create("figures", showWarnings = FALSE)

raw_data <- read.table(data_path, header = FALSE, fill = TRUE)

colnames(raw_data) <- c("id", "chol", "stab.glu", "hdl", "glyhb", "location",
                     "age", "gender", "height", "weight", "frame",
                     "bp.1s", "bp.1d", "waist", "hip")
head(raw_data)

diabetes <- raw_data %>%
  mutate(
    ID = as.character(id),
    GHB = as.numeric(glyhb),
    SGLU = as.numeric(stab.glu),
    CHOL = as.numeric(chol),
    HDL = as.numeric(hdl),
    AGE = as.numeric(age),
    GENDER = as.factor(gender),
    LOCATION = as.factor(location),
    HHT = as.numeric(height),
    WHT = as.numeric(weight),
    FRAME = as.factor(frame),
    SBP = as.numeric(bp.1s),
    DSP = as.numeric(bp.1d),
    WAIST = as.numeric(waist),
    HIP = as.numeric(hip)
  ) %>%

  select(ID, GHB, SGLU, CHOL, HDL, AGE, GENDER, LOCATION, HHT, WHT, FRAME, SBP, DSP, WAIST, HIP) %>%
  na.omit()

print(paste("Complete observations:", nrow(diabetes)))
print(head(diabetes))

#--------------Validation of the Gamma Hypothesis (Exploratory Analysis)--------

library(ggplot2)


p1 <- ggplot(diabetes, aes(x = GHB)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "#69b3a2", color = "black", alpha = 0.7) +
  geom_density(color = "red", lwd = 1.2) +
  labs(
    title = "Distribution of Glycosolated Hemoglobin (GHB)",
    subtitle = "Right-skewness supports the choice of a Gamma Model",
    x = "Glycosolated Hemoglobin (GHB)",
    y = "Density"
  ) +
  theme_minimal() +
  theme(

    plot.title = element_text(size = 18, face = "bold"),

    plot.subtitle = element_text(size = 14),

    axis.title = element_text(size = 15, face = "bold"),

    axis.text = element_text(size = 12)
  )

print(p1)


p2 <- ggplot(diabetes, aes(x = SGLU, y = GHB, color = GENDER)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_smooth(method = "glm", method.args = list(family = Gamma(link="log")), se = FALSE, lwd = 1.5) +
  labs(
    title = "Relationship: GHB vs. SGLU",
    subtitle = "Positive correlation observed with Gamma fit",
    x = "Stabilized Glucose (SGLU)",
    y = "Glycosolated Hemoglobin (GHB)",
    color = "Gender" # Nome da legenda
  ) +
  scale_color_manual(values = c("male" = "blue", "female" = "red"),
                     labels = c("Male", "Female")) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 15, face = "bold"),
    axis.text = element_text(size = 12),

    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12)
  )

print(p2)


ggsave("figures/histogram_ghb.png",
       plot = p1,
       width = 8,
       height = 4,
       dpi = 300)

ggsave("figures/scatterplot_glucose.png",
       plot = p2,
       width = 8,
       height = 6,
       dpi = 300)


#------------Selection of Risk Factors (Correlations)---------------------------

vars_to_plot <- diabetes %>%
  select(GHB, SGLU, AGE, CHOL, WAIST)

panel.cor_large <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)


  final_cex <- 2 + (r * 1.5)

  text(0.5, 0.5, txt, cex = final_cex, font = 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#69b3a2", border="white")
}

panel.dots <- function(x, y, ...) {
  points(x, y, pch = 19, col = rgb(0, 0, 1, 0.2), cex = 0.8)
  abline(lm(y~x), col="red", lwd=1.5)
}


png("figures/correlation_matrix_final_large.png", width = 1200, height = 1000, res = 150)

pairs(vars_to_plot,
      lower.panel = panel.dots,
      upper.panel = panel.cor_large,
      diag.panel = panel.hist,
      main = "Correlation Matrix of Key Variables",
      cex.main = 1.5,
      gap = 0.5)

dev.off()

pairs(vars_to_plot,
      lower.panel = panel.dots,
      upper.panel = panel.cor_large,
      diag.panel = panel.hist,
      main = "Correlation Matrix of Key Variables")

while (!is.null(dev.list())) dev.off()

vars_all_numeric <- diabetes %>%
  select(GHB, SGLU, CHOL, HDL, AGE, WAIST, HIP, HHT, WHT, SBP, DSP)

panel.cor_large_full <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {

  par(usr = c(0, 1, 0, 1))

  r <- abs(cor(x, y, use="complete.obs"))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)

  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)

  final_cex <- 1.5 + (r * 2)

  text(0.5, 0.5, txt, cex = final_cex, font = 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr"); on.exit(par(usr = usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks; nB <- length(breaks)
  y <- h$counts; y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "#69b3a2", border="white")
}


panel.dots <- function(x, y, ...) {
  points(x, y, pch = 19, col = rgb(0, 0, 1, 0.2), cex=0.5)
  abline(lm(y~x), col="red", lwd=1)
}

png("figures/correlation_matrix_full_large.png", width = 2000, height = 2000, res = 200)

pairs(vars_all_numeric,
      lower.panel = panel.dots,
      upper.panel = panel.cor_large_full,
      diag.panel = panel.hist,
      main = "Full Correlation Matrix (All Numeric Variables)",
      cex.main = 2,
      gap = 0)

dev.off()

cat("Plot saved: figures/correlation_matrix_full_large.png\n")


#------------------Analysis of Demographic Groups ------------------------------

library(ggplot2)
library(patchwork)

theme_large <- theme_minimal() +
  theme(

    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),


    axis.title = element_text(size = 20, face = "bold"),


    axis.text = element_text(size = 18, color = "black"),


    panel.grid.minor = element_blank()
  )

p_gender <- ggplot(diabetes, aes(x = GENDER, y = GHB, fill = GENDER)) +
  geom_boxplot(alpha = 0.7, lwd = 1, outlier.size = 3) +
  scale_fill_manual(values = c("female"="#E69F00", "male"="#56B4E9")) +
  labs(title = "GHB by Gender", x = "", y = "GHB Levels") +
  theme_large +
  theme(legend.position = "none")

p_loc <- ggplot(diabetes, aes(x = LOCATION, y = GHB, fill = LOCATION)) +
  geom_boxplot(alpha = 0.7, lwd = 1, outlier.size = 3) +
  scale_fill_manual(values = c("Buckingham"="#009E73", "Louisa"="#CC79A7")) +
  labs(title = "GHB by Location", x = "", y = "") +
  theme_large +
  theme(legend.position = "none")

if(require(patchwork)) {
  p_combined <- p_gender + p_loc

  ggsave("figures/boxplots_categorical.png", plot = p_combined, width = 14, height = 7, dpi = 300)

} else {
  ggsave("figures/boxplot_gender.png", plot = p_gender, width = 8, height = 7, dpi = 300)
  ggsave("figures/boxplot_location.png", plot = p_loc, width = 8, height = 7, dpi = 300)
}

print(p_gender + p_loc)

#----------------MCMC Modeling (The Gamma Model) -------------------------------

library(rjags)

jags_data <- list(
  N = nrow(diabetes),
  y = diabetes$GHB,
  sglu = diabetes$SGLU,
  age = diabetes$AGE,
  chol = diabetes$CHOL,
  waist = diabetes$WAIST,
  sbp = diabetes$SBP
)

model_string <- "model {
    for (i in 1:N) {
        # Likelihood: Gamma distribution
        y[i] ~ dgamma(shape, rate[i])

        # Link Function: Log
        log(mu[i]) <- beta0 +
                      beta_sglu * sglu[i] +
                      beta_age * age[i] +
                      beta_chol * chol[i] +
                      beta_waist * waist[i] +
                      beta_sbp * sbp[i]

        # Rate definition (JAGS uses the shape/rate parameterisation)
        rate[i] <- shape / mu[i]
    }

    # --- PRIORS (as described in the report) ---

    # 1. Coefficients (vague Normal: mean 0, precision 0.001)
    beta0 ~ dnorm(0, 0.001)
    beta_sglu ~ dnorm(0, 0.001)
    beta_age ~ dnorm(0, 0.001)
    beta_chol ~ dnorm(0, 0.001)
    beta_waist ~ dnorm(0, 0.001)
    beta_sbp ~ dnorm(0, 0.001)

    # 2. Shape parameter (Uniform on 0-100)
    # Note: called 'r' in the report.
    shape ~ dunif(0, 100)
}"


model <- jags.model(textConnection(model_string),
                    data = jags_data,
                    n.chains = 3)

update(model, 2000)

samples <- coda.samples(model,
                        variable.names = c("beta0", "beta_sglu", "beta_age",
                                           "beta_chol", "beta_waist", "beta_sbp",
                                           "shape"),
                        n.iter = 5000)

summary(samples)
gelman.diag(samples)
plot(samples)

dic_res <- dic.samples(model, n.iter = 5000, type = "pD")

print(dic_res)

#------------Robust Validation (Bootstrap)--------------------------------------

library(boot)

glm_coef_function <- function(data, indices) {

  d <- data[indices, ]

  fit <- glm(GHB ~ SGLU + AGE + CHOL + WAIST + SBP,
             family = Gamma(link = "log"),
             data = d)

  return(coef(fit))
}


set.seed(123)
boot_results <- boot(data = diabetes, statistic = glm_coef_function, R = 1000)


print(boot_results)


cat("\n--- Bootstrap CI: Glucose (Index 2) ---\n")
boot.ci(boot_results, type = "perc", index = 2)

cat("\n--- Bootstrap CI: Age (Index 3) ---\n")
boot.ci(boot_results, type = "perc", index = 3)

modelo_glm_final <- glm(GHB ~ SGLU + AGE + CHOL + WAIST + SBP,
                        family = Gamma(link = "log"),
                        data = diabetes)

cat("\n--- AIC of the frequentist GLM ---\n")
print(AIC(modelo_glm_final))


#---------------Model Selection Criteria (DIC): full 14-variable model----------

library(rjags)

jags_data_full <- list(
  N = nrow(diabetes),
  y = diabetes$GHB,

  sglu = diabetes$SGLU,
  chol = diabetes$CHOL,
  hdl = diabetes$HDL,
  age = diabetes$AGE,
  height = diabetes$HHT,
  weight = diabetes$WHT,
  sbp = diabetes$SBP,
  dsp = diabetes$DSP,
  waist = diabetes$WAIST,
  hip = diabetes$HIP,

  gender = as.numeric(diabetes$GENDER),
  location = as.numeric(diabetes$LOCATION),
  frame = as.numeric(diabetes$FRAME)
)

model_string_full <- "model {
    for (i in 1:N) {
        # Likelihood
        y[i] ~ dgamma(shape, rate[i])

        # Link function: log (all 14 candidate predictors)
        log(mu[i]) <- beta0 +
                      beta_sglu * sglu[i] +
                      beta_chol * chol[i] +
                      beta_hdl * hdl[i] +
                      beta_age * age[i] +
                      beta_height * height[i] +
                      beta_weight * weight[i] +
                      beta_sbp * sbp[i] +
                      beta_dsp * dsp[i] +
                      beta_waist * waist[i] +
                      beta_hip * hip[i] +
                      beta_gender * gender[i] +
                      beta_location * location[i] +
                      beta_frame * frame[i]

        rate[i] <- shape / mu[i]
    }

    # --- PRIORS (14 betas + intercept + shape) ---
    beta0 ~ dnorm(0, 0.001)

    beta_sglu ~ dnorm(0, 0.001)
    beta_chol ~ dnorm(0, 0.001)
    beta_hdl ~ dnorm(0, 0.001)
    beta_age ~ dnorm(0, 0.001)
    beta_height ~ dnorm(0, 0.001)
    beta_weight ~ dnorm(0, 0.001)
    beta_sbp ~ dnorm(0, 0.001)
    beta_dsp ~ dnorm(0, 0.001)
    beta_waist ~ dnorm(0, 0.001)
    beta_hip ~ dnorm(0, 0.001)
    beta_gender ~ dnorm(0, 0.001)
    beta_location ~ dnorm(0, 0.001)
    beta_frame ~ dnorm(0, 0.001)

    shape ~ dunif(0, 100)
}"

model_full <- jags.model(textConnection(model_string_full),
                         data = jags_data_full,
                         n.chains = 3)

update(model_full, 2000)

dic_full <- dic.samples(model_full, n.iter = 5000, type = "pD")

print("--- DIC (FULL MODEL) ---")
print(dic_full)
