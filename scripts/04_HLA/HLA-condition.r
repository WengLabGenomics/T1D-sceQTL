res_df = read.csv('./HLA_res/disease_sig_hla_loci.csv',row.names = 1)
res_df = res_df[res_df$p < 5e-8,]
hla_vars <- res_df$name

merge_df <- read.csv('./HLA_matrix/Discovery_cohor_HLA_dosage.csv',row.names = 1)
hla_df = merge_df[c(hla_vars,'sex','province')]
hla_df$Disease = merge_df2$disease_status

candidate <- hla_vars
selected <- character(0)

step_results <- list()      
all_step_results <- list()   
step <- 1

p_threshold <- 0.05

while(length(candidate) > 0){

  message("Step ", step)

  test_results <- list()

  for(a in candidate){

    full_vars <- c(
      selected,
      a,
      "sex",
      "province"
    )

    full_formula <- as.formula(
      paste(
        "Disease ~",
        paste(full_vars, collapse = "+")
      )
    )

    null_vars <- c(
      selected,
      "sex",
      "province"
    )

    null_formula <- as.formula(
      paste(
        "Disease ~",
        paste(null_vars, collapse = "+")
      )
    )

    vars_for_model <- unique(
      c(
        "Disease",
        selected,
        a,
        "sex",
        "province"
      )
    )

    temp_df <- hla_df[
      complete.cases(hla_df[, vars_for_model]),
    ]


    fit <- glm(
      full_formula,
      data = temp_df,
      family = binomial(link = "logit")
    )

    null_fit <- glm(
      null_formula,
      data = temp_df,
      family = binomial(link = "logit")
    )


    lrt <- anova(
      null_fit,
      fit,
      test = "Chisq"
    )

    p_value <- lrt$`Pr(>Chi)`[2]


    coef_table <- summary(fit)$coefficients

    res <- coef_table[a, ]


    test_results[[a]] <- data.frame(

      allele = a,

      beta = res["Estimate"],

      SE = res["Std. Error"],

      P = p_value,

      OR = exp(res["Estimate"]),

      lower95 = exp(
        res["Estimate"] -
          1.96 * res["Std. Error"]
      ),

      upper95 = exp(
        res["Estimate"] +
          1.96 * res["Std. Error"]
      )

    )

  }


  test_results <- do.call(
    rbind,
    test_results
  )


  best <- test_results[
    which.min(test_results$P),
  ]


  if(best$P > 0.05){
    break
  }


  step_results[[step]] <- data.frame(

    Step = step,

    Selected_allele = best$allele,

    Conditioned_on = paste(
      selected,
      collapse = ", "
    ),

    beta = best$beta,

    SE = best$SE,

    OR = best$OR,

    CI_lower = best$lower95,

    CI_upper = best$upper95,

    P = best$P

  )

  selected <- c(
    selected,
    best$allele
  )

  candidate <- setdiff(
    candidate,
    best$allele
  )

  step <- step + 1
}

HLA_conditional_res <- do.call(
  rbind,
  step_results
)

rownames(HLA_conditional_res) <- NULL

write.csv(HLA_conditional_res,'./res/HLA_condition_res.csv')

library(ggplot2)
library(dplyr)
library(forcats)

step_df = read.csv('./res/HLA_condition_res.csv')


step_df$neglog10P <- -log10(step_df$P)
step_df$Step2 = paste0('step',step_df$Step,'-',step_df$Selected_allele)

step_df$Step2 = factor(step_df$Step2,levels = step_df$Step2)

library(ggplot2)

p1 <- ggplot(
  step_df,
  aes(
    x = Step2,
    y = neglog10P
  )
) +
  geom_segment(
    aes(
      x = Step2,
      xend = Step2,
      y = 0,
      yend = neglog10P
    ),
    linewidth = 0.8
  ) +
  geom_point(
    size = 3
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "-log10(P)"
  ) +
  theme(
    aspect.ratio = 0.7,
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

p1