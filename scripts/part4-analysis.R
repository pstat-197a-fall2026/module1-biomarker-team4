
# Tools and data
library(tidyverse)
library(tidymodels)
library(randomForest)
library(modelr)

load("data/biomarker-clean.RData")

set.seed(101422)

part4_split <- initial_split(
  biomarker_clean,
  prop = 0.8,
  strata = group
)

part4_train <- training(part4_split)
part4_test <- testing(part4_split)

source("scripts/scratch/part4baseline.R", local = TRUE)

# Rank proteins on training data
protein_names <- setdiff(
  names(part4_train),
  c("group", "ados")
)

part4_ranking <- purrr::map_dfr(protein_names, function(protein) {

  result <- t.test(
    part4_train[[protein]][part4_train$group == "ASD"],
    part4_train[[protein]][part4_train$group == "TD"]
  )

  tibble(
    protein = protein,
    p_value = result$p.value
  )
}) %>%
  arrange(p_value)

set.seed(101422)

part4_folds <- vfold_cv(
  part4_train,
  v = 5,
  strata = group
)

evaluate_fold <- function(split) {

  fold_train <- analysis(split)
  fold_valid <- assessment(split)

  ranking <- purrr::map_dfr(protein_names, function(protein) {

    result <- t.test(
      fold_train[[protein]][fold_train$group == "ASD"],
      fold_train[[protein]][fold_train$group == "TD"]
    )

    tibble(
      protein = protein,
      p_value = result$p.value
    )
  }) %>%
    arrange(p_value)

  purrr::map_dfr(1:3, function(k) {

    selected <- head(ranking$protein, k)

    model_data <- fold_train %>%
      select(all_of(selected)) %>%
      mutate(class = fold_train$group == "ASD")

    model <- glm(
      class ~ .,
      data = model_data,
      family = binomial()
    )

    probabilities <- predict(
      model,
      newdata = fold_valid,
      type = "response"
    )

    tibble(
      panel_size = k,
      accuracy = mean(
        (probabilities > 0.5) == (fold_valid$group == "ASD")
      )
    )
  })
}

part4_cv_results <- part4_folds %>%
  mutate(scores = purrr::map(splits, evaluate_fold)) %>%
  select(id, scores) %>%
  unnest(scores)

part4_cv_summary <- part4_cv_results %>%
  group_by(panel_size) %>%
  summarise(
    mean_accuracy = mean(accuracy),
    sd_accuracy = sd(accuracy),
    .groups = "drop"
  )

# Choose size and proteins
part4_best <- part4_cv_summary %>%
  arrange(desc(mean_accuracy), panel_size) %>%
  slice(1)


best_k <- part4_best$panel_size

best_proteins <- part4_ranking %>%
  slice_head(n = best_k) %>%
  pull(protein)

# Fit final model
final_train <- part4_train %>%
  select(all_of(best_proteins)) %>%
  mutate(
    class = factor(part4_train$group, levels = c("TD", "ASD"))
  )

final_model <- glm(
  class ~ .,
  data = final_train,
  family = binomial()
)

test_probabilities <- predict(
  final_model,
  newdata = part4_test %>% select(all_of(best_proteins)),
  type = "response"
)

test_predictions <- factor(
  ifelse(test_probabilities > 0.5, "ASD", "TD"),
  levels = c("TD", "ASD")
)

part4_test_results <- tibble(
  truth = factor(part4_test$group, levels = c("TD", "ASD")),
  estimate = test_predictions
)

alternative_metrics <- part4_test_results %>%
  mutate(prob_asd = test_probabilities) %>%
  class_metrics(
    truth = truth,
    estimate = estimate,
    prob_asd,
    event_level = "second"
  )

part4_comparison <- bind_rows(
  baseline %>%
    mutate(
      model = "Baseline: t-test + random forest intersection",
      panel_size = length(proteins_sstar_question3)
    ),
  alternative_metrics %>%
    mutate(
      model = "Alternative: t-test only",
      panel_size = length(best_proteins)
    )
) %>%
  select(model, panel_size, .metric, .estimate) %>%
  pivot_wider(names_from = .metric, values_from = .estimate)


baseline_accuracy <- baseline %>%
  filter(.metric == "accuracy") %>%
  pull(.estimate)

alternative_accuracy <- alternative_metrics %>%
  filter(.metric == "accuracy") %>%
  pull(.estimate)

accuracy_difference_pp <- 100 * (alternative_accuracy - baseline_accuracy)
protein_reduction <- length(proteins_sstar_question3) - length(best_proteins)
