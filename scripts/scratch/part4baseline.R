# function to compute tests
test_fn_question3 <- function(.df){
  t_test(.df, 
         formula = level ~ group,
         order = c('ASD', 'TD'),
         alternative = 'two-sided',
         var.equal = F)
}

ttests_out_question3 <- part4_train %>%
  # drop ADOS score
  select(-ados) %>%
  # arrange in long format
  pivot_longer(-group, 
               names_to = 'protein', 
               values_to = 'level') %>%
  # nest by protein
  nest(data = c(level, group)) %>% 
  # compute t tests
  mutate(ttest = map(data, test_fn_question3)) %>%
  unnest(ttest) %>%
  # sort by p-value
  arrange(p_value) %>%
  # multiple testing correction
  mutate(m = n(),
         hm = log(m) + 1/(2*m) - digamma(1),
         rank = row_number(),
         p.adj = m*hm*p_value/rank)

# select significant proteins
proteins_s1_question3 <- ttests_out_question3 %>%
  slice_min(p.adj, n = 10) %>%
  pull(protein)

predictors_question3 <- part4_train %>%
  select(-c(group, ados))

response_question3 <- part4_train %>% pull(group) %>% factor()

# fit RF
rf_out_question3 <- randomForest(x = predictors_question3, 
                                 y = response_question3, 
                                 ntree = 1000, 
                                 importance = T)

# check errors
rf_out_question3$confusion

# compute importance scores
proteins_s2_question3 <- rf_out_question3$importance %>% 
  as_tibble() %>%
  mutate(protein = rownames(rf_out_question3$importance)) %>%
  slice_max(MeanDecreaseGini, n = 10) %>%
  pull(protein)
proteins_sstar_question3 <- intersect(proteins_s1_question3, proteins_s2_question3)

part4_sstar_train <- part4_train %>%
  select(group, any_of(proteins_sstar_question3)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

# fit logistic regression model to training set
fit_question3 <- glm(class ~ ., 
                     data = part4_sstar_train, 
                     family = 'binomial')

# evaluate errors on test set
class_metrics <- metric_set(sensitivity, 
                            specificity, 
                            accuracy,
                            roc_auc)

# create dataset with same columns for test set
part4_sstar_test <- part4_test %>%
  select(group, any_of(proteins_sstar_question3)) %>%
  mutate(class = (group == 'ASD')) %>%
  select(-group)

baseline <- part4_sstar_test %>%
  add_predictions(fit_question3, type = "response") %>%
  mutate(
    actual_class = factor(class, levels = c(FALSE, TRUE)),
    predicted_class = factor(pred > 0.5, levels = c(FALSE, TRUE))
  ) %>%
  class_metrics(
    truth = actual_class,
    estimate = predicted_class, 
    pred,
    event_level = "second"
  )
print(baseline)