library(dplyr)
library(tidyr)
library(readxl)
library(forecast)
library(zoo)
library(TTR)
library(quadprog)
library(ROI)
library(ROI.plugin.quadprog)
library(ggplot2)
library(dplyr)
library(lpSolve)
library(nloptr)
library(gt)
library(webshot2)
library(CVXR)
library(writexl)
library(patchwork)
library(reshape2)


set.seed(123)

# Read the dataset
df <- read_excel('Final Dataset.xlsx')


# Main function for solving optimization problem and resampling
optimize_funding <- function(df,
                             perf_cols,          
                             a0 = 0.5,
                             a1 = 0.5,
                             lambda = 0.1,
                             bootstrap = FALSE,
                             n_boot = 100,
                             bootstrap_type = c("block", "jackknife_uni")
) {
  F <- as.matrix(df %>% select(all_of(perf_cols)))
  n_metrics <- ncol(F)
  
  # Single experiment
  if(!bootstrap) {
    
    # Mean enrollment for modified model 
    mean_E_by_year <- df %>%
      group_by(Year) %>%
      summarise(mean_E = mean(EnrollmentPart), .groups = "drop")
    
    df <- df %>% left_join(mean_E_by_year, by = "Year")
    
    # Weights
    w <- Variable(n_metrics, nonneg = TRUE)
    
    y_pred <- a0 * df$EnrollmentPart + a1 * (F %*% w) * (df$EnrollmentPart / df$mean_E)
    years <- unique(df$Year)
    
    # Constraints
    year_constraints <- lapply(years, function(yr) {
      idx <- which(df$Year == yr)
      sum(y_pred[idx]) <= sum(df$Funding[idx])
    })
    constraints <- c(
      year_constraints,
      list(w >= 0)
    )
    
    # Objective function
    objective <- Minimize(sum_squares(y_pred - df$Funding) + lambda * sum_squares(w))
    
    problem <- Problem(objective, constraints)
    
    result <- solve(problem)
    
    # Results
    w_opt <- result$getValue(w)
    w_opt <- as.vector(result$getValue(w))  
    w_norm <- w_opt / sum(w_opt)
    
    # Predicted funding
    df$Predicted_Funding <- a0 * df$EnrollmentPart + a1 * (F %*% w_opt) * (df$EnrollmentPart / df$mean_E)
    
    return(list(
      w_opt = w_opt,
      w_norm = w_norm,
      result_df = df
    ))
    
    # Resampling
  } else if(bootstrap) {
    
    metric_names <- perf_cols
    w_boot <- matrix(NA, nrow = n_boot, ncol = n_metrics)
    colnames(w_boot) <- metric_names
    fail_count <- 0
  
      if (bootstrap_type == "block") {
        for (i in 1:n_boot) {
        years <- sort(unique(df$Year))
        n_years <- length(years)
        # Preparing blocks
        b <- 2  # block length
        n_blocks <- ceiling(n_years / b)
        n_full_blocks <- floor(n_years / b)
        
        blocks <- split(years[1:(n_full_blocks * b)], rep(1:n_full_blocks, each = b))
        
        if (n_years %% b != 0) {
          blocks[[n_full_blocks]] <- c(blocks[[n_full_blocks]], years[(n_full_blocks * b + 1):n_years])
        }
        
        sampled_blocks <- sample(seq_along(blocks), length(blocks), replace = TRUE)
        sampled_years <- unlist(blocks[sampled_blocks])
        
        df_sample <- do.call(rbind, lapply(sampled_years, function(yr) df[df$Year == yr, ]))
      
        # Mean enrollment for modified model
        mean_E_by_year <- df_sample %>%
          group_by(Year) %>%
          summarise(mean_E = mean(EnrollmentPart), .groups = "drop")
        
        df_sample <- df_sample %>% left_join(mean_E_by_year, by = "Year")
        
        w_sample <- Variable(n_metrics, nonneg = TRUE)
      
        F_sample <- as.matrix(df_sample %>% select(all_of(metric_names)))
        y_pred_sample <- a0 * df_sample$EnrollmentPart + a1 * (F_sample %*% w_sample) * (df_sample$EnrollmentPart / df_sample$mean_E)
        
        year_counts <- table(df_sample$Year)  
        # Constraints
        year_constraints_sample <- lapply(unique(df_sample$Year), function(yr) {
          idx <- which(df_sample$Year == yr)
          sum(y_pred_sample[idx]) <= sum(df_sample$Funding[df_sample$Year == yr])
        })
        constraints <- c(
          year_constraints_sample,
          list(w_sample >= 0)
        )
        
        objective_sample <- Minimize(
          sum_squares(y_pred_sample - df_sample$Funding) +
            lambda * sum_squares(w_sample)
        )
        problem_sample <- Problem(objective_sample, constraints = constraints)
        
        result_sample <- tryCatch(
          solve(problem_sample),
          error = function(e) { fail_count <<- fail_count + 1; NULL },
          warning = function(w) { fail_count <<- fail_count + 1; NULL }
        )
        
        # Results
        if (!is.null(result_sample)) {
          w_boot[i, ] <- as.vector(result_sample$getValue(w_sample))
          
        }
        } 
        w_mean <- colMeans(w_boot, na.rm = TRUE)
        w_mean_norm <- w_mean / sum(w_mean)
        
        w_boot_norm <- t(apply(w_boot, 1, function(x) x / sum(x, na.rm = TRUE)))
        
        w_sd   <- apply(w_boot, 2, sd, na.rm = TRUE)
        w_sd_norm <- apply(w_boot_norm, 2, sd, na.rm = TRUE)
        
        w_ci <- apply(w_boot_norm, 2, function(x) {
          quantile(x, probs = c(0.025, 0.975), na.rm = TRUE)
        })
      }
      else if (bootstrap_type == "jackknife_uni") {
        unis <- unique(df$University_Code)
        n_unis <- length(unis)
        w_jack <- matrix(NA, nrow = n_unis, ncol = n_metrics)
        colnames(w_jack) <- metric_names
        
        mean_E_by_year <- df %>% group_by(Year) %>% summarise(mean_E = mean(EnrollmentPart), .groups="drop")
        
        # Full dataset optimization for jacckinfe bias calculation later
        df_full <- df %>% left_join(mean_E_by_year, by = "Year")
        w_full_var <- Variable(n_metrics, nonneg = TRUE)
        F_full <- as.matrix(df_full %>% select(all_of(metric_names)))
        y_pred_full <- a0 * df_full$EnrollmentPart + a1 * (F_full %*% w_full_var) * (df_full$EnrollmentPart / df_full$mean_E)
        year_constraints_full <- lapply(unique(df_full$Year), function(yr) {
          idx <- which(df_full$Year == yr)
          sum(y_pred_full[idx]) <= sum(df_full$Funding[idx])
        })
        problem_full <- Problem(Minimize(sum_squares(y_pred_full - df_full$Funding) + lambda * sum_squares(w_full_var)),
                                constraints = c(year_constraints_full, list(w_full_var >= 0)))
        result_full <- solve(problem_full)
        w_full <- as.vector(result_full$getValue(w_full_var))
        w_full_norm <- w_full / sum(w_full)

        for (i in seq_along(unis)) {
          df_jack <- df[df$University_Code != unis[i], ]
          mean_E_by_year <- df_jack %>% group_by(Year) %>% summarise(mean_E = mean(EnrollmentPart), .groups="drop")
          df_jack <- df_jack %>% left_join(mean_E_by_year, by="Year")
          
          w_var <- Variable(n_metrics, nonneg = TRUE)
          F_jack <- as.matrix(df_jack %>% select(all_of(metric_names)))
          y_pred <- a0 * df_jack$EnrollmentPart + a1 * (F_jack %*% w_var) * (df_jack$EnrollmentPart / df_jack$mean_E)
          year_constraints <- lapply(unique(df_jack$Year), function(yr) {
            idx <- which(df_jack$Year == yr)
            sum(y_pred[idx]) <= sum(df_jack$Funding[idx])
          })
          constraints <- c(
            year_constraints,
            list(w_var >= 0)
          )
          problem <- Problem(Minimize(sum_squares(y_pred - df_jack$Funding) + lambda * sum_squares(w_var)),
                             constraints = constraints)
          result <- solve(problem)
          w_jack[i, ] <- as.vector(result$getValue(w_var))
        }
        # Results
        w_jack_norm <- t(apply(w_jack, 1, function(x) x / sum(x, na.rm = TRUE)))
        w_mean <- colMeans(w_jack)
        w_sd <- sqrt((n_unis - 1) / n_unis * colSums((w_jack - w_mean)^2))
        w_mean_norm <- w_mean / sum(w_mean)
        w_bias <- (n_unis - 1) * (w_mean_norm - w_full_norm)
        w_sd_norm <- apply(w_jack_norm, 2, sd, na.rm = TRUE)
        w_ci <- apply(w_jack_norm, 2, function(x) quantile(x, probs=c(0.025, 0.975), na.rm=TRUE))
      
    }
    cat("Bootstraping simulations failed:", fail_count, "\n")
    
    mean_E_by_year <- df %>%
      group_by(Year) %>%
      summarise(mean_E = mean(EnrollmentPart), .groups = "drop")
    df <- df %>% left_join(mean_E_by_year, by = "Year")
    
    df$Predicted_Funding <- a0 * df$EnrollmentPart + a1 * (F %*% w_mean) * (df$EnrollmentPart / df$mean_E)
    return(list(
      w_boot = w_boot,
      w_mean = w_mean,
      w_mean_norm = w_mean_norm,
      w_sd = w_sd,
      w_sd_norm = w_sd_norm,
      w_bias = if(bootstrap_type=="jackknife_uni") w_bias else NULL,
      w_ci = t(w_ci),
      result_df = df
    ))
  
}
}

# Function for results analysis
analyze_funding_results <- function(df, plot, save_excel,
                                    excel_data_path, image_data_path) {
  
  year_summary <- df %>%
    group_by(Year) %>%
    summarise(
      Actual = sum(Funding),
      Predicted = sum(Predicted_Funding),
      Slack = Actual - Predicted
    )
  
  total_actual <- sum(df$Funding)
  total_predicted <- sum(df$Predicted_Funding)
  slack_total <- total_actual - total_predicted
  slack_percent <- 100 * slack_total / total_actual
  
  cat("\nTotal unallocated funds across all years:\n")
  cat("Unallocated:", slack_total, "\n")
  cat("Percent of total funding not allocated:", slack_percent, "%\n")
  
  
  print("Yearly budget check:")
  print(year_summary)
  
  if(any(year_summary$Slack < 0)) {
    warning("Predicted funding exceeds actual budget in some years!")
  }
  
  uni_year_summary <- df %>%
    mutate(
      Change = Predicted_Funding - Funding,
      Percent_Change = 100 * Change / Funding
    ) %>%
    arrange(University_Code, Year)
  
  print("University-level funding change per year:")
  print(uni_year_summary %>% select(University_Code, Year, Funding, Predicted_Funding, Change, Percent_Change))
  
  uni_total_summary <- uni_year_summary %>%
    group_by(University_Code) %>%
    summarise(
      Actual_Total = sum(Funding),
      Predicted_Total = sum(Predicted_Funding),
      Change_Total = sum(Change),
      Percent_Change_Total = 100 * Change_Total / Actual_Total
    ) %>%
    arrange(desc(Percent_Change_Total))
  
  print("University-level total funding change:")
  print(uni_total_summary)
  
  df <- df %>%
    mutate(
      Perf_Part_Amount = Predicted_Funding - EnrollmentPart  # remaining money from performance
    )
  
  
  perf_summary <- df %>%
    group_by(University) %>%
    summarise(
      Actual_Total = sum(Funding),
      Perf_Total = sum(Perf_Part_Amount),
      Enrollment_Model_Total = sum(Predicted_Funding - Perf_Part_Amount),
      Extra_From_Performance = Perf_Total - 0,          # baseline actual has no performance component
      Percent_Extra = 100 * Perf_Total / Actual_Total
    ) %>%
    arrange(desc(Percent_Extra))
  
  cat("\nPerformance-based allocation summary per university:\n")
  print(perf_summary)
  
  
  if(plot) {
    df_list <- df %>% split(.$University_Code)
    plot_uni <- function(data, uni_name) {
      df_long <- data %>%
        select(Year, Funding, Predicted_Funding) %>%
        tidyr::pivot_longer(
          cols = c(Funding, Predicted_Funding),
          names_to = "Type",
          values_to = "Amount"
        )
      
      ggplot(df_long, aes(x = factor(Year), y = Amount, fill = Type)) +
        geom_col(position = "dodge") +
        labs(title = uni_name, x = "Year", y = "Funding") +
        theme_minimal() +
        theme(
          legend.position = "bottom",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", hjust = 0.5)
        )
    }
    plots <- lapply(names(df_list), function(u) plot_uni(df_list[[u]], u))
    final_plot <- wrap_plots(plots, ncol = 4)
    
    print(final_plot)
    
    ggsave(image_data_path,
           final_plot,
           width = 14, height = 10, dpi = 300)
    
  }
  if(save_excel) {
    write_xlsx(
      list(
        "Yearly_Budget" = year_summary,
        "University_Yearly" = uni_year_summary,
        "University_Totals" = uni_total_summary,
        "Performance_Summary" = perf_summary,
        "Full_Data" = df
      ),
      path = excel_data_path
    )
    
    message("\nSaved Excel file to: ", excel_data_path)
  }
  
  return(list(
    year_summary = year_summary,
    uni_year_summary = uni_year_summary,
    uni_total_summary = uni_total_summary,
    plot = if(plot) final_plot else NULL
  ))
}

# Function for blocked cross validation to get optimal lambda
cv_lambda <- function(df, metric_names, lambda_values, a0, a1) {
  years <- sort(unique(df$Year))
  cv_results <- data.frame(lambda = lambda_values, MAE = NA)
  all_weights <- list()
  for (i in seq_along(lambda_values)) {
    lam <- lambda_values[i]
    fold_errors <- c()
    fold_weights <- matrix(NA, nrow = length(metric_names), ncol = 0)
    
    for (yr in years) {
      test_years <- c(yr, yr + 1)
      if (!all(test_years %in% years)) next
      
      # Train on all years except the test block
      train_df <- df %>% filter(!Year %in% test_years)
      
      # Test on the two-year block
      test_df  <- df %>% filter(Year %in% test_years)
      
      fit <- optimize_funding(
        train_df,
        perf_cols = metric_names,
        a0 = a0,
        a1 = a1,
        lambda = lam,
        bootstrap = FALSE
      )
      
      w <- fit$w_opt
      w_norm <- fit$w_norm
      fold_weights <- cbind(fold_weights, w_norm)

      mean_E_by_year <- test_df %>%
        group_by(Year) %>%
        summarise(mean_E = mean(EnrollmentPart), .groups = "drop")
      
      test_df <- test_df %>% left_join(mean_E_by_year, by = "Year")
      
      test_pred <- test_df %>%
        mutate(Predicted = a0 * EnrollmentPart + a1 * rowSums(as.matrix(select(., all_of(metric_names))) * w) * (EnrollmentPart / mean_E))

      fold_mae <- mean(abs(test_pred$Funding - test_pred$Predicted))
      fold_errors <- c(fold_errors, fold_mae)
    }
    
    cv_results$MAE[i] <- mean(fold_errors)
    all_weights[[i]] <- rowMeans(fold_weights)
  }
  best_lambda <- cv_results$lambda[which.min(cv_results$MAE)]
  cv_plot <- ggplot(cv_results, aes(x = lambda, y = MAE)) +
    geom_line(size = 1, color = "steelblue") +
    geom_point(size = 2, color = "darkred") +
    geom_vline(xintercept = best_lambda, linetype = "dashed", color = "darkgreen", size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste0("Cross-Validation MAE vs Lambda (a0=", a0, ", a1=", a1, ")"),
      subtitle = paste0("Best lambda = ", best_lambda),
      x = expression(lambda),
      y = "Mean Absolute Error (MAE)"
    )
  
  print(cv_plot)
  
  ggsave(
    filename = paste0("lambda_CV_a0_", a0*100, "_a1_", a1*100, ".png"),
    plot = cv_plot,
    width = 8,
    height = 4,
    dpi = 300
  )
  
  weight_df <- data.frame(
    metric = rep(metric_names, times = length(lambda_values)),
    lambda = rep(lambda_values, each = length(metric_names)),
    weight = unlist(all_weights)
  )
  
  weight_plot <- ggplot(weight_df, aes(x = lambda, y = weight, color = metric)) +
    geom_line(size = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = "Metric Weights vs Lambda",
      x = expression(lambda),
      y = "Weight"
    )
  
  print(weight_plot)
  ggsave(
    filename = paste0("weights_CV_a0_", a0*100, "_a1_", a1*100, ".png"),
    plot = weight_plot,
    width = 8,
    height = 4,
    dpi = 300
  )
  return(list(
    results = cv_results,
    best_lambda = best_lambda)
  )
}


##########################
#   Perform Simulations  #
##########################

# Calculate performance metrics
df <- df %>%
  mutate(
    StudentEmployeeRate   = `Number of Employees` / `Number current students`,
    GradRate = (`Bachelor graduates`+`Master graduates`+ `Continious graduates`) / `Number current students`,
    ExchangeRate = (`Bachelors exchange students`+`Master exchange students`+`Continuous exchange students`)/ `Number current students`
  )


metric_names <- c(
  "GradRate",
  "ExchangeRate",
  "StudentEmployeeRate"
)

######################
##### Experiments ####
######################

lambda_values <- seq(0, 1, by = 0.1)


# a0=0.85, a1=0.15

cv_result_85_15 <- cv_lambda(
  df = df,
  metric_names = metric_names,
  lambda_values = lambda_values,
  a0 = 0.85,
  a1 = 0.15
)

lambda_85_15 <- cv_result_85_15$best_lambda
lambda_85_15
lambda_85_15 = 0.1

a0_85_a1_15_result <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.85,
  a1 = 0.15,
  lambda = lambda_85_15,
)

a0_85_a1_15_result$w_norm
a0_85_a1_15_result$w_opt


a0_85_a1_15_analysis <- analyze_funding_results(a0_85_a1_15_result$result_df,
                                                save_excel=TRUE,
                                                excel_data_path="../rezultatai/a0_85_a1_a15_modified_results.xlsx",
                                                plot=TRUE,
                                                image_data_path="../images/results/a0_85_a1_a15_modified_results.png")




# a0=0.75, a1=0.25

optimize_lambda(df, metric_names, a0=0.75, a1=0.25)

lambda_75_25=0.1



cv_result_75_25 <- cv_lambda(
  df = df,
  metric_names = metric_names,
  lambda_values = lambda_values,
  a0 = 0.75,
  a1 = 0.25
)
lambda_75_25 <- cv_result_75_25$best_lambda

lambda_75_25
lambda_75_25 = 0.1
a0_75_a1_25_result <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.75,
  a1 = 0.25,
  lambda = lambda_75_25
)

a0_75_a1_25_result$w_norm

a0_75_a1_25_analysis <- analyze_funding_results(a0_75_a1_25_result$result_df,
                                                save_excel=TRUE,
                                                excel_data_path="../rezultatai/a0_75_a1_a25_modified_results.xlsx",
                                                plot=TRUE,
                                                image_data_path="../images/results/a0_75_a1_a25_modified_results.png")




# a0=0.5, a1=0.5
lambda_values <- seq(0, 1, by = 0.05)


cv_result_5_5 <- cv_lambda(
  df = df,
  metric_names = metric_names,
  lambda_values = lambda_values,
  a0 = 0.5,
  a1 = 0.5
)
lambda_5_5 <- cv_result_5_5$best_lambda
lambda_5_5
lambda_5_5 = 0.1

a0_5_a1_5_result <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.5,
  a1 = 0.5,
  lambda = lambda_5_5
)

a0_5_a1_5_result$w_norm

a0_5_a1_5_analysis <- analyze_funding_results(a0_5_a1_5_result$result_df,
                                              save_excel=TRUE,
                                              excel_data_path="../rezultatai/a0_5_a1_a5_modified_results.xlsx",
                                              plot=TRUE,
                                              image_data_path="../images/results/a0_5_a1_a5_modified_results.png")




# BOOTSTRAPING


# BLOCK

# a0=0.85, a1=0.15

a0_85_a1_15_result_block <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.85,
  a1 = 0.15,
  lambda = lambda_85_15,
  bootstrap=TRUE,
  n_boot=100,
  bootstrap_type="block"
)

a0_85_a1_15_result_block$w_mean_norm

a0_85_a1_15_result_block$w_ci

a0_85_a1_15_result_block$w_sd_norm
sum(a0_85_a1_15_result_boot$w_mean_norm)

a0_85_a1_15_analysis_block <- analyze_funding_results(a0_85_a1_15_result_block$result_df,
                                                      plot=TRUE)

# JACK UNI

a0_85_a1_15_result_jack_uni <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.85,
  a1 = 0.15,
  lambda = lambda_85_15,
  bootstrap=TRUE,
  bootstrap_type="jackknife_uni"
)

a0_85_a1_15_result_jack_uni$w_mean_norm

a0_85_a1_15_result_jack_uni$w_ci

a0_85_a1_15_result_jack_uni$w_sd_norm
a0_85_a1_15_result_jack_uni$w_bias

# a0=0.75, a1=0.25

# BLOCK

a0_75_a1_25_result_block <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.75,
  a1 = 0.25,
  lambda = lambda_75_25,
  bootstrap=TRUE,
  bootstrap_type="block"
)

a0_75_a1_25_result_block$w_mean_norm
a0_75_a1_25_result_block$w_ci

a0_75_a1_25_result_block$w_sd_norm

a0_75_a1_25_analysis_boot <- analyze_funding_results(a0_75_a1_25_result_boot$result_df)

# JACK UNI


a0_75_a1_25_result_jack_uni <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.75,
  a1 = 0.25,
  lambda = lambda_75_25,
  bootstrap=TRUE,
  bootstrap_type="jackknife_uni"
)

a0_75_a1_25_result_jack_uni$w_mean_norm
a0_75_a1_25_result_jack_uni$w_ci
a0_75_a1_25_result_jack_uni$w_sd_norm
a0_75_a1_25_result_jack_uni$w_bias

# a0=0.5, a1=0.5

# BLOCK

a0_5_a1_5_result_block <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.5,
  a1 = 0.5,
  lambda = lambda_5_5,
  bootstrap=TRUE,
  n_boot=100,
  bootstrap_type="block"
)

a0_5_a1_5_result_block$w_mean_norm

a0_5_a1_5_result_block$w_ci
a0_5_a1_5_result_block$w_sd_norm
a0_5_a1_5_analysis <- analyze_funding_results(a0_5_a1_5_result$result_df)

# JACK UNI


a0_5_a1_5_result_jack_uni <- optimize_funding(
  df,
  perf_cols = metric_names,
  a0 = 0.5,
  a1 = 0.5,
  lambda = lambda_5_5,
  bootstrap=TRUE,
  bootstrap_type="jackknife_uni"
)

a0_5_a1_5_result_jack_uni$w_mean_norm

a0_5_a1_5_result_jack_uni$w_ci
a0_5_a1_5_result_jack_uni$w_sd_norm
a0_5_a1_5_result_jack_uni$w_bias
