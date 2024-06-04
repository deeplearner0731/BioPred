#' XGBoostSub_sur: Function for Training XGBoost Model with Customized Loss Function for survival outcomes
#'
#' This function trains an XGBoost model using a customized loss function based on the A-learning and weight-learning.
#'
#' @title XGBoost Model with Modified Loss Function for Subgroup Identification with Survival Outcomes
#' @description Function for training XGBoost model with customized loss function for survival outcomes
#' @param X_data The input features matrix.
#' @param y_data The input y matrix.
#' @param trt The treatment indicator vector. Should take values of 1 or -1, where 1 represents the treatment group and -1 represents the control group.
#' @param pi The propensity scores vector, which should range from 0 to 1, representing the probability of assignment to treatment.
#' @param censor The censor status vector. Should take values of 1 or 0, where 1 represents censoring and 0 represents an observed event.
#' @param Loss_type Type of loss function to use: "A_learning" or "Weight_learning".
#' @param params A list of additional parameters for the xgb.train function.
#' @param nrounds Number of boosting rounds. Default is 50.
#' @param disable_default_eval_metric If 1, default evaluation metric will be disabled.
#' @param verbose Logical. If TRUE, training progress will be printed; if FALSE, no progress will be printed.
#' @return Trained XGBoostSub_sur model.
#' @details
#' This function requires the 'xgboost' library. Make sure to install and load the 'xgboost' library before using this function.
#' @import xgboost
#' @export
XGBoostSub_sur <- function(X_data, y_data, trt, pi,censor, Loss_type = "Weight_learning", params = list(), nrounds = 50, disable_default_eval_metric = 1, verbose = TRUE) {

  risk_set_matrix <- function(time_to_event) {
    dim <- length(time_to_event)
    risk_set <- matrix(0, nrow = dim, ncol = dim)
    for (k in 1:dim) {
      arr <- numeric(dim)
      dummy <- which(time_to_event >= time_to_event[k])
      arr[dummy] <- 1
      risk_set[, k] <- arr
    }
    return(risk_set)
  }
  d1 <- risk_set_matrix(y_data)
  if (Loss_type == "A_learning") {
    gradient_sur_A<- function(preds, dmatrix, X_trt, pi_trt, delta) {
      c <- (X_trt + 1.0) / 2.0 - pi_trt
      dummy1 <- exp(c * preds)
      dummy2 <- dummy1 * c
      dummy4 <- as.vector(dummy1 %*% d1)
      value <- c - dummy2 / dummy4
      return(delta * value)
    }
    hessian_sur_A <- function(preds, dmatrix, X_trt, pi_trt, delta) {
      c<- X_trt
      dummy1 <- c * delta
      dummy2 <- exp(preds * c)
      dummy3 <- as.vector(dummy2 %*% d1)
      dummy4 <- dummy3^2
      dummy5 <- c * dummy2 * dummy3
      dummy6 <- (dummy2^2) * c
      value <- dummy1 * (dummy5 - dummy6) / dummy4
      return(value)
    }
  partial_log_sur <- function(X_trt, pi_trt, delta) {
      function(preds, dmatrix) {
        grad <- gradient_sur_A(preds, dmatrix, X_trt, pi_trt, delta)
        hess <- hessian_sur_A(preds, dmatrix, X_trt, pi_trt, delta)
        return(list(grad = grad, hess = hess))
      }
  }

  partial_sur_loss <- function(X_trt, pi_trt, delta) {
    function(preds, dmatrix) {
      c <- (X_trt + 1.0) / 2.0 - pi_trt
      dummy1 <- exp(c * preds)
      dummy2 <- as.vector(dummy1 %*% d1)
      dummy3 <- log(dummy2)
      dummy4 <- c * preds - dummy3
      value <-  delta * dummy4
      return(list(metric = "A_loss", value = sum(value) / length(value)))
    }
  }
  }

  if (Loss_type == "Weight_learning") {
    partial_log_sur <- function(X_trt, pi_trt, delta) {
      function(preds, dmatrix) {
        c <- X_trt
        const <- 1.0 / (X_trt * pi_trt + (1 - X_trt) / 2.0)
        dummy1 <- exp(c * preds)
        dummy2 <- dummy1 * c
        dummy4 <- as.vector(dummy1 %*% d1)
        value <- (c - dummy2 / dummy4)
        return(list(grad =  delta * value * const, hess = hessian_sur_weight(preds, dmatrix, X_trt, pi_trt, delta)))
      }
    }
  hessian_sur_weight <- function(preds, dmatrix, X_trt, pi_trt, delta) {
      c<- X_trt
      const <- 1.0 / (X_trt * pi_trt + (1 - X_trt) / 2.0)
      dummy1 <- c * delta * const
      dummy2 <- exp(preds * c)
      dummy3 <- as.vector(dummy2 %*% d1)
      dummy4 <- dummy3^2
      dummy5 <- c * dummy2 * dummy3
      dummy6 <- (dummy2^2) * c
      value <-  dummy1 * (dummy5 - dummy6) / dummy4
      return(value)
    }
  partial_sur_loss<- function(X_trt, pi_trt, delta) {
    function(preds, dmatrix) {
      c <- X_trt
      c1=(X_trt + 1.0) / 2.0 - pi_trt
      const <- 1.0 / (X_trt * pi_trt + (1 - X_trt) / 2.0)
      dummy1 <- exp(c * preds)
      dummy2 <- as.vector(dummy1 %*% d1)
      dummy3 <- log(dummy2)
      dummy4 <- c * preds - dummy3
      loss <-  delta * dummy4 * const
      loss= sum(loss) /length(preds)
      return(list(metric = "Weight_loss", value = loss))
    }
  }
  }

  # Create training matrix
  dtrain <- xgb.DMatrix(data = as.matrix(X_data), label = y_data)
  # Set additional parameters for training
  X_train_trt <- trt
  pi_train <- pi
  delta<-censor
  # Define objective and evaluation metric
  objective <- partial_log_sur(X_train_trt, pi_train,delta)
  eval_metric <- partial_sur_loss(X_train_trt, pi_train,delta)
  # Merge parameters
  all_params <- c(list(objective = objective), params)

  # Train the model
  model <- xgb.train(data = dtrain,
                     params = all_params,
                     watchlist = list(train = dtrain),
                     nrounds = nrounds, verbose = verbose,
                     disable_default_eval_metric = disable_default_eval_metric,
                     eval_metric = eval_metric)
  if (verbose) {
    cat("XGBoost model training finished.\n")
  }
  return(model)
}





#' eval_metric: Function for Evaluating XGBoostSub_con Model Performance
#'
#' This function evaluates the performance of an XGBoostSub_con model using a A-learning or weight-learning function.
#'
#' @title Evaluation Metrics for XGBoostSub_sur Model
#' @description Function for evaluating XGBoostSub_sur model performance.
#' @param model The trained XGBoostSub_sur model object.
#' @param X_feature The input features matrix.
#' @param y_label The input y matrix.
#' @param trt The treatment indicator vector. Should take values of 1 or -1, where 1 represents the treatment group and -1 represents the control group.
#' @param pi The propensity scores vector, which should range from 0 to 1, representing the probability of assignment to treatment.
#' @param censor The censor status vector. Should take values of 1 or 0, where 1 represents censoring and 0 represents an observed event.
#' @param Loss_type Type of loss function to use: "A_learning" or "Weight_learning".
#' @return Evaluation result of the XGBoostSub_sur model.
#' @import xgboost
#' @export
eval_metric_sur <- function(model, X_feature, y_label, pi, trt, censor, Loss_type = "A_learning") {
    risk_set_matrix <- function(time_to_event) {
    dim <- length(time_to_event)
    risk_set <- matrix(0, nrow = dim, ncol = dim)
    for (k in 1:dim) {
      arr <- numeric(dim)
      dummy <- which(time_to_event >= time_to_event[k])
      arr[dummy] <- 1
      risk_set[, k] <- arr
    }
    return(risk_set)
  }


if (Loss_type == "A_learning") {
  partial_loss <- function(X_trt, pi_trt, delta) {
      function(preds, dmatrix) {
        c <- (X_trt + 1.0) / 2.0 - pi_trt
        dummy1 <- exp(c * preds)
        dummy2 <- as.vector(dummy1 %*% d1)
        dummy3 <- log(dummy2)
        dummy4 <- c * preds - dummy3
        value <-  delta * dummy4
        return(list(metric = "A_loss", value = sum(value) / length(value)))
      }
  }
   }

  if (Loss_type == "Weight_learning") {
    partial_loss <- function(X_trt, pi_trt, delta) {
      function(preds, dmatrix) {
        c <- X_trt
        c1=(X_trt + 1.0) / 2.0 - pi_trt
        const <- 1.0 / (X_trt * pi_trt + (1 - X_trt) / 2.0)
        dummy1 <- exp(c * preds)
        dummy2 <- as.vector(dummy1 %*% d1)
        dummy3 <- log(dummy2)
        dummy4 <- c * preds - dummy3
        loss <-  delta * dummy4 * const
        loss= sum(loss) /length(preds)
        return(list(metric = "lossE", value = loss))
      }
    }
  }
    d1 <- risk_set_matrix(y_label)
    delta<- censor
    dtest <- xgb.DMatrix(data = as.matrix(X_feature), label = y_label)
    eval_metric_test <- partial_loss(trt, pi,delta)
    eval_result_test <- eval_metric_test(stats::predict(model, dtest), dtest)
  return(eval_result_test)
}









