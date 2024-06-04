#' XGBoostSub_bin: Function for Training XGBoost Model with Customized Loss Function for binary outcomes
#'
#' This function trains an XGBoost model using a customized loss function based on the A-learning and weight-learning.
#'
#' @title XGBoost Model with Modified Loss Function for Subgroup Identification with Binary Outcomes
#' @description Function for training XGBoost model with customized loss function for binary outcomes
#' @param X_data The input features matrix.
#' @param y_data The input y matrix.
#' @param trt The treatment indicator vector. Should take values of 1 or -1, where 1 represents the treatment group and -1 represents the control group.
#' @param pi The propensity scores vector, which should range from 0 to 1, representing the probability of assignment to treatment.
#' @param Loss_type Type of loss function to use: "A_learning" or "Weight_learning".
#' @param params A list of additional parameters for the xgb.train function.
#' @param nrounds Number of boosting rounds. Default is 50.
#' @param disable_default_eval_metric If 1, default evaluation metric will be disabled.
#' @param verbose Logical. If TRUE, training progress will be printed; if FALSE, no progress will be printed.
#' @return Trained XGBoostSub_bin model.
#' @details
#' This function requires the 'xgboost' library. Make sure to install and load the 'xgboost' library before using this function.
#'
#' After running this function, the returned model can be used like a regular xgboost model.
#' @import xgboost
#' @export
XGBoostSub_bin <- function(X_data, y_data, trt, pi, Loss_type = "A_learning", params = list(), nrounds = 50, disable_default_eval_metric = 1, verbose = TRUE) {

  if (Loss_type == "A_learning") {
    squared_log_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        preds <- 1.0 / (1.0 + exp(-preds))

        # Define gradient function
        gradient <- function(preds, dmatrix, X_trt, pi_trt) {
          y <- getinfo(dmatrix, "label")
          c <- (X_trt + 1.0) / 2.0 - pi_trt
          return(-y * c - c * exp(-c * preds) / (exp(-c * preds) + 1))
        }

        # Define hessian function
        hessian <- function(preds, dmatrix, X_trt, pi_trt) {
          y <- getinfo(dmatrix, "label")
          c <- (X_trt + 1.0) / 2.0 - pi_trt
          return(c^2 * exp(-c * preds) / (exp(-c * preds) + 1)^2)
        }

        grad <- gradient(preds, dmatrix, X_trt, pi_trt)
        hess <- hessian(preds, dmatrix, X_trt, pi_trt)
        return(list(grad = grad, hess = hess))
      }
    }
    rmsle_version_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        y <- getinfo(dmatrix, "label")
        preds <- 1.0 / (1 + exp(-preds))
        c <- (X_trt + 1.0) / 2.0 - pi_trt
        elements <- -y * c * preds + log(1 + exp(-c * preds))
        loss <- sqrt(sum(elements) / length(y))
        return(list(metric = "A_loss", value = loss))
      }
    }
  }

  if (Loss_type == "Weight_learning") {
    squared_log_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        preds <- 1.0 / (1.0 + exp(-preds))

        # Define gradient function
        gradient <- function(preds, dmatrix, X_trt, pi_trt) {
          y <- getinfo(dmatrix, "label")
          c <- (1.0- X_trt) / 2.0 + pi_trt*X_trt

          return (1 / c * (-y * X_trt - X_trt * exp(-X_trt * preds) / (exp(-X_trt  * preds) + 1)))
        }

        # Define hessian function
        hessian <- function(preds, dmatrix, X_trt, pi_trt) {
          y <- getinfo(dmatrix, "label")
          c <- (1.0- X_trt) / 2.0 + pi_trt*X_trt
          return (1 / c * (X_trt^2 * exp(-X_trt * preds) / (exp(-X_trt * preds) + 1)^2))
        }

        grad <- gradient(preds, dmatrix, X_trt, pi_trt)
        hess <- hessian(preds, dmatrix, X_trt, pi_trt)
        return(list(grad = grad, hess = hess))
      }
    }
    rmsle_version_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        y <- getinfo(dmatrix, "label")
        preds <- 1.0 / (1 + exp(-preds))

        c <- (1.0- X_trt) / 2.0 + pi_trt*X_trt

        elements <- 1 / c * (-y * X_trt * preds + log(1 + exp(-X_trt * preds)))
        loss <- sqrt(sum(elements) / length(y))
        return(list(metric = "Weight_loss", value = loss))
      }
    }
  }

  # Create training matrix
  dtrain <- xgb.DMatrix(data = as.matrix(X_data), label = y_data)

  # Set additional parameters for training
  X_train_trt <- trt
  pi_train <- pi

  # Define objective and evaluation metric
  objective <- squared_log_binary(X_train_trt, pi_train)
  eval_metric <- rmsle_version_binary(X_train_trt, pi_train)

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











#' eval_metric: Function for Evaluating XGBoostSub_bin Model Performance
#'
#' This function evaluates the performance of an XGBoostSub_bin model using a A-learning or weight-learning function.
#'
#' @title Evaluation Metrics for XGBoostSub_bin Model
#' @description Function for evaluating XGBoostSub_bin model performance.
#' @param model The trained XGBoostSub_bin model object.
#' @param X_feature The input features matrix.
#' @param y_label The input y matrix.
#' @param trt The treatment indicator vector. Should take values of 1 or -1, where 1 represents the treatment group and -1 represents the control group.
#' @param pi The propensity scores vector, which should range from 0 to 1, representing the probability of assignment to treatment.
#' @param Loss_type Type of loss function to use: "A_learning" or "Weight_learning".
#' @return Evaluation result of the XGBoostSub_bin model.
#' @import xgboost
#' @export
eval_metric_bin <- function(model, X_feature, y_label, pi,trt, Loss_type = "A_learning") {
  if (Loss_type == "A_learning") {
    rmsle_version_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        y <- getinfo(dmatrix, "label")
        preds <- 1.0 / (1 + exp(-preds))
        c <- (X_trt + 1.0) / 2.0 - pi_trt
        elements <- -y * c * preds + log(1 + exp(-c * preds))
        loss <- sqrt(sum(elements) / length(y))
        return(list(metric = "A_loss", value = loss))
      }
    }
    X_train_trt<-trt
    dtest <- xgb.DMatrix(data = as.matrix(X_feature), label = y_label)
    eval_metric_test <- rmsle_version_binary(trt, pi)
    eval_result_test <- eval_metric_test(stats::predict(model, dtest), dtest)
  } else if (Loss_type == "Weight_learning") {
    rmsle_version_binary <- function(X_trt, pi_trt) {
      function(preds, dmatrix) {
        y <- getinfo(dmatrix, "label")
        preds <- 1.0 / (1 + exp(-preds))

        c <- (1.0- X_trt) / 2.0 + pi_trt*X_trt

        elements <- 1 / c * (-y * X_trt * preds + log(1 + exp(-X_trt * preds)))
        loss <- sqrt(sum(elements) / length(y))
        return(list(metric = "Weight_loss", value = loss))
      }
    }
    X_train_trt<-trt
    dtest <- xgb.DMatrix(data = as.matrix(X_feature), label = y_label)
    eval_metric_test <- rmsle_version_binary(X_train_trt, pi)
    eval_result_test <- eval_metric_test(stats::predict(model, dtest), dtest)
  }
  return(eval_result_test)
}


