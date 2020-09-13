#' Weighted log-rank test
#'
#' @param df A data frame. Assume standard structure for time-to-event data.
#' @param trt_colname A character string. The name of the treatment column in \code{df}.
#' @param time_colname A character string. The name of the column in \code{df} with the survival times.
#' @param event_colname A character string. The name of the column in \code{df} with the event status. Assumes 1 means event; 0 means censored.
#' @param wlr The type of weighted log-rank test. Either the default "lr" for a standard log-rank test, "mw" for a modestly-weighted log-rank test, or "fh" for the Fleming-Harrington rho-gamma family.
#' @param s_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), \code{s_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param t_star This is a parameter for the "mw"-test. Either \code{s_star} or \code{t_star} must be specified. The weights are defined as w(t) = 1 / max(S(t-), S(\code{t_star}), where S is the Kaplan-Meier estimate in the pooled data.
#' @param rho Rho parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @param gamma Gamma parameter in the "fh"-test. The weights are defined as w(t) = S(t-)^\code{rho} (1 - S(t-))^\code{gamma}, where S is the Kaplan-Meier estimate in the pooled data.
#' @return A data frame. The outcome of the weighted log-rank test. There is a column \code{o_minus_e_trt} to indicate which treatment the "obs - exp" refers to.
#' @export

wlrt <- function(df,
                 trt_colname,
                 time_colname,
                 event_colname,
                 wlr = "lr",
                 s_star = NULL,
                 t_star = NULL,
                 rho = NULL,
                 gamma = NULL){

  if(!all(is.character(c(trt_colname,
                    time_colname,
                    event_colname)))) stop("trt_colname, time_colname and event_colname must all be character strings")

  if (any(is.na(df[,c(trt_colname,
                      time_colname,
                      event_colname)]))) stop("NA's in data set. wlrt doesn't have a default for missing data.")

  if (!(wlr %in% c("lr", "fh", "mw"))) stop("wlr must be one of: 'lr', 'fh', 'mw'")

  #### get summary
  s_sum <- summary(survival::survfit(survival::Surv(eval(as.name(time_colname)),
                                                    eval(as.name(event_colname))) ~ 1,
                                    data= df))

  ### get survival probabilities for the pooled data
  s_pool <- s_sum$surv

  if (wlr == "lr"){
    ### standard log-rank test weights
    w <- rep(1, length(s_pool))
  }
  else if (wlr == "fh"){
    if (is.null(rho) || is.null(gamma)) stop("must specify rho and gamma")
    if (rho < 0 || gamma < 0) stop("rho and gamma must be non-negative")

    w <- c(1, s_pool[-length(s_pool)]) ^ rho * (1 - c(1, s_pool[-length(s_pool)])) ^ gamma

  }
  else{
    if (is.null(t_star) && is.null(s_star)) stop("must specify either t_star or s_star")
    if (!is.null(t_star) && !is.null(s_star)) stop("must specify either t_star or s_star (not both)")
    ### modest weights
    if (!is.null(t_star)){
      w <- pmin(1 / c(1, s_pool[-length(s_pool)]),
                1 / s_pool[max(which(s_sum$time < t_star))])
    }
    else {
      w <- pmin(1 / c(1, s_pool[-length(s_pool)]),
                s_star)
    }
  }

  ### produce risk table
  rt_row <- function(t,
                     df,
                     trt_colname,
                     time_colname,
                     event_colname){

    df1 = df[df[[trt_colname]] == unique(df[[trt_colname]])[1],]
    df2 = df[df[[trt_colname]] == unique(df[[trt_colname]])[2],]

    data.frame(time = t,
               r1 = sum(df1[[time_colname]] >= t),
               r2 = sum(df2[[time_colname]] >= t),
               e1 = sum(df1[[time_colname]] == t & df1[[event_colname]] == 1),
               e2 = sum(df2[[time_colname]] == t & df2[[event_colname]] == 1))
  }

  #########################################
  rt <- purrr::map_df(unique(df[[time_colname]]),
                      rt_row,
                      df = df,
                      trt_colname = trt_colname,
                      time_colname = time_colname,
                      event_colname = event_colname)

  rt <- dplyr::filter(dplyr::arrange(rt, time),
                      e1 > 0 | e2 > 0)


  ## test statistics
  u = sum(w * (rt$e2 - (rt$e1+rt$e2) * rt$r2 / (rt$r1+rt$r2)))
  v_u = sum(w^2 * rt$r1 * rt$r2 * (rt$e1+rt$e2) * (rt$r1+rt$r2 - (rt$e1+rt$e2)) / ((rt$r1+rt$r2)^2 * (rt$r1+rt$r2-1)), na.rm = TRUE)

  z = u / sqrt(v_u)

  data.frame(u = u,
             v_u = v_u,
             z = z,
             o_minus_e_trt = unique(df[[trt_colname]])[2])

}
