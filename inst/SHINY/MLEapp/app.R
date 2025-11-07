#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)

# app.R
library(shiny)
library(ggplot2)
library(reshape2)

# ---------------------------
# Helper: negative log-likelihoods (on transformed params)
# transforms:
#  - unconstrained params returned by optim: we convert to actual params
# ---------------------------
unpack_params <- function(par, dist, extra = list()){
  # returns a named list of actual params and a function to compute Jacobian if needed
  if(dist == "normal"){
    mu <- par[1]
    sigma <- exp(par[2])   # par[2] is log(sigma)
    list(params = list(mu = mu, sigma = sigma),
         inv_transform = function(p) c(p$mu, log(p$sigma)))
  } else if(dist == "exponential"){
    rate <- exp(par[1])
    list(params = list(rate = rate),
         inv_transform = function(p) log(p$rate))
  } else if(dist == "poisson"){
    lambda <- exp(par[1])
    list(params = list(lambda = lambda),
         inv_transform = function(p) log(p$lambda))
  } else if(dist == "binomial"){
    # par is logit(p)
    p <- plogis(par[1])
    size <- extra$size
    list(params = list(p = p, size = size),
         inv_transform = function(p) qlogis(p$p))
  } else if(dist == "gamma"){
    shape <- exp(par[1])
    rate  <- exp(par[2])
    list(params = list(shape = shape, rate = rate),
         inv_transform = function(p) c(log(p$shape), log(p$rate)))
  } else if(dist == "beta"){
    a <- exp(par[1])
    b <- exp(par[2])
    list(params = list(a = a, b = b),
         inv_transform = function(p) c(log(p$a), log(p$b)))
  } else stop("Unknown dist")
}

nll_fun <- function(par, x, dist, extra = list()){
  up <- unpack_params(par, dist, extra)$params
  if(dist == "normal"){
    -sum(dnorm(x, mean = up$mu, sd = up$sigma, log = TRUE))
  } else if(dist == "exponential"){
    -sum(dexp(x, rate = up$rate, log = TRUE))
  } else if(dist == "poisson"){
    -sum(dpois(x, lambda = up$lambda, log = TRUE))
  } else if(dist == "binomial"){
    -sum(dbinom(x, size = extra$size, prob = up$p, log = TRUE))
  } else if(dist == "gamma"){
    -sum(dgamma(x, shape = up$shape, rate = up$rate, log = TRUE))
  } else if(dist == "beta"){
    -sum(dbeta(x, shape1 = up$a, shape2 = up$b, log = TRUE))
  } else stop("Unknown dist")
}

# compute MLE via optim, and transform back plus SE via Hessian and delta method
compute_mle <- function(x, dist, extra = list()){
  # initial guesses by method-of-moments where possible, then transform to unconstrained space
  if(dist == "normal"){
    mu0 <- mean(x)
    sd0 <- sd(x)
    par0 <- c(mu0, log(sd0))
    lower <- -Inf; upper <- Inf
  } else if(dist == "exponential"){
    rate0 <- 1/mean(x)
    par0 <- log(rate0)
  } else if(dist == "poisson"){
    lam0 <- mean(x)
    par0 <- log(lam0)
  } else if(dist == "binomial"){
    # x are counts; size must be in extra$size
    size <- extra$size
    p0 <- mean(x)/size
    p0 <- min(max(p0, 1e-6), 1-1e-6)
    par0 <- qlogis(p0)
  } else if(dist == "gamma"){
    m <- mean(x); v <- var(x)
    shape0 <- max(m^2 / v, 1e-6)
    rate0  <- max(m / v, 1e-6)
    par0 <- c(log(shape0), log(rate0))
  } else if(dist == "beta"){
    m <- mean(x); v <- var(x)
    # method of moments for alpha, beta
    tmp <- (m*(1 - m) / v) - 1
    a0 <- max(m * tmp, 1e-6)
    b0 <- max((1-m) * tmp, 1e-6)
    par0 <- c(log(a0), log(b0))
  } else stop("Unknown dist")
  
  res <- tryCatch({
    optim(par = par0,
          fn = function(p) nll_fun(p, x, dist, extra),
          method = "BFGS",
          hessian = TRUE,
          control = list(maxit = 10000))
  }, error = function(e) list(convergence = 1, message = e$message))
  
  if(is.list(res) && res$convergence != 0){
    return(list(success = FALSE, res = res, message = "optim failed or did not converge"))
  }
  
  # transform back
  up <- unpack_params(res$par, dist, extra)$params
  # get covariance on transformed (optim's) scale
  cov_trans <- tryCatch(solve(res$hessian), error = function(e) {
    # near-singular; return NA
    matrix(NA, nrow = length(res$par), ncol = length(res$par))
  })
  
  # delta method: compute Jacobian of mapping from transformed par to actual params
  # For each dist define jacobian J (actual_params x transformed_params)
  J <- switch(dist,
              normal = {
                mu <- 1; # d(mu)/d(mu) = 1
                # d(sigma)/d(log_sigma) = exp(log_sigma) = sigma
                sigma <- up$sigma
                matrix(c(1, 0,
                         0, sigma), nrow = 2, byrow = TRUE)
              },
              exponential = {
                rate <- up$rate
                matrix(rate, nrow = 1, ncol = 1)
              },
              poisson = {
                lambda <- up$lambda
                matrix(lambda, nrow = 1, ncol = 1)
              },
              binomial = {
                # param: p = plogis(q); dq -> dp = p(1-p)
                p <- up$p
                matrix(p * (1 - p), nrow = 1, ncol = 1)
              },
              gamma = {
                shape <- up$shape; rate <- up$rate
                matrix(c(shape, 0,
                         0, rate), nrow = 2, byrow = TRUE)
              },
              beta = {
                a <- up$a; b <- up$b
                matrix(c(a, 0,
                         0, b), nrow = 2, byrow = TRUE)
              }
  )
  
  # covariance on actual parameter scale
  cov_actual <- tryCatch({
    J %*% cov_trans %*% t(J)
  }, error = function(e) {
    matrix(NA, nrow = nrow(J), ncol = nrow(J))
  })
  
  ses <- sqrt(diag(cov_actual))
  ests <- unlist(up)
  names(ests) <- names(up)
  
  # return list
  list(success = TRUE,
       estimates = ests,
       ses = ses,
       cov_actual = cov_actual,
       loglik = -res$value,
       optim_res = res)
}

# UI
ui <- fluidPage(
  titlePanel("MLE demo for several univariate distributions"),
  sidebarLayout(
    sidebarPanel(
      selectInput("dist", "Choose distribution:",
                  choices = c("Normal" = "normal",
                              "Exponential" = "exponential",
                              "Poisson" = "poisson",
                              "Binomial" = "binomial",
                              "Gamma" = "gamma",
                              "Beta" = "beta")),
      numericInput("n", "Sample size (N):", value = 100, min = 10, step = 10),
      # parameter inputs - dynamic
      uiOutput("param_ui"),
      numericInput("seed", "Random seed (for reproducibility, empty = random):", value = NA),
      actionButton("gen", "Generate sample & fit MLE"),
      hr(),
      helpText("Notes: MLE computed with optim(); standard errors via Hessian + delta method.")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Sample & Fit",
                 br(),
                 fluidRow(column(6, plotOutput("mainPlot")),
                          column(6,
                                 tableOutput("estTable"),
                                 verbatimTextOutput("loglikText")
                          )),
                 hr(),
                 uiOutput("likUI")
        ),
        tabPanel("Diagnostics",
                 verbatimTextOutput("optimSummary"),
                 plotOutput("likProfile", height = "500px"))
      )
    )
  )
)

# Server
server <- function(input, output, session){
  
  # dynamic parameter UI
  output$param_ui <- renderUI({
    switch(input$dist,
           normal = tagList(
             numericInput("true_mu", "True mu:", value = 0),
             numericInput("true_sigma", "True sigma (>0):", value = 1, min = 1e-6)
           ),
           exponential = tagList(
             numericInput("true_rate", "True rate (>0):", value = 1, min = 1e-6)
           ),
           poisson = tagList(
             numericInput("true_lambda", "True lambda (>0):", value = 3, min = 1e-6)
           ),
           binomial = tagList(
             numericInput("true_p", "True p (0-1):", value = 0.3, min = 1e-6, max = 1-1e-6),
             numericInput("bin_size", "Binomial size (trials per obs):", value = 10, min = 1, step = 1)
           ),
           gamma = tagList(
             numericInput("true_shape", "True shape (>0):", value = 2, min = 1e-6),
             numericInput("true_rate_g", "True rate (>0):", value = 1, min = 1e-6)
           ),
           beta = tagList(
             numericInput("true_a", "True alpha (>0):", value = 2, min = 1e-6),
             numericInput("true_b", "True beta (>0):", value = 5, min = 1e-6)
           )
    )
  })
  
  data_and_fit <- eventReactive(input$gen, {
    if(!is.na(input$seed)) set.seed(input$seed)
    N <- input$n
    dist <- input$dist
    extra <- list()
    x <- NULL
    true_params <- list()
    if(dist == "normal"){
      x <- rnorm(N, mean = input$true_mu, sd = input$true_sigma)
      true_params <- list(mu = input$true_mu, sigma = input$true_sigma)
    } else if(dist == "exponential"){
      x <- rexp(N, rate = input$true_rate)
      true_params <- list(rate = input$true_rate)
    } else if(dist == "poisson"){
      x <- rpois(N, lambda = input$true_lambda)
      true_params <- list(lambda = input$true_lambda)
    } else if(dist == "binomial"){
      size <- input$bin_size
      x <- rbinom(N, size = size, prob = input$true_p)
      true_params <- list(p = input$true_p, size = size)
      extra$size <- size
    } else if(dist == "gamma"){
      x <- rgamma(N, shape = input$true_shape, rate = input$true_rate_g)
      true_params <- list(shape = input$true_shape, rate = input$true_rate_g)
    } else if(dist == "beta"){
      x <- rbeta(N, shape1 = input$true_a, shape2 = input$true_b)
      true_params <- list(a = input$true_a, b = input$true_b)
    }
    
    fit <- compute_mle(x = x, dist = dist, extra = extra)
    list(x = x, fit = fit, true = true_params, extra = extra, dist = dist)
  })
  
  output$mainPlot <- renderPlot({
    dat <- data_and_fit()
    x <- dat$x
    fit <- dat$fit
    dist <- dat$dist
    df <- data.frame(x = x)
    if(dist %in% c("normal", "exponential", "gamma", "beta")){
      # continuous: histogram + fitted density
      p <- ggplot(df, aes(x = x)) + geom_histogram(aes(y = ..density..), bins = 30, fill = NA, color = "black")
      if(fit$success){
        params <- fit$estimates
        if(dist == "normal"){
          p <- p + stat_function(fun = function(z) dnorm(z, mean = params["mu"], sd = params["sigma"]), size = 1)
        } else if(dist == "exponential"){
          p <- p + stat_function(fun = function(z) dexp(z, rate = params["rate"]), size = 1)
        } else if(dist == "gamma"){
          p <- p + stat_function(fun = function(z) dgamma(z, shape = params["shape"], rate = params["rate"]), size = 1)
        } else if(dist == "beta"){
          p <- p + stat_function(fun = function(z) dbeta(z, shape1 = params["a"], shape2 = params["b"]), size = 1)
        }
      }
      p + theme_minimal() + labs(title = "Sample histogram with fitted density")
    } else {
      # discrete: barplot of counts and fitted pmf
      counts <- as.data.frame(table(x))
      counts$xnum <- as.numeric(as.character(counts$x))
      counts$Freq <- counts$Freq / sum(counts$Freq)
      p <- ggplot(counts, aes(x = xnum, y = Freq)) + geom_bar(stat = "identity", fill = NA, color = "black") +
        theme_minimal() + labs(x = "Value", y = "Relative frequency", title = "Sample frequencies and fitted pmf")
      if(fit$success){
        params <- fit$estimates
        if(dist == "poisson"){
          lambda <- params["lambda"]
          support <- min(counts$xnum):max(counts$xnum)
          pmf <- dpois(support, lambda)
          dfp <- data.frame(xnum = support, pmf = pmf)
          p <- p + geom_point(data = dfp, aes(x = xnum, y = pmf), color = "red", size = 2)
        } else if(dist == "binomial"){
          size <- dat$true$size
          p_est <- params["p"]
          support <- 0:size
          pmf <- dbinom(support, size = size, prob = p_est)
          dfp <- data.frame(xnum = support, pmf = pmf)
          p <- p + geom_point(data = dfp, aes(x = xnum, y = pmf), color = "red", size = 2)
        }
      }
      p
    }
  })
  
  output$estTable <- renderTable({
    dat <- data_and_fit()
    fit <- dat$fit
    if(!fit$success){
      return(data.frame(Message = fit$message))
    }
    ests <- fit$estimates
    ses <- fit$ses
    param_names <- names(ests)
    ci_low <- ests - 1.96 * ses
    ci_high <- ests + 1.96 * ses
    df <- data.frame(Parameter = param_names,
                     Estimate = round(ests, 6),
                     SE = round(ses, 6),
                     CI_lower = round(ci_low, 6),
                     CI_upper = round(ci_high, 6),
                     row.names = NULL)
    # show true params if available
    true <- dat$true
    for(i in seq_len(nrow(df))){
      pn <- df$Parameter[i]
      if(pn %in% names(true)) df$True[i] <- true[[pn]]
      else df$True[i] <- NA
    }
    df <- df[, c("Parameter", "True", "Estimate", "SE", "CI_lower", "CI_upper")]
    df
  })
  
  output$loglikText <- renderText({
    dat <- data_and_fit()
    fit <- dat$fit
    if(!fit$success) return("No fit")
    paste("Log-likelihood at MLE:", round(fit$loglik, 6))
  })
  
  output$optimSummary <- renderPrint({
    dat <- data_and_fit()
    fit <- dat$fit
    if(!fit$success){
      print(fit$res)
    } else {
      print(fit$optim_res)
    }
  })
  
  # Likelihood surface or profile
  output$likUI <- renderUI({
    dat <- data_and_fit()
    fit <- dat$fit
    if(!fit$success) return(NULL)
    if(dat$dist %in% c("normal", "gamma", "beta")){
      helpText("Below: likelihood surface (contours) over a grid of two parameters")
    } else {
      helpText("Below: profile likelihood over a range for the single parameter")
    }
  })
  
  output$likProfile <- renderPlot({
    dat <- data_and_fit()
    fit <- dat$fit
    x <- dat$x
    dist <- dat$dist
    extra <- dat$extra
    if(!fit$success) return(NULL)
    
    if(dist %in% c("normal", "gamma", "beta")){
      # 2D grid: vary first param +/- 40% and second param +/- 40% of estimates
      ests <- fit$estimates
      if(dist == "normal"){
        mu0 <- ests["mu"]; sigma0 <- ests["sigma"]
        mu_seq <- seq(mu0 - 1.5 * sd(x), mu0 + 1.5 * sd(x), length.out = 80)
        sigma_seq <- exp(seq(log(max(1e-6, sigma0/3)), log(sigma0*3), length.out = 80))
        grid <- expand.grid(mu = mu_seq, sigma = sigma_seq)
        grid$nll <- mapply(function(mu, sigma){
          -sum(dnorm(x, mean = mu, sd = sigma, log = TRUE))
        }, grid$mu, grid$sigma)
        ggplot(grid, aes(x = mu, y = sigma, z = -nll)) +
          geom_contour_filled() + scale_y_log10() +
          geom_point(aes(x = mu0, y = sigma0), color = "red", size = 3) +
          labs(title = "Log-likelihood surface (Normal)", x = "mu", y = "sigma (log scale)") +
          theme_minimal()
      } else if(dist == "gamma"){
        shape0 <- ests["shape"]; rate0 <- ests["rate"]
        shape_seq <- exp(seq(log(max(1e-4, shape0/3)), log(shape0*3), length.out = 80))
        rate_seq <- exp(seq(log(max(1e-4, rate0/3)), log(rate0*3), length.out = 80))
        grid <- expand.grid(shape = shape_seq, rate = rate_seq)
        grid$nll <- mapply(function(sh, ra){
          -sum(dgamma(x, shape = sh, rate = ra, log = TRUE))
        }, grid$shape, grid$rate)
        ggplot(grid, aes(x = shape, y = rate, z = -nll)) +
          geom_contour_filled() + scale_x_log10() + scale_y_log10() +
          geom_point(aes(x = shape0, y = rate0), color = "red", size = 3) +
          labs(title = "Log-likelihood surface (Gamma)", x = "shape (log)", y = "rate (log)") +
          theme_minimal()
      } else if(dist == "beta"){
        a0 <- ests["a"]; b0 <- ests["b"]
        a_seq <- exp(seq(log(max(1e-4, a0/3)), log(a0*3), length.out = 80))
        b_seq <- exp(seq(log(max(1e-4, b0/3)), log(b0*3), length.out = 80))
        grid <- expand.grid(a = a_seq, b = b_seq)
        grid$nll <- mapply(function(a, b) -sum(dbeta(x, a, b, log = TRUE)), grid$a, grid$b)
        ggplot(grid, aes(x = a, y = b, z = -nll)) +
          geom_contour_filled() + scale_x_log10() + scale_y_log10() +
          geom_point(aes(x = a0, y = b0), color = "red", size = 3) +
          labs(title = "Log-likelihood surface (Beta)", x = "alpha", y = "beta") +
          theme_minimal()
      }
    } else {
      # 1D profile plot
      ests <- fit$estimates
      if(dist == "exponential"){
        rate0 <- ests["rate"]
        seqr <- seq(max(1e-6, rate0/5), rate0*5, length.out = 200)
        ll <- sapply(seqr, function(r) sum(dexp(x, rate = r, log = TRUE)))
        df <- data.frame(param = seqr, ll = ll)
        ggplot(df, aes(x = param, y = ll)) + geom_line() +
          geom_vline(xintercept = rate0, color = "red") + labs(x = "rate", y = "log-likelihood") + theme_minimal()
      } else if(dist == "poisson"){
        lam0 <- ests["lambda"]
        seqr <- seq(max(1e-6, lam0/4), lam0*4, length.out = 200)
        ll <- sapply(seqr, function(l) sum(dpois(x, lambda = l, log = TRUE)))
        df <- data.frame(param = seqr, ll = ll)
        ggplot(df, aes(x = param, y = ll)) + geom_line() +
          geom_vline(xintercept = lam0, color = "red") + labs(x = "lambda", y = "log-likelihood") + theme_minimal()
      } else if(dist == "binomial"){
        p0 <- ests["p"]
        seqr <- seq(1e-5, 1-1e-5, length.out = 200)
        ll <- sapply(seqr, function(p) sum(dbinom(x, size = dat$true$size, prob = p, log = TRUE)))
        df <- data.frame(param = seqr, ll = ll)
        ggplot(df, aes(x = param, y = ll)) + geom_line() +
          geom_vline(xintercept = p0, color = "red") + labs(x = "p", y = "log-likelihood") + theme_minimal()
      }
    }
  })
  
}

shinyApp(ui, server)


# Run the application 
shinyApp(ui = ui, server = server)
