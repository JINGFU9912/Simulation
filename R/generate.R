library(lme4) 
library(dplyr) 
####Function belows ####
# Enumerate feasible (G, R) combinations under budget constraint
get_feasible_combinations <- function(B, c1, c2, max_G = 1000, max_R = 1000, min_G = 10, min_R = 2) {
    feasible <- expand.grid(G = min_G:max_G, R = min_R:max_R) %>%
        mutate(cost = G * c1 + G * (R - 1) * c2) %>%
        filter(cost == B)
    return(feasible)
}

#Simulate
generate_data <- function(alpha, beta, gamma2, sigma2, G, R, distribution = "normal") {
    # Generate cluster-level random effects/ epsilon
    cluster_effect <- rnorm(G, mean = 0, sd = sqrt(gamma2)) #epsilon
    treatment <- rbinom(G, 1, 0.5)
    if(length(unique(treatment)) == 1){
        ## avoid all same treatment
        k <- unique(treatment)
        treatment[sample(G,1,prob = rep(1, G))] <- 1 - k
    }
    data <- do.call(rbind, lapply(1:G, function(i) {
        mu <- if (distribution == "normal") {
            alpha + beta * treatment[i] + cluster_effect[i]
        } else {
            exp(alpha + beta * treatment[i] + cluster_effect[i]) 
        }
        y <- if (distribution == "normal") {
            rnorm(R, mean = mu, sd = sqrt(sigma2))
        } else {
            rpois(R, lambda = mu)
        }
        data.frame(cluster = i, treatment = treatment[i], y = y)
    }))
    
    return(data)
}

# Calculate Estimate Parameters
fit_model <- function(data, distribution = "normal") {
    if (distribution == "normal") {
        model <- lmer(y ~ treatment + (1 | cluster), data = data)
    } else {
        # glm for Poisson
        model <- glmer(y ~ treatment + (1 | cluster), data = data, family = poisson(link = "log"))
    }
    
    est_beta <- fixef(model)["treatment"]  # Extract treatment effect estimate
    se_beta <- sqrt(vcov(model)["treatment", "treatment"])  # Extract standard error
    return(c(beta_hat = est_beta, se_beta = se_beta))
}

#Calculate Performance Measures
calculate_performance <- function(sim_results, true_beta) {
    bias <- mean(sim_results$beta_hat) - true_beta
    se <- sd(sim_results$beta_hat)
    coverage <- mean((sim_results$beta_hat - 1.96 * sim_results$se_beta <= true_beta) &
                         (sim_results$beta_hat + 1.96 * sim_results$se_beta >= true_beta))
    return(data.frame(bias = bias, se = se, coverage = coverage))
}

#Run Simulation Across Parameter Combinations
run_simulation <- function(sim_params, c1, c2,  n_sim = 1000, distribution = "normal") {
    feasible_combinations <- get_feasible_combinations(sim_params$B, c1, c2)
    results <- list()
    param_grid <- expand.grid(sim_params$alpha, sim_params$beta, sim_params$gamma2, sim_params$sigma2)
    param_grid <- do.call(rbind, lapply(1:nrow(feasible_combinations), function(i) {
        cbind(param_grid, G = feasible_combinations$G[i], R = feasible_combinations$R[i])
    }))
    colnames(param_grid) <- c("alpha", "beta", "gamma2", "sigma2", "G","R")
    
    for (i in 1:nrow(param_grid)) {
        params <- param_grid[i, ]
        sim_results <- replicate(n_sim, {
            data <- generate_data(params$alpha, params$beta, params$gamma2,
                                  params$sigma2, params$G, params$R, distribution = distribution)
            fit_model(data, distribution = distribution)
        }, simplify = FALSE)
        sim_results <- do.call(rbind, sim_results)
        colnames(sim_results) <- c("beta_hat","se_beta")
        sim_results <- data.frame(sim_results)
        performance <- calculate_performance(sim_results, params$beta)
        results[[i]] <- c(params, performance)
    }
    
                            
    return(do.call(rbind, results))
}


####Run in each senario ####
sim_params <- list(
    alpha = c(5, 10),  # Baseline mean
    beta = c(0.1, 0.5),  # Treatment effect
    gamma2 = c(0.1, 0.5),  # Cluster random effect variance
    sigma2 = c(1, 2),  # Within-cluster variance
    B = 1000
)

c1_c2_ratios <- list(
    c(20, 5),  
    c(30, 10)
)


# For Normal Distribution
set.seed(202412)
all_results <- lapply(seq_along(c1_c2_ratios), function(i) {
    cost <- c1_c2_ratios[[i]]
    c1 <- cost[1]
    c2 <- cost[2]
    result <- run_simulation(sim_params, c1 = c1, c2 = c2, n_sim = 1000, distribution = "normal")
    result <- data.frame(result)
    result$c1 <- c1
    result$c2 <- c2
    result$ratio <- c1 / c2
    return(result)
})
all_results <- do.call(rbind, all_results)
all_results<- as.data.frame(lapply(all_results, unlist))
write.csv(all_results,file = "/Users/fusei/Desktop/24FALL/PHP2550/Project3/normal_simulation.csv")

# For Poisson Distribution
set.seed(202412)
poi_results <- lapply(seq_along(c1_c2_ratios), function(i) {
    cost <- c1_c2_ratios[[i]]
    c1 <- cost[1]
    c2 <- cost[2]
    result <- run_simulation(sim_params, c1 = c1, c2 = c2, n_sim = 1000, distribution = "poisson")
    result <- data.frame(result)
    result$c1 <- c1
    result$c2 <- c2
    result$ratio <- c1 / c2
    return(result)
})
poi_results <- do.call(rbind, poi_results)
poi_results<- as.data.frame(lapply(poi_results, unlist))
write.csv(poi_results,file = "/Users/fusei/Desktop/24FALL/PHP2550/Project3/poi_simulation.csv")

