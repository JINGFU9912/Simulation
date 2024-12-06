source('../R/generate.R')
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
