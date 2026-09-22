cat("\014")   # Clear the console
library(httr)
library(jsonlite)


# initialize variables 
testID <- "example1"
administeredItems <- data.frame()

init_sessionData <- list(
  bayesian_theta = list(
    # individualized prior for Bayesian CAT
    vector = c(0.1001, 0.2512, 2.806, -1.924, -2.653, -0.053, 1.741, 1.421, -0.535, 1.696, 1.789, 2.083, -0.695, 1.795, 0.941, -1.279, -0.392, 3.241, 4.403, -4.081, -5.018, -0.07, 3.849, -1.145, -3.644, -1.38, 2.974, -3.531, 5.054, -2.21, -1.212, 3.234, 2.733, -0.172, 2.406, -4.526, 1.461, -0.898, 2.343, 5.058, -0.664, 5.679, 0.675, -1.914, 0.326, -4.554, 0.917, 1.964, 1.614, -2.718, -1.304, -0.629, 0.24, 0.862, 4.979, -3.659, 2.677, 6.963, 2.253, -3.721, -1.664, -1.32, 3.974, -0.142, -2.825, 2.912, -4.246, -3.784, 0.31, 1.686, -1.642, 5.319, -4.069, -0.372, 4.397, -1.508, -1.427, -0.046, -2.191, -0.779, -4.296, -2.707, 0.804, 1.323, 1.63, -4.616, -1.521, -1.137, 2.47, 5.162, 0.598, 1.104, 0.635, -1.036, -3.135, -4.069, 3.17, 0.516, 3.182, 3.148, 2.74, 3.008, -0.109, 0.678, 3.659, 0.575, -0.066, 7.729, 3.049, -5.675, 5.154, -1.795, -3.947, -1.32, -1.338, 0.599, 2.485, -3.905, 1.644, 3.392, -0.51, 0.779, 0.748, 0.324, -7.357, 4.553, -2.585, -4.838, -1.01, -2.277, 2.012, -2.66, 0.796, -0.017, 5.66, 1.383, -1.396, 7.173, -2.694, -1.149, -5.679, -2.105, -3.918, -6.304, 1.115, 6.834, 2.552, 2.57, 0.439, -0.503, -1.995, -4.598, -0.126, 1.271, 8.051, 1.799, 6.543, -1.038, -3.78, 0.418, -3.453, 0.937, 0.093, 4.34, -1.446, 4.909, -0.612, -5.95, 7.042, 1.022, 5.707, 5.882, 0.885, -1.596, 2.398, -1.215, -3.884, -1.654, 0.308, -4.805, 1.265, 2.918, -4.312, 1.365, -1.963, 2.716, -5.862, -1.14, -4.224, 1.039, -1.371, 0.973, -0.357, 1.723, 0.076, -0.435, 5.151, 0.908, -3.494, 1.056, 2.313, -0.065, 3.777, -0.563, 7.913, -1.979, -1.439, 2.045, 4.875, -1.478, -0.285, -5.727, -1.005, -1.076, -1.542, 1.363, 0.855, 3.646, -6.094, -4.433, -3.227, -0.122, -4.03, 2.486, 4.239, 3.461, -0.927, 3.886, 1.248, -0.779, 1.881, 0.164, 0.968, 1.996, 1.485, 2.573, 5.976, -1.12, -7.942, 0.747, -2.038, -0.879, -1.788, 2.287, -5.973, 1.243, -1.461, -0.053, 0.382, -0.29, 0.809, 1.164, 0.98, -2.413, 2.795, -0.588, -0.795, 4.312, 1.553, 1.796, -1.162, -5.624, 1.047, 6.345, 0.224, 1.547, -0.601, -4.434, 2.389, 5.615, 3.814, -3.905, -0.663, 0.079, -3.484, -2.205, -3.648, -3.012, -1.677, -4.613, 2.2, 2.263, 1.37, -3.005, -5.231, -3.912, -3.546, -2.501, -3.895, 4.036, -0.713, 5.962, 4.893, 2.451, -5.299, -4.272, -0.499, 2.109, 3.53, -4.255, -1.208, 0.039, 3.974, 1.287, -2.468, -2.896, 2.568, -2.679, -4.41, 0.324, -3.898, 3.381, -0.306, 0.874, 2.531, -0.39, 1.771, -2.809, 1.786, -4.307, -3.197, 5.215, 0.032, -3.689, -1.177, -1.222, 3.862, 0.437, -3.577, -6.366, 2.645, -2.215, -3.305, 3.605, 4.922, 1.869, -1.268, 2.04, 2.681, 0.18, 2.788, 4.698, -1.603, -2.512, 2.005, -7.186, 3.146, 4.147, -2.594, -2.571, 1.106, 1.715, -2.03, 1.496, -2.487, -1.381, -1.909, 0.143, 0.253, 5.361, 8.337, -2.726, 1.966, -4.685, 3.029, -5.708, -2.969, 1.437, -2.446, 3.68, -0.464, 2.599, -1.155, 2.16, 1.784, -4.445, 1.544, -3.937, -5.685, 1.782, -1.201, -4.268, -2.578, 0.552, 2.623, 0.997, 4.601, -5.381, -4.145, -2.177, 1.514, 3.596, 2.761, -1.422, 2.448, -1.08, -2.526, 3.05, -0.268, -3.032, -0.641, -2.227, -3.064, 1.316, -2.094, 1.714, -1.239, -2.561, 0.812, 4.827, 2.325, -1.257, 2.244, 1.703, -3.346, 2.135, -2.515, -3.99, -3.87, -4.659, 2.428, 1.812, 6.52, -2.794, 3.064, 0.209, -0.478, 2.114, -1.943, 0.837, -1.453, 3.421, -2.342, 11.137, 2.177, 0.657, -2.281, -4.117, 0.096, -0.592, -0.99, 3.952, 5.517, -3.83, 5.703, -4.874, -1.727, -4.457, 5.257, 0.567, -3.094, 4.765, 3.215, 3.397, 2.271, -2.519, 0.396, -2.498, 1.824, -3.433, 3.114, 0.395, -2.961, 1.712, 0.844, 0.91, 3.135, 8.395, 0.415, -0.445, -2.248, -2.463, 2.329, -1.611, 2.997, -1.148, -0.744, -3.182, 3.233, -1.343, 7.655, -1.902, -5.135, -4.136, -4.039, -2.348, 3.195, -3.33, 1.645, -1.271, -1.923, -2.717, 5.677, 4.019, 4.096, 0.997, 0.931, 0.023, -1.378, 1.305)        
  )
)


# create json body to call /CAT end-point
json_body <- toJSON(list(testID = testID, administeredItems = administeredItems, sessionData = init_sessionData), auto_unbox = TRUE)
# cat(prettify(json_body))


# make first call to /CAT end-point
r <- POST("http://localhost:8080/cat", body = json_body, encode = "raw", add_headers(accept="application/json"), content_type("application/json"))
d <- fromJSON(content(r, as="text", encoding="UTF-8"))
# d




# repeat call to /CAT until test complete d$testComplelted=T
call_counter = 1

while (d$status=="OK" & d$testCompleted==F){                                                  
  
  cat("-------------------------------------------------------------------------------------\n")
  cat(paste("call_counter = ", call_counter), "\n")
  cat(paste("bayesian vector=", paste(head(d$sessionData$bayesian_theta$vector), collapse = " "), "...", "\n"))
  cat(paste("bayesian mean=", mean(d$sessionData$bayesian_theta$vector), "\n")) 
  cat(paste("bayesian sd=", sd(d$sessionData$bayesian_theta$vector), "\n"))
  cat("\n")
  
  n_items <- length(d$nextItemsToAdminister[[1]])                                             # number of the next items to administer
  d$nextItemsToAdminister$itemScore = round(runif(n_items))                                   # assign random 0 or 1 scores

  cat("next items to adminster and their responses:", "\n")
  print(d$nextItemsToAdminister)                                                                
  
  administeredItems <- rbind(administeredItems, as.data.frame(d$nextItemsToAdminister))       # append to administeredItems
  
  json_body <- toJSON(list(testID = testID, administeredItems=administeredItems, sessionData = d$sessionData), auto_unbox = TRUE)    # always post back sessionData
  r <- POST("http://localhost:8080/cat", body = json_body, encode = "raw", add_headers(accept="application/json"), content_type("application/json"))
  d <- fromJSON(content(r, as="text", encoding="UTF-8"))
  
  call_counter = call_counter + 1
  
}


# final test score
d$status
d$testCompleted
d$sessionData$bayesian_theta
d$score$overall




# To retrieve item bank
# r <- GET(url = paste("http://localhost:8080/", testID, sep=""))
# d <- fromJSON(content(r, as="text", encoding="UTF-8"))
# head(d$item_bank)




#rescore
json_body <- toJSON(list(testID = testID, administeredItems=administeredItems, sessionData = init_sessionData), auto_unbox = TRUE)    # set sessionData to init_sessionData
r <- POST("http://localhost:8080/rescore", body = json_body, encode = "raw", add_headers(accept="application/json"), content_type("application/json"))
d <- fromJSON(content(r, as="text", encoding="UTF-8"))
d$score$bayesian_theta

