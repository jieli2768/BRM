library(httr)
library(jsonlite)

testID <- "example1"


# first call to /CAT
administeredItems <- data.frame()

json_body <- toJSON(list(testID = testID, administeredItems = administeredItems, sessionData = NULL), auto_unbox = TRUE)
r <- POST("http://localhost:8080/cat", body = json_body, encode = "raw", add_headers(accept="application/json"), content_type("application/json"))
d <- fromJSON(content(r, as="text", encoding="UTF-8"))
d


# repeat call to /CAT until test complete d$testComplelted=T
while (d$status=="OK" & d$testCompleted==F){                                                  
  
  n_items <- length(d$nextItemsToAdminister[[1]])                                             # number of items to administer
  d$nextItemsToAdminister$itemScore = round(runif(n_items))                                   # assign random 0 or 1 scores
  print(d$nextItemsToAdminister)
  
  administeredItems <- rbind(administeredItems, as.data.frame(d$nextItemsToAdminister))       # append to administeredItems
  
  json_body <- toJSON(list(testID = testID, administeredItems=administeredItems, sessionData = d$sessionData), auto_unbox = TRUE)    # always post back sessionData
  r <- POST("http://localhost:8080/cat", body = json_body, encode = "raw", add_headers(accept="application/json"), content_type("application/json"))
  d <- fromJSON(content(r, as="text", encoding="UTF-8"))
  
}


# score
d$score$overall





# GET /{testID}, retrieve item bank
# r <- GET(url = paste("http://localhost:8080/", testID, sep=""))
# d <- fromJSON(content(r, as="text", encoding="UTF-8"))
# head(d$item_bank)


