
categorize_ship_types <- function(type_names){
  
  type_names <- as.numeric(type_names)
  
  group_tags <- list(
    "Cargo" = c(70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 1003, 1004, 1016),
    #"Tanker" = "Tanker|tanker",
    "Tug_Tow" = c(21,22,31,32,1023,1025),
    "Fishing" = c(30,1001,1002),
    #"Military" = "Military|military",
    "Passenger" = c(60,61,62,63,64,65,66,67,68,69,1012,1013,1014,1015),
    "Pleasure" = c(36,37,1019)
  )
  
  out_categories <- rep("Other", length(type_names))
  
  for(n in seq(length(group_tags))){
    
    out_categories[type_names %in% group_tags[[n]]] <- names(group_tags)[n]
    
  }
  
  return(out_categories)
}



