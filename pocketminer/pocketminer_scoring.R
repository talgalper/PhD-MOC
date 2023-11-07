### Takes pocketminer results as input and converts into a formatted data frame ###

pocketminer_results <- list.files("results/txt_results/")

# Function to calculate the average of a residue and its 10 neighboring residues
calculate_average <- function(data, residue_index) {
  # Define the window size
  window_size <- 10
  
  # Calculate the starting and ending indices of the window
  start_index <- max(1, residue_index - window_size)
  end_index <- min(nrow(data), residue_index + window_size)
  
  # Extract the values within the window
  window_values <- data$scores[start_index:end_index]
  
  # Calculate the average of the window values
  average_value <- mean(window_values)
  
  return(average_value)
}


# loop the function over all the files 
for (i in seq_along(pocketminer_results)) {
  file <- pocketminer_results[i]
  
  print(paste0("Scoring file ", file, " : ", i, " of ", length(pocketminer_results)))
  
  data <- read.table(file.path("results/txt_results/", file), header = F)
  colnames(data)[1] <- "scores"
  
  # Create an empty vector to store the calculated averages
  averages <- numeric(nrow(data))
  
  # Calculate the averages for each residue and store them in the 'averages' vector
  for (x in 1:nrow(data)) {
    averages[x] <- calculate_average(data, x)
  }
  
  # Add the calculated averages to the original data frame
  data$averages <- averages
  
  write.table(data, file = file.path("results/pocket_results/", file), quote = F, row.names = F)
}


# identify the highest scoring pocket and the number of pockets above 0.7
pocketminer_averages <- list.files("results/pocket_results/")

# initialise lists
IDs <- list()
largest_averages <- list()
number_of_hits <- list()
uniprot_IDs <- list()


# create list of IDs and count scores
for (i in seq_along(pocketminer_averages)) {
  print(paste0("Progress: ", i, "/", length(pocketminer_averages)))
  file <- pocketminer_averages[i]
  
  # Get uniprot/subunit ID
  split_id <- strsplit(file, split = "-")
  split_id <- unlist(split_id)
  id <- paste0(split_id[1], "-", split_id[2])
  IDs <- append(IDs, id)
  
  uniprot_id <- split_id[1]
  uniprot_IDs <- append(uniprot_IDs, uniprot_id)
  
  # read in the data
  data <- read.table(file.path("results/pocket_results", file), header = T)
  
  # Find the largest value in the "averages" column
  largest_average <- max(data$averages)
  largest_averages <- append(largest_averages, largest_average)
  
  # Count the number of values above 0.7 in the "averages" column
  above_0.7_count <- sum(data$averages > 0.7)
  number_of_hits <- append(number_of_hits, above_0.7_count)
}


# combine data into data frame
pocketminer_master <- data.frame(ID = unlist(IDs),
                                 uniprot_id = unlist(uniprot_IDs),
                                 max_hit = unlist(largest_averages),
                                 num_hits = unlist(number_of_hits))

# save results
write.csv(pocketminer_master, "results/pocketminer_results.csv", row.names = F)






