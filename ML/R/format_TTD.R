

# Read the lines from the file
lines <- readLines("data/TTD/P1-01-TTD_target_download.txt")

# Initialize variables
current_target_id <- NA
current_target_info <- list()
targets_list <- list()
drug_info_list <- list()

# Loop over each line in the file
for (line in lines) {
  # Skip empty lines
  if (line == "") next
  
  # Skip lines without tabs (e.g., headers or separators)
  if (!grepl("\t", line)) next
  
  # Split the line by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Skip lines with fewer than 3 fields
  if (length(fields) < 3) next
  
  # Extract the target ID, field name, and field values
  target_id <- fields[1]
  field_name <- fields[2]
  field_values <- fields[3:length(fields)]
  
  # Check if we're starting a new target
  if (is.na(current_target_id) || target_id != current_target_id) {
    # Save the previous target information
    if (!is.na(current_target_id)) {
      # Add the current target info to the list
      targets_list[[current_target_id]] <- current_target_info
    }
    # Initialize new target information
    current_target_id <- target_id
    current_target_info <- list()
    current_target_info$TARGETID <- target_id  # Store the TARGETID
  }
  
  if (field_name == "DRUGINFO") {
    # Handle DRUGINFO entries
    # Expected fields: TTD Drug ID, Drug Name, Highest Clinical Status
    if (length(field_values) >= 3) {
      drug_entry <- data.frame(
        TARGETID = target_id,
        TTDDrugID = field_values[1],
        DrugName = field_values[2],
        HighestClinicalStatus = field_values[3],
        stringsAsFactors = FALSE
      )
    } else {
      # Handle cases with missing fields
      warning(paste("DRUGINFO line has fewer than expected fields:", line))
      next
    }
    # Append the drug entry to the drug info list
    drug_info_list[[length(drug_info_list) + 1]] <- drug_entry
  } else {
    # Handle other fields
    field_value <- paste(field_values, collapse = "\t")
    # If the field already exists, append the new value
    if (is.null(current_target_info[[field_name]])) {
      current_target_info[[field_name]] <- field_value
    } else {
      current_target_info[[field_name]] <- paste(current_target_info[[field_name]], field_value, sep = "; ")
    }
  }
}

# Save the last target information
if (!is.na(current_target_id)) {
  targets_list[[current_target_id]] <- current_target_info
}

# Get all unique field names
all_field_names <- unique(unlist(lapply(targets_list, names)))

# Create a dataframe for target information
target_info_list <- lapply(targets_list, function(x) {
  # Add missing fields with NA
  missing_fields <- setdiff(all_field_names, names(x))
  for (mf in missing_fields) {
    x[[mf]] <- NA
  }
  # Order the fields
  x <- x[all_field_names]
  data.frame(x, stringsAsFactors = FALSE)
})

# Combine all target information into a single dataframe
target_info_df <- do.call(rbind, target_info_list)

# Combine all drug information into a single dataframe
drug_info_df <- do.call(rbind, drug_info_list)






# Read the lines from the file
lines <- readLines("data/TTD/P1-05-Drug_disease.txt")

# Initialize empty vectors to store the data
ttd_drug_ids <- c()
drug_names <- c()
indications <- c()
disease_entries <- c()
icd11_codes <- c()
clinical_statuses <- c()

# Variables to hold current drug ID and name
current_ttd_drug_id <- NA
current_drug_name <- NA

# Loop over each line in the file
for (line in lines) {
  # Skip empty lines
  if (line == "") next
  
  # Split the line by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Check the first field to determine the type of line
  if (fields[1] == "TTDDRUID") {
    # Update the current drug ID
    current_ttd_drug_id <- fields[2]
  } else if (fields[1] == "DRUGNAME") {
    # Update the current drug name
    current_drug_name <- fields[2]
  } else if (fields[1] == "INDICATI") {
    # Process the indication line
    # Handle cases where the 'Disease entry' might be missing
    if (length(fields) == 5) {
      indication <- fields[2]
      disease_entry <- fields[3]
      icd11 <- fields[4]
      clinical_status <- fields[5]
    } else if (length(fields) == 4) {
      indication <- fields[2]
      disease_entry <- NA
      icd11 <- fields[3]
      clinical_status <- fields[4]
    } else {
      # If the number of fields is unexpected, issue a warning
      warning(paste("Unexpected number of fields in line:", line))
      next
    }
    
    # Append the extracted data to the vectors
    ttd_drug_ids <- c(ttd_drug_ids, current_ttd_drug_id)
    drug_names <- c(drug_names, current_drug_name)
    indications <- c(indications, indication)
    disease_entries <- c(disease_entries, disease_entry)
    icd11_codes <- c(icd11_codes, icd11)
    clinical_statuses <- c(clinical_statuses, clinical_status)
  }
}

# Create a data frame from the vectors
drugDisease_data <- data.frame(
  TTDDRUID = ttd_drug_ids,
  DRUGNAME = drug_names,
  INDICATION = indications,
  DISEASE_ENTRY = disease_entries,
  ICD11 = icd11_codes,
  CLINICAL_STATUS = clinical_statuses,
  stringsAsFactors = FALSE
)

drugDisease_data <- drugDisease_data[-1, ]



# Read the lines from the file
lines <- readLines("data/TTD/P1-06-Target_disease.txt")

# Initialize empty vectors to store the data
target_ids <- c()
target_names <- c()
clinical_statuses <- c()
disease_entries <- c()
icd11_codes <- c()

# Variables to hold current target ID and name
current_target_id <- NA
current_target_name <- NA

# Loop over each line in the file
for (line in lines) {
  # Skip empty lines
  if (line == "") next
  
  # Split the line by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Check if there are at least 2 fields
  if (length(fields) < 2) next
  
  # The first field is the TargetID
  # The second field is the field name ('TARGETID', 'TARGNAME', 'INDICATI')
  
  # Check the field name to determine the type of line
  if (fields[2] == "TARGETID") {
    # Update the current target ID
    current_target_id <- fields[3]
  } else if (fields[2] == "TARGNAME") {
    # Update the current target name
    current_target_name <- fields[3]
  } else if (fields[2] == "INDICATI") {
    # Process the INDICATI line
    # Extract Clinical Status, Disease Entry, ICD-11 Code
    # Check if there are at least 5 fields
    if (length(fields) >= 5) {
      clinical_status <- fields[3]
      disease_entry <- fields[4]
      icd11_code <- fields[5]
    } else {
      # If missing fields, set to NA
      clinical_status <- fields[3]
      disease_entry <- fields[4]
      icd11_code <- NA
    }
    
    # Append the data to the vectors
    target_ids <- c(target_ids, current_target_id)
    target_names <- c(target_names, current_target_name)
    clinical_statuses <- c(clinical_statuses, clinical_status)
    disease_entries <- c(disease_entries, disease_entry)
    icd11_codes <- c(icd11_codes, icd11_code)
  }
}

# Create a data frame from the vectors
targetDisease_data <- data.frame(
  TARGETID = target_ids,
  TARGNAME = target_names,
  CLINICAL_STATUS = clinical_statuses,
  DISEASE_ENTRY = disease_entries,
  ICD11_CODE = icd11_codes,
  stringsAsFactors = FALSE
)





# Read the lines from your text file
lines <- readLines("data/TTD/P2-02-TTD_uniprot_successful.txt")

# Find the starting index of the data (after the header and abbreviations)
data_start_index <- which(grepl("^TARGETID\t", lines))[2]

# Extract the data lines starting from the data_start_index
data_lines <- lines[data_start_index:length(lines)]

# Initialize an empty list to store each record
record_list <- list()
record_index <- 1
i <- 1

# Loop through the data lines to separate records based on blank lines
while (i <= length(data_lines)) {
  if (data_lines[i] == "") {
    i <- i + 1  # Skip blank lines
  } else {
    record <- c()
    # Collect all lines belonging to a single record
    while (i <= length(data_lines) && data_lines[i] != "") {
      record <- c(record, data_lines[i])
      i <- i + 1
    }
    record_list[[record_index]] <- record
    record_index <- record_index + 1
  }
}

# Process each record to extract key-value pairs
records_list <- lapply(record_list, function(record) {
  key_values <- setNames(
    sapply(record, function(line) unlist(strsplit(line, "\t"))[2]),
    sapply(record, function(line) unlist(strsplit(line, "\t"))[1])
  )
  return(key_values)
})

# Get all unique keys (column names)
all_keys <- unique(unlist(lapply(records_list, names)))

# Convert the list of records into a data frame
TTD_approved <- do.call(rbind, lapply(records_list, function(x) {
  x <- x[all_keys]
  x[sapply(x, is.null)] <- NA  # Fill missing keys with NA
  return(as.data.frame(as.list(x), stringsAsFactors = FALSE))
}))

# Ensure columns are in the desired order
TTD_approved <- TTD_approved[, all_keys]



# Read the lines from your text file
lines <- readLines("data/TTD/P2-03-TTD_uniprot_clinical.txt")

# Find the starting index of the data (after the header and abbreviations)
data_start_index <- which(grepl("^TARGETID\t", lines))[2]

# Extract the data lines starting from the data_start_index
data_lines <- lines[data_start_index:length(lines)]

# Initialize an empty list to store each record
record_list <- list()
record_index <- 1
i <- 1

# Loop through the data lines to separate records based on blank lines
while (i <= length(data_lines)) {
  if (data_lines[i] == "") {
    i <- i + 1  # Skip blank lines
  } else {
    record <- c()
    # Collect all lines belonging to a single record
    while (i <= length(data_lines) && data_lines[i] != "") {
      record <- c(record, data_lines[i])
      i <- i + 1
    }
    record_list[[record_index]] <- record
    record_index <- record_index + 1
  }
}

# Process each record to extract key-value pairs
records_list <- lapply(record_list, function(record) {
  key_values <- setNames(
    sapply(record, function(line) unlist(strsplit(line, "\t"))[2]),
    sapply(record, function(line) unlist(strsplit(line, "\t"))[1])
  )
  return(key_values)
})

# Get all unique keys (column names)
all_keys <- unique(unlist(lapply(records_list, names)))

# Convert the list of records into a data frame
TTD_clinical <- do.call(rbind, lapply(records_list, function(x) {
  x <- x[all_keys]
  x[sapply(x, is.null)] <- NA  # Fill missing keys with NA
  return(as.data.frame(as.list(x), stringsAsFactors = FALSE))
}))

# Ensure columns are in the desired order
TTD_clinical <- TTD_clinical[, all_keys]




lines <- readLines("data/TTD/P1-03-TTD_crossmatching.txt")

data_start <- which(grepl("^D\\d", lines))[1]

# Extract only the data section
data_lines <- lines[data_start:length(lines)]

# Initialize variables
current_drug_id <- NA
current_drug_info <- list()
drugs_list <- list()

# Loop over each line in the data section
for (line in data_lines) {
  # Skip empty lines
  if (line == "") next
  
  # Skip lines without tabs (unlikely in the data section, but just in case)
  if (!grepl("\t", line)) next
  
  # Split the line by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Skip lines with fewer than 3 fields
  if (length(fields) < 3) next
  
  # Extract the drug ID, field name, and field values
  drug_id <- fields[1]
  field_name <- fields[2]
  field_values <- fields[3:length(fields)]
  
  # Check if we're starting a new drug
  if (is.na(current_drug_id) || drug_id != current_drug_id) {
    # If we have a previously stored drug, save its info
    if (!is.na(current_drug_id)) {
      drugs_list[[current_drug_id]] <- current_drug_info
    }
    # Initialize new drug information
    current_drug_id <- drug_id
    current_drug_info <- list()
    current_drug_info$TTDDRUID <- drug_id  # Store the TTD Drug ID
  }
  
  # Join multiple field values into a single string
  field_value <- paste(field_values, collapse = "; ")
  
  # If the field already exists, append the new value (in case of multiple lines for the same field)
  if (is.null(current_drug_info[[field_name]])) {
    current_drug_info[[field_name]] <- field_value
  } else {
    current_drug_info[[field_name]] <- paste(current_drug_info[[field_name]], field_value, sep = "; ")
  }
}

# Save the last drug information
if (!is.na(current_drug_id)) {
  drugs_list[[current_drug_id]] <- current_drug_info
}

# Get all unique field names
all_field_names <- unique(unlist(lapply(drugs_list, names)))

# Create a dataframe for drug information
drug_info_list <- lapply(drugs_list, function(x) {
  # Add missing fields with NA
  missing_fields <- setdiff(all_field_names, names(x))
  for (mf in missing_fields) {
    x[[mf]] <- NA
  }
  # Order the fields
  x <- x[all_field_names]
  data.frame(x, stringsAsFactors = FALSE)
})

# Combine all drug information into a single dataframe
drug_crossID_df <- do.call(rbind, drug_info_list)


# get drug types
file_path <- "data/TTD/P1-02-TTD_drug_download.txt"

# Read all lines
lines <- readLines(file_path)

# Identify where data begins (i.e., the first line starting with 'D')
data_start <- which(grepl("^D\\d", lines))[1]

# Extract only the data section
data_lines <- lines[data_start:length(lines)]

# Initialise variables
current_drug_id <- NA
current_drug_info <- list()
drugs_list <- list()

# Loop over each line in the data section
for (line in data_lines) {
  # Skip empty lines
  if (line == "") next
  
  # Skip lines without tabs (unlikely in the data, but just in case)
  if (!grepl("\t", line)) next
  
  # Split by tabs
  fields <- strsplit(line, "\t")[[1]]
  
  # Skip lines with fewer than 3 fields
  if (length(fields) < 3) next
  
  # Extract the Drug ID, field name, and field value(s)
  drug_id <- fields[1]
  field_name <- fields[2]
  field_values <- fields[3:length(fields)]
  
  # If we detect a new Drug ID, store the old one and start a new record
  if (is.na(current_drug_id) || drug_id != current_drug_id) {
    if (!is.na(current_drug_id)) {
      # Save the record for the previous drug
      drugs_list[[current_drug_id]] <- current_drug_info
    }
    # Begin new drug record
    current_drug_id <- drug_id
    current_drug_info <- list()
    current_drug_info$DRUG__ID <- drug_id  # or "TTDDRUID", whichever label you prefer
  }
  
  # Combine field values (in case there are multiple fields to store)
  field_value <- paste(field_values, collapse = "; ")
  
  # If the field already exists, append to it (some entries have multiple lines for the same field)
  if (is.null(current_drug_info[[field_name]])) {
    current_drug_info[[field_name]] <- field_value
  } else {
    current_drug_info[[field_name]] <- paste(current_drug_info[[field_name]], field_value, sep = "; ")
  }
}

# Save the last drug record
if (!is.na(current_drug_id)) {
  drugs_list[[current_drug_id]] <- current_drug_info
}

# Extract all unique field names
all_field_names <- unique(unlist(lapply(drugs_list, names)))

# Create a data frame from the list of drugs
drug_info_list <- lapply(drugs_list, function(x) {
  # Fill missing fields
  missing_fields <- setdiff(all_field_names, names(x))
  for (mf in missing_fields) {
    x[[mf]] <- NA
  }
  # Order columns consistently
  x <- x[all_field_names]
  data.frame(x, stringsAsFactors = FALSE)
})

# Combine all drug entries into a single data frame
drug_rawinfo_df <- do.call(rbind, drug_info_list)
drug_rawinfo_df$DRUG__ID <- rownames(drug_rawinfo_df)
# Inspect your final data frame
# View(drug_rawinfo_df)

# Optionally save it to RData or CSV
# save(drug_rawinfo_df, file = "data_general/TTD_drugs_rawinfo_data.RData")
# write.csv(drug_rawinfo_df, "data_general/TTD_drugs_rawinfo_data.csv", row.names = FALSE)


save(TTD_approved, TTD_clinical, target_info_df, drug_info_df, targetDisease_data, drugDisease_data, drug_rawinfo_df,
     file = "RData/TTD_data.RData")

