library(httr)
library(jsonlite)

# Set disease_id variable for breast cancer
disease_id <- "MONDO_0007254"

query_string = "
query KnownDrugsQuery(
  $efoId: String!
  $cursor: String
  $freeTextQuery: String
  $size: Int = 100
) {
  disease(efoId: $efoId) {
    id
    knownDrugs(cursor: $cursor, freeTextQuery: $freeTextQuery, size: $size) {
      count
      cursor
      rows {
        phase
        status
        urls {
          name
          url
        }
        disease {
          id
          name
        }
        drug {
          id
          name
          mechanismsOfAction {
            rows {
              actionType
              targets {
                id
              }
            }
          }
        }
        urls {
          url
          name
        }
        drugType
        mechanismOfAction
        target {
          id
          approvedName
          approvedSymbol
        }
      }
    }
  }
}
"


# Set base URL of GraphQL API endpoint
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"

# Set initial variables object of arguments to be passed to endpoint
variables <- list("efoId" = disease_id, "cursor" = NULL, "size" = 100)

fetch_known_drugs <- function(variables, base_url, query_string) {
  all_known_drugs <- list()
  total_fetched <- 0
  cursor <- NULL
  
  repeat {
    # Update cursor in variables
    variables$cursor <- cursor
    
    # Construct POST request body object with query string and variables
    post_body <- list(query = query_string, variables = variables)
    
    # Perform POST request
    r <- POST(url=base_url, body=post_body, encode='json')
    
    # Extract content and convert to JSON
    content_data <- content(r, as = "text", encoding = "UTF-8")
    json_data <- fromJSON(content_data)
    
    # Debugging information
    if (!is.null(json_data$errors)) {
      print(json_data$errors)
      stop("GraphQL query error")
    }
    
    # Extract known drugs data
    known_drugs_data <- json_data$data$disease$knownDrugs$rows
    
    # Append to all_known_drugs
    all_known_drugs <- c(all_known_drugs, list(known_drugs_data))
    
    # Update progress tracker
    count_fetched <- length(known_drugs_data)
    total_fetched <- total_fetched + count_fetched
    cat("Fetched", total_fetched, "entries...\n")
    
    # Get the next cursor
    new_cursor <- json_data$data$disease$knownDrugs$cursor
    
    # Debugging information for cursor
    cat("Current cursor:", if (is.null(cursor)) "NULL" else cursor, "\n")
    cat("New cursor:", if (is.null(new_cursor)) "NULL" else new_cursor, "\n")
    cat("New cursor type:", typeof(new_cursor), "\n")
    
    # Break the loop if the cursor does not change, is NULL, or unexpected value
    if (is.null(new_cursor) || length(new_cursor) == 0 || (!is.null(cursor) && new_cursor == cursor)) {
      cat("Breaking the loop due to cursor condition.\n")
      break
    }
    
    cursor <- new_cursor
  }
  
  return(all_known_drugs)
}

# Fetch all known drugs
all_known_drugs <- fetch_known_drugs(variables, base_url, query_string)

# Initialize an empty data frame with the correct columns
known_drugs_df <- data.frame(
  phase = character(),
  status = character(),
  disease_id = character(),
  disease_name = character(),
  drug_id = character(),
  drug_name = character(),
  drugType = character(),
  mechanismOfAction = character(),
  target_id = character(),
  target_approvedName = character(),
  target_approvedSymbol = character(),
  stringsAsFactors = FALSE
)

# Populate the data frame
for (x in all_known_drugs) {
  temp_df <- data.frame(
    phase = x$phase,
    status = x$status,
    disease_id = x$disease$id,
    disease_name = x$disease$name,
    drug_id = x$drug$id,
    drug_name = x$drug$name,
    drugType = x$drugType,
    mechanismOfAction = ifelse(is.null(x$mechanismOfAction), NA, x$mechanismOfAction),
    target_id = ifelse(is.null(x$target$id), NA, x$target$id),
    target_approvedName = ifelse(is.null(x$target$approvedName), NA, x$target$approvedName),
    target_approvedSymbol = ifelse(is.null(x$target$approvedSymbol), NA, x$target$approvedSymbol),
    stringsAsFactors = FALSE
  )
  
  known_drugs_df <- rbind(known_drugs_df, temp_df)
  
  rm(x)
}



