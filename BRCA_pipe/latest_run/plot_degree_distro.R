# Sort data frame by degree in descending order
df <- df[order(-df$degree), ]
df$rank <- seq_len(nrow(df))

# Identify ranks of interest mentioned by reviewer
ranks_of_interest <- c(1:5, 9, 10, 30, 50, 100, 1000, nrow(df))

# Load required packages
library(ggplot2)
library(ggrepel)
library(patchwork) # For combining plots easily

# Main plot: full distribution of degree centrality
ggplot(df, aes(x = rank, y = degree_centrality)) +
  geom_line(colour = "steelblue") +
  # Highlight and label points of interest
  geom_point(
    data = subset(df, rank %in% ranks_of_interest), 
    aes(x = rank, y = degree_centrality), 
    colour = "red",
    size = 3
  ) +
  geom_text_repel(
    data = subset(df, rank %in% ranks_of_interest),
    aes(label = paste0("Rank ", rank, ": ", external_gene_name)),
    size = 6, 
    colour = "black",
    hjust = 0,            # Left-align text, ensuring it appears to the right
    nudge_x = 100,
    nudge_y = 2,
    direction = "y",       # Encourages vertical stacking
    segment.size = 0.3,    # Adjust as needed for the connecting line thickness
    segment.color = "grey20", # Make the segment line clearly visible
    box.padding = 0.35,
    point.padding = 0.3,
  ) +
  labs(
    title = NULL,
    x = "Rank",
    y = "Degree Centrality"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15, colour = "black"),
    axis.title = element_text(size = 20),
    axis.text.x = element_text(margin = margin(t=-15)),
    axis.text.y = element_text(margin = margin(r=-15)),
    axis.title.x = element_text(margin = margin(t=10)),
    axis.title.y = element_text(margin = margin(r=10))
  )

# Zoomed-in plot: focus on top 200 nodes to show granularity
ggplot(subset(df, rank <= 200), aes(x = rank, y = degree_centrality)) +
  geom_line(colour = "steelblue") +
  # Highlight points of interest within top 200
  geom_point(
    data = subset(df, rank %in% ranks_of_interest & rank <= 200), 
    aes(x = rank, y = degree_centrality), 
    colour = "red"
  ) +
  geom_text_repel(
    data = subset(df, rank %in% ranks_of_interest & rank <= 200),
    aes(label = paste0("Rank ", rank, ": ", external_gene_name)),
    size = 3, 
    colour = "black",
    hjust = 0,         # Left-align again
    nudge_x = 100, 
    nudge_y = 2,
    direction = "y",    # Encourages vertical spreading
    segment.size = 0.3,
    segment.color = "grey20",
    box.padding = 0.35,
    point.padding = 0.3
  ) +
  labs(
    title = "Zoomed-In View of Top 200 Genes",
    x = "Rank",
    y = "Degree Centrality"
  ) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 20)
  ) 

# Combine the two plots side by side
combined_plot <- p_full + p_zoom

# Print the combined plot
print(combined_plot)
print(p_full)
print(p_zoom)
