addMAR <- function(data, beta = 0.4)
{
  # Add missing values to a dataset under missingness at random.
  #
  # Args:
  #   data, original dataset;
  #   beta, the percentage of missing values.
  #
  # Return: the dataset with added missing values.
  #
  
  #
  if (beta>0.5)
    stop('beta should be smaller or equal to 0.5.')
  #
  n <- dim(data)[1]
  p <- dim(data)[2]
  #
  for (j in 1:round(p/2 - 0.1)){
    temp <- data[,2*j-1]
    ind <- (1:n)[ temp < quantile(temp, runif(1, 0, 2*beta))]
    ind <- sample(ind, 0.99*length(ind))
    data[ind, 2*j] <- NA
  }
  #
  data
}