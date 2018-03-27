gaussCItestLocal <- function (x, y, S, suffStat) 
{
  # Using Fisher's z-transformation of the partial correlation, test for zero partial correlation. 
  #
  # Args:
  #   see function gaussCItest() in Package 'pcalg';
  #   'suffStat' contains an addtional element 'ESS.Mat', that is, the matrix of effective sample size. 
  #
  # Returns:
  #   the p-value of the current test.
  
  subMat <- suffStat$ESS.Mat[c(x,y,S),c(x,y,S)]
  ne <- mean(subMat[upper.tri(subMat)], na.rm = T)
  z <- zStat(x, y, S, C = suffStat$C, n = ne)
  2 * pnorm(abs(z), lower.tail = FALSE)
}