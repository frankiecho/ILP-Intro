plotSpatialGrid <- function(mat, label = '', colorScheme = 'viridis') ggplot(melt(mat, c("x", "y"), value.name = "z"), aes(x=x,y=y,fill=z)) + 
  geom_tile() + scale_fill_viridis_c('', option = colorScheme, direction = -1) + ggtitle(label) + theme_void() + theme(legend.position="bottom")

meshgrid <- function(x, y = x) {
  ## meshgrid convenience function from pracma
  if (!is.numeric(x) || !is.numeric(y))
    stop("Arguments 'x' and 'y' must be numeric vectors.")
  
  x <- c(x); y <- c(y)
  n <- length(x)
  m <- length(y)
  
  X <- matrix(rep(x, each = m),  nrow = m, ncol = n)
  Y <- matrix(rep(y, times = n), nrow = m, ncol = n)
  
  return(list(X = X, Y = Y))
}

spatialWeightsMatrix <- function(nx, ny, order = 1, rowStandardise = T) {
  grid <- meshgrid(1:nx,1:ny)
  X    <- as.vector(grid$X)
  Y    <- as.vector(grid$Y)
  P    <- matrix(nrow = nx*ny, ncol = nx*ny)
  W    <- P
  for (i in 1:length(X)) {
    manhattanDist <- abs(X - X[i]) + abs(Y - Y[i])
    P[,i] <- ifelse(manhattanDist <= order, 1, 0)
    P[i,i] <- 0
  }
  if (rowStandardise) {
    for (i in 1:length(X)) {
      W[i,] = P[i,] / sum(P[i,])
    }
    return(W)
  } else {
    return(P)
  }
}

correlatedSpatialGrid <- function(g, nx, ny) {
  W <- spatialWeightsMatrix(nx,ny)
  return(W %*% runif(nx*ny) * g)
}
