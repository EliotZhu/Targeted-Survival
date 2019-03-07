check.overlap <- function(x,
                          trt,
                          propensity.func,
                          type  = c("histogram", "density", "both"),
                          bins  = 50L,
                          alpha = ifelse(type == "both", 0.35, 0.5))
{
  require(ggthemes)
  require(ggplot2)
  
  type <- match.arg(type)
  bins <- as.integer(bins[1])
  
  # compute propensity scores
  pi.x <- drop(propensity.func(x = x, trt = trt))
  
  # make sure the resulting propensity scores are in the
  # acceptable range (ie 0-1)
  rng.pi <- range(pi.x)
  
  if (rng.pi[1] <= 0 | rng.pi[2] >= 1) stop("propensity.func() should return values between 0 and 1")
  
  # should be FALSE for treatment/control scenario,
  # TRUE for multiple treatment scenario
  multiplot <- FALSE
  
  dim.pi.x <- dim(pi.x)
  if (!is.null(dim.pi.x))
  {
    if (length(dim.pi.x) == 1)
    {
      pi.x <- as.vector(pi.x)
      prop.scores <- data.frame(Treatment = as.factor(trt), prop.score = pi.x)
      
    } else if (length(dim.pi.x) > 2)
    {
      stop("propensity.func() returns a multidimensional array; it can only return a vector or matrix.")
    }
    
    
    trt.names <- colnames(pi.x)
    
    if (is.null(trt.names))
    {
      if (is.factor(trt))
      {
        # drop any unused levels of trt
        trt         <- droplevels(trt)
        trt.names   <- levels(trt)
      } else
      {
        trt.names   <- sort(unique(trt))
      }
      
    }
    
    prop.scores <- data.frame(Treatment_Received = as.factor(rep(trt, NCOL(pi.x))),
                              Treatment = rep(trt.names, NROW(trt)),
                              prop.score = as.vector(t(pi.x)))
    
    levels(prop.scores$Treatment_Received) <- paste(levels(prop.scores$Treatment_Received), "Group")
    
    multiplot <- TRUE
    
  } else
  {
    prop.scores <- data.frame(Treatment = as.factor(trt), prop.score = pi.x)
  }
  
  
  Treatment <- prop.score <- NULL
  
  if (type == "density")
  {
    pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
      geom_density(alpha = alpha, colour = "grey50") +
      geom_rug(aes(colour = Treatment)) +
      theme(legend.position = "bottom") +
      ggtitle("Densities of propensity scores by treatment group") +
      xlab("Propensity Score")+theme_hc()
  } else if (type == "histogram")
  {
    pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
      geom_histogram(bins = bins, alpha = alpha, position = "identity") +
      geom_rug(aes(colour = Treatment)) +
      theme(legend.position = "bottom") +
      ggtitle("Histograms of propensity scores by treatment group") +
      xlab("Propensity Score")+theme_hc()
  } else
  {
    pl.obj <- ggplot(prop.scores, aes(x = prop.score, fill = Treatment)) +
      geom_histogram(aes(y = ..density..), bins = bins, alpha = alpha, position = "identity") +
      geom_rug(aes(colour = Treatment)) +
      geom_density(alpha = alpha) +
      theme(legend.position = "bottom") +
      ggtitle("Densities and histograms of propensity scores by treatment group") +
      xlab("Propensity Score")+theme_hc()
  }
  
  if (multiplot)
  {
    pl.obj <- pl.obj + facet_grid(Treatment_Received ~ .)
  }
  pl.obj
}