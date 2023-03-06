bal.plot <- function(obj, var.name, ..., which, which.sub = NULL, cluster = NULL, which.cluster = NULL, imp = NULL, which.imp = NULL, which.treat = NULL, which.time = NULL, size.weight = FALSE, mirror = FALSE, type = c("density", "histogram"), colors = NULL) {
  
  tryCatch(identity(obj), error = function(e) stop(conditionMessage(e), call. = FALSE))
  
  #Replace .all and .none with NULL and NA respectively
  .call <- match.call(expand.dots = TRUE)
  if (any(sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all") || identical(as.character(.call[[x]]), ".none")))) {
    .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".all"))] <- expression(NULL)
    .call[sapply(seq_along(.call), function(x) identical(as.character(.call[[x]]), ".none"))] <- expression(NA)
    return(eval(.call))
  }
  
  args <- list(...)
  
  obj <- is.designmatch(obj)
  obj <- is.time.list(obj)
  
  X <- x2base(obj, ..., cluster = cluster, imp = imp, s.d.denom = "treated") #s.d.denom to avoid x2base warning
  
  if (is_not_null(X$covs.list)) {
    if (missing(var.name)) {
      var.name <- names(X$covs.list[[1]])[1]
      message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
    }
    var.list <- vector("list", length(X$covs.list))
    appears.in.time <- rep(TRUE, length(X$covs.list))
    for (i in seq_along(X$covs.list)) {
      if (var.name %in% names(X$covs.list[[i]])) var.list[[i]] <- X$covs.list[[i]][[var.name]]
      else if (!is.null(X$addl.list) && var.name %in% names(X$addl.list[[i]])) var.list[[i]] <- X$addl[[var.name]]
      else if (!is.null(X$distance.list) && var.name %in% names(X$distance.list[[i]])) var.list[[i]] <- X$distance.list[[i]][[var.name]]
      else appears.in.time[i] <- FALSE
    }
    if (all(sapply(var.list, is_null))) stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
    X$var <- unlist(var.list[appears.in.time])
    X$time <- rep(seq_along(X$covs.list)[appears.in.time], each = NROW(X$covs.list[[1]]))
    X$treat <- unlist(X$treat.list[appears.in.time])
    if (is_not_null(names(X$treat.list))) treat.names <- names(X$treat.list)
    else treat.names <- seq_along(X$treat.list)
    if (is_not_null(X$weights)) X$weights <- do.call("rbind", replicate(sum(appears.in.time), list(X$weights)))
  }
  else {
    if (missing(var.name)) {
      var.name <- names(X$covs)[1]
      message(paste0("No var.name was provided. Dispalying balance for ", var.name, "."))
    }
    if (var.name %in% names(X$covs)) X$var <- X$covs[[var.name]]
    else if (!is.null(X$addl) && var.name %in% names(X$addl)) X$var <- X$addl[[var.name]]
    else if (!is.null(args$addl) && var.name %in% names(args$addl)) X$var <- args$addl[[var.name]]
    else if (!is.null(X$distance) && var.name %in% names(X$distance)) X$var <- X$distance[[var.name]]
    else if (!is.null(args$distance) && var.name %in% names(args$distance)) X$var <- args$distance[[var.name]]
    else stop(paste0("\"", var.name, "\" is not the name of a variable in any available data set input."), call. = FALSE)
  }
  
  #Density arguments supplied through ...
  if (is_not_null(args$bw)) bw <- args$bw else bw <- "nrd0"
  if (is_not_null(args$adjust)) adjust <- args$adjust else adjust <- 1
  if (is_not_null(args$kernel)) kernel <- args$kernel else kernel <- "gaussian"
  if (is_not_null(args$n)) n <- args$n else n <- 512
  
  if (missing(which)) {
    if (is_not_null(args$un)) {
      message("Note: \'un\' is deprecated; please use \'which\' for the same and added functionality.")
      if (args$un) which <- "unadjusted"
      else which <- "adjusted"
    }
    else {
      which <- "adjusted"
    }
  }
  else {
    if (is_null(X$weights) && is_null(X$subclass)) which <- "unadjusted"
    else {
      which <- match_arg(tolower(which), c("adjusted", "unadjusted", "both"))
    }
  }
  
  title <- paste0("Distributional Balance for \"", var.name, "\"")
  subtitle <- NULL
  
  facet <- NULL
  is.categorical.var <- is_binary(X$var) || is.factor(X$var) || is.character(X$var)
  
  if (is_not_null(X$subclass)) {
    if (which %in% c("adjusted", "both")) {
      if (is_not_null(X$cluster)) stop("Subclasses are not supported with clusters.", call. = FALSE)
      if (is_not_null(X$imp)) stop("Subclasses are not supported with multiple imputations.", call. = FALSE)
      if (is_null(which.sub)) { 
        which.sub <- levels(X$subclass)
      }
      if (is.numeric(which.sub) && any(which.sub %in% levels(X$subclass))) {
        if (any(!which.sub %in% levels(X$subclass))) {
          w.l <- word_list(which.sub[!which.sub %in% levels(X$subclass)])
          warning(paste(w.l, ifelse(attr(w.l, "plural"), "do", "does"), "not correspond to any subclass in the object and will be ignored."), call. = FALSE)
          which.sub <- which.sub[which.sub %in% levels(X$subclass)]
        }
        in.sub <- !is.na(X$subclass) & X$subclass %in% which.sub
        D <- setNames(as.data.frame(matrix(0, nrow = sum(in.sub), ncol = 4)),
                      c("weights", "treat", "var", "subclass"))
        D$weights <- rep(1, NROW(D))
        D$treat <- X$treat[in.sub]
        D$var <- X$var[in.sub]
        D$subclass <- paste("Subclass", X$subclass[in.sub])
        #title <- paste0(title, "\nin subclass ", which.sub)
        subtitle <- NULL
      }
      else stop("The argument to which.sub must be a number vector corresponding to the subclass for which distributions are to be displayed.", call. = FALSE)
      facet <- "subclass"
    }
    #D$weights <- rep(1, length(treat))
    if (which == "both") {
      D2 <- setNames(as.data.frame(matrix(0, nrow = length(X$treat), ncol = 4)),
                     c("weights", "treat", "var", "subclass"))
      D2$weights <- rep(1, NROW(D2))
      D2$treat <- X$treat
      D2$var <- X$var
      D2$subclass <- rep("Unadjusted Sample", length(X$treat))
      D <- rbind(D2, D, stringsAsFactors = TRUE)
      D$subclass <- relevel(factor(D$subclass), "Unadjusted Sample")
    }
    
  }
  else if (is_null(X$subclass) && is_not_null(which.sub)) {
    warning("which.sub was specified but no subclasses were supplied. Ignoring which.sub.", call. = FALSE)
  }
  else if (which == "unadjusted" && is_not_null(which.sub)) {
    warning("which.sub was specified but the unadjusted sample was requested. Ignoring which.sub.", call. = FALSE)
  }
  
  if ("subclass" %nin% facet) {
    
    facet <- "facet.which"
    
    if (is_null(X$weights) || which == "unadjusted") {
      X$weights <- data.frame(rep(1, length(X$treat)))
      names(X$weights) <- "Unadjusted Sample"
      #nweights <- 1
    }
    else {
      if (ncol(X$weights) == 1) {
        names(X$weights) <- "Adjusted Sample"
      }
      if (which == "both") {
        X$weights <- setNames(cbind(rep(1, length(X$treat)), X$weights),
                              c("Unadjusted Sample", names(X$weights)))
      }
    }
    
    nweights <- ncol(X$weights)
    weight.names <- names(X$weights)
    
    #NULL: all; NA: none
    in.imp <- rep(TRUE, length(X$var))
    if (is_not_null(X$imp)) {
      if (is_null(which.imp) || all(is.na(which.imp))) {
        in.imp <- !is.na(X$imp)
      }
      else {
        if (is.numeric(which.imp)) {
          if (all(which.imp %in% seq_len(nlevels(X$imp)))) {
            in.imp <- !is.na(X$imp) & sapply(X$imp, function(x) !is.na(match(x, levels(X$imp)[which.imp])))
          }
          else {
            stop(paste0("The following inputs to which.imp do not correspond to given imputations:\n\t", word_list(which.imp[!which.imp %in% seq_len(nlevels(X$imp))])), call. = FALSE)
          }
        }
        else stop("The argument to which.imp must be the indices corresponding to the imputations for which distributions are to be displayed.", call. = FALSE)
      }
      facet <- c("imp", facet)
    }
    else if (is_not_null(which.imp)) {
      warning("which.imp was specified but no imp values were supplied. Ignoring which.imp.", call. = FALSE)
    }
    
    in.cluster <- rep(TRUE, length(X$var))
    if (is_not_null(X$cluster)) {
      if (is_null(which.cluster)|| all(is.na(which.cluster))) {
        in.cluster <- !is.na(X$cluster)
      }
      else {
        if (is.numeric(which.cluster)) {
          if (all(which.cluster %in% seq_len(nlevels(X$cluster)))) {
            in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, levels(X$cluster)[which.cluster])))
          }
          else {
            stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[!which.cluster %in% seq_len(nlevels(X$cluster))])), call. = FALSE)
          }
        }
        else if (is.character(which.cluster)) {
          if (all(!is.na(match(which.cluster, levels(X$cluster))))) {
            in.cluster <- !is.na(X$cluster) & sapply(X$cluster, function(x) !is.na(match(x, which.cluster)))
          }
          else {
            stop(paste0("The following inputs to which.cluster do not correspond to given clusters:\n\t", word_list(which.cluster[is.na(match(which.cluster, levels(X$cluster)))])), call. = FALSE)
          }
        }
        else stop("The argument to which.cluster must be the names or indices corresponding to the clusters for which distributions are to be displayed.", call. = FALSE)
      }
      facet <- c("cluster", facet)
    }
    else if (is_not_null(which.cluster)) {
      warning("which.cluster was specified but no cluster values were supplied. Ignoring which.cluster.", call. = FALSE)
    }
    
    in.time <- rep(TRUE, length(X$var))
    if (is_not_null(X$time)) {
      if (is_null(which.time) || all(is.na(which.time))) {
        in.time <- !is.na(X$time)
      }
      else {
        if (is.numeric(which.time)) {
          if (all(which.time %in% seq_along(X$covs.list))) {
            if (all(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
              #nothing; which.time is good
            }
            else if (any(which.time %in% seq_along(X$covs.list)[appears.in.time])) {
              warning(paste0(var.name, " does not appear in time period ", word_list(which.time[!which.time %in% seq_along(X$covs.list)[appears.in.time]], "or"), "."), call. = FALSE)
              which.time <- which.time[which.time %in% seq_along(X$covs.list)[appears.in.time]]
            }
            else {
              stop(paste0(var.name, " does not appear in time period ", word_list(which.time, "or"), "."), call. = FALSE)
            }
            in.time <- !is.na(X$time) & X$time %in% which.time
          }
          else {
            stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% seq_along(X$covs.list)])), call. = FALSE)
          }
        }
        else if (is.character(which.time)) {
          if (all(which.time %in% treat.names)) {
            if (all(which.time %in% treat.names[appears.in.time])) {
              #nothing; which.time is good
            }
            else if (any(which.time %in% treat.names[appears.in.time])) {
              time.periods <- word_list(which.time[!which.time %in% treat.names[appears.in.time]], "and")
              warning(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                             "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                             time.periods, "."), call. = FALSE)
              which.time <- which.time[which.time %in% treat.names[appears.in.time]]
            }
            else {
              time.periods <- word_list(which.time, "and")
              stop(paste0(var.name, " does not appear in the time period", ifelse(attr(time.periods, "plural"), "s ", " "),
                          "corresponding to treatment", ifelse(attr(time.periods, "plural"), "s ", " "),
                          time.periods, "."), call. = FALSE)
            }
            in.time <- !is.na(X$time) & treat.names[X$time] %in% which.time
            
          }
          else {
            stop(paste0("The following inputs to which.time do not correspond to given time periods:\n\t", word_list(which.time[!which.time %in% treat.names])), call. = FALSE)
          }
        }
        else stop("The argument to which.time must be the names or indices corresponding to the time periods for which distributions are to be displayed.", call. = FALSE)
      }
      facet <- c("time", facet)
    }
    else if (is_not_null(which.time)) {
      warning("which.time was specified but a point treatment was supplied. Ignoring which.time.", call. = FALSE)
    }
    
    nobs <- sum(in.imp & in.cluster & in.time)
    if (nobs == 0) stop("No observations to display.", call. = FALSE)
    
    Q <- setNames(vector("list", nweights), weight.names)
    for (i in weight.names) {
      Q[[i]] <- setNames(as.data.frame(matrix(0, ncol = 7, nrow = nobs)),
                         c("imp", "cluster", "time", "treat", "var", "weights", "facet.which"))
      Q[[i]]$imp <- Q[[i]]$cluster <- Q[[i]]$time <- character(nobs)
      Q[[i]]$treat <- X$treat[in.imp & in.cluster & in.time]
      Q[[i]]$var <- X$var[in.imp & in.cluster & in.time]
      Q[[i]]$weights <-  X$weights[in.imp & in.cluster & in.time, i]
      Q[[i]]$facet.which <- rep(i, nobs)
      
      if ("imp" %in% facet) Q[[i]]$imp <- paste("Imputation", X$imp[in.imp & in.cluster & in.time])
      if ("cluster" %in% facet) Q[[i]]$cluster <- factor(X$cluster[in.imp & in.cluster & in.time])
      if ("time" %in% facet) Q[[i]]$time <- paste("Time", X$time[in.imp & in.cluster & in.time])
    }
    D <- do.call("rbind", Q)
    D$facet.which <- factor(D$facet.which, levels = c(weight.names[weight.names == "Unadjusted Sample"],
                                                      weight.names[weight.names != "Unadjusted Sample"]))
    
  }
  
  if ("facet.which" %in% facet) {
    if (nlevels(D$facet.which) == 1) {
      subtitle <- levels(D$facet.which)[1]
      facet <- facet[!facet %in% "facet.which"]
    }
  }
  
  if (!is_binary(D$treat) && is.numeric(D$treat)) { #Continuous treatments
    if ("subclass" %in% facet) {
      if (is.categorical.var) {
        weights <- with(D, ave(weights, subclass, var, FUN = function(x) x/sum(x)))
      }
      else {
        weights <- with(D, ave(weights, subclass, FUN = function(x) x/sum(x)))
      }
      d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
    }
    else {
      if (is.categorical.var) {
        weights <- with(D, ave(weights, cluster, imp, time, facet.which, var, FUN = function(x) x/sum(x)))
      }
      else {
        weights <- with(D, ave(weights, cluster, imp, time, facet.which, FUN = function(x) x/sum(x)))
      }
      d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
      
    }
    
    if (is.categorical.var) { #Categorical vars
      d$var <- factor(d$var)
      cat.sizes <- tapply(rep(1, NROW(d)), d$var, sum)
      smallest.cat <- names(cat.sizes)[which.min(cat.sizes)]
      if (is.character(bw)) {
        if (is.function(get0(paste0("bw.", bw)))) {
          bw <- get0(paste0("bw.", bw))(d$treat[d$var == smallest.cat])
        }
        else {
          stop(paste(bw, "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
        }
      }
      
      #Color
      ntypes <- length(cat.sizes)
      if (is_not_null(args$colours)) colors <- args$colours
      
      if (is_null(colors)) {
        colors <- gg_color_hue(ntypes)
      }
      else {
        if (length(colors) > ntypes) {
          colors <- colors[seq_len(ntypes)]
          warning(paste("Only using first", ntypes, "values in colors."), call. = FALSE)
        }
        else if (length(colors) < ntypes) {
          warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
          colors <- gg_color_hue(ntypes)
        }
        
        if (!all(sapply(colors, isColor))) {
          warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
          colors <- gg_color_hue(ntypes)
        }
        
      }
      
      bp <- ggplot(d, mapping = aes(treat, fill = var, weight = weights)) + 
        geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n, trim = TRUE) + 
        labs(fill = var.name, y = "Density", x = "Treat", title = title, subtitle = subtitle) +
        scale_fill_manual(values = colors) + geom_hline(yintercept = 0)
    }
    else { #Continuous vars
      bp <- ggplot(d, mapping = aes(x = var, y = treat, weight = weights))
      if (which == "unadjusted" || !isTRUE(size.weight)) bp <- bp + geom_point(alpha = .9)
      else bp <- bp + geom_point(aes(size = weights), alpha = .9)
      bp <- bp + geom_smooth(method = "loess", se = FALSE, alpha = .1) + 
        geom_smooth(method = "lm", se = FALSE, linetype = 2, alpha = .4) + 
        geom_hline(yintercept = w.m(d$treat, d$weights), linetype = 1, alpha = .9) + 
        labs(y = "Treat", x = var.name, title = title, subtitle = subtitle)
    }
  }
  else { #Categorical treatments (multinomial supported)
    D$treat <- factor(D$treat)
    
    if (is_null(which.treat)) 
      which.treat <- character(0)
    else if (is.numeric(which.treat)) {
      which.treat <- levels(D$treat)[seq_along(levels(D$treat)) %in% which.treat]
      if (is_null(which.treat)) {
        warning("No numbers in which.treat correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
        which.treat <- character(0)
      }
    }
    else if (is.character(which.treat)) {
      which.treat <- levels(D$treat)[levels(D$treat) %in% which.treat]
      if (is_null(which.treat)) {
        warning("No names in which.treat correspond to treatment values. All treatment groups will be displayed.", call. = FALSE)
        which.treat <- character(0)
      }
    }
    else if (is.na(which.treat)) {
      which.treat <- character(0)
    }
    else {
      warning("The argument to which.treat must be NA, NULL, or a vector of treatment names or indices. All treatment groups will be displayed.", call. = FALSE)
      which.treat <- character(0)
    }
    if (is_not_null(which.treat) && all(!is.na(which.treat))) D <- D[D$treat %in% which.treat,]
    
    for (i in names(D)[sapply(D, is.factor)]) D[[i]] <- factor(D[[i]])
    
    if (is_not_null(facet)) {
      if ("subclass" %in% facet) {
        weights <- with(D, ave(weights, treat, subclass, FUN = function(x) x/sum(x)))
        d <- data.frame(weights = weights, treat = D$treat, var = D$var, subclass = D$subclass)
      }
      else {
        weights <- with(D, ave(weights, treat, cluster, imp, time, facet.which, FUN = function(x) x/sum(x)))
        d <- data.frame(weights = weights, treat = D$treat, var = D$var, cluster = D$cluster, imp = D$imp, time = D$time, facet.which = D$facet.which)
      }
    }
    else {
      weights <- with(D, ave(weights, treat, FUN = function(x) x/sum(x)))
      d <- data.frame(weights = weights, treat = D$treat, var = D$var)
    }
    
    #Color
    ntypes <- nlevels.treat <- nlevels(d$treat)
    if (is_not_null(args$colours)) colors <- args$colours
    
    if (is_null(colors)) {
      colors <- gg_color_hue(ntypes)
    }
    else {
      if (length(colors) > ntypes) {
        colors <- colors[seq_len(ntypes)]
        warning(paste("Only using first", ntypes, "values in colors."), call. = FALSE)
      }
      else if (length(colors) < ntypes) {
        warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
        colors <- gg_color_hue(ntypes)
      }
      
      if (!all(sapply(colors, isColor))) {
        warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
        colors <- gg_color_hue(ntypes)
      }
    }
    names(colors) <- levels(d$treat)
    
    if (is_binary(d$var) || is.factor(d$var) || is.character(d$var)) { #Categorical vars
      d$var <- factor(d$var)
      bp <- ggplot(d, mapping = aes(x = var, fill = treat, weight = weights)) + 
        geom_bar(position = "dodge", alpha = .8, color = "black") + 
        labs(x = var.name, y = "Proportion", fill = "Treat", title = title, subtitle = subtitle) + 
        scale_x_discrete(drop=FALSE) + scale_fill_manual(drop=FALSE, values = colors) +
        geom_hline(yintercept = 0)
    }
    else { #Continuous vars
      t.sizes <- tapply(rep(1, NROW(d)), d$treat, sum)
      smallest.t <- names(t.sizes)[which.min(t.sizes)]
      if (is.character(bw)) {
        if (is.function(get0(paste0("bw.", bw)))) {
          bw <- get0(paste0("bw.", bw))(d$var[d$treat == smallest.t])
        }
        else {
          stop(paste(bw, "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
        }
      }
      
      # bp <- ggplot(d, mapping = aes(var, fill = treat, weight = weights)) +
      #     geom_density(alpha = .4, bw = bw, adjust = adjust, kernel = kernel, n = n, trim = TRUE) +
      #     labs(x = var.name, y = "Density", fill = "Treat", title = title, subtitle = subtitle)
      
      #colors <- setNames(gg_color_hue(nlevels.treat), levels(d$treat))
      
      if (nlevels.treat <= 2 && mirror == TRUE) {
        posneg <- c(1, -1)
        alpha <- .8
      }
      else {
        posneg <- rep(1, nlevels.treat)
        alpha <- .4
      }
      type <- match_arg(type)
      if (type == "histogram" && nlevels.treat <= 2) {
        if (is_null(args$bins) || !is.numeric(args$bins)) args$bins <- 12
        geom_fun <- function(t) geom_histogram(data = d[d$treat == levels(d$treat)[t],],
                                               mapping = aes(x = var, y = posneg[t]*stat(count), weight = weights,
                                                             fill = names(colors)[t]),
                                               alpha = alpha, bins = args$bins, color = "black")
        ylab <- "Proportion"
        legend.alpha <- alpha
      }
      else {
        geom_fun <- function(t) geom_density(data = d[d$treat == levels(d$treat)[t],],
                                             mapping = aes(x = var, y = posneg[t]*stat(density), weight = weights,
                                                           fill = names(colors)[t]),
                                             alpha = alpha, bw = bw, adjust = adjust,
                                             kernel = kernel, n = n, trim = TRUE)
        ylab <- "Density"
        legend.alpha <- alpha/nlevels.treat
      }
      
      bp <- Reduce("+", c(list(ggplot(),
                               labs(x = var.name, y = ylab, fill = "Treat", title = title, subtitle = subtitle)),
                          lapply(seq_len(nlevels.treat), geom_fun))) +
        scale_fill_manual(values = colors, guide = guide_legend(override.aes = list(alpha = legend.alpha, size = .25))) +
        geom_hline(yintercept = 0)
    }
  }
  
  if (is_not_null(facet)) {
    if (length(facet) >= 2) {
      if ("facet.which" %in% facet) {
        facet.formula <- f.build("facet.which", facet[!facet %in% "facet.which"])
      }
      else if ("imp" %in% facet) {
        facet.formula <- f.build("imp", facet[!facet %in% "imp"])
      }
      else {
        facets <- data.frame(facet = facet, length = sapply(facet, function(x) nlevels(d[[x]])),
                             stringsAsFactors = FALSE)
        facets <- facets[with(facets, order(length, facet, decreasing = c(FALSE, TRUE))), ]
        facet.formula <- formula(facets)
      }
    }
    else facet.formula <- f.build(".", facet)
    bp <- bp + facet_grid(facet.formula, drop = FALSE, scales = ifelse("subclass" %in% facet, "free_x", "fixed"))
  }
  return(bp)
}

#Balance functions for use in programming

col_w_mean <- function(mat, weights = NULL, s.weights = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  possibly.supplied <- c("mat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  return(col.w.m(mat[subset, , drop = FALSE], w = weights[subset], na.rm = na.rm))
}
col_w_sd <- function(mat, weights = NULL, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  
  possibly.supplied <- c("mat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  return(sqrt(col.w.v(mat[subset, , drop = FALSE], w = weights[subset], 
                      bin.vars = bin.vars, na.rm = na.rm)))
}

col_w_smd <- function(mat, treat, weights = NULL, std = TRUE, s.d.denom = "pooled", abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
  allowable.s.d.denoms <- c("pooled", "all", "treated", "control")
  unique.treats <- as.character(unique(treat, nmax = 2))
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
  if (!is.atomic(std) || any(is.na(as.logical(std))) ||
      length(std) %nin% c(1L, NCOL(mat))) {
    stop("std must be a logical vector with length equal to 1 or the number of columns of mat.")
  }
  
  
  if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
    stop("abs must be a logical of length 1.")
  }
  
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  
  possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  if (abs || !(length(s.d.denom) == 1L && as.character(s.d.denom) %in% c("treated", "control"))) tval1 <- treat[1]
  else tval1 <- treat[binarize(treat)==1][1]
  
  if (length(std) == 1L) std <- rep(std, NCOL(mat))
  
  m1 <- col.w.m(mat[treat==tval1 & subset, , drop = FALSE], weights[treat==tval1 & subset], na.rm = na.rm)
  m0 <- col.w.m(mat[treat!=tval1 & subset, , drop = FALSE], weights[treat!=tval1 & subset], na.rm = na.rm)
  diffs <- m1 - m0
  zeros <- check_if_zero(diffs)
  
  if (any(to.sd <- std & !zeros)) {
    if (is.atomic(s.d.denom) && length(s.d.denom) == 1L) {
      s.d.denom <- match_arg(as.character(s.d.denom), c(allowable.s.d.denoms, unique.treats))
      
      if (s.d.denom %in% unique.treats) s.d.denom <- sqrt(col.w.v(mat[treat == s.d.denom, to.sd, drop = FALSE], 
                                                                  w = s.weights[treat == s.d.denom],
                                                                  bin.vars = bin.vars[to.sd], na.rm = na.rm))
      else if (s.d.denom == "pooled") s.d.denom <- sqrt(.5*(col.w.v(mat[treat == tval1, to.sd, drop = FALSE], 
                                                                    w = s.weights[treat == tval1],
                                                                    bin.vars = bin.vars[to.sd], na.rm = na.rm) +
                                                              col.w.v(mat[treat != tval1, to.sd, drop = FALSE], 
                                                                      w = s.weights[treat != tval1],
                                                                      bin.vars = bin.vars[to.sd], na.rm = na.rm)))
      else if (s.d.denom == "all") s.d.denom <- sqrt(col.w.v(mat[, to.sd, drop = FALSE], 
                                                             w = s.weights,
                                                             bin.vars = bin.vars[to.sd], na.rm = na.rm))
      else if (s.d.denom == "treated") s.d.denom <- sqrt(col.w.v(mat[treat == tval1, to.sd, drop = FALSE], 
                                                                 w = s.weights[treat == tval1],
                                                                 bin.vars = bin.vars[to.sd], na.rm = na.rm))
      else if (s.d.denom == "control") s.d.denom <- sqrt(col.w.v(mat[treat != tval1, to.sd, drop = FALSE], 
                                                                 w = s.weights[treat != tval1],
                                                                 bin.vars = bin.vars[to.sd], na.rm = na.rm))
    }
    else {
      if (!is.numeric(s.d.denom) || length(s.d.denom) %nin% c(sum(to.sd), length(diffs))) {
        stop(paste0("s.d.denom must be one of ", word_list(c(allowable.s.d.denoms, unique.treats), and.or = "or", quotes = TRUE), 
                    " or a numeric vector of with length equal to the number of columns of mat"))
      }
      else if (length(s.d.denom) == length(diffs)) s.d.denom <- s.d.denom[to.sd]
    }
    
    diffs[to.sd] <- diffs[to.sd]/s.d.denom
  }
  
  if (abs) diffs <- abs(diffs)
  
  return(setNames(diffs, colnames(mat)))
  
}
col_w_vr <- function(mat, treat, weights = NULL, abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
  
  if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
    stop("abs must be a logical of length 1.")
  }
  
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  
  possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  if (abs) tval1 <- treat[1]
  else tval1 <- treat[binarize(treat)==1][1]
  
  v1 <- col.w.v(mat[treat==tval1, , drop = FALSE], weights[treat==tval1], bin.vars = bin.vars, na.rm = na.rm)
  v0 <- col.w.v(mat[treat!=tval1, , drop = FALSE], weights[treat!=tval1], bin.vars = bin.vars, na.rm = na.rm)
  
  v.ratios = v1/v0
  
  if (abs) v.ratios <- abs_(v.ratios, ratio = TRUE)
  
  return(setNames(v.ratios, colnames(mat)))
  
}
col_w_ks <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  if (!(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
  
  possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  tval1 <- treat[1]
  ks <- rep(NA_real_, NCOL(mat))
  
  if (any(!bin.vars)) {
    weights_ <- weights
    weights_[treat == tval1] <-  weights[treat == tval1]/sum(weights[treat == tval1])
    weights_[treat != tval1] <- -weights[treat != tval1]/sum(weights[treat != tval1])
    ks[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(x) {
      if (na.rm) x <- x[!is.na(x)]
      if (!na.rm && anyNA(x)) return(NA_real_)
      else {
        ordered.index <- order(x)
        cumv <- abs(cumsum(weights_[ordered.index]))[diff(x[ordered.index]) != 0]
        return(if (is_null(cumv)) 0 else max(cumv))
      }
    })
  }
  if (any(bin.vars)) {
    ks[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1], na.rm = na.rm) - 
                          col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1], na.rm = na.rm))
  }
  return(setNames(ks, colnames(mat)))
  
}
col_w_ovl <- function(mat, treat, weights = NULL, s.weights = NULL, bin.vars = NULL, integrate = FALSE, subset = NULL, na.rm = TRUE, ...) {
  
  A <- list(...)
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  if (!(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
  
  possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  weights <- weights * s.weights
  
  weights <- weights[subset]
  treat <- treat[subset]
  mat <- mat[subset, , drop = FALSE]
  
  tval1 <- treat[1]
  
  t.sizes <- setNames(vapply(unique(treat, nmax = 2), function(x) sum(treat == x), numeric(1L)),
                      unique(treat, nmax = 2))
  smallest.t <- names(t.sizes)[which.min(t.sizes)]
  ovl <- setNames(numeric(ncol(mat)), colnames(mat))
  if (any(!bin.vars)) {
    if (is_null(A[["bw"]])) A[["bw"]] <- "nrd" else A[["bw"]] <- A[["bw"]]
    A[names(A) %nin% names(formals(density.default))] <- NULL
    
    ovl[!bin.vars] <- apply(mat[, !bin.vars, drop = FALSE], 2, function(cov) {
      if (na.rm) cov <- cov[!is.na(cov)]
      if (!na.rm && anyNA(cov)) return(NA_real_)
      else {
        if (is.function(get0(paste0("bw.", A[["bw"]])))) {
          A[["bw"]] <- get0(paste0("bw.", A[["bw"]]))(cov[treat == smallest.t])
        }
        else {
          stop(paste(A[["bw"]], "is not an acceptable entry to bw. See ?stats::density for allowable options."), call. = FALSE)
        }
        
        f1 <- approxfun(do.call(density.default, c(list(cov[treat==tval1], 
                                                        weights = weights[treat==tval1]/sum(weights[treat==tval1])), A)))
        f0 <- approxfun(do.call(density.default, c(list(cov[treat!=tval1], 
                                                        weights = weights[treat!=tval1]/sum(weights[treat!=tval1])), A)))
        fn <- function(x) {
          f1_ <- f1(x)
          f1_[is.na(f1_)] <- 0
          f0_ <- f0(x)
          f0_[is.na(f0_)] <- 0
          pmin(f1_, f0_)
        }
        min.c <- min(cov)
        max.c <- max(cov)
        range <- max.c - min.c
        # min.c.ext <- min.c - .01 * range
        # max.c.ext <- max.c + .01 * range
        if (isTRUE(integrate)) {
          
          s <- try(integrate(fn, lower = min.c,
                             upper = max.c)$value,
                   silent = TRUE)
        }
        else {
          seg <- seq(min.c, max.c, length = 1001)
          mids <- .5 * (seg[2:length(seg)] + seg[1:(length(seg)-1)])
          s <- sum(fn(mids))*(seg[2]-seg[1])
        }
        
        if (inherits(s, "try-error"))  return(NA_real_)
        else return(1 - s) #Reverse: measure imbalance
      }
    })
  }
  if (any(bin.vars)) {
    ovl[bin.vars] <- abs(col.w.m(mat[treat == tval1, bin.vars, drop = FALSE], weights[treat == tval1]) - 
                           col.w.m(mat[treat != tval1, bin.vars, drop = FALSE], weights[treat != tval1]))
  }
  
  return(ovl)
  
}
col_w_corr <- function(mat, treat, weights = NULL, type = "pearson", abs = FALSE, s.weights = NULL, bin.vars = NULL, subset = NULL, na.rm = TRUE, ...) {
  
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        A <- list(...)
        A <- A[names(A) %in% names(formals(splitfactor)) & 
                 names(A) %nin% c("data", "var.name", "drop.first",
                                  "drop.level")]
        mat <- do.call(splitfactor, c(list(mat, drop.first ="if2"),
                                      A))
      }
      else mat <- as.matrix(mat)
    }
    else if (is.numeric(mat)) mat <- matrix(mat, ncol = 1)
    else stop("mat must be a data.frame or numeric matrix.")
  }
  else if (!is.numeric(mat)) stop("mat must be a data.frame or numeric matrix.")
  
  if (missing(treat) || !is.atomic(treat) || !is.numeric(treat)) stop("treat must be a numeric variable.")
  if (!(is.atomic(abs) && length(abs) == 1L && !is.na(as.logical(abs)))) {
    stop("abs must be a logical of length 1.")
  }
  if (is_null(bin.vars)) bin.vars <- apply(mat, 2, is_binary)
  else if (!is.atomic(bin.vars) || any(is.na(as.logical(bin.vars))) ||
           length(bin.vars) != NCOL(mat)) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.")
  }
  
  possibly.supplied <- c("mat", "treat", "weights", "s.weights", "subset")
  lengths <- setNames(vapply(mget(possibly.supplied), len, integer(1L)),
                      possibly.supplied)
  supplied <- lengths > 0
  if (!all_the_same(lengths[supplied])) {
    stop(paste(word_list(possibly.supplied[supplied]), "must have the same number of units."))
  }
  
  if (lengths["weights"] == 0) weights <- rep(1, NROW(mat))
  if (lengths["s.weights"] == 0) s.weights <- rep(1, NROW(mat))
  
  if (lengths["subset"] == 0) subset <- rep(TRUE, NROW(mat))
  else if (anyNA(as.logical(subset))) stop("subset must be a logical vector.")
  
  type <- match_arg(tolower(type), c("pearson", "spearman"))
  if (type == "spearman") {
    for (i in 1:ncol(mat)) mat[,i] <- rank(mat[,i], na.last = "keep")
    treat <- rank(treat, na.last = "keep")
  }
  
  weights <- weights * s.weights
  
  if (all(subset)) {
    corrs <- col.w.r(mat, treat, w = weights, s.weights = s.weights, bin.vars = bin.vars,
                     na.rm = na.rm)
  }
  else {
    cov <- col.w.cov(mat[subset, , drop = FALSE], y = treat[subset], w = weights[subset], na.rm = na.rm)
    den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) * 
      sqrt(col.w.v(treat, w = s.weights, na.rm = na.rm))
  }
  
  if (abs) corrs <- abs(corrs)
  
  return(setNames(corrs, colnames(mat)))
  
}

#Scalar balance functions 
bal.sum <- function(mat, treat, weights = NULL, type, s.weights = NULL, check = TRUE, ...) {
  uni.type <- c("smd", "ks", "ovl")
  agg.funs <- c("max", "mean", "rms")
  sample.type <- c("mahalanobis", "gwd", "cstat", "wr2", "design.effect")
  uni.type.expanded <- expand.grid_string(uni.type, agg.funs, collapse = ".")
  shortcuts <- c("all", "rec")
  allowable.type <- c(uni.type.expanded, sample.type, shortcuts)
  if (missing(type)) stop("type must be specified.", call. = FALSE)
  else type <- match_arg(type, allowable.type, several.ok = TRUE)
  
  if ("all" %in% type) type <- unique(c(type[type != "all"], uni.type.expanded, sample.type))
  if ("rec" %in% type) type <- unique(c(type[type != "rec"], 
                                        "smd.mean", "smd.rms",
                                        "ks.mean", "ks.rms",
                                        "mahalanobis", "gwd", "wr2"))
  
  A <- list(...)
  
  if (check) {
    bad.mat <- FALSE
    if (missing(mat)) bad.mat <- TRUE
    else {
      if (is.data.frame(mat)) {
        if (any(vapply(mat, function(x) is_(x, c("character", "factor")), logical(1L)))) 
          mat <- splitfactor(mat)
        mat <- as.matrix.data.frame(mat)
      }
      else if (is.vector(mat, "numeric")) mat <- matrix(mat, ncol = 1)
      else if (!is.matrix(mat) || !is.numeric(mat)) bad.mat <- TRUE
    }
    if (bad.mat) stop("mat must be a numeric matrix.")
    
    if (missing(treat) || !(is.factor(treat) || is.atomic(treat)) || !is_binary(treat)) stop("treat must be a binary variable.")
    
    if (is_null(weights)) weights <- rep(1, NROW(mat))
    if (is_null(s.weights)) s.weights <- rep(1, NROW(mat))
    if (!all_the_same(c(NROW(mat), length(treat), length(weights), length(s.weights)))) {
      stop("mat, treat, weights, and s.weights (if supplied) must have the same number of units.")
    }
  }
  
  if (any(paste.("smd", agg.funs) %in% type)) {
    if (!exists("s.d.denom")) s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = A[["estimand"]], weights = data.frame(weights), treat = treat)
    smd <- col_w_smd(mat, treat, weights, std = TRUE, s.d.denom, abs = TRUE, s.weights = s.weights, check = FALSE)
    if (is_null(A[["smd.weights"]])) smd.weights <- rep(1, ncol(mat))
    else if (!is.vector(A[["smd.weights"]], "numeric")) {
      warning("smd.weights is not numeric. Ignoring smd.weights", 
              call. = FALSE, immediate. = TRUE)
      smd.weights <- rep(1, ncol(mat))
    }
    else if (length(A[["smd.weights"]]) == ncol(mat)) {
      smd.weights <- A[["smd.weights"]]
    }
    else {
      warning("smd.weights should be of length ncol(mat). Ignoring smd.weights", 
              call. = FALSE, immediate. = TRUE)
      smd.weights <- rep(1, ncol(mat))
    }
    smd.weights <- smd.weights/mean(smd.weights) #Make sum to ncol(mat)
    smd <- smd.weights*smd
  }
  if (any(paste.("ks", agg.funs) %in% type)) {
    ks <- col_w_ks(mat, treat, weights, bin.vars = A[["bin.vars"]], check = is_null(A[["bin.vars"]]))
  }
  if (any(paste.("ovl", agg.funs) %in% type)) {
    ovl <- do.call(col_w_ovl, c(list(mat = mat, treat = treat, weights = weights, check = is_null(A[["bin.vars"]])), A))
  }
  if ("gwd" %in% type) {
    if (!exists("s.d.denom")) s.d.denom <- get.s.d.denom(A[["s.d.denom"]], estimand = A[["estimand"]], weights = data.frame(weights), treat = treat)
  }
  
  bal <- setNames(vapply(type, function(m) {
    if (endsWith(m, ".mean")) {
      agg <- mean
      m <- substr(m, 1, nchar(m) - nchar(".mean"))
    }
    else if (endsWith(m, ".max")) {
      agg <- max
      m <- substr(m, 1, nchar(m) - nchar(".max"))
    }
    else if (endsWith(m, ".rms")) {
      agg <- function(x, ...) {sqrt(mean(x^2, ...))}
      m <- substr(m, 1, nchar(m) - nchar(".rms"))
    }
    
    if (m %in% uni.type) {
      return(agg(get0(m), na.rm = TRUE))
    }
    else if (m == "mahalanobis") {
      if (is_null(s.weights)) s.weights <- rep(1, nrow(mat))
      mdiff <- matrix(col_w_smd(mat, treat, weights, std = FALSE, abs = FALSE, check = FALSE), ncol = 1)
      wcov <- cov.wt(mat, s.weights)$cov
      mahal <- crossprod(mdiff, solve(wcov)) %*% mdiff
      return(mahal)
    }
    else if (m == "cstat") {
      tval1 <- treat[1]
      d <- data.frame(treat, mat)
      f <- formula(d)
      pred <- glm(f, data = d, family = quasibinomial(),
                  weights = weights)$fitted
      wi <- wilcox.test(pred ~ treat)
      cstat <- wi$statistic/(sum(treat==tval1)*sum(treat!=tval1))
      cstat <- 2*max(cstat, 1-cstat)-1
      return(cstat)
    }
    else if (m == "wr2") {
      tval1 <- treat[1]
      d <- data.frame(treat, mat)
      f <- formula(d)
      fit <- glm(f, data = d, family = quasibinomial(),
                 weights = weights)
      r2 <- 1 - (pi^2/3)/(3*var(fit$linear.predictors) + pi^2/3)
      return(r2)
    }
    else if (m == "gwd") {
      co.names <- setNames(lapply(colnames(mat), function(x) setNames(list(x, "base"), c("component", "type"))), colnames(mat))
      new <- int.poly.f(mat, int = TRUE, poly = 2, center = isTRUE(A[["center"]]), 
                        sep = getOption("cobalt_int_sep", default = " * "), co.names = co.names)
      mat_ <- cbind(mat, new)
      
      smd <- col_w_smd(mat_, treat, weights, std = TRUE, s.d.denom, abs = TRUE, check = FALSE)
      if (is_null(A[["gwd.weights"]])) gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
      else if (!is.vector(A[["gwd.weights"]], "numeric")) {
        warning("gwd.weights is not numeric. Ignoring gwd.weights.", 
                call. = FALSE, immediate. = TRUE)
        gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
      }
      else if (length(A[["gwd.weights"]]) == 1L) {
        gwd.weights <- rep(A[["gwd.weights"]], ncol(mat_))
      }
      else if (length(A[["gwd.weights"]]) == 2L) {
        gwd.weights <- c(rep(A[["gwd.weights"]][1], ncol(mat)), rep(A[["gwd.weights"]][2], ncol(new)))
      }
      else {
        warning("gwd.weights should be of length 1 or 2. Ignoring gwd.weights.", 
                call. = FALSE, immediate. = TRUE)
        gwd.weights <- c(rep(1, ncol(mat)), rep(.5, ncol(new)))
      }
      
      gwd.weights <- gwd.weights/sum(gwd.weights) #Make sum to 1
      return(sum(gwd.weights*smd, na.rm = TRUE))
    }
    else if (m == "design.effect") {
      tval1 <- treat[1]
      q <- sum(treat == tval1)/length(treat)
      des.eff <- function(w) length(w)*sum(w^2)/sum(w)^2
      des <- c(des.eff(weights[treat == tval1]), des.eff(weights[treat != tval1]))
      de <- des[1]*(1-q) + des[2]*q
      return(de)
    }
  }, numeric(1L)), type)
  
  return(bal)
}

bal.tab <- function(...) {
  
  if (...length() == 0L) stop("No arguments were supplied.", call. = FALSE)
  .obj <- ...elt(1)
  
  .obj <- is.designmatch(.obj)
  .obj <- is.time.list(.obj)
  
  #Replace .all and .none with NULL and NA respectively
  .call <- match.call(expand.dots = TRUE)
  .alls <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.all)), logical(1L))
  .nones <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.none)), logical(1L))
  if (any(c(.alls, .nones))) {
    .call[.alls] <- expression(NULL)
    .call[.nones] <- expression(NA)
    return(eval(.call))
  }
  
  UseMethod("bal.tab", .obj)
}

#Point treatments
bal.tab.matchit <- function(m, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.matchit", c(list(m), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.ps <- function(ps, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.ps", c(list(ps), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.mnps <- function(mnps, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.mnps", c(list(mnps), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  
  return(out)
}
bal.tab.ps.cont <- function(ps.cont, stop.method, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, r.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.ps.cont", c(list(ps.cont), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.Match <- function(M, formula = NULL, data = NULL, treat = NULL, covs = NULL, int = FALSE, poly = 1, distance = NULL, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base", c(list(M), args), quote = TRUE) 
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.formula <- function(formula, data = NULL, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.formula", c(list(formula = formula), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  
  return(out)
}
bal.tab.data.frame <- function(covs, treat, data = NULL, weights = NULL, distance = NULL, subclass = NULL, match.strata = NULL, method, int = FALSE, poly = 1, addl = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, estimand = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-c(1, length(formals()))])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  X <- do.call("x2base.data.frame", c(list(covs), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  
  return(out)
}
bal.tab.numeric <- function(treat, covs, ...) {
  bal.tab(covs, treat = treat, ...)
}
bal.tab.factor <- bal.tab.character <- bal.tab.logical <- bal.tab.numeric
bal.tab.CBPS <- function(cbps, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, pairwise = TRUE, focal = NULL, s.weights = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base", c(list(cbps), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.weightit <- function(weightit, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL,  continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, cluster = NULL, imp = NULL, pairwise = TRUE, focal = NULL, abs = FALSE, subset = NULL, quick = TRUE, ... ) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base", c(list(weightit), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.ebalance <- function(ebal, ...) {
  bal.tab.Match(ebal, ...)
}
bal.tab.ebalance.trim <- bal.tab.ebalance
bal.tab.optmatch <- function(optmatch, ...) {
  bal.tab.Match(optmatch, ...)
}
bal.tab.designmatch <- function(dm, ...) {
  class(dm) <- "designmatch"
  bal.tab.Match(dm, ...)
}
bal.tab.mimids <- function(mimids, int = FALSE, poly = 1, distance = NULL, addl = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, cluster = NULL, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.mimids", c(list(mimids), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.wimids <- bal.tab.mimids

#MSMs wth multiple time points
bal.tab.formula.list <- function(formula.list, data = NULL, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  X <- do.call("x2base.formula.list", c(list(formula.list = formula.list), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.data.frame.list <- function(covs.list, treat.list = NULL, data = NULL, weights = NULL, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, method, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, pairwise = TRUE, s.weights = NULL, estimand = "ATE", abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  X <- do.call("x2base.data.frame.list", c(list(covs.list = covs.list), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.iptw <- function(iptw, stop.method, int = FALSE, poly = 1, distance.list = NULL, addl.list = NULL, data = NULL, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, pairwise = TRUE, abs = FALSE, subset = NULL, quick = TRUE, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  if (is_not_null(args$cluster)) stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
  if (is_not_null(args$imp)) stop("Multiply imputed data sets are not yet supported with longitudinal treatments.", call. = FALSE)
  
  #Initializing variables
  X <- do.call("x2base.iptw", c(list(iptw), args), quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  return(out)
}
bal.tab.CBMSM <- bal.tab.CBPS

#default method
bal.tab.default <- function(obj, ...) {
  
  args <- c(as.list(environment()), list(...))[-1]
  
  #Adjustments to arguments
  args.with.choices <- names(formals()[-1])[vapply(formals()[-c(1, length(formals()))], function(x) length(x)>1, logical(1L))]
  for (i in args.with.choices) args[[i]] <- eval(parse(text=paste0("match_arg(", i, ")")))
  
  blank.args <- vapply(formals()[-c(1, length(formals()))], function(x) identical(x, quote(expr =)), logical(1L))
  if (any(blank.args)) {
    for (arg.name in names(blank.args)[blank.args]) {
      if (identical(args[[arg.name]], quote(expr = ))) {
        args[[arg.name]] <- NULL
      }
    }
  }
  
  #Initializing variables
  X <- do.call("x2base.default", c(list(obj = obj), args),
               quote = TRUE)
  
  args <- args[names(args) %nin% names(X)]
  
  out <- do.call(paste.("base.bal.tab", class(X)), c(X, args),
                 quote = TRUE)
  
  return(out)
}

base.bal.tab.binary <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, addl.sds = NULL, ...) {
  #Preparations
  args <- list(...)
  
  if (nunique(treat) != 2) {
    stop("Treatment indicator must be a binary (0, 1) variable---i.e., treatment (1) or control (0)", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    if (is.factor(treat)) treat.names <- unique.treat <- levels(treat)
    else treat.names <- unique.treat <- unique(treat, nmax = 2)
  }
  else {
    treat.names <- c("Control", "Treated")
    unique.treat <- sort(unique(treat, nmax = 2))
  }
  names(treat.names) <- unique.treat
  
  check_if_zero_weights(weights, treat, unique.treat, treat.type = "cat")
  
  treat <- binarize(treat)
  if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
  if (is_not_null(v.threshold)) {
    v.threshold <- max(v.threshold, 1/v.threshold)
    disp.v.ratio <- TRUE
  }
  if (is_null(ks.threshold) && is_null(args$k.threshold)) {
    ks.threshold <- args$k.threshold
  }
  if (is_not_null(ks.threshold)) {
    if (ks.threshold > 1) {
      warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
      ks.threshold <- NULL
    }
    else disp.ks <- TRUE
  }
  if (is_null(weights) && is_null(subclass)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat))
  }
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
  
  #Actions
  if (nunique.gt(cluster, 1)) {
    out.names <- c("Cluster.Balance", 
                   "Cluster.Balance.Across.Subclass", 
                   "Cluster.Summary", "Observations",
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    cluster <- factor(cluster)
    
    C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
    C.list <- setNames(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), 
                       levels(cluster))
    types <- get.types(C)
    co.names <- attr(C, "co.names")
    
    if (length(method) == 1 && method == "subclassification") {
      stop("Subclassification with clusters is not yet supported.", call. = FALSE)
      #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab")
    }
    else {
      out[["Cluster.Balance"]] <- setNames(lapply(levels(cluster), function(c) setNames(list(balance.table(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights[cluster == c], abs = abs, no.adj = no.adj, types = types, quick = quick),
                                                                                             samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                                        c("Balance", "Observations"))),
                                           levels(cluster))
      balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
      observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
      
      if (cluster.summary || !quick) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                               weight.names = names(weights),
                                                                                               no.adj = no.adj,
                                                                                               abs = abs,
                                                                                               quick = quick,
                                                                                               types = types)
      
      for (i in names(attr(balance.tables[[1]], "disp"))) {
        if (all(vapply(balance.tables, function(x) !attr(x, "disp")[i], logical(1L)))) assign(paste0("disp.", i), FALSE)
      }
      for (i in names(attr(balance.tables[[1]], "disp"))) {
        if (all(vapply(balance.tables, function(x) is_null(attr(x, "threshold")[i]), logical(1L)))) assign(paste0(i, ".threshold"), NULL)
      }
      
      # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["means"], logical(1L)))) {disp.means <- FALSE}
      # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["sds"], logical(1L)))) {disp.sds <- FALSE}
      # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["v.ratio"], logical(1L)))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
      # if (all(vapply(balance.tables, function(x) !attr(x, "disp")["ks"], logical(1L)))) {disp.ks <- FALSE; ks.threshold <- NULL}
      
      out <- out[names(out) %nin% "Cluster.Balance.Across.Subclass"]
      out[["Observations"]] <- samplesize.across.clusters(observations)
      out[["call"]] <- call
      attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                         v.threshold=v.threshold,
                                         ks.threshold=ks.threshold,
                                         imbalanced.only = imbalanced.only,
                                         un=un, 
                                         disp.means=disp.means, 
                                         disp.sds=disp.sds,
                                         disp.v.ratio=disp.v.ratio, 
                                         disp.ks=disp.ks, 
                                         disp.adj=!no.adj, 
                                         disp.subclass=disp.subclass,
                                         disp.bal.tab = disp.bal.tab, 
                                         which.cluster=which.cluster,
                                         cluster.summary=cluster.summary,
                                         cluster.fun = cluster.fun,
                                         abs = abs,
                                         continuous = continuous,
                                         binary = binary,
                                         quick = quick,
                                         nweights = ifelse(no.adj, 0, ncol(weights)),
                                         weight.names = names(weights),
                                         treat.names = treat.names,
                                         co.names = co.names)
      class(out) <- c("bal.tab.cluster", "bal.tab")
    }
    
  }
  else {
    if (length(method) == 1 && method == "subclassification") {
      if (is_not_null(subclass)) {
        out.names <- c("Subclass.Balance", "Balance.Across.Subclass", 
                       expand.grid_string(c("Balanced", "Max.Imbalance"),
                                          c("Means", "Variances", "KS"),
                                          "Subclass", collapse = "."), 
                       "Subclass.Observations", "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
        co.names <- attr(C, "co.names")
        
        if (is_not_null(list(...)$sub.by)) sub.by <- list(...)$sub.by
        else sub.by <- call$sub.by
        out[["Subclass.Balance"]] <- balance.table.subclass(C, weights=weights[[1]], treat=treat, subclass=subclass, continuous=continuous, binary=binary, s.d.denom=s.d.denom[1], m.threshold=m.threshold, v.threshold=v.threshold, ks.threshold = ks.threshold, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, abs = abs, quick = quick)
        out[["Subclass.Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
        out[["Balance.Across.Subclass"]] <- balance.table.across.subclass(balance.table = balance.table(C, 
                                                                                                        weights = weights[[1]], 
                                                                                                        treat = treat, 
                                                                                                        continuous = continuous, 
                                                                                                        binary = binary, 
                                                                                                        s.d.denom = s.d.denom[1], 
                                                                                                        m.threshold = m.threshold, 
                                                                                                        v.threshold = v.threshold, 
                                                                                                        ks.threshold = ks.threshold,
                                                                                                        un = un, 
                                                                                                        disp.means = disp.means, 
                                                                                                        disp.sds = disp.sds,
                                                                                                        disp.v.ratio = disp.v.ratio, 
                                                                                                        disp.ks = disp.ks,
                                                                                                        abs = abs, 
                                                                                                        no.adj = TRUE, quick = quick), 
                                                                          balance.table.subclass.list=out[["Subclass.Balance"]], 
                                                                          subclass.obs=out[["Subclass.Observations"]], 
                                                                          sub.by=sub.by, 
                                                                          m.threshold=m.threshold, 
                                                                          v.threshold=v.threshold, 
                                                                          ks.threshold=ks.threshold,
                                                                          s.d.denom = s.d.denom[1])
        #Reassign disp... and ...threshold based on balance table output
        for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
          assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
        }
        for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
          assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
        }
        
        S <- list(diff = list(threshold = m.threshold,
                              Names = "Means",
                              Threshold = "M.Threshold",
                              Stat = "Diff"),
                  v.ratio = list(threshold = v.threshold,
                                 Names = "Variances",
                                 Threshold = "V.Threshold",
                                 Stat = "V.Ratio"),
                  ks = list(threshold = ks.threshold,
                            Names = "KS",
                            Threshold = "KS.Threshold",
                            Stat = "KS"))
        
        for (s in S) {
          if (is_not_null(s[["threshold"]])) {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
            names(out[[paste.("Balanced", s[["Names"]], "Subclass")]]) <- paste("Subclass", levels(subclass))
            max.imbal.list <- lapply(levels(subclass), function(x) {
              return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
            } )
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- data.frame(do.call("rbind", max.imbal.list), 
                                                                                   row.names = paste("Subclass", levels(subclass)))
          }
          else {
            out[[paste.("Balanced", s[["Names"]], "Subclass")]] <- NULL
            out[[paste.("Max.Imbalance", s[["Names"]], "Subclass")]] <- NULL
          }
        }
        
        out[["call"]] <- call
        attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                           v.threshold=v.threshold, 
                                           ks.threshold=ks.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un, 
                                           disp.means=disp.means, 
                                           disp.sds=disp.sds,
                                           disp.v.ratio=disp.v.ratio, 
                                           disp.ks=disp.ks, 
                                           disp.adj=!no.adj, 
                                           disp.subclass=disp.subclass,
                                           disp.bal.tab = disp.bal.tab, 
                                           abs = abs,
                                           continuous = continuous,
                                           binary = binary,
                                           quick = quick,
                                           treat.names = treat.names,
                                           co.names = co.names)
        class(out) <- c("bal.tab.subclass", "bal.tab")
      }
      else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
    }
    else {
      out.names <- c("Balance", 
                     expand.grid_string(c("Balanced", "Max.Imbalance"),
                                        c("Means", "Variances", "KS"),
                                        collapse = "."), 
                     "Observations", "call")
      out <- vector("list", length(out.names))
      names(out) <- out.names
      
      C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
      co.names <- attr(C, "co.names")
      
      out[["Balance"]] <- balance.table(C, weights, treat, continuous, binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, s.weights = s.weights, abs = abs, no.adj = no.adj, quick = quick, addl.sds = addl.sds)
      
      #Reassign disp... and ...threshold based on balance table output
      for (i in names(attr(out[["Balance"]], "disp"))) {
        assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
      }
      for (i in names(attr(out[["Balance"]], "threshold"))) {
        assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
      }
      
      S <- list(diff = list(threshold = m.threshold,
                            Names = "Means",
                            Threshold = "M.Threshold",
                            Stat = "Diff"),
                v.ratio = list(threshold = v.threshold,
                               Names = "Variances",
                               Threshold = "V.Threshold",
                               Stat = "V.Ratio"),
                ks = list(threshold = ks.threshold,
                          Names = "KS",
                          Threshold = "KS.Threshold",
                          Stat = "KS"))
      
      for (s in S) {
        if (is_not_null(s[["threshold"]])) {
          if (no.adj) {
            out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
          }
          else if (ncol(weights) == 1) {
            out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
          }
          else if (ncol(weights) > 1) {
            out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                names(weights))
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(weights),
                                                                  do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                               c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                  stringsAsFactors = FALSE)
          }
        }
        else {
          out[[paste.("Balanced", s[["Names"]])]] <- NULL
          out[[paste.("Max.Imbalance", s[["Names"]])]] <- NULL
        }
      }
      
      out[["Observations"]] <- samplesize(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded, treat.names = treat.names)
      out[["call"]] <- call
      attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                         v.threshold=v.threshold, 
                                         ks.threshold=ks.threshold, 
                                         imbalanced.only = imbalanced.only,
                                         un=un, 
                                         disp.means=disp.means, 
                                         disp.sds=disp.sds,
                                         disp.v.ratio=disp.v.ratio, 
                                         disp.ks=disp.ks, 
                                         disp.adj=!no.adj,
                                         disp.bal.tab = disp.bal.tab, 
                                         abs = abs,
                                         continuous = continuous,
                                         binary = binary,
                                         quick = quick,
                                         nweights = ifelse(no.adj, 0, ncol(weights)),
                                         weight.names = names(weights),
                                         treat.names = treat.names,
                                         co.names = co.names)
      class(out) <- "bal.tab"
    }
  }
  
  #attr(out, "int") <- int
  return(out)
}
base.bal.tab.cont <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
  
  #Preparations
  args <- list(...)
  
  check_if_zero_weights(weights, treat.type = "cont")
  
  if (is_not_null(r.threshold)) {
    r.threshold <- abs(r.threshold)
    if (r.threshold > 1) {
      warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
      r.threshold <- NULL
    }
  }
  if (is_null(weights) && is_null(subclass)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat))
  }    
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
  
  #Actions
  if (nlevels(cluster) > 0) {
    out.names <- c("Cluster.Balance", 
                   "Cluster.Balance.Across.Subclass", 
                   "Cluster.Summary", "Observations",
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    cluster <- factor(cluster)
    
    out[["Cluster.Balance"]] <- vector("list", length(levels(cluster)))
    names(out[["Cluster.Balance"]]) <- levels(cluster)
    
    C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
    co.names <- attr(C, "co.names")
    C.list <- structure(lapply(levels(cluster), function(x) C[cluster == x, , drop = FALSE]), names = levels(cluster))
    types <- get.types(C)
    
    if (length(method) == 1 && method == "subclassification") {
      stop("Subclassification with clusters is not yet supported.", call. = FALSE)
      #class(out) <- c("bal.tab.cluster", "bal.tab.subclass", "bal.tab") #add more for subclasses
    }
    else {
      out[["Cluster.Balance"]] <- lapply(levels(cluster), function(c) setNames(list(balance.table.cont(C = C.list[[c]], weights = weights[cluster == c, , drop = FALSE], treat = treat[cluster == c], r.threshold = r.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, s.weights = s.weights[cluster == c], abs = abs, no.adj = no.adj, types = types, quick = quick),
                                                                                    samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, cluster = cluster, which.cluster = c, discarded = discarded)), 
                                                                               c("Balance", "Observations")))
      names(out[["Cluster.Balance"]]) <- levels(cluster)
      
      balance.tables <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Balance"]])
      observations <- lapply(levels(cluster), function(c) out[["Cluster.Balance"]][[c]][["Observations"]])
      
      if (!(!cluster.summary && quick)) out[["Cluster.Summary"]] <- balance.table.cluster.summary(balance.tables,
                                                                                                  weight.names = names(weights),
                                                                                                  no.adj = no.adj,
                                                                                                  abs = abs,
                                                                                                  quick = quick,
                                                                                                  types = types)
      out <- out[names(out) %nin% "Cluster.Balance.Across.Subclass"]
      out[["Observations"]] <- samplesize.across.clusters(observations)
      
      out[["call"]] <- call
      attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                         imbalanced.only = imbalanced.only,
                                         un=un, 
                                         disp.means=disp.means, 
                                         disp.sds = disp.sds,
                                         disp.adj=!no.adj, 
                                         disp.bal.tab = disp.bal.tab,
                                         which.cluster=which.cluster,
                                         cluster.summary=cluster.summary,
                                         cluster.fun=cluster.fun,
                                         abs = abs,
                                         quick = quick,
                                         nweights = ifelse(no.adj, 0, ncol(weights)),
                                         weight.names = names(weights),
                                         co.names = co.names)
      class(out) <- c("bal.tab.cont.cluster", "bal.tab.cluster", "bal.tab.cont", "bal.tab")
    }
    
  }
  else {
    if (length(method) == 1 && method == "subclassification") {
      if (is_not_null(subclass)) {
        out.names <- c("Subclass.Balance", 
                       "Balanced.Corr.Subclass", "Max.Imbalance.Corr.Subclass", 
                       "Subclass.Observations", "call")
        out <- vector("list", length(out.names))
        names(out) <- out.names
        
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
        co.names <- attr(C, "co.names")
        
        # if (length(list(...)$sub.by > 0)) sub.by <- list(...)$sub.by
        # else sub.by <- call$sub.by
        
        out[["Subclass.Balance"]] <- balance.table.subclass.cont(C, weights=weights[[1]], treat=treat, subclass=subclass, r.threshold=r.threshold, disp.means = disp.means, disp.sds = disp.sds, quick = quick)
        out[["Subclass.Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
        
        #Reassign disp... and ...threshold based on balance table output
        for (i in names(attr(out[["Subclass.Balance"]], "disp"))) {
          assign(paste.("disp", i), attr(out[["Subclass.Balance"]], "disp")[i])
        }
        for (i in names(attr(out[["Subclass.Balance"]], "threshold"))) {
          assign(paste.(i, "threshold"), attr(out[["Subclass.Balance"]], "threshold")[i])
        }
        
        S <- list(corr = list(threshold = r.threshold,
                              Names = "Corr",
                              Threshold = "R.Threshold",
                              Stat = "Corr"))
        
        for (s in S) {
          if (is_not_null(s[["threshold"]])) {
            out[[paste.("Balanced", s[["Stat"]], "Subclass")]] <- as.data.frame(lapply(levels(subclass), function(x) baltal(out[["Subclass.Balance"]][[x]][[s[["Threshold"]]]])))
            names(out[[paste.("Balanced", s[["Stat"]], "Subclass")]]) <- paste("Subclass", levels(subclass))
            mi.list <- lapply(levels(subclass), function(x) {
              return(max.imbal(out[["Subclass.Balance"]][[x]][out[["Subclass.Balance"]][[x]][["Type"]] != "Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio"))
            } )
            mi <- do.call("rbind", mi.list)
            out[[paste.("Max.Imbalance", s[["Stat"]], "Subclass")]] <- data.frame(mi, row.names = paste("Subclass", levels(subclass)))
          }
        }
        
        out[["call"]] <- call
        attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                           imbalanced.only = imbalanced.only,
                                           un=un,
                                           disp.means=disp.means, 
                                           disp.sds=disp.sds,
                                           disp.adj=!no.adj, 
                                           disp.subclass=disp.subclass,
                                           disp.bal.tab = disp.bal.tab,
                                           abs = abs,
                                           quick = quick,
                                           co.names = co.names)
        class(out) <- c("bal.tab.subclass.cont", "bal.tab.subclass", "bal.tab.cont", "bal.tab")
      }
      else stop("Method specified as subclassification, but no subclasses were specified.", call. = FALSE)
      
    }
    else {
      out.names <- c("Balance", "Balanced.Corr", 
                     "Max.Imbalance.Corr", 
                     "Observations", 
                     "call")
      out <- vector("list", length(out.names))
      names(out) <- out.names
      
      C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, cluster = cluster, ...)
      co.names <- attr(C, "co.names")
      
      out[["Balance"]] <- balance.table.cont(C, weights, treat, r.threshold, un = un, disp.means = disp.means, disp.sds = disp.sds, s.weights = s.weights, abs = abs, no.adj = no.adj, quick = quick)
      
      #Reassign disp... and ...threshold based on balance table output
      for (i in names(attr(out[["Balance"]], "disp"))) {
        assign(paste0("disp.", i), attr(out[["Balance"]], "disp")[i])
      }
      for (i in names(attr(out[["Balance"]], "threshold"))) {
        assign(paste0(i, ".threshold"), attr(out[["Balance"]], "threshold")[i])
      }
      
      S <- list(corr = list(threshold = r.threshold,
                            Names = "Corr",
                            Threshold = "R.Threshold",
                            Stat = "Corr"))
      
      for (s in S) {
        if (is_not_null(s[["threshold"]])) {
          if (no.adj) {
            out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[paste.(s[["Threshold"]], "Un")]])
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Un"), paste.(s[["Threshold"]], "Un"), ratio = s$Stat == "V.Ratio")
          }
          else if (ncol(weights) == 1) {
            out[[paste.("Balanced", s[["Names"]])]] <- baltal(out[["Balance"]][[s[["Threshold"]]]])
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], "Adj"), s[["Threshold"]], ratio = s$Stat == "V.Ratio")
          }
          else if (ncol(weights) > 1) {
            out[[paste.("Balanced", s[["Names"]])]] <- setNames(do.call("cbind", lapply(names(weights), function(x) baltal(out[["Balance"]][[paste.(s[["Threshold"]], x)]]))),
                                                                names(weights))
            out[[paste.("Max.Imbalance", s[["Names"]])]] <- cbind(Weights = names(weights),
                                                                  do.call("rbind", lapply(names(weights), function(x) setNames(max.imbal(out[["Balance"]][out[["Balance"]][["Type"]]!="Distance", , drop = FALSE], paste.(s[["Stat"]], x), paste.(s[["Threshold"]], x), ratio = s$Stat == "V.Ratio"),
                                                                                                                               c("Variable", s[["Stat"]], s[["Threshold"]])))),
                                                                  stringsAsFactors = FALSE)
          }
        }
      }
      
      out[["Observations"]] <- samplesize.cont(treat = treat, weights = weights, subclass = subclass, s.weights = s.weights, method = method, discarded = discarded)
      out[["call"]] <- call
      attr(out, "print.options") <- list(r.threshold=r.threshold, 
                                         imbalanced.only = imbalanced.only,
                                         un=un, 
                                         disp.means=disp.means, 
                                         disp.sds = disp.sds,
                                         disp.adj=!no.adj,
                                         disp.bal.tab = disp.bal.tab,
                                         abs = abs,
                                         quick = quick,
                                         nweights = ifelse(no.adj, 0, ncol(weights)),
                                         weight.names = names(weights),
                                         co.names = co.names)
      class(out) <- c("bal.tab.cont", "bal.tab")
    }
  }
  
  attr(out, "int") <- int
  return(out)
}
base.bal.tab.imp <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), imp = NULL, which.imp = NA, imp.summary = getOption("cobalt_imp.summary", TRUE), imp.fun = getOption("cobalt_imp.fun", NULL), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
  args <- list(...)
  #Preparations
  if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
  if (is_not_null(v.threshold)) {
    v.threshold <- max(v.threshold, 1/v.threshold)
    disp.v.ratio <- TRUE
  }
  if (is_not_null(ks.threshold)) {
    if (ks.threshold > 1) {
      warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
      ks.threshold <- NULL
    }
    else disp.ks <- TRUE
  }
  if (is_not_null(r.threshold)) {
    r.threshold <- abs(r.threshold)
    if (r.threshold > 1) {
      warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
      r.threshold <- NULL
    }
  }
  if (is_null(weights)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat))
  }    
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- imp.fun <- args[["agg.fun"]]
  
  #Setup output object
  out.names <- c("Imputation.Balance", 
                 "Cluster.Balance.Across.Imputations",
                 "Balance.Across.Imputations", 
                 "Observations", 
                 "call")
  out <- vector("list", length(out.names))
  names(out) <- out.names
  
  #Get list of bal.tabs for each imputation
  if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
    args <- args[names(args) %nin% names(formals(base.bal.tab.cont))]
    out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) do.call("base.bal.tab.cont", 
                                                                           c(list(weights = weights[imp==i, , drop  = FALSE], 
                                                                                  treat = treat[imp==i], 
                                                                                  distance = distance[imp==i, , drop = FALSE], 
                                                                                  subclass = subclass[imp==i], 
                                                                                  covs = covs[imp == i, , drop = FALSE], 
                                                                                  call = call, 
                                                                                  int = int, 
                                                                                  poly = poly, 
                                                                                  addl = addl[imp = i, , drop = FALSE], 
                                                                                  r.threshold = r.threshold, 
                                                                                  imbalanced.only = imbalanced.only, 
                                                                                  un = un, 
                                                                                  disp.means = disp.means, 
                                                                                  disp.sds = disp.sds, 
                                                                                  disp.bal.tab = disp.bal.tab, 
                                                                                  method = method, 
                                                                                  cluster = cluster[imp==i], 
                                                                                  which.cluster = which.cluster, 
                                                                                  cluster.summary = cluster.summary, 
                                                                                  s.weights = s.weights[imp==i], 
                                                                                  discarded = discarded[imp==i], 
                                                                                  abs = abs,
                                                                                  quick = quick), 
                                                                             args), quote = TRUE))
  }
  else if (isTRUE(get.treat.type(treat) == "multinomial") || (is_(treat, c("factor", "character")) && !is_binary(treat))) {
    args <- args[names(args) %nin% names(formals(base.bal.tab.multi))]
    stop("Multinomial treatments are not yet supported with multiply imputed data.", call. = FALSE)
  }
  else {#if binary treatment
    args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
    out[["Imputation.Balance"]] <- lapply(levels(imp), function(i) do.call("base.bal.tab.binary", 
                                                                           c(list(weights = weights[imp==i, , drop = FALSE], 
                                                                                  treat = treat[imp==i], 
                                                                                  distance = distance[imp==i, , drop = FALSE], 
                                                                                  subclass = subclass[imp==i], 
                                                                                  covs = covs[imp==i, , drop = FALSE], 
                                                                                  call = call, 
                                                                                  int = int, 
                                                                                  poly = poly, 
                                                                                  addl = addl[imp==i, , drop = FALSE], 
                                                                                  continuous = continuous, 
                                                                                  binary = binary, 
                                                                                  s.d.denom = s.d.denom, 
                                                                                  m.threshold = m.threshold, 
                                                                                  v.threshold = v.threshold, 
                                                                                  ks.threshold = ks.threshold, 
                                                                                  imbalanced.only = imbalanced.only, 
                                                                                  un = un, 
                                                                                  disp.means = disp.means, 
                                                                                  disp.sds = disp.sds, 
                                                                                  disp.v.ratio = disp.v.ratio, 
                                                                                  disp.ks = disp.ks, 
                                                                                  disp.subclass = disp.subclass, 
                                                                                  disp.bal.tab = disp.bal.tab,
                                                                                  method = method,
                                                                                  cluster = cluster[imp==i],
                                                                                  which.cluster = which.cluster,
                                                                                  cluster.summary = cluster.summary,
                                                                                  s.weights = s.weights[imp==i],
                                                                                  discarded = discarded[imp==i],
                                                                                  abs = abs,
                                                                                  quick = quick), 
                                                                             args), quote = TRUE))
  }
  
  names(out[["Imputation.Balance"]]) <- levels(imp)
  
  #Create summary of lists
  
  if ("bal.tab.cluster" %in% class(out[["Imputation.Balance"]][[1]])) {
    if (imp.summary || !quick) {
      out[["Cluster.Balance.Across.Imputations"]] <- lapply(levels(cluster), 
                                                            function(c) setNames(list(balance.table.imp.summary(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Balance"]]), 
                                                                                                                weight.names = names(weights),
                                                                                                                no.adj = no.adj,
                                                                                                                abs = abs, quick = quick),
                                                                                      samplesize.across.imps(lapply(out[["Imputation.Balance"]], function(i) i[["Cluster.Balance"]][[c]][["Observations"]]))), 
                                                                                 c("Cluster.Balance", "Cluster.Observations")))
      names(out[["Cluster.Balance.Across.Imputations"]]) <- levels(cluster)
      balance.tables <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Balance"]])
      observations <- lapply(out[["Cluster.Balance.Across.Imputations"]], function(c) c[["Cluster.Observations"]])
      
      out[["Balance.Across.Imputations"]] <- balance.table.clust.imp.summary(balance.tables,
                                                                             weight.names = names(weights),
                                                                             no.adj = no.adj,
                                                                             abs = abs,
                                                                             quick = quick,
                                                                             types = NULL)
      out[["Observations"]] <- samplesize.across.clusters(observations)
    }
    
    classes <- c("bal.tab.imp.cluster", "bal.tab.imp")
  }
  else {
    if ("bal.tab.subclass" %in% class(out[["Imputation.Balance"]][[1]])) {
      #Put something here
      stop("Subclassification cannot be used with multiply imputed data.", call. = FALSE)
    }
    else {
      if (imp.summary || !quick) out[["Balance.Across.Imputations"]] <- balance.table.imp.summary(bal.tab.imp.list = out[["Imputation.Balance"]], 
                                                                                                  weight.names = names(weights),
                                                                                                  no.adj = no.adj,
                                                                                                  abs = abs,
                                                                                                  quick = quick,
                                                                                                  types = NULL)
      observations <- lapply(out[["Imputation.Balance"]], function(x) x[["Observations"]])
      
      out[["Observations"]] <- samplesize.across.imps(observations)
      classes <- "bal.tab.imp"
    }
  }
  
  out[["call"]] <- call
  attr(out, "print.options") <- list(m.threshold=m.threshold,
                                     v.threshold=v.threshold,
                                     ks.threshold=ks.threshold,
                                     r.threshold=r.threshold,
                                     imbalanced.only = imbalanced.only,
                                     un=un, 
                                     disp.adj=!no.adj, 
                                     which.cluster=which.cluster,
                                     cluster.summary=cluster.summary,
                                     cluster.fun = cluster.fun,
                                     which.imp=which.imp,
                                     imp.summary=imp.summary,
                                     imp.fun = imp.fun,
                                     abs = abs,
                                     continuous = continuous,
                                     binary = binary,
                                     quick = quick,
                                     disp.means=disp.means, 
                                     disp.sds = disp.sds,
                                     disp.v.ratio=disp.v.ratio, 
                                     disp.ks=disp.ks,
                                     disp.bal.tab = disp.bal.tab,
                                     nweights = ifelse(no.adj, 0, ncol(weights)),
                                     weight.names = names(weights),
                                     co.names = attr(out[["Imputation.Balance"]][[1]], "print.options")[["co.names"]])
  class(out) <- unique(c(classes, sapply(out[["Imputation.Balance"]], class)))
  
  return(out)
}
base.bal.tab.multi <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = getOption("cobalt_multi.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
  #Preparations
  args <- list(...)
  if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
  if (is_not_null(v.threshold)) {
    v.threshold <- max(v.threshold, 1/v.threshold)
    disp.v.ratio <- TRUE
  }
  if (is_not_null(ks.threshold)) {
    if (ks.threshold > 1) {
      warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
      ks.threshold <- NULL
    }
    else disp.ks <- TRUE
  }
  if (is_null(weights) && is_null(subclass)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat))
  }
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
  
  #Treat is a factor variable of 3+ levels
  if (is_null(focal)) {
    if (pairwise) treat.combinations <- combn(levels(treat), 2, list)
    else treat.combinations <- lapply(levels(treat), function(x) c(x, "Others"))
  }
  else if (length(focal) == 1) {
    if (is.numeric(focal)) {
      focal <- levels(treat)[focal]
    }
    if (is.character(focal)) {
      treat <- relevel(treat, focal)
    }
    else {
      stop("focal must be the name or index of the focal treatment group.", call. = FALSE)
    }
    treat.combinations <- lapply(levels(treat)[levels(treat) != focal], function(x) rev(c(focal, x)))
    pairwise <- TRUE
  }
  else stop("focal must be a vector of length 1 containing the name or index of the focal treatment group.", call. = FALSE)
  treat.names <- levels(treat)
  names(treat.names) <- treat.names
  
  if (is_not_null(cluster)) {
    stop("Clusters are not yet supported with multiple categorical treatments.", call. = FALSE)
  }
  else {
    #Setup output object
    out.names <- c("Pair.Balance", 
                   "Balance.Across.Pairs", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    addl.sds <- list()
    if (any(s.d.denom %in% c("pooled", "all"))) {
      if (any(s.d.denom == "pooled")) {
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
        addl.sds[["pooled"]] <- sqrt(rowMeans(do.call("cbind", lapply(levels(treat), function(t) col.w.v(C[treat == t, , drop = FALSE], 
                                                                                                         s.weights[treat == t])))))
      }
      if (any(s.d.denom == "all")) {
        C <- get.C(covs = covs, int = int, poly = poly, addl = addl, distance = distance, ...)
        addl.sds[["all"]] <- sqrt(col.w.v(C, s.weights))
      }
    }
    
    if (pairwise || is_not_null(focal)) {
      args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
      balance.tables <- lapply(treat.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick, addl.sds = addl.sds), args), quote = TRUE))
    }
    else {
      if (any(treat.names == "Others")) stop ("\"Others\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
      args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
      balance.tables <- lapply(treat.combinations, function(t) {
        treat_ <- factor(treat, levels = c(levels(treat), "Others"))
        treat_[treat_ != t[1]] <- "Others"
        treat_ <- factor(treat_, rev(t))
        do.call("base.bal.tab.binary", c(list(weights = weights, treat = treat_, distance = distance, subclass = subclass, covs = covs, call = NULL, int = int, poly = poly, addl = addl, continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick, addl.sds = addl.sds), args), quote = TRUE)
      })
    }
    for (i in seq_along(balance.tables)) {
      names(balance.tables)[i] <- paste(rev(treat.combinations[[i]]), collapse = " vs. ")
    }
    
    out[["Pair.Balance"]] <- balance.tables
    
    out[["Observations"]] <- samplesize.multi(balance.tables, treat.names, focal)
    
    if (multi.summary || !quick) {
      out[["Balance.Across.Pairs"]] <- balance.table.multi.summary(balance.tables, 
                                                                   weight.names = names(weights),
                                                                   m.threshold = m.threshold,
                                                                   v.threshold = v.threshold,
                                                                   ks.threshold = ks.threshold,
                                                                   no.adj = no.adj,
                                                                   quick = quick,
                                                                   types = NULL)
    }
    
    out[["call"]] <- call
    
    attr(out, "print.options") <- list(m.threshold=m.threshold,
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.adj=!no.adj, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks,
                                       disp.bal.tab = disp.bal.tab,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       treat.names = treat.names,
                                       which.treat = which.treat,
                                       multi.summary = multi.summary,
                                       pairwise = pairwise,
                                       co.names = attr(out[["Pair.Balance"]][[1]], "print.options")[["co.names"]])
    
    class(out) <- c("bal.tab.multi", "bal.tab")
  }
  return(out)
  
}
base.bal.tab.msm <- function(weights, treat.list, distance.list = NULL, subclass = NULL, covs.list, call = NULL, int = FALSE, poly = 1, addl.list = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), pairwise = TRUE, focal = NULL, which.treat = NA, multi.summary = getOption("cobalt_multi.summary", TRUE), which.time = NULL, msm.summary = getOption("cobalt_msm.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
  #One vector of weights
  #treat.list should be a df/list of treatment vectors, one for each time period
  #cov.list should be a list of covariate data.frames, one for each time period; 
  #   should include all covs from previous time points, but no treatment statuses
  
  #Preparations
  args <- list(...)
  
  if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
  if (is_not_null(v.threshold)) {
    v.threshold <- max(v.threshold, 1/v.threshold)
    disp.v.ratio <- TRUE
  }
  if (is_not_null(ks.threshold) && is_not_null(args$k.threshold)) {
    ks.threshold <- args$k.threshold
  }
  if (is_not_null(ks.threshold)) {
    if (ks.threshold > 1) {
      warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
      ks.threshold <- NULL
    }
    else disp.ks <- TRUE
  }
  if (is_not_null(r.threshold)) {
    r.threshold <- abs(r.threshold)
    if (r.threshold > 1) {
      warning("r.threshold must be between 0 and 1; ignoring r.threshold.", call. = FALSE)
      r.threshold <- NULL
    }
  }
  if (is_null(weights) && is_null(subclass)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat.list[[1]]))
  }
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
  
  if (nunique.gt(cluster, 1)) {
    stop("Clusters are not yet supported with longitudinal treatments.", call. = FALSE)
  }
  else {
    #Setup output object
    out.names <- c("Time.Balance", 
                   "Balance.Across.Times", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    out[["Time.Balance"]] <- vector("list", length(covs.list))
    
    treat.type <- vapply(treat.list, function(x) get.treat.type(assign.treat.type(x)), character(1L))
    
    #Get list of bal.tabs for each time period
    out[["Time.Balance"]] <- lapply(seq_along(treat.list), function(ti) {
      if (treat.type[ti] == "continuous") {
        args <- args[names(args) %nin% names(formals(base.bal.tab.cont))]
        out_ <- do.call("base.bal.tab.cont", c(list(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, poly = poly, addl = addl.list[[ti]], r.threshold = r.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick), args),
                        quote = TRUE)
      }
      else if (treat.type[ti] == "multinomial") {
        args <- args[names(args) %nin% names(formals(base.bal.tab.multi))]
        out_ <- do.call("base.bal.tab.multi", c(list(weights = weights, treat = treat.list[[ti]], distance=distance.list[[ti]], 
                                                     covs=covs.list[[ti]], call=NULL, int=int, addl=addl.list[[ti]], 
                                                     continuous=continuous, binary=binary, s.d.denom=s.d.denom, 
                                                     m.threshold=m.threshold, v.threshold=v.threshold, 
                                                     ks.threshold=ks.threshold, 
                                                     imbalanced.only = imbalanced.only,
                                                     un=un, 
                                                     disp.bal.tab = disp.bal.tab, 
                                                     disp.means=disp.means,
                                                     disp.sds=disp.sds,
                                                     disp.v.ratio=disp.v.ratio, 
                                                     disp.ks=disp.ks, 
                                                     method=method, 
                                                     cluster = cluster, which.cluster = which.cluster, 
                                                     cluster.summary = cluster.summary, pairwise = pairwise, focal = NULL,
                                                     which.treat = which.treat, multi.summary = multi.summary,
                                                     s.weights = s.weights, quick = quick),
                                                args),
                        quote = TRUE)
      }
      else if (treat.type[ti] == "binary") {
        args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
        out_ <- do.call("base.bal.tab.binary", c(list(weights = weights, treat = treat.list[[ti]], distance = distance.list[[ti]], subclass = NULL, covs = covs.list[[ti]], call = NULL, int = int, poly = poly, addl = addl.list[[ti]], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = FALSE, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster, which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights, discarded = discarded, quick = quick), args),
                        quote = TRUE)
      }
      else stop("Each treatment must be binary, multinomial, or continuous.", call. = FALSE)
      
      return(out_)
    })
    
    if (length(names(treat.list)) == length(treat.list)) {
      names(out[["Time.Balance"]]) <- names(treat.list)
    }
    else names(out[["Time.Balance"]]) <- seq_along(treat.list)
    
    out[["Observations"]] <- lapply(out[["Time.Balance"]], function(x) x$Observations)
    
    if (!(quick && !msm.summary) && all_the_same(treat.type) && !any(treat.type == "multinomial")) {
      out[["Balance.Across.Times"]] <- balance.table.msm.summary(out[["Time.Balance"]],
                                                                 weight.names = names(weights),
                                                                 no.adj = no.adj,
                                                                 m.threshold = m.threshold, 
                                                                 v.threshold = v.threshold, 
                                                                 ks.threshold = ks.threshold, 
                                                                 r.threshold = r.threshold, 
                                                                 quick = quick, 
                                                                 types = NULL)
    }
    
    out[["call"]] <- call
    
    attr(out, "print.options") <- list(m.threshold=m.threshold, 
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       r.threshold = r.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.means=disp.means, 
                                       disp.sds=disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks, 
                                       disp.adj=!no.adj, 
                                       disp.bal.tab = disp.bal.tab, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       which.time = which.time,
                                       msm.summary = msm.summary,
                                       co.names = attr(out[["Time.Balance"]][[1]], "print.options")[["co.names"]])
    
    class(out) <- c("bal.tab.msm", "bal.tab")
  }
  
  return(out)
}
base.bal.tab.target <- function(weights, treat, distance = NULL, subclass = NULL, covs, call = NULL, int = FALSE, poly = 1, addl = NULL, continuous = getOption("cobalt_continuous", "std"), binary = getOption("cobalt_binary", "raw"), m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, imbalanced.only = getOption("cobalt_imbalanced.only", FALSE), un = getOption("cobalt_un", FALSE), disp.means = getOption("cobalt_disp.means", FALSE), disp.sds = getOption("cobalt_disp.sds", FALSE), disp.v.ratio = getOption("cobalt_disp.v.ratio", FALSE), disp.ks = getOption("cobalt_disp.ks", FALSE), disp.subclass = getOption("cobalt_disp.subclass", FALSE), disp.bal.tab = getOption("cobalt_disp.bal.tab", TRUE), method, cluster = NULL, which.cluster = NULL, cluster.summary = getOption("cobalt_cluster.summary", TRUE), cluster.fun = getOption("cobalt_cluster.fun", NULL), which.treat = NA, target.summary = getOption("cobalt_target.summary", TRUE), s.weights = NULL, discarded = NULL, abs = FALSE, quick = TRUE, ...) {
  #Preparations
  args <- list(...)
  
  if (is_not_null(m.threshold)) m.threshold <- abs(m.threshold)
  if (is_not_null(v.threshold)) {
    v.threshold <- max(v.threshold, 1/v.threshold)
    disp.v.ratio <- TRUE
  }
  if (is_not_null(ks.threshold)) {
    if (ks.threshold > 1) {
      warning("ks.threshold must be between 0 and 1; ignoring ks.threshold.", call. = FALSE)
      ks.threshold <- NULL
    }
    else disp.ks <- TRUE
  }
  if (is_null(weights) && is_null(subclass)) {
    un <- TRUE
    no.adj <- TRUE
  }
  else {
    no.adj <- FALSE
    if (is_not_null(weights) && ncol(weights) == 1) names(weights) <- "Adj"
  }
  if (is_null(s.weights)) {
    s.weights <- rep(1, length(treat))
  }
  if (is_not_null(args[["agg.fun"]])) cluster.fun <- args[["agg.fun"]]
  
  #Create new Target group
  target.name <- "Target"
  n <- length(treat)
  
  if (isTRUE(get.treat.type(treat) == "continuous") || (is.numeric(treat) && !is_binary(treat))) {#if continuous treatment
    covs <- data.frame(treat = treat, covs)
    treat <- factor(rep(c("All", target.name), each = n))
    target.summary <- FALSE
    which.treat <- NULL
    needs.summary <- FALSE
    treat.names <- unique.treat <- "All"
  }
  else {
    if (is.factor(treat) || is.character(treat)) {
      if (is.factor(treat)) treat.names <- unique.treat <- levels(treat)
      else treat.names <- unique.treat <- unique(treat, nmax = n - 1)
    }
    else {
      treat.names <- c("Control", "Treated")
      unique.treat <- sort(unique(treat, nmax = 2))
    }
    names(treat.names) <- unique.treat
    
    treat <- factor(c(treat.names[as.character(treat)], rep(target.name, n)))
    needs.summary <- TRUE
  }
  
  covs <- rbind(covs, covs)
  if (is_not_null(weights)) weights <- rbind(weights, as.data.frame(array(1, dim = dim(weights), 
                                                                          dimnames = dimnames(weights))))
  distance <- rbind(distance, distance)
  addl <- rbind(addl, addl)
  s.weights <- c(s.weights, s.weights)
  if (is_not_null(discarded)) discarded <- c(discarded, rep(FALSE, length(discarded)))
  s.d.denom <- "treated"
  
  treat.target.combinations <- lapply(treat.names, function(x) c(x, target.name))
  
  if (is_not_null(cluster)) {
    stop("Clusters are not yet supported with target balance assessment.", call. = FALSE)
  }
  else if (is_not_null(subclass)) {
    stop("Subclassification is not yet supported with target balance assessment.", call. = FALSE)
  }
  else {
    #Setup output object
    out.names <- c("Target.Balance", 
                   "Balance.Across.Treatments", 
                   "Observations", 
                   "call")
    out <- vector("list", length(out.names))
    names(out) <- out.names
    
    
    if (any(treat.names == "Target")) stop ("\"Target\" cannot be the name of a treatment level. Please rename your treatments.", call. = FALSE)
    args <- args[names(args) %nin% names(formals(base.bal.tab.binary))]
    balance.tables <- lapply(treat.target.combinations, function(t) do.call("base.bal.tab.binary", c(list(weights = weights[treat %in% t, , drop = FALSE], treat = factor(treat[treat %in% t], t), distance = distance[treat %in% t, , drop = FALSE], subclass = subclass[treat %in% t], covs = covs[treat %in% t, , drop = FALSE], call = NULL, int = int, poly = poly, addl = addl[treat %in% t, , drop = FALSE], continuous = continuous, binary = binary, s.d.denom = s.d.denom, m.threshold = m.threshold, v.threshold = v.threshold, ks.threshold = ks.threshold, imbalanced.only = imbalanced.only, un = un, disp.means = disp.means, disp.sds = disp.sds, disp.v.ratio = disp.v.ratio, disp.ks = disp.ks, disp.subclass = disp.subclass, disp.bal.tab = disp.bal.tab, method = method, cluster = cluster[treat %in% t], which.cluster = which.cluster, cluster.summary = cluster.summary, s.weights = s.weights[treat %in% t], discarded = discarded[treat %in% t], quick = quick), args), quote = TRUE))
    
    for (i in seq_along(balance.tables)) {
      names(balance.tables)[i] <- paste(treat.target.combinations[[i]], collapse = " vs. ")
      balance.tables[[i]][["Observations"]][[2]] <- NULL
    }
    
    out[["Target.Balance"]] <- balance.tables
    
    out[["Observations"]] <- samplesize.target(balance.tables, treat.names, target.name) 
    
    if (needs.summary && (target.summary || !quick)) {
      out[["Balance.Across.Treatments"]] <- balance.table.target.summary(balance.tables, 
                                                                         weight.names = names(weights),
                                                                         m.threshold = m.threshold,
                                                                         v.threshold = v.threshold,
                                                                         ks.threshold = ks.threshold,
                                                                         no.adj = no.adj,
                                                                         quick = quick,
                                                                         types = NULL)
    }
    
    out[["call"]] <- call
    
    attr(out, "print.options") <- list(m.threshold=m.threshold,
                                       v.threshold=v.threshold,
                                       ks.threshold=ks.threshold,
                                       imbalanced.only = imbalanced.only,
                                       un=un, 
                                       disp.adj=!no.adj, 
                                       which.cluster=which.cluster,
                                       cluster.summary=cluster.summary,
                                       cluster.fun = cluster.fun,
                                       abs = abs,
                                       continuous = continuous,
                                       binary = binary,
                                       quick = quick,
                                       disp.means=disp.means, 
                                       disp.sds = disp.sds,
                                       disp.v.ratio=disp.v.ratio, 
                                       disp.ks=disp.ks,
                                       disp.bal.tab = disp.bal.tab,
                                       nweights = ifelse(no.adj, 0, ncol(weights)),
                                       weight.names = names(weights),
                                       treat.names = treat.names,
                                       target.name = target.name,
                                       which.treat = which.treat,
                                       target.summary = target.summary,
                                       co.names = attr(out[["Target.Balance"]][[1]], "print.options")[["co.names"]])
    
    class(out) <- c("bal.tab.target", "bal.tab")
  }
  return(out)
  
}

#bal.tab
is.designmatch <- function(x) {
  dm.b.names <- c("obj_total", "obj_dist_mat", "t_id", 
                  "c_id", "group_id", "time")
  dm.n.names <- c("obj_total", "obj_dist_mat", "id_1", 
                  "id_2", "group_id", "time")
  if (length(x) >= min(length(dm.b.names), length(dm.n.names)) && 
      (all(dm.b.names %in% names(x)) || all(dm.n.names %in% names(x)))) {
    class(x) <- c("designmatch")
  }
  return(x)
}
is.time.list <- function(x) {
  if (is.vector(x, mode = "list")) {
    if (all(vapply(x, is.formula, logical(1)))) {
      class(x) <- c("formula.list", "time.list", class(x))
    }
    else if (all(vapply(x, is.data.frame, logical(1)))) {
      class(x) <- c("data.frame.list", "time.list", class(x))
    }
  }
  return(x)
}

#x2base
match.strata2weights <- function(match.strata, treat, covs = NULL) {
  #Process match.strata into weights (similar to weight.subclass from MatchIt)
  if (is_null(covs)) names(treat) <- seq_along(treat)
  else names(treat) <- row.names(covs)
  matched <- !is.na(match.strata); unmatched <- !matched
  treat.matched <- treat[matched]
  match.strata.matched <- match.strata[matched]
  
  labels <- names(treat.matched)
  tlabels <- labels[treat.matched == 1]
  clabels <- labels[treat.matched == 0]
  
  weights.matched <- rep(0, length(treat.matched))
  names(weights.matched) <- labels
  weights.matched[tlabels] <- 1
  
  for (j in unique(match.strata.matched)){
    qn0 <- sum(treat.matched==0 & match.strata.matched==j)
    qn1 <- sum(treat.matched==1 & match.strata.matched==j)
    weights.matched[treat.matched==0 & match.strata.matched==j] <- qn1/qn0
  }
  if (all(check_if_zero(weights.matched[clabels]))) #all control weights are 0
    weights.matched[clabels] <- rep(0, length(weights.matched[clabels]))
  else {
    ## Number of C units that were matched to at least 1 T
    num.cs <- sum(!check_if_zero(weights.matched[clabels]))
    weights.matched[clabels] <- weights.matched[clabels]*num.cs/sum(weights.matched[clabels])
  }
  
  if (any(unmatched)) {
    weights.unmatched <- rep(0, sum(unmatched))
    names(weights.unmatched) <- names(treat[unmatched])
    weights <- c(weights.matched, weights.unmatched)[names(treat)]
  }
  else {
    weights <- weights.matched
  }
  
  if (all(check_if_zero(weights))) 
    stop("No units were matched", call. = FALSE)
  else if (all(check_if_zero(weights[tlabels])))
    stop("No treated units were matched", call. = FALSE)
  else if (all(check_if_zero(weights[clabels])))
    stop("No control units were matched", call. = FALSE)
  return(weights)
}
weights.same.as.strata <- function(weights, match.strata, treat) {
  weights == match.strata2weights(match.strata, treat)
}
use.tc.fd <- function(formula = NULL, data = NULL, treat = NULL, covs = NULL, needs.treat = TRUE, needs.covs = TRUE) {
  if (is_not_null(formula) && class(formula) == "formula") {
    D <- NULL
    if (is_not_null(data)) D <- data
    if (is_not_null(covs)) if (is_not_null(D)) D <- cbind(D, covs) else D <- covs
    t.c <- get.covs.and.treat.from.formula(formula, D, treat = treat)
    t.c <- list(treat = t.c[["treat"]], covs = t.c[["reported.covs"]], treat.name = t.c[["treat.name"]])
    attr(t.c, "which") <- "fd"
  }
  else {
    if (is.matrix(covs)) covs <- as.data.frame(covs)
    else if (!is.data.frame(covs)) stop("covs must be a data.frame of covariates.", call. = FALSE)
    if (!is.atomic(treat)) stop("treat must be an atomic vector of treatment statuses.", call. = FALSE)
    t.c <- list(treat = treat, covs = covs)
    attr(t.c, "which") <- "tc"
  }
  
  if (needs.covs && is_null(t.c[["covs"]])) stop("No covariates were specified.", call. = FALSE)
  if (needs.treat && is_null(t.c[["treat"]])) stop("No treatment variable was specified.", call. = FALSE)
  
  return(t.c)
}
process.val <- function(val, i, treat, covs, ...) {
  if (is.numeric(val)) {
    val.df <- setNames(data.frame(val), i)
  }
  else if (is.character(val)) {
    data.sets <- list(...)
    data.sets <- data.sets[!vapply(data.sets, is_null, logical(1))]
    if ((is_not_null(data.sets) && length(val) > max(vapply(data.sets, ncol, numeric(1)))) || length(val) == NROW(covs) || length(val) == length(treat)){
      val.df <- setNames(data.frame(val), i)
    }
    else {
      if (is_not_null(data.sets)) {
        val <- unique(val)
        val.df <- setNames(as.data.frame(matrix(NA, ncol = length(val), nrow = max(vapply(data.sets, nrow, numeric(1))))),
                           val)
        not.found <- setNames(rep(FALSE, length(val)), val)
        for (v in val) {
          found <- FALSE
          k <- 1
          while (found == FALSE && k <= length(data.sets)) {
            if (v %in% names(data.sets[[k]])) {
              val.df[[v]] <- data.sets[[k]][[v]]
              found <- TRUE
            }
            else k <- k + 1
          }
          if (!found) not.found[v] <- TRUE
        }
        if (any(not.found)) {
          warning(paste("The following variable(s) named in", i, "are not in any available data sets and will be ignored: ",
                        paste(val[not.found])), call. = FALSE)
          val.df <- val.df[!not.found]
        }
      }
      else {
        val.df <- NULL
        warning(paste0("Names were provided to ", i, ", but no argument to data was provided. Ignoring ", i,"."), 
                call. = FALSE)
      }
    }
  }
  else if (is.data.frame(val)) {
    val.df <- val
  }
  else stop(paste("The argument supplied to", i, "must be a vector, a data.frame, or the names of variables in an available data set."), call. = FALSE)
  
  return(val.df)
}
data.frame.process <- function(i, df, treat, covs, ...) {
  val <- df
  val.df <- NULL
  if (is_not_null(val)) {
    if (is.vector(val, mode = "list")) {
      val.list <- lapply(val, function(x) process.val(x, i, treat, covs, ...))
      if (is_null(names(val.list)) || "" %in% names(val.list)) {
        stop(paste("All entries in", i, "must have names."), call. = FALSE)
      }
      val.list <- lapply(seq_along(val.list), function(x) {
        if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
        return(val.list[[x]])})
      if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
        stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
      }
      
      val.df <- setNames(do.call("cbind", val.list),
                         c(sapply(val.list, names)))
    }
    else {
      val.df <- process.val(val, i, treat, covs, ...)
    }
    if (is_not_null(val.df)) { if (sum(is.na(val.df)) > 0) {
      stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
    }
  }
  return(val.df)
}
list.process <- function(i, List, ntimes, call.phrase, treat.list, covs.list, ...) {
  val.List <- List
  if (is_not_null(val.List)) {
    if (class(val.List)[1] != "list") {
      val.List <- list(val.List)
    }
    if (length(val.List) == 1) {
      val.List <- replicate(ntimes, val.List)
    }
    else if (length(val.List) == ntimes) {
      
    }
    else {
      stop(paste0("The argument to ", i, " must be a list of the same length as the number of time points in ",  call.phrase, "."), call. = FALSE)
    }
    for (ti in seq_along(val.List)) {
      val <- val.List[[ti]]
      val.df <- NULL
      if (is_not_null(val)) {
        if (is.vector(val, mode = "list")) {
          val.list <- lapply(val, function(x) process.val(x, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], ...))
          val.list <- lapply(seq_along(val.list), function(x) {
            if (NCOL(val.list[[x]]) == 1) names(val.list[[x]]) <- names(val.list)[x]
            val.list[[x]]})
          if (!all_the_same(vapply(val.list, nrow, numeric(1)))) {
            stop(paste("Not all items in", i, "have the same length."), call. = FALSE)
          }
          
          val.df <- setNames(do.call("cbind", val.list),
                             c(vapply(val.list, names, character(1))))
        }
        else {
          val.df <- process.val(val, strsplit(i, ".list", fixed = TRUE)[[1]], treat.list[[ti]], covs.list[[ti]], ...)
        }
        if (is_not_null(val.df)) { if (sum(is.na(val.df)) > 0) {
          stop(paste0("Missing values exist in ", i, "."), call. = FALSE)}
        }
        val.List[[ti]] <- val.df
      }
      
    }
    val.df.lengths <- vapply(val.List[lengths(val.List) > 0], nrow, numeric(1))
    if (max(val.df.lengths) != min(val.df.lengths)) {
      stop(paste("All columns in", i, "need to have the same number of rows."), call. = FALSE)
    }
  }
  return(val.List)
}
get.s.d.denom <- function(s.d.denom, estimand = NULL, weights = NULL, treat = NULL, focal = NULL, method = NULL) {
  check.estimand <- check.weights <- check.focal <- bad.s.d.denom <- bad.estimand <- FALSE
  s.d.denom.specified <- is_not_null(s.d.denom)
  estimand.specified <- is_not_null(estimand)
  
  if (s.d.denom.specified) {
    try.s.d.denom <- tryCatch(match_arg(s.d.denom, c("treated", "control", "pooled", "all"), several.ok = TRUE),
                              error = function(cond) FALSE)
    if (any(try.s.d.denom == FALSE)) {
      check.estimand <- TRUE
      bad.s.d.denom <- TRUE
    }
    else {
      if (length(try.s.d.denom) > 1 && length(try.s.d.denom) != ncol(weights)) {
        stop("s.d.denom must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
      }
      else s.d.denom <- try.s.d.denom
    }
  }
  else {
    check.estimand <- TRUE
  }
  
  if (check.estimand == TRUE) {
    if (estimand.specified) {
      try.estimand <- tryCatch(match_arg(tolower(estimand), c("att", "atc", "ate"), several.ok = TRUE),
                               error = function(cond) FALSE)
      if (any(try.estimand == FALSE)) {
        check.focal <- TRUE
        bad.estimand <- TRUE
      }
      else {
        if (length(try.estimand) > 1 && length(try.estimand) != ncol(weights)) {
          stop("estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
        }
        else s.d.denom <- vapply(try.estimand, switch, character(1L), att = "treated", atc = "control", ate = "pooled")
      }
    }
    else {
      check.focal <- TRUE
    }
  }
  if (check.focal == TRUE) {
    if (is_not_null(focal)) {
      s.d.denom <- "treated"
      estimand <- "att"
    }
    else check.weights <- TRUE
  }
  if (check.weights == TRUE) {
    if (is_null(weights)) {
      s.d.denom <- "pooled"
      estimand <- "ate"
    }
    else {
      s.d.denom <- estimand <- character(ncol(weights))
      for (i in seq_len(ncol(weights))) {
        if (is_null(method) || method[i] == "weighting") {
          if (is_binary(treat)) {
            if (all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])]) &&
                !all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])])
            ) { #if treated weights are the same and control weights differ; ATT
              estimand[i] <- "att"
              s.d.denom[i] <- "treated"
            }
            else if (all_the_same(weights[[i]][treat==0 & !check_if_zero(weights[[i]])]) &&
                     !all_the_same(weights[[i]][treat==1 & !check_if_zero(weights[[i]])])
            ) { #if control weights are the same and treated weights differ; ATC
              estimand[i] <- "atc"
              s.d.denom[i] <- "control"
            }
            else {
              estimand[i] <- "ate"
              s.d.denom[i] <- "pooled"
            }
          }
          else {
            if (length(focal) == 1) {
              estimand[i] <- "att"
              s.d.denom[i] <- "treated"
            }
            else {
              estimand[i] <- "ate"
              s.d.denom[i] <- "pooled"
            }
          }
        }
        else {
          estimand[i] <- "att"
          s.d.denom[i] <- "treated"
        }
      }
    }
  }
  if (is_not_null(weights) && length(s.d.denom) == 1) s.d.denom <- rep(s.d.denom, ncol(weights))
  
  if (s.d.denom.specified && bad.s.d.denom && (!estimand.specified || bad.estimand)) {
    message("Warning: s.d.denom should be one of \"treated\", \"control\", \"pooled\", or \"all\".\n         Using \"", word_list(s.d.denom), "\" instead.")
  }
  else if (estimand.specified && bad.estimand) {
    message("Warning: estimand should be one of \"ATT\", \"ATC\", or \"ATE\". Using \"", ifelse(all_the_same(estimand), toupper(estimand)[1], word_list(toupper(estimand))), "\" instead.")
  }
  else if (check.focal || check.weights) {
    message("Note: estimand and s.d.denom not specified; assuming ", ifelse(all_the_same(toupper(estimand)), toupper(unique(estimand)), word_list(toupper(estimand))), " and ", ifelse(all_the_same(s.d.denom), unique(s.d.denom), word_list(s.d.denom)), ".")
  }
  
  if (all(method %in% c("weighting", "matching"))) {
    if (is_not_null(weights) && length(s.d.denom) != ncol(weights)) {
      stop("Valid inputs to s.d.denom or estimand must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
  }
  return(s.d.denom)
}
get.X.class <- function(X) {
  if (is_not_null(X[["imp"]])) {
    if (is_not_null(X[["treat.list"]])) stop("Multiply imputed data is not yet supported with longitudinal treatments.", call. = FALSE)
    else if (is_(X[["treat"]], c("factor", "character")) && !is_binary(X[["treat"]])) stop("Multiply imputed data is not yet supported with multinomial treatments.", call. = FALSE)
    else X.class <- "imp"
  }
  else if (is_not_null(X[["treat.list"]])) X.class <- "msm"
  else if (is_binary(X[["treat"]])) X.class <- "binary"
  else if (is_(X[["treat"]], c("factor", "character"))) X.class <- "multi"
  else if (is.numeric(X[["treat"]])) X.class <- "cont"
  else probably.a.bug()
  
  return(X.class)
}
subset_X <- function(X, subset = NULL) {
  if (is_not_null(subset)) {
    if (!any(subset)) {
      stop("All subset set to FALSE.", call. = FALSE)
    }
    else if (!all(subset)) {
      n <- length(subset)
      subset_X_internal <- function(x, subset) {
        if (is_not_null(x)) {
          if (is.factor(x) && length(x) == n) factor(x[subset])
          else if (is.atomic(x) && length(x) == n) x[subset]
          else if ((is.matrix(x) || is.data.frame(x)) && nrow(x) == n) x[subset, , drop = FALSE]
          else if (is.list(x)) lapply(x, subset_X_internal, subset = subset)
          else x
        }
        else x
      }
      lapply(X["weights"], subset_X_internal, subset)
    }
    else X
  }
  else X
}
imp.complete <- function(data) {
  if (!is_(data, "mids")) stop("'data' not of class 'mids'")
  
  single.complete <- function (data, where, imp, ell) {
    if (is.null(where)) 
      where <- is.na(data)
    idx <- seq_len(ncol(data))[apply(where, 2, any)]
    for (j in idx) {
      if (is_null(imp[[j]])) data[where[, j], j] <- NA
      else data[where[, j], j] <- imp[[j]][, ell]
    }
    return(data)
  }
  
  m <- as.integer(data$m)
  idx <- 1L:m
  
  mylist <- lapply(idx, function(i) single.complete(data$data, data$where, data$imp, i))
  
  cmp <- data.frame(.imp = rep(idx, each = nrow(data$data)), 
                    .id = rep.int(1L:nrow(data$data), length(idx)), 
                    do.call("rbind", mylist))
  
  if (is.integer(attr(data$data, "row.names"))) 
    row.names(cmp) <- seq_len(nrow(cmp))
  else row.names(cmp) <- as.character(seq_len(nrow(cmp)))
  
  return(cmp)
}

#get.C
#Functions to turn input covariates into usable form
#int.poly.f creates interactions and polynomials
#splitfactor splits factor variable into indicators (now in utilities)
#binarize transforms 2-value variable into binary (0,1)
#get.C controls flow and handles redunancy
#get.types gets variables types (contin./binary)

int.poly.f <- function(mat, ex=NULL, int=FALSE, poly=1, center = FALSE, sep, co.names) {
  #Adds to data frame interactions and polynomial terms; interaction terms will be named "v1_v2" and polynomials will be named "v1_2"
  #Only to be used in base.bal.tab; for general use see int.poly()
  #mat=matrix input
  #ex=matrix of variables to exclude in interactions and polynomials; a subset of df
  #int=whether to include interactions or not; currently only 2-way are supported
  #poly=degree of polynomials to include; will also include all below poly. If 1, no polynomial will be included
  #nunder=number of underscores between variables
  
  if (is_not_null(ex)) d <- mat[, colnames(mat) %nin% colnames(ex), drop = FALSE]
  else d <- mat
  binary.vars <- apply(d, 2, is_binary)
  if (center) {
    d[,!binary.vars] <- center(d[, !binary.vars, drop = FALSE])
  }
  nd <- NCOL(d)
  nrd <- NROW(d)
  no.poly <- binary.vars
  npol <- nd - sum(no.poly)
  new <- matrix(0, ncol = (poly-1)*npol + int*(.5*(nd)*(nd-1)), nrow = nrd)
  nc <- NCOL(new)
  new.co.names <- vector("list", (nc))
  if (poly > 1 && npol != 0) {
    for (i in 2:poly) {
      new[, (1 + npol*(i - 2)):(npol*(i - 1))] <- apply(d[, !no.poly, drop = FALSE], 2, function(x) x^i)
      new.co.names[(1 + npol*(i - 2)):(npol*(i - 1))] <- lapply(colnames(d)[!no.poly], function(x) setNames(list(c(co.names[[x]][["component"]], num_to_superscript(i)), c(co.names[[x]][["type"]], "power")), c("component", "type")))
      
    }
  }
  if (int && nd > 1) {
    new[,(nc - .5*nd*(nd-1) + 1):nc] <- matrix(t(apply(d, 1, combn, 2, prod)), nrow = nrd)
    new.co.names[(nc - .5*nd*(nd-1) + 1):nc] <- lapply(as.data.frame(combn(colnames(d), 2), stringsAsFactors = FALSE), 
                                                       function(x) setNames(list(c(co.names[[x[1]]][["component"]], sep, co.names[[x[2]]][["component"]]),
                                                                                 c(co.names[[x[1]]][["type"]], "isep", co.names[[x[2]]][["type"]])),
                                                                            c("component", "type")))
  }
  
  colnames(new) <- vapply(new.co.names, function(x) paste0(x[["component"]], collapse = ""), character(1))
  names(new.co.names) <- colnames(new)
  new <- new[, !apply(new, 2, all_the_same), drop = FALSE]
  attr(new, "co.names") <- new.co.names
  return(new)
}
get.C <- function(covs, int = FALSE, poly = 1, addl = NULL, distance = NULL, cluster = NULL, ...) {
  #gets C data.frame, which contains all variables for which balance is to be assessed. Used in balance.table.
  A <- list(...)
  if (is_null(A[["int_sep"]])) A[["int_sep"]] <- getOption("cobalt_int_sep", default = " * ")
  if (is_null(A[["factor_sep"]])) A[["factor_sep"]] <- getOption("cobalt_factor_sep", default = "_")
  if (is_null(A[["center"]]) || A[["center"]] %nin% c(TRUE, FALSE)) A[["center"]] <- getOption("cobalt_center", default = FALSE)
  
  C <- covs
  if (!is.null(addl)) {
    if (!is.data.frame(addl)) {
      if (is.character(addl)) stop("The argument to addl must be a data.frame containing the the values of the additional variables you want to include in the balance assessment.", call. = FALSE)
      else stop("The argument to addl must be a data.frame. Wrap data.frame() around the argument if it is a matrix or vector.", call. = FALSE)
    }
    else {
      repeat.name.indices <- vapply(names(addl), function(x) x %in% names(C), logical(1))
      if (any(repeat.name.indices)) {
        warning(paste("The following variables in addl have the same name as covariates and will be ignored:\n",
                      paste(names(addl)[repeat.name.indices], collapse = " ")), call. = FALSE)
        addl <- addl[!repeat.name.indices]
      }
      C <- cbind(C, addl)
    }
  } 
  
  covs.with.inf <- vapply(C, function(x) !is.character(x) && any(!is.finite(x) & !is.na(x)), logical(1L))
  if (any(covs.with.inf)) {
    s <- if (sum(covs.with.inf) == 1) c("", "s") else c("s", "")
    stop(paste0("The variable", s[1], " ", word_list(names(C)[covs.with.inf], quotes = TRUE), 
                " contain", s[2], " non-finite values, which are not allowed."), call. = FALSE)
  }
  
  vars.w.missing <- data.frame(placed.after = names(C),
                               has.missing = FALSE, 
                               has.Inf = FALSE,
                               row.names = names(C),
                               stringsAsFactors = FALSE)
  co.names <- setNames(lapply(names(C), function(x) setNames(list(x, "base"), c("component", "type"))), names(C))
  #component types: base, fsep, isep, power, na, level
  is.0.1.cov <- setNames(rep(FALSE, ncol(C)), names(C))
  for (i in names(C)) {
    if (nunique(C[[i]]) == 2) {
      #if (is.logical(C[[i]])) C[[i]] <- as.numeric(C[[i]])
      if (((is.numeric(C[[i]]) || is.logical(C[[i]])) && 
           all(check_if_zero(as.numeric(C[[i]]) - binarize(C[[i]])), na.rm = TRUE)) ||
          all(as.character(C[[i]]) %in% c("0", "1"))) {
        is.0.1.cov[i] <- TRUE
      }
      
      C[[i]] <- factor(C[[i]])
      C[[i]] <- relevel(C[[i]], levels(C[[i]])[2])
    }
    else if (is.character(C[[i]]) || is.factor(C[[i]])) C[[i]] <- factor(C[[i]])
    if (is_not_null(cluster) && !nunique.gt(C[[i]], nunique(cluster)) && 
        equivalent.factors(C[[i]], cluster)) {
      C <- C[names(C) != i] #Remove variable if it is the same (linear combo) as cluster variable
    }
    else {
      if (anyNA(C[[i]])) vars.w.missing[i, "has.missing"] <- TRUE
      if (!is.numeric(C[[i]])) {
        old.C.names <- names(C)
        C <- splitfactor(C, i, replace = TRUE, sep = A[["factor_sep"]], drop.first = FALSE, 
                         drop.singleton = FALSE)
        newly.added.names <- names(C)[names(C) %nin% old.C.names]
        vars.w.missing[i, "placed.after"] <- newly.added.names[length(newly.added.names)]
        co.names <- c(co.names, setNames(lapply(newly.added.names, function(x) {
          split.points <- c(nchar(i), nchar(i) + nchar(A[["factor_sep"]]))
          split.names <- substring(x,
                                   c(1, split.points[1] + 1, split.points[2] + 1),
                                   c(split.points[1], split.points[2], nchar(x))
          )
          setNames(list(split.names, c("base", "fsep", "level")), 
                   c("component", "type"))
        }), newly.added.names))
      }
    }
  }
  
  if (NCOL(C) == 0) stop("There are no variables for which to display balance.", call. = FALSE)
  
  #Make sure categorical variable have missingness indicators done correctly
  
  C <- C[!vapply(C, all_the_same, logical(1L))]
  C <- as.matrix(C)
  
  #Process int and poly
  if (length(int) != 1L || !is.finite(int) || !(is.logical(int) || is.numeric(int))) {
    stop("int must be TRUE, FALSE, or a numeric value of length 1.", call. = FALSE)
  }
  if (int < 0 || !check_if_zero(abs(int - round(int)))) {
    stop("int must be TRUE, FALSE, or a numeric (integer) value greater than 1.", call. = FALSE)
  }
  int <- as.integer(round(int))
  
  if (length(poly) != 1L || !is.finite(poly) || !is.numeric(poly)) {
    stop("poly must be a numeric value of length 1.", call. = FALSE)
  }
  if (poly < 0 || !check_if_zero(abs(poly - round(poly)))) {
    stop("poly must be a numeric (integer) value greater than 1.", call. = FALSE)
  }
  poly <- as.integer(round(poly))
  
  if (int || poly) {
    if (int) { 
      #Prevent duplicate var names with `sep`s
      nsep <- 1
      repeat {
        all.possible.names <- outer(colnames(C), colnames(C), paste, sep = paste0(rep(A[["int_sep"]], nsep), collapse = ""))
        if (!any(colnames(C) %in% all.possible.names)) break
        else nsep <- nsep + 1
      }
      
      if (poly < int) poly <- int
      
      int <- TRUE
    }
    
    new <- int.poly.f(C, int = int, poly = poly, center = A[["center"]], 
                      sep = rep(A[["int_sep"]], nsep), co.names = co.names)
    C <- cbind(C, new)
    co.names <- c(co.names, attr(new, "co.names"))
  }
  
  #Add missingness indicators
  vars.w.missing <- vars.w.missing[vars.w.missing$placed.after %in% colnames(C) & vars.w.missing$has.missing, , drop = FALSE]
  if (NROW(vars.w.missing) > 0) {
    missing.ind <- apply(C[,colnames(C) %in% vars.w.missing$placed.after, drop = FALSE], 2, function(x) as.numeric(is.na(x)))
    colnames(missing.ind) <- rownames(vars.w.missing)
    vars.w.missing <- vars.w.missing[colnames(missing.ind), , drop = FALSE]
    colnames(missing.ind) <- paste0(colnames(missing.ind), ":<NA>")
    original.var.order <- setNames(seq_len(NCOL(C)), colnames(C))
    new.var.order <- original.var.order + cumsum(c(0,(colnames(C) %in% vars.w.missing$placed.after)[-NCOL(C)]))
    new.C <- matrix(NA, nrow = NROW(C), ncol = NCOL(C) + NCOL(missing.ind),
                    dimnames = list(rownames(C), seq_len(NCOL(C) + NCOL(missing.ind))))
    new.C[, new.var.order] <- C
    new.C[, -new.var.order] <- missing.ind
    colnames(new.C)[new.var.order] <- colnames(C)
    colnames(new.C)[-new.var.order] <- colnames(missing.ind)
    miss.co.names <- setNames(lapply(rownames(vars.w.missing), function(x) setNames(list(c(x, ":<NA>"),
                                                                                         c("base", "na")), c("component", "type"))),
                              colnames(missing.ind))
    C <- new.C
    missing.ind.C <- colnames(missing.ind)
    if (int) {
      new <- int.poly.f(missing.ind, int = TRUE, poly = 1, sep = rep(A[["int_sep"]], nsep), co.names = miss.co.names)
      C <- cbind(C, new)
      missing.ind.C <- c(missing.ind.C, colnames(new))
      miss.co.names <- c(miss.co.names, attr(new, "co.names"))
    }
    co.names <- c(co.names, miss.co.names)
    
  }
  
  #Remove duplicate & redundant variables
  C <- remove.perfect.col(C) 
  
  if (is_not_null(distance)) {
    if (any(names(distance) %in% colnames(C))) stop("distance variable(s) share the same name as a covariate. Please ensure each variable name is unique.", call. = FALSE)
    if (any(apply(distance, 2, function(x) anyNA(x)))) stop("Missing values are not allowed in the distance measure.", call. = FALSE)
    C <- cbind(as.matrix(distance), C, row.names = NULL)
    dist.co.names <- setNames(lapply(names(distance), function(x) setNames(list(x, "base"), c("component", "type"))), names(distance))
    co.names <- c(co.names, dist.co.names)
  }
  
  co.names <- co.names[colnames(C)]
  
  #Get rid of _1 for binary covs
  for (i in colnames(C)) {
    in.is.0.1.cov <- vapply(names(is.0.1.cov)[is.0.1.cov], 
                            function(i1) length(co.names[[i]][["component"]]) == 3 && 
                              co.names[[i]][["component"]][1] == i1 && 
                              co.names[[i]][["component"]][2] == A[["factor_sep"]] &&
                              co.names[[i]][["component"]][3] %in% c("1", "TRUE"), 
                            logical(1L))
    
    if (any(in.is.0.1.cov)) {
      name.index <- which(colnames(C) == i)
      new.name <- names(is.0.1.cov)[is.0.1.cov][in.is.0.1.cov][1]
      colnames(C)[name.index] <- new.name
      names(co.names)[name.index] <- new.name
      co.names[[name.index]][["component"]] <- new.name
      co.names[[name.index]][["type"]] <- "base"
    }
  }
  
  attr(co.names, "seps") <- c(factor = A[["factor_sep"]], int = A[["int_sep"]])
  attr(C, "co.names") <- co.names
  if (is_not_null(distance)) attr(C, "distance.names") <- names(distance)
  if (NROW(vars.w.missing) > 0) attr(C, "missing.ind") <- missing.ind.C
  
  return(C)
  
}
get.types <- function(C) {
  vapply(colnames(C), function(x) {
    if (any(attr(C, "distance.names") == x)) "Distance"
    else if (is_binary(C[,x]))  "Binary"
    else "Contin."
  }, character(1))
}
remove.perfect.col <- function(C) {
  C.no.miss <- C[,colnames(C) %nin% attr(C, "missing.ind"), drop = FALSE]
  #If many rows, select subset to test redundancy
  if (NROW(C.no.miss) > 1500) {
    repeat {
      mini.C.no.miss <- C.no.miss[sample(seq_len(NROW(C.no.miss)), 1000),,drop=FALSE]
      single.value <- apply(mini.C.no.miss, 2, all_the_same)
      if (all(!single.value)) break
    }
    suppressWarnings(C.cor <- cor(mini.C.no.miss, use = "pairwise.complete.obs"))
  }
  else suppressWarnings(C.cor <- cor(C.no.miss, use = "pairwise.complete.obs"))
  
  s <- !lower.tri(C.cor, diag=TRUE) & !is.na(C.cor) & check_if_zero(1 - abs(C.cor))
  redundant.vars <- colnames(C.no.miss)[apply(s, 2, any)]
  C <- C[, colnames(C) %nin% redundant.vars, drop = FALSE] 
  return(C)
}

#base.bal.tab
check_if_zero_weights <- function(weights.df, treat, unique.treat = NULL, treat.type = "cat") {
  #Checks if all weights are zero in each treat group for each set of weights
  if (treat.type == "cat") {
    if (is_null(unique.treat)) unique.treat <- unique(treat)
    w.t.mat <- expand.grid(colnames(weights.df), unique.treat, stringsAsFactors = FALSE)
    if (NROW(w.t.mat) > 0) {
      problems <- vapply(seq_len(NROW(w.t.mat)), function(x) all(check_if_zero(weights.df[treat == w.t.mat[x,2], w.t.mat[x, 1]])), logical(1L))
      if (any(problems)) {
        prob.w.t.mat <- droplevels(w.t.mat[problems,])
        if (NCOL(weights.df) == 1) {
          error <- paste0("All weights are zero when ", word_list(paste("treat =", prob.w.t.mat[, 2]), "or"), ".")
        }
        else {
          errors <- setNames(character(nlevels(prob.w.t.mat[,1])), levels(prob.w.t.mat[,1]))
          
          for (i in levels(prob.w.t.mat[,1])) {
            errors[i] <- paste0("\"", i, "\" weights are zero when ", word_list(paste("treat =", prob.w.t.mat[prob.w.t.mat[,1] == i, 2]), "or"))
          }
          errors <- paste(c("All", rep("all", length(errors)-1)), errors)
          error <- paste0(word_list(errors, "and"), ".")
        }
        stop(error, call. = FALSE)
      }
    }
  }
  else if (treat.type == "cont") {
    if (length(colnames(weights.df)) > 0) {
      problems <- vapply(colnames(weights.df), function(wn) all(check_if_zero(weights.df[, wn])), logical(1L))
      if (any(problems)) {
        prob.wts <- colnames(weights.df)[problems]
        if (NCOL(weights.df) == 1) {
          error <- paste0("All weights are zero.")
        }
        else {
          errors <- setNames(character(length(prob.wts)), prob.wts)
          
          for (i in prob.wts) {
            errors[i] <- paste0("\"", i, "\" weights are zero")
          }
          errors <- paste(c("All", rep("all", length(errors)-1)), errors)
          error <- paste0(word_list(errors, "and"), ".")
        }
        stop(error, call. = FALSE)
      }
    }
  }
  else stop("treat.type must be either \"cat\" or \"cont\".")
  
}
.col.std.diff <- function(mat, treat, weights, subclass = NULL, which.sub = NULL, x.types, continuous, binary, s.d.denom, no.weights = FALSE, s.weights = rep(1, length(treat)), pooled.sds = NULL) {
  if (no.weights) weights <- rep(1, NROW(mat))
  w <- weights*s.weights
  sw <- s.weights
  
  #Check continuous and binary
  if (missing(continuous) || is_null(continuous)) {
    continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
  }
  else continuous <- match_arg(continuous, c("std", "raw"))
  if (missing(binary) || is_null(binary)) {
    binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
  }
  else binary <- match_arg(binary, c("raw", "std"))
  
  no.sub <- is_null(which.sub)
  if (no.sub) ss <- sw > 0
  else {
    ss <- (!is.na(subclass) & subclass == which.sub & sw > 0)
    
    if (sum(treat==0 & ss) == 0) {
      warning(paste0("There are no control units in subclass ", which.sub, "."), call. = FALSE)
      return(rep(NA_real_, NCOL(mat)))
    }
    if (sum(treat==1 & ss) == 0) {
      warning(paste0("There are no treated units in subclass ", which.sub, "."), call. = FALSE)
      return(rep(NA_real_, NCOL(mat)))
    }
  }
  
  diffs <- col.w.m(mat[treat == 1 & ss, , drop = FALSE], w[treat == 1 & ss]) - 
    col.w.m(mat[treat == 0 & ss, , drop = FALSE], w[treat == 0 & ss])
  diffs[check_if_zero(diffs)] <- 0
  denoms <- rep(1, NCOL(mat))
  denoms.to.std <- ifelse(x.types == "Binary", binary == "std", continuous == "std")
  
  if (any(denoms.to.std)) {
    if (s.d.denom == "control") {
      denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]))
    }
    else if (s.d.denom == "treated") {
      denoms[denoms.to.std] <- sqrt(col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss]))
    }
    else if (s.d.denom == "pooled") {
      if (is_not_null(pooled.sds)) {
        denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
      }
      else {
        denoms[denoms.to.std] <-  sqrt(.5*(col.w.v(mat[treat == 0 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 0 & ss]) +
                                             col.w.v(mat[treat == 1 & ss, denoms.to.std, drop = FALSE], s.weights[treat == 1 & ss])))
      }
    }
  }
  
  std.diffs <- ifelse(check_if_zero(diffs), 0, diffs/denoms)
  if (any(!is.finite(std.diffs))) {
    warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
    std.diffs[!is.finite(std.diffs)] <- NA_real_
  }
  
  return(std.diffs)
}
std.diffs <- function(m0, s0, m1, s1, x.types, continuous, binary, s.d.denom, pooled.sds = NULL) {
  #Check continuous and binary
  if (missing(continuous) || is_null(continuous)) {
    continuous <- match_arg(getOption("cobalt_continuous", "std"), c("std", "raw"))
  }
  else continuous <- match_arg(continuous, c("std", "raw"))
  if (missing(binary) || is_null(binary)) {
    binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
  }
  else binary <- match_arg(binary, c("raw", "std"))
  
  diffs <- m1 - m0
  
  diffs[check_if_zero(diffs)] <- 0
  denoms <- rep(1, length(diffs))
  denoms.to.std <- ifelse(x.types == "Binary", binary == "std", continuous == "std")
  
  if (any(denoms.to.std)) {
    if (s.d.denom == "control") {
      denoms[denoms.to.std] <- s0[denoms.to.std]
    }
    else if (s.d.denom == "treated") {
      denoms[denoms.to.std] <- s1[denoms.to.std]
    }
    else if (s.d.denom == "pooled") {
      if (is_not_null(pooled.sds)) {
        denoms[denoms.to.std] <- pooled.sds[denoms.to.std]
      }
      else {
        denoms[denoms.to.std] <-  sqrt(.5*(s0[denoms.to.std]^2 + s1[denoms.to.std]^2))
      }
    }
  }
  
  std.diffs <- ifelse(check_if_zero(diffs), 0, diffs/denoms)
  if (any(!is.finite(std.diffs))) {
    warning("Some standardized mean differences were not finite. This can result from no variation in one of the treatment groups.", call. = FALSE)
    std.diffs[!is.finite(std.diffs)] <- NA_real_
  }
  
  return(std.diffs)
}
col.ks <- function(mat, treat, weights, x.types, no.weights = FALSE) {
  ks <- rep(NA_integer_, NCOL(mat))
  if (no.weights) weights <- rep(1, NROW(mat))
  weights_ <- weights
  weights_[treat == 1] <- weights[treat==1]/sum(weights[treat==1])
  weights_[treat == 0] <- -weights[treat==0]/sum(weights[treat==0])
  binary <- x.types == "Binary"
  if (any(!binary)) {
    ks[!binary] <- apply(mat[, !binary, drop = FALSE], 2, function(x_) {
      x <- x_[!is.na(x_)]
      ordered.index <- order(x)
      cumv <- abs(cumsum(weights_[ordered.index]))[diff(x[ordered.index]) != 0]
      return(if (is_null(cumv)) 0 else max(cumv))
    })
  }
  if (any(binary)) {
    ks[binary] <- abs(col.w.m(mat[treat == 1, binary, drop = FALSE], weights[treat == 1]) - 
                        col.w.m(mat[treat == 0, binary, drop = FALSE], weights[treat == 0]))
  }
  return(ks)
}
.col.var.ratio <- function(mat, treat, weights, x.types, no.weights = FALSE) {
  if (no.weights) weights <- rep(1, NROW(mat))
  ratios <- rep(NA_real_, NCOL(mat))
  non.binary <- x.types != "Binary"
  ratios[non.binary] <- col.w.v(mat[treat == 1, non.binary, drop = FALSE], w = weights[treat == 1]) / col.w.v(mat[treat == 0, non.binary, drop = FALSE], w = weights[treat == 0])
  return(pmax(ratios, 1/ratios))
}
var.ratios <- function(s0, s1, x.types) {
  ratios <- rep(NA_real_, length(s0))
  non.binary <- x.types != "Binary"
  ratios[non.binary] <- (s1[non.binary]^2) / (s0[non.binary]^2)
  return(ratios)
  # return(pmax(ratios, 1/ratios))
}
baltal <- function(threshold) {
  #threshold: vector of threshold values (i.e., "Balanced"/"Not Balanced")
  threshnames <- names(table(threshold))
  balstring <- threshnames[nchar(threshnames) > 0][1]
  thresh.val <- substring(balstring, 1 + regexpr("[><]", balstring), nchar(balstring))
  b <- data.frame(count=c(sum(threshold==paste0("Balanced, <", thresh.val)), 
                          sum(threshold==paste0("Not Balanced, >", thresh.val))))
  rownames(b) <- c(paste0("Balanced, <", thresh.val), paste0("Not Balanced, >", thresh.val))
  return(b)
}
samplesize <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL, treat.names = c("Control", "Treated")) {
  #Computes sample size info. for unadjusted and adjusted samples.
  # method is what method the weights are to be used for. 
  # method="subclassification" is for subclass sample sizes only.
  
  if (is_not_null(cluster) && is_not_null(which.cluster)) in.cluster <- cluster == which.cluster
  else in.cluster <- rep(TRUE, length(treat))
  if (is_null(s.weights)) s.weights <- rep(1, length(treat))
  if (is_null(discarded)) discarded <- rep(FALSE, length(treat))
  
  if (length(method) == 1 && method == "subclassification") {
    if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
    qbins <- nlevels(subclass)
    
    nn <- as.data.frame(matrix(0, 3, qbins))
    
    dimnames(nn) <- list(c(treat.names[1], treat.names[2], "Total"), 
                         paste("Subclass", levels(subclass)))
    
    matched <- !is.na(subclass)
    k <- 0
    for (i in levels(subclass)) {
      qi <- subclass[matched]==i
      qt <- treat[matched][qi]
      if (sum(qt==1)<2|(sum(qt==0)<2)){
        if (sum(qt==1)<2)
          warning("Not enough treatment units in subclass ", i, call. = FALSE)
        else if (sum(qt==0)<2)
          warning("Not enough control units in subclass ", i, call. = FALSE)
      }
      k <- k + 1
      nn[, k] <- c(sum(qt==0), sum(qt==1), length(qt))
    }
    attr(nn, "tag") <- "Sample sizes by subclass"
  }
  else if (is_null(weights)) {
    
    t <- treat[in.cluster]
    sw <- s.weights[in.cluster]
    
    nn <- as.data.frame(matrix(0, ncol = 2, nrow = 1))
    nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
    dimnames(nn) <- list(c("All"), 
                         c(treat.names[1], treat.names[2]))
    if (nunique.gt(s.weights, 2) || !any(s.weights==1) || any(s.weights %nin% c(0,1))) {
      attr(nn, "ss.type") <- c("ess")
    }
    else {
      attr(nn, "ss.type") <- c("ss")
    }
    
  }
  else if (NCOL(weights) == 1) {
    if (method=="matching") {
      nn <- as.data.frame(matrix(0, ncol=2, nrow=5))
      nn[1, ] <- c(sum(in.cluster & treat==0), 
                   sum(in.cluster & treat==1))
      nn[2, ] <- c(ESS(weights[in.cluster & treat==0, 1]), 
                   ESS(weights[in.cluster & treat==1, 1]))
      nn[3, ] <- c(sum(in.cluster & treat==0 & weights[,1]>0), 
                   sum(in.cluster & treat==1 & weights[,1]>0))
      nn[4, ] <- c(sum(in.cluster & treat==0 & weights[,1]==0 & discarded==0), 
                   sum(in.cluster & treat==1 & weights[,1]==0 & discarded==0))
      nn[5, ] <- c(sum(in.cluster & treat==0 & weights[,1]==0 & discarded==1), 
                   sum(in.cluster & treat==1 & weights[,1]==0 & discarded==1))
      dimnames(nn) <- list(c("All", "Matched (ESS)", "Matched (Unweighted)", "Unmatched", "Discarded"), 
                           c(treat.names[1], treat.names[2]))
      
      attr(nn, "ss.type") <- rep("ss", NROW(nn))
      
      if (!any(discarded)) {
        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
      }
    }
    else if (method == "weighting") {
      
      t <- treat[in.cluster]
      w <- weights[in.cluster, 1]
      sw <- s.weights[in.cluster]
      dc <- discarded[in.cluster]
      
      nn <- as.data.frame(matrix(0, ncol = 2, nrow = 3))
      nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
      nn[2, ] <- c(ESS(w[t==0]*sw[t==0]), ESS(w[t==1]*sw[t==1]))
      nn[3, ] <- c(sum(t==0 & dc==1), 
                   sum(t==1 & dc==1))
      dimnames(nn) <- list(c("Unadjusted", "Adjusted", "Discarded"), 
                           c(treat.names[1], treat.names[2]))
      attr(nn, "ss.type") <- c("ss", "ess", "ss")
      
      if (!any(dc)) {
        attr(nn, "ss.type") <- attr(nn, "ss.type")[rownames(nn) != "Discarded"]
        nn <- nn[rownames(nn) != "Discarded", ,drop = FALSE]
      }
    }
  }
  else {
    t <- treat[in.cluster]
    sw <- s.weights[in.cluster]
    
    nn <- as.data.frame(matrix(0, ncol=2, nrow=1+NCOL(weights)))
    nn[1, ] <- c(ESS(sw[t==0]), ESS(sw[t==1]))
    for (i in seq_len(NCOL(weights))) {
      if (method[i] == "matching") {
        w <- weights[in.cluster, i]
        nn[1+i,] <- c(ESS(w[t==0]), ESS(w[t==1]))
      }
      else if (method[i] == "weighting") {
        w <- weights[in.cluster, i]
        nn[1+i,] <- c(ESS(w[t==0]*sw[t==0]), ESS(w[t==1]*sw[t==1]))
      }
      
    }
    dimnames(nn) <- list(c("All", names(weights)), 
                         c(treat.names[1], treat.names[2]))
    attr(nn, "ss.type") <- c("ss", rep("ess", length(method)))
    
  }
  if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
    attr(nn, "tag") <- "Effective sample sizes"
  }
  else attr(nn, "tag") <- "Sample sizes"
  return(nn)
}
samplesize.across.clusters <- function(samplesize.list) {
  samplesize.list <- clear_null(samplesize.list)
  obs <- Reduce("+", samplesize.list)
  attr(obs, "tag") <- paste0("Total ", tolower(attr(samplesize.list[[1]], "tag")), " across clusters")
  return(obs)
}
max.imbal <- function(balance.table, col.name, thresh.col.name, ratio = FALSE) {
  balance.table.clean <- balance.table[balance.table$Type != "Distance" & is.finite(balance.table[, col.name]),]
  maxed <- balance.table.clean[which.max(abs_(balance.table.clean[, col.name], ratio = ratio)), match(c(col.name, thresh.col.name), names(balance.table.clean))]
  maxed <- data.frame(Variable = rownames(maxed), maxed)
  return(maxed)
}
balance.table <- function(C, weights, treat, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, un = FALSE, disp.means = FALSE, disp.sds = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, 
                          s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, types = NULL, addl.sds = NULL, disp.pop = FALSE, pop.means = NULL, pop.sds = NULL, quick = TRUE) {
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  if (no.adj) weight.names <- "Adj"
  else weight.names <- names(weights)
  names(s.d.denom) <- weight.names
  
  #B=Balance frame
  Bnames <- c("Type", 
              apply(expand.grid(c("M.0", "SD.0", "M.1", "SD.1", 
                                  # "M.Pop", "SD.Pop", 
                                  "Diff", "M.Threshold", 
                                  # "Diff.0.Pop", "Diff.1.Pop", 
                                  "V.Ratio", "V.Threshold", "KS", "KS.Threshold"),
                                c("Un", weight.names)), 1, paste, collapse = "."))
  B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
  colnames(B) <- Bnames
  rownames(B) <- colnames(C)
  
  #Set var type (binary/continuous)
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- get.types(C)
  bin.vars <- B[["Type"]] == "Binary"
  
  #Means for each group
  # if (!((!un || !disp.means) && quick)) {
  for (t in c(0, 1)) {
    B[[paste.("M", t, "Un")]] <- col_w_mean(C[treat == t, , drop = FALSE], weights = NULL,
                                            s.weights = s.weights[treat==t])
  }
  # if (!(!disp.pop && quick)) {
  #     B[["M.Pop.Un"]] <- col.w.m(C, w = s.weights)
  # }
  # }
  if (!no.adj) {
    # if (!no.adj && !(!disp.means && quick)) {
    for (i in weight.names) {
      for (t in c(0, 1)) {
        B[[paste.("M", t, i)]] <- col_w_mean(C[treat == t, , drop = FALSE], weights = weights[[i]][treat==t],
                                             s.weights = s.weights[treat==t])
      }
      # if (!(!disp.pop && quick)) {
      #     B[[paste.("M.Pop.Un", i)]] <- col.w.m(C, w = weights[[i]]*s.weights)
      # }
    }
  }
  
  #SDs for each group
  if (missing(binary) || is_null(binary)) {
    binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
  }
  else binary <- match_arg(binary, c("raw", "std"))
  
  sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
  # if (!((!un || !disp.sds) && quick)) {
  for (t in c(0, 1)) {
    sds <- rep(NA_real_, NCOL(C))
    if (any(sd.computable)) {
      sds[sd.computable] <- col_w_sd(C[treat == t, sd.computable, drop = FALSE],
                                     weights = NULL, s.weights = s.weights[treat==t],
                                     bin.vars = bin.vars[sd.computable])
    }
    B[[paste.("SD", t, "Un")]] <- sds
  }
  # if (!(!disp.pop && quick)) {
  #     sds <- rep(NA_real_, NCOL(C))
  #     sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = s.weights))
  #     B[["SD.Pop.Un"]] <- sds
  # }
  # }
  # if (!no.adj && !(!disp.sds && quick)) {
  if (!no.adj) {
    for (i in weight.names) {
      for (t in c(0, 1)) {
        sds <- rep(NA_real_, NCOL(C))
        if (any(sd.computable)) {
          sds[sd.computable] <- col_w_sd(C[treat == t, sd.computable, drop = FALSE],
                                         weights = weights[[i]][treat==t], s.weights = s.weights[treat==t],
                                         bin.vars = bin.vars[sd.computable])
        }
        B[[paste.("SD", t, i)]] <- sds
      }
      # if (!(!disp.pop && quick)) {
      #     sds <- rep(NA_real_, NCOL(C))
      #     sds[non.binary] <- sqrt(col.w.v(C[, non.binary, drop = FALSE], w = weights[[i]]*s.weights))
      #     B[[paste.("SD.Pop", i)]] <- sds
      # }
    }
  }
  if (!any(sapply(B[startsWith(names(B), "SD.")], is.finite))) {disp.sds <- FALSE}
  
  #Mean differences
  # if (!(!un && quick)) #Always compute unadjusted diffs
  B[["Diff.Un"]] <- col_w_smd(C, treat = treat, weights = NULL,
                              std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                              s.d.denom = if (is_null(addl.sds)) s.d.denom[1] else addl.sds[[s.d.denom[1]]],
                              abs = abs, s.weights = s.weights, bin.vars = bin.vars)
  
  if (!no.adj) {
    for (i in weight.names) {
      B[[paste.("Diff", i)]] <- col_w_smd(C, treat = treat, weights = weights[[i]],
                                          std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                          s.d.denom = if (is_null(addl.sds[[s.d.denom[i]]])) s.d.denom[i] 
                                          else addl.sds[[s.d.denom[i]]],
                                          abs = abs, s.weights = s.weights, bin.vars = bin.vars)
    }
  }
  
  #Variance ratios
  if (!(!disp.v.ratio && quick)) {
    if (!(!un && quick)) {
      vrs <- rep(NA_real_, NCOL(C))
      if (any(!bin.vars)) {
        vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = NULL, abs = abs, 
                                   s.weights = s.weights, bin.vars = bin.vars[!bin.vars])
      }
      B[["V.Ratio.Un"]] <- vrs
    }
    if (!no.adj) {
      for (i in weight.names) {
        vrs <- rep(NA_real_, NCOL(C))
        if (any(!bin.vars)) {
          vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = weights[[i]], abs = abs, 
                                     s.weights = s.weights, bin.vars = bin.vars[!bin.vars])
        }
        B[[paste.("V.Ratio", i)]] <- vrs
        
      }
    }
  }
  if (!any(sapply(B[startsWith(names(B), "V.Ratio.")], is.finite))) {disp.v.ratio <- FALSE; v.threshold <- NULL}
  
  #KS Statistics
  if (!(!disp.ks && quick)) {
    if (!(!un && quick)) {
      B[["KS.Un"]] <- col_w_ks(C, treat = treat, weights = NULL, s.weights = s.weights, bin.vars = bin.vars)
    }
    if (!no.adj) {
      for (i in weight.names) {
        B[[paste.("KS", i)]] <- col_w_ks(C, treat = treat, weights = weights[[i]], s.weights = s.weights, bin.vars = bin.vars)
      }
    }
  }
  if (!any(sapply(B[startsWith(names(B), "KS.")], is.finite))) {disp.ks <- FALSE; ks.threshold <- NULL}
  
  
  if (is_not_null(m.threshold)) {
    if (no.adj) {
      B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Diff.Un"]]), paste0(ifelse(abs_(B[["Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)), "")
      }
    }
    
  }
  if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
  
  if (is_not_null(v.threshold)) {
    if (no.adj) {
      B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["V.Ratio.Un"]]), paste0(ifelse(abs_(B[["V.Ratio.Un"]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("V.Ratio", i)]]), paste0(ifelse(abs_(B[[paste.("V.Ratio", i)]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
      }
    }
    
  }
  if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
  
  if (is_not_null(ks.threshold)) {
    if (no.adj) {
      B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["KS.Un"]]), paste0(ifelse(B[["KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("KS", i)]]), paste0(ifelse(B[[paste.("KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
      }
    }
    
  }
  if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
  
  attr(B, "thresholds") <- c(m = m.threshold,
                             v = v.threshold,
                             ks = ks.threshold)
  attr(B, "disp") <- c(means = disp.means,
                       sds = disp.sds,
                       v.ratio = disp.v.ratio,
                       ks = disp.ks)
  
  return(B)
}
balance.table.subclass <- function(C, weights = NULL, treat, subclass, continuous, binary, s.d.denom, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, disp.means = FALSE, disp.sds = FALSE, disp.v.ratio = FALSE, disp.ks = FALSE, s.weights = rep(1, length(treat)), types = NULL, abs = FALSE, quick = TRUE) {
  #Creates list SB of balance tables for each subclass
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  #B=Balance frame
  Bnames <- c("Type", "M.0.Adj", "SD.0.Adj", "M.1.Adj", "SD.1.Adj", "Diff.Adj", "M.Threshold", "V.Ratio.Adj", "V.Threshold", "KS.Adj", "KS.Threshold")
  B <- as.data.frame(matrix(NA_real_, nrow = NCOL(C), ncol = length(Bnames)))
  colnames(B) <- Bnames
  rownames(B) <- colnames(C)
  #Set var type (binary/continuous)
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- get.types(C)
  bin.vars <- B[["Type"]] == "Binary"
  
  SB <- vector("list", nlevels(subclass))
  names(SB) <- levels(subclass)
  
  if (missing(binary) || is_null(binary)) {
    binary <- match_arg(getOption("cobalt_binary", "raw"), c("raw", "std"))
  }
  else binary <- match_arg(binary, c("raw", "std"))
  
  #-------------------------------------
  for (i in levels(subclass)) {
    
    SB[[i]] <- B
    in.subclass <- !is.na(subclass) & subclass==i
    
    #Means for each group
    # if (!(!disp.means && quick)) {
    for (t in c(0,1)) {
      SB[[i]][[paste.("M", t, "Adj")]] <- col_w_mean(C, subset = treat==t & in.subclass)
    }
    # }
    
    #SDs for each group
    # if (!(!disp.sds && quick)) {
    sd.computable <- if (binary == "std") rep(TRUE, nrow(B)) else !bin.vars
    for (t in c(0, 1)) {
      sds <- rep(NA_real_, NCOL(C))
      sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE], subset = treat == t & in.subclass)
      SB[[i]][[paste.("SD", t, "Adj")]] <- sds
    }
    
    # }
    
    #Mean differences
    SB[[i]][["Diff.Adj"]] <- col_w_smd(C, treat = treat, weights = NULL,
                                       std = (bin.vars & binary == "std") | (!bin.vars & continuous == "std"),
                                       s.d.denom = s.d.denom, abs = FALSE, s.weights = NULL, 
                                       bin.vars = bin.vars, subset = in.subclass)
    
    #Variance ratios
    if (!(!disp.v.ratio && quick)) {
      vrs <- rep(NA_real_, NCOL(C))
      if (any(!bin.vars)) {
        vrs[!bin.vars] <- col_w_vr(C[, !bin.vars, drop = FALSE], treat, weights = NULL, abs = abs, 
                                   s.weights = NULL, bin.vars = bin.vars[!bin.vars],
                                   subset = in.subclass)
      }
      SB[[i]][["V.Ratio.Adj"]] <- vrs
    }
    
    #KS Statistics
    if (!(!disp.ks && quick)) {
      SB[[i]][["KS.Adj"]] <- col_w_ks(C, treat = treat, weights = NULL, s.weights = NULL, bin.vars = bin.vars,
                                      subset = in.subclass)
    }
  }
  
  if (all(sapply(SB, function(x) !any(is.finite(c(x[["SD.0.Adj"]], x[["SD.1.Adj"]])))))) {
    attr(SB, "dont.disp.sds") <- TRUE
    disp.sds <- FALSE
  }
  
  if (is_not_null(m.threshold)) {
    for (i in levels(subclass)) {
      SB[[i]][["M.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                         paste0(ifelse(is.finite(SB[[i]][["Diff.Adj"]]) & abs_(SB[[i]][["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), round(m.threshold, 3)))
    }
  }
  
  if (all(sapply(SB, function(x) !any(is.finite(x[["V.Ratio.Adj"]]))))) {
    attr(SB, "dont.disp.v.ratio") <- TRUE; v.threshold <- NULL
    disp.v.ratio <- FALSE
  }
  if (is_not_null(v.threshold)) {
    for (i in levels(subclass)) {
      SB[[i]][["V.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["V.Ratio.Adj"]]), 
                                         paste0(ifelse(abs_(SB[[i]][["V.Ratio.Adj"]], ratio = TRUE) < v.threshold, "Balanced, <", "Not Balanced, >"), round(v.threshold, 3)), "")
    }
  }
  if (all(sapply(SB, function(x) !any(is.finite(x[["KS.Adj"]]))))) {
    attr(SB, "dont.disp.ks") <- TRUE
    disp.ks <- FALSE
  }
  if (is_not_null(ks.threshold)) {
    for (i in levels(subclass)) {
      SB[[i]][["KS.Threshold"]] <- ifelse(SB[[i]][["Type"]]!="Distance" & is.finite(SB[[i]][["KS.Adj"]]), 
                                          paste0(ifelse(SB[[i]][["KS.Adj"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), round(ks.threshold, 3)), "")
    }
  }
  
  attr(SB, "thresholds") <- c(m = m.threshold,
                              v = v.threshold,
                              ks = ks.threshold)
  attr(SB, "disp") <- c(means = disp.means,
                        sds = disp.sds,
                        v.ratio = disp.v.ratio,
                        ks = disp.ks)
  
  return(SB)
}
balance.table.across.subclass <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, s.d.denom = NULL) {
  #Variance ratio, v.threshold, and KS not yet supported
  if (is_not_null(s.d.denom)){
    sub.by <- switch(s.d.denom, treated = "treat",
                     pooled = "all", control = "control")
  }
  if (sub.by=="treat") {
    wsub <- "Treated"
  } else if (sub.by=="control") {
    wsub <- "Control"
  } else if (sub.by=="all") {
    wsub <- "Total"
  }
  
  B.A <- balance.table.subclass.list[[1]][c("M.0.Adj", "M.1.Adj", "Diff.Adj")]
  
  for(i in rownames(B.A)) {
    for(j in colnames(B.A)) {
      B.A[[i, j]] <- sum(vapply(seq_along(balance.table.subclass.list),
                                function(s) subclass.obs[[wsub, s]]/sum(subclass.obs[wsub, ]) * (balance.table.subclass.list[[s]][[i, j]]), numeric(1)))
    }
  }
  B.A.df <- data.frame(balance.table[c("Type", "M.0.Un", "SD.0.Un", "M.1.Un", "SD.1.Un", "Diff.Un", "V.Ratio.Un", "KS.Un")], 
                       B.A, M.Threshold = NA_character_)
  if (is_not_null(m.threshold)) B.A.df[["M.Threshold"]] <- ifelse(B.A.df[["Type"]]=="Distance", "", paste0(ifelse(is.finite(B.A.df[["Diff.Adj"]]) & abs_(B.A.df[["Diff.Adj"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold))
  return(B.A.df)
}
balance.table.cluster.summary <- function(balance.table.clusters.list, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
  
  balance.table.clusters.list <- clear_null(balance.table.clusters.list)
  cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(balance.table.clusters.list, names)))
  if (no.adj) weight.names <- "Adj"
  
  Brownames <- unique(do.call("c", lapply(balance.table.clusters.list, rownames)))
  #cluster.functions <- c("Min", "Mean", "Median", "Max")
  cluster.functions <- c("Min", "Mean", "Max")
  stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
  Bcolnames <- c("Type", apply(expand.grid(cluster.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
  B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
  names(B) <- Bcolnames
  
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(balance.table.clusters.list, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
  
  abs0 <- function(x) {if (abs) abs(x) else (x)}
  funs <- vfuns <- structure(vector("list", length(cluster.functions)), names = cluster.functions)
  for (Fun in cluster.functions) {
    funs[[Fun]] <- function(x, ...) {
      if (!any(is.finite(x))) NA_real_
      else get(tolower(Fun))(x, ...)
    }
    vfuns[[Fun]] <- function(x, ...) {
      if (!any(is.finite(x))) NA_real_
      else if (Fun == "Mean") geom.mean(x, ...)
      else get(tolower(Fun))(x, ...)
    }
    for (sample in c("Un", weight.names)) {
      if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
        if (cont.treat) {
          B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs0(y[[x, paste.("Corr", sample)]])), na.rm = TRUE), numeric(1))
        }
        else {
          B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) abs0(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
          B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
          B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(balance.table.clusters.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
        }            
      }
    }
  }
  
  return(B)
}

#base.bal.tab.cont
samplesize.cont <- function(treat, weights = NULL, subclass = NULL, s.weights = NULL, method=c("matching", "weighting", "subclassification"), cluster = NULL, which.cluster = NULL, discarded = NULL) {
  #Computes sample size info. for unadjusted and adjusted samples.
  # method is what method the weights are to be used for. 
  # method="subclassification" is for subclass sample sizes only.
  #method <- match_arg(method)
  if (nlevels(cluster) > 0 && is_not_null(which.cluster)) in.cluster <- cluster == which.cluster
  else in.cluster <- rep(TRUE, length(treat))
  if (is_null(discarded)) discarded <- rep(0, length(treat))
  
  if (length(method) == 1 && method == "subclassification") {
    #stop("Subclassification is not yet surpported with continuous treatments.", call. = FALSE)
    if (is_null(subclass)) stop("subclass must be a vector of subclasses.")
    qbins <- nlevels(subclass)
    
    nn <- as.data.frame(matrix(0, nrow = 1, ncol = qbins))
    
    dimnames(nn) <- list(c("Total"), 
                         paste("Subclass", levels(subclass)))
    
    matched <- !is.na(subclass)
    k <- 0
    for (i in levels(subclass)) {
      qi <- subclass[matched]==i
      qt <- treat[matched][qi]
      if (length(qt)<2){
        if (sum(qt==1)<2)
          warning("Not enough units in subclass ", i, call. = FALSE)
      }
      k <- k + 1
      nn[, k] <- c(length(qt))
    }
    attr(nn, "tag") <- "Sample sizes by subclass"
  }
  else if (is_null(weights)) {
    nn <- as.data.frame(matrix(0, ncol = 1, nrow = 1))
    if (nunique.gt(s.weights, 2) || !any(s.weights==1) || !all(s.weights %in% c(0,1))) {
      sw <- s.weights[in.cluster]
      
      nn[1, ] <- ESS(sw)
    }
    else {
      nn[1, ] <- sum(in.cluster)
      
    }
    dimnames(nn) <- list(c("All"), 
                         c("Total"))
    attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
  }
  else if (length(weights) == 1) {
    if (method=="matching") {
      
      nn <- as.data.frame(matrix(0, ncol = 1, nrow = 3))
      nn[1, ] <- c(length(treat[in.cluster]))
      nn[2, ] <- c(sum(in.cluster & weights[,1] > 0))
      nn[3, ] <- c(sum(in.cluster & weights[,1] == 0))
      dimnames(nn) <- list(c("All", "Matched", "Unmatched"), 
                           c("Total"))
      attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
      
      #attr(nn, "tag") <- "Sample sizes"
    }
    else if (method == "weighting") {
      w <- weights[in.cluster, 1]
      sw <- s.weights[in.cluster]
      
      nn <- as.data.frame(matrix(0, ncol = 1, nrow = 2))
      nn[1, ] <- ESS(sw)
      nn[2, ] <- ESS(w*sw)
      dimnames(nn) <- list(c("Unadjusted", "Adjusted"), 
                           c("Total"))
      attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
      #attr(nn, "tag") <- "Effective sample sizes"
    }
  }
  else {
    #t <- treat[in.cluster]
    sw <- s.weights[in.cluster]
    nn <- as.data.frame(matrix(0, ncol=1, nrow=1+NCOL(weights)))
    nn[1, ] <- ESS(sw)
    for (i in seq_len(NCOL(weights))) {
      if (method[i] == "matching") {
        nn[1+i,] <- c(sum(in.cluster & weights[,i] > 0))
      }
      else if (method[i] == "weighting") {
        w <- weights[in.cluster, i]
        nn[1+i,] <- ESS(w*sw)
      }
      
    }
    dimnames(nn) <- list(c("Unadjusted", names(weights)), 
                         c("Total"))
    attr(nn, "ss.type") <- c("ss", ifelse(method == "weighting", "ess", "ss"))
    # if (all(obs$ss.type == "ess")) attr(obs, "tag") <- "Effective sample sizes"
    # else attr(obs, "tag") <- "Sample sizes"
    
  }
  if (length(attr(nn, "ss.type")) > 1 && all(attr(nn, "ss.type")[-1] == "ess")) {
    attr(nn, "tag") <- "Effective sample sizes"
  }
  else attr(nn, "tag") <- "Sample sizes"
  
  return(nn)
}
balance.table.cont <- function(C, weights, treat, r.threshold = NULL, un = FALSE, disp.means = FALSE, disp.sds = FALSE, s.weights = rep(1, length(treat)), abs = FALSE, no.adj = FALSE, types = NULL, target.means = NULL, target.sds = NULL, quick = TRUE) {
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  if (no.adj) weight.names <- "Adj"
  else weight.names <- names(weights)
  
  #B=Balance frame
  Bnames <- c("Type", 
              apply(expand.grid(c("M", "SD", "Corr", "R.Threshold"),
                                c("Un", weight.names)), 1, paste, collapse = "."))
  B <- as.data.frame(matrix(nrow = NCOL(C), ncol = length(Bnames)))
  colnames(B) <- Bnames
  rownames(B) <- colnames(C)
  
  #Set var type (binary/continuous)
  if (is_not_null(types)) B[,"Type"] <- types
  else B[["Type"]] <- get.types(C)
  bin.vars <- B[["Type"]] == "Binary"
  
  #Means
  if (!((!un || !disp.means) && quick)) {
    B[["M.Un"]] <- col_w_mean(C, weights = NULL, s.weights = s.weights)
  }
  if (!no.adj && !(!disp.means && quick)) {
    for (i in weight.names) {
      B[[paste.("M", i)]] <- col_w_mean(C, weights = weights[[i]], s.weights = s.weights)
    }
  }
  
  #SDs
  sd.computable <- !bin.vars
  if (!((!un || !disp.sds) && quick)) {
    sds <- rep(NA_real_, NCOL(C))
    if (any(sd.computable)) {
      sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                     weights = NULL, s.weights = s.weights,
                                     bin.vars = bin.vars[sd.computable])
    }
    B[["SD.Un"]] <- sds
  }
  if (!no.adj && !(!disp.sds && quick)) {
    for (i in weight.names) {
      sds <- rep(NA_real_, NCOL(C))
      if (any(sd.computable)) {
        sds[sd.computable] <- col_w_sd(C[, sd.computable, drop = FALSE],
                                       weights = weights[[i]], s.weights = s.weights,
                                       bin.vars = bin.vars[sd.computable])
      }
      B[[paste.("SD", i)]] <- sds
    }
  }
  if (!any(sapply(B[startsWith(names(B), "SD.")], is.finite))) {disp.sds <- FALSE}
  
  #Correlations
  # if (!(!un && quick)) #Always calculate unadjusted corrs
  B[["Corr.Un"]] <- col_w_corr(C, treat, weights = NULL, abs = abs, s.weights = s.weights, 
                               bin.vars = bin.vars, na.rm = TRUE)
  if (!no.adj) {
    for (i in weight.names) {
      B[[paste.("Corr", i)]]  <- col_w_corr(C, treat, weights = weights[[i]], abs = abs, s.weights = s.weights, 
                                            bin.vars = bin.vars, na.rm = TRUE)
    }
  }
  
  if (is_not_null(r.threshold)) {
    if (no.adj) {
      B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Corr.Un"]]), paste0(ifelse(abs_(B[["Corr.Un"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("R.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Corr", i)]]), paste0(ifelse(abs_(B[[paste.("Corr", i)]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)), "")
      }
    }
  }
  if (no.adj || NCOL(weights) <= 1) names(B)[names(B) == "R.Threshold.Adj"] <- "R.Threshold"
  
  
  attr(B, "thresholds") <- c(r = r.threshold)
  attr(B, "disp") <- c(means = disp.means,
                       sds = disp.sds)
  return(B)
  
}
balance.table.subclass.cont <- function(C, weights = NULL, treat, subclass, r.threshold = NULL, disp.means = FALSE, disp.sds = FALSE, s.weights = rep(1, length(treat)), types = NULL, quick = TRUE) {
  #Creates list SB of balance tables for each subclass
  #C=frame of variables, including distance; distance name (if any) stores in attr(C, "distance.name")
  
  #B=Balance frame
  Bnames <- c("Type", "M.Adj", "SD.Adj", "Corr.Adj", "R.Threshold")
  B <- as.data.frame(matrix(nrow=NCOL(C), ncol=length(Bnames)))
  colnames(B) <- Bnames
  rownames(B) <- colnames(C)
  #Set var type (binary/continuous)
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- get.types(C)
  
  SB <- vector("list", nlevels(subclass))
  names(SB) <- levels(subclass)
  
  #-------------------------------------
  for (i in levels(subclass)) {
    
    SB[[i]] <- B
    in.subclass <- !is.na(subclass) & subclass==i
    
    if (!(!disp.means && quick)) {
      SB[[i]][["M.Adj"]] <- colMeans(C[in.subclass, , drop = FALSE])
    }
    if (!(!disp.sds && quick)) {
      non.binary <- B[["Type"]] != "Binary"
      sds <- rep(NA_real_, NCOL(C))
      sds[non.binary] <- apply(C[in.subclass, non.binary, drop = FALSE], 2, sd)
      SB[[i]][["SD.Adj"]] <- sds
    }
    
    #Correlations
    # SB[[i]][["Corr.Adj"]] <- apply(C, 2, function(x) w.r(x[in.subclass], y = treat[in.subclass]))
    SB[[i]][["Corr.Adj"]] <- col.w.r(C[in.subclass,], y = treat[in.subclass])
    
  }
  
  if (all(sapply(SB, function(x) !any(is.finite(x[["SD.Adj"]]))))) {
    attr(SB, "dont.disp.sds") <- TRUE
    disp.sds <- FALSE
  }
  
  if (is_not_null(r.threshold)) {
    for (i in levels(subclass)) {
      SB[[i]][["R.Threshold"]] <- ifelse(SB[[i]][["Type"]]=="Distance", "", 
                                         paste0(ifelse(is.finite(SB[[i]][["Corr.Adj"]]) & abs_(SB[[i]][["Corr.Adj"]]) < r.threshold, "Balanced, <", "Not Balanced, >"), round(r.threshold, 3)))
    }
  }
  
  attr(SB, "thresholds") <- c(r = r.threshold)
  attr(SB, "disp") <- c(means = disp.means,
                        sds = disp.sds)
  
  return(SB)
}
balance.table.across.subclass.cont <- function(balance.table, balance.table.subclass.list, subclass.obs, sub.by = NULL, r.threshold = NULL) {
  #Not specified
}

#base.bal.tab.imp
balance.table.imp.summary <- function(bal.tab.imp.list, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
  if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.imp.list, class)))) {
    bal.tab.imp.list <- lapply(bal.tab.imp.list, function(x) x[["Balance"]])}
  cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(bal.tab.imp.list, names)))
  if (length(weight.names) <= 1) weight.names <- "Adj"
  bal.tab.imp.list <- clear_null(bal.tab.imp.list)
  
  Brownames <- unique(do.call("c", lapply(bal.tab.imp.list, rownames)))
  #imp.functions <- c("Min", "Mean", "Median", "Max")
  imp.functions <- c("Min", "Mean", "Max")
  stats <- if (cont.treat) "Corr" else c("Diff", "V.Ratio", "KS")
  Bcolnames <- c("Type", apply(expand.grid(imp.functions, stats, c("Un", weight.names)), 1, paste, collapse = "."))
  B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
  names(B) <- Bcolnames
  
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(sapply(bal.tab.imp.list, function(y) y[[x, "Type"]])); return(u[!is.na(u)])}), use.names = FALSE)
  
  abs0 <- function(x) {if (is_null(x)) NA_real_ else if (abs) abs(x) else (x)}
  funs <- vfuns <- structure(vector("list", length(imp.functions)), names = imp.functions)
  for (Fun in imp.functions) {
    funs[[Fun]] <- function(x, ...) {
      if (!any(is.finite(x))) NA_real_
      else get(tolower(Fun))(x, ...)
    }
    vfuns[[Fun]] <- function(x, ...) {
      if (!any(is.finite(x))) NA_real_
      else if (Fun == "Mean") geom.mean(x, ...)
      else get(tolower(Fun))(x, ...)
    }
    for (sample in c("Un", weight.names)) {
      if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
        if (cont.treat) {
          B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs0(y[x, paste.("Corr", sample)])), na.rm = TRUE), numeric(1))
        }
        else {
          B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) abs0(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
          B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
          B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(bal.tab.imp.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
        }
      }
    }
  }
  return(B)
}
balance.table.clust.imp.summary <- function(summary.tables, weight.names = NULL, no.adj = FALSE, abs = FALSE, quick = TRUE, types = NULL) {
  #cont.treat <- !is.na(match("bal.tab.cont", unique(do.call("c", lapply(bal.tab.imp.list, class)))))
  #clusters <- unique(do.call("c", lapply(bal.tab.imp.list, function(x) names(x[["Cluster.Balance"]]))))
  #cluster.tables <- lapply(clusters, function(x) lapply(bal.tab.imp.list, function(y) y[["Cluster.Balance"]][[x]]))
  #cluster.balance.across.imps <- lapply(cluster.tables, balance.table.imp.summary, no.adj, quick, types)
  #names(cluster.balance.across.imps) <- clusters
  
  if (!all(vapply(summary.tables, is_null, logical(1)))) {
    Brownames <- unique(do.call("c", lapply(summary.tables, rownames)))
    Bcolnames <- unique(do.call("c", lapply(summary.tables, colnames)))
    cont.treat <- !is.na(charmatch("Mean.Corr.Un", Bcolnames))
    if (length(weight.names) <= 1) weight.names <- "Adj"
    #imp.functions <- c("Min", "Mean", "Median", "Max")
    imp.functions <- c("Min", "Mean", "Max")
    
    B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)))
    dimnames(B) <- list(Brownames, Bcolnames)
    
    if (is_not_null(types)) B[["Type"]] <- types
    else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(summary.tables, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
    
    abs0 <- function(x) {if (abs) abs(x) else (x)}
    funs <- vfuns <- structure(vector("list", length(imp.functions)), names = imp.functions)
    for (Fun in imp.functions) {
      funs[[Fun]] <- function(x, ...) {
        if (!any(is.finite(x))) NA_real_
        else get(tolower(Fun))(x, ...)
      }
      vfuns[[Fun]] <- function(x, ...) {
        if (!any(is.finite(x))) NA_real_
        else if (Fun == "Mean") geom.mean(x, ...)
        else get(tolower(Fun))(x, ...)
      }
      for (sample in c("Un", weight.names)) {
        if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
          if (cont.treat) {
            B[[paste.(Fun, "Corr", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs0(y[[x, paste.(Fun, "Corr", sample)]])), na.rm = TRUE), numeric(1))
          }
          else {
            B[[paste.(Fun, "Diff", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) abs0(y[[x, paste.(Fun, "Diff", sample)]])), na.rm = TRUE), numeric(1))
            B[[paste.(Fun, "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else vfuns[[Fun]](sapply(summary.tables, function(y) y[[x, paste.(Fun, "V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
            B[[paste.(Fun, "KS", sample)]] <- vapply(Brownames, function(x) funs[[Fun]](sapply(summary.tables, function(y) y[[x, paste.(Fun, "KS", sample)]]), na.rm = TRUE), numeric(1))
          }
        }
      }
    }
  }
  else B <- NULL
  
  return(B)
}
samplesize.across.imps <- function(obs.list) {
  #obs.list <- lapply(bal.tab.imp.list, function(x) x[["Observations"]])
  obs.list <- clear_null(obs.list)
  
  obs <- Reduce("+", obs.list)/length(obs.list)
  attr(obs, "tag") <- paste0("Average ", tolower(attr(obs.list[[1]], "tag")), " across imputations")
  return(obs)
}

#base.bal.tab.multi
balance.table.multi.summary <- function(bal.tab.multi.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, quick = TRUE, types = NULL) {
  if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.multi.list, class)))) {
    bal.tab.multi.list <- lapply(bal.tab.multi.list, function(x) x[["Balance"]])}
  if (length(weight.names) <= 1) weight.names <- "Adj"
  bal.tab.multi.list <- clear_null(bal.tab.multi.list)
  
  Brownames <- unique(do.call("c", lapply(bal.tab.multi.list, rownames)))
  Bcolnames <- c("Type", expand.grid_string(c("Max.Diff", "M.Threshold", "Max.V.Ratio", "V.Threshold", "Max.KS", "KS.Threshold"), 
                                            c("Un", weight.names), collapse = "."))
  B <- as.data.frame(matrix(nrow = length(Brownames), ncol = length(Bcolnames)), row.names = Brownames)
  names(B) <- Bcolnames
  
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(bal.tab.multi.list, function(y) y[[x, "Type"]], character(1))); return(u[!is.na(u)])}), use.names = FALSE)
  
  max_ <- function(x, na.rm = TRUE) {
    if (!any(is.finite(x))) NA_real_
    else max(x, na.rm = na.rm)
  }
  for (sample in c("Un", weight.names)) {
    if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
      B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) abs(y[[x, paste.("Diff", sample)]])), na.rm = TRUE), numeric(1))
      B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("V.Ratio", sample)]]), na.rm = TRUE), numeric(1))
      B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) max_(sapply(bal.tab.multi.list, function(y) y[[x, paste.("KS", sample)]]), na.rm = TRUE), numeric(1))
    }
  }
  
  if (is_not_null(m.threshold)) {
    if (no.adj) {
      B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs_(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
  
  if (is_not_null(v.threshold)) {
    if (no.adj) {
      B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.V.Ratio", i)]]), paste0(ifelse(B[[paste.("Max.V.Ratio", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
  
  if (is_not_null(ks.threshold)) {
    if (no.adj) {
      B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.KS", i)]]), paste0(ifelse(B[[paste.("Max.KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
  
  return(B)
}
samplesize.multi <- function(bal.tab.multi.list, treat.names, focal) {
  if (is_not_null(focal)) which <- c(treat.names[treat.names != focal], focal)
  else which <- treat.names
  bal.tab.multi.list <- clear_null(bal.tab.multi.list)
  obs <- do.call("cbind", unname(lapply(bal.tab.multi.list, function(x) x[["Observations"]])))[, which]
  attr(obs, "tag") <- attr(bal.tab.multi.list[[1]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.multi.list[[1]][["Observations"]], "ss.type")
  return(obs)
}

#base.bal.tab.msm
balance.table.msm.summary <- function(bal.tab.msm.list, weight.names = NULL, no.adj = FALSE, m.threshold = NULL, v.threshold = NULL, ks.threshold = NULL, r.threshold = NULL, quick = TRUE, types = NULL) {
  if ("bal.tab" %in% unique(do.call("c", lapply(bal.tab.msm.list, class)))) {
    bal.tab.msm.list <- lapply(bal.tab.msm.list, function(x) x[["Balance"]])}
  cont.treat <- "Corr.Un" %in% unique(do.call("c", lapply(bal.tab.msm.list, names)))
  if (length(weight.names) <= 1) weight.names <- "Adj"
  
  Brownames <- unique(do.call("c", lapply(bal.tab.msm.list, rownames)))
  Brownames.appear <- vapply(Brownames, function(x) paste(seq_along(bal.tab.msm.list)[sapply(bal.tab.msm.list, function(y) x %in% rownames(y))], collapse = ", "), character(1))
  if (cont.treat) {
    Bcolnames <- c("Type", expand.grid_string(c("Max.Corr", "R.Threshold"), 
                                              c("Un", weight.names), collapse = "."))
  }
  else {
    Bcolnames <- c("Type", expand.grid_string(c("Max.Diff", "M.Threshold", "Max.V.Ratio", "V.Threshold", "Max.KS", "KS.Threshold"), 
                                              c("Un", weight.names), collapse = "."))
  }
  
  B <- as.data.frame(matrix(NA, nrow = length(Brownames), ncol = 1 + length(Bcolnames)), row.names = Brownames)
  names(B) <- c("Times", Bcolnames)
  
  if (is_not_null(types)) B[["Type"]] <- types
  else B[["Type"]] <- unlist(lapply(Brownames, function(x) {u <- unique(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, "Type"]] else NA_character_, character(1))); return(u[!is.na(u)])}), use.names = FALSE)
  
  B[["Times"]] <- Brownames.appear[Brownames]
  
  max_ <- function(x, na.rm = TRUE) {
    if (!any(is.finite(x))) NA_real_
    else max(x, na.rm = na.rm)
  }
  for (sample in c("Un", weight.names)) {
    if (sample == "Un" || !no.adj) { #Only fill in "stat".Adj if no.adj = FALSE
      if (cont.treat) {
        B[[paste.("Max", "Corr", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Corr", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
      }
      else {
        B[[paste.("Max", "Diff", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) abs(y[[x, paste.("Diff", sample)]]) else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
        B[[paste.("Max", "V.Ratio", sample)]] <- vapply(Brownames, function(x) if (B[[x, "Type"]]!="Contin.") NA_real_ else max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("V.Ratio", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
        B[[paste.("Max", "KS", sample)]] <- vapply(Brownames, function(x) max_(vapply(bal.tab.msm.list, function(y) if (x %in% rownames(y)) y[[x, paste.("KS", sample)]] else NA_real_, numeric(1)), na.rm = TRUE), numeric(1))
      }
    }
  }
  
  if (is_not_null(m.threshold)) {
    if (no.adj) {
      B[["M.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Diff.Un"]]), paste0(ifelse(abs_(B[["Max.Diff.Un"]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("M.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Diff", i)]]), paste0(ifelse(abs_(B[[paste.("Max.Diff", i)]]) < m.threshold, "Balanced, <", "Not Balanced, >"), m.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "M.Threshold.Adj"] <- "M.Threshold"
  
  if (is_not_null(v.threshold)) {
    if (no.adj) {
      B[["V.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.V.Ratio.Un"]]), paste0(ifelse(B[, "Max.V.Ratio.Un"] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("V.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.V.Ratio", i)]]), paste0(ifelse(B[[paste.("Max.V.Ratio", i)]] < v.threshold, "Balanced, <", "Not Balanced, >"), v.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "V.Threshold.Adj"] <- "V.Threshold"
  
  if (is_not_null(ks.threshold)) {
    if (no.adj) {
      B[["KS.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.KS.Un"]]), paste0(ifelse(B[["Max.KS.Un"]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("KS.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.KS", i)]]), paste0(ifelse(B[[paste.("Max.KS", i)]] < ks.threshold, "Balanced, <", "Not Balanced, >"), ks.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "KS.Threshold.Adj"] <- "KS.Threshold"
  
  if (is_not_null(r.threshold)) {
    if (no.adj) {
      B[["R.Threshold.Un"]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[["Max.Corr.Un"]]), paste0(ifelse(B[["Max.Corr.Un"]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
    }
    else {
      for (i in weight.names) {
        B[[paste.("R.Threshold", i)]] <- ifelse(B[["Type"]]!="Distance" & is.finite(B[[paste.("Max.Corr", i)]]), paste0(ifelse(B[[paste.("Max.Corr", i)]] < r.threshold, "Balanced, <", "Not Balanced, >"), r.threshold), "")
      }
    }
  }
  if (no.adj || length(weight.names) <= 1) names(B)[names(B) == "R.Threshold.Adj"] <- "R.Threshold"
  
  
  return(B)
}
samplesize.msm <- function(bal.tab.msm.list) {
  obs <- do.call("cbind", lapply(bal.tab.msm.list, function(x) x[["Observations"]]))
  attr(obs, "tag") <- attr(bal.tab.msm.list[[1]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.msm.list[[1]][["Observations"]], "ss.type")
  return(obs)
}

#base.bal.tab.target
balance.table.target.summary <- balance.table.multi.summary
samplesize.target <- function(bal.tab.target.list, treat.names, target.name) {
  which <- treat.names[treat.names != target.name]
  obs <- do.call("cbind", unname(lapply(bal.tab.target.list, function(x) x[["Observations"]])))[, which]
  attr(obs, "tag") <- attr(bal.tab.target.list[[1]][["Observations"]], "tag")
  attr(obs, "ss.type") <- attr(bal.tab.target.list[[1]][["Observations"]], "ss.type")
  return(obs)
}

#love.plot
isColor <- function(x) {
  tryCatch(is.matrix(col2rgb(x)), 
           error = function(e) FALSE)
}
f.recode <- function(f, ...) {
  #Simplified version of forcats::fct_recode
  f <- factor(f)
  new_levels <- unlist(list(...), use.names = TRUE)
  old_levels <- levels(f)
  idx <- match(new_levels, old_levels)
  
  old_levels[idx] <- names(new_levels)
  
  levels(f) <- old_levels
  return(f)
}
seq_int_cycle <- function(begin, end, max) {
  seq(begin, end, by = 1) - max*(seq(begin-1, end-1, by = 1) %/% max)
}
assign.shapes <- function(colors, default.shape = "circle") {
  if (nunique(colors) < length(colors)) {
    shapes <- seq_int_cycle(19, 19 + length(colors) - 1, max = 25)
  }
  else shapes <- rep(default.shape, length(colors))
  return(shapes)
}
shapes.ok <- function(shapes, nshapes) {
  shape_names <- c(
    "circle", paste("circle", c("open", "filled", "cross", "plus", "small")), "bullet",
    "square", paste("square", c("open", "filled", "cross", "plus", "triangle")),
    "diamond", paste("diamond", c("open", "filled", "plus")),
    "triangle", paste("triangle", c("open", "filled", "square")),
    paste("triangle down", c("open", "filled")),
    "plus", "cross", "asterisk"
  )
  shape_nums <- 1:25
  return((length(shapes) == 1 || length(shapes) == nshapes) && ((is.numeric(shapes) && all(shapes %in% shape_nums)) || (is.character(shapes) && all(shapes %in% shape_names))))
}
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
ggarrange_simple <- function (plots, nrow = NULL, ncol = NULL) {
  #A thin version of egg:ggarrange
  
  gtable_frame <- function (g, width = unit(1, "null"), height = unit(1, "null")) {
    panels <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
    pargins <- g[["layout"]][grepl("panel", g[["layout"]][["name"]]),]
    ll <- unique(panels$l)
    margins <- if (length(ll) == 1) unit(0, "pt") else g$widths[ll[-length(ll)] + 2]
    tt <- unique(panels$t)
    fixed_ar <- g$respect
    if (fixed_ar) {
      ar <- as.numeric(g$heights[tt[1]])/as.numeric(g$widths[ll[1]])
      height <- width * (ar/length(ll))
      g$respect <- FALSE
    }
    core <- g[seq(min(tt), max(tt)), seq(min(ll), max(ll))]
    top <- g[seq(1, min(tt) - 1), seq(min(ll), max(ll))]
    bottom <- g[seq(max(tt) + 1, nrow(g)), seq(min(ll), max(ll))]
    left <- g[seq(min(tt), max(tt)), seq(1, min(ll) - 1)]
    right <- g[seq(min(tt), max(tt)), seq(max(ll) + 1, ncol(g))]
    fg <- grid::nullGrob()
    if (length(left)) {
      lg <- gtable::gtable_add_cols(left, unit(1, "null"), 0)
      lg <- gtable::gtable_add_grob(lg, fg, 1, l = 1)
    }
    else {
      lg <- fg
    }
    if (length(right)) {
      rg <- gtable::gtable_add_cols(right, unit(1, "null"))
      rg <- gtable::gtable_add_grob(rg, fg, 1, l = ncol(rg))
    }
    else {
      rg <- fg
    }
    if (length(top)) {
      tg <- gtable::gtable_add_rows(top, unit(1, "null"), 0)
      tg <- gtable::gtable_add_grob(tg, fg, t = 1, l = 1)
    }
    else {
      tg <- fg
    }
    if (length(bottom)) {
      bg <- gtable::gtable_add_rows(bottom, unit(1, "null"), 
                                    -1)
      bg <- gtable::gtable_add_grob(bg, fg, t = nrow(bg), l = 1)
    }
    else {
      bg <- fg
    }
    grobs <- list(fg, tg, fg, lg, core, rg, fg, bg, fg)
    widths <- grid::unit.c(sum(left$widths), width, sum(right$widths))
    heights <- grid::unit.c(sum(top$heights), height, sum(bottom$heights))
    all <- gtable::gtable_matrix("all", grobs = matrix(grobs, ncol = 3, nrow = 3, byrow = TRUE), 
                                 widths = widths, heights = heights)
    
    all[["layout"]][5, "name"] <- "panel"
    if (fixed_ar) 
      all$respect <- TRUE
    all
  }
  
  n <- length(plots)
  
  grobs <- lapply(plots, ggplot2::ggplotGrob)
  
  if (is_null(nrow) && is_null(ncol)) {
    nm <- grDevices::n2mfrow(n)
    nrow <- nm[1]
    ncol <- nm[2]
  }
  
  hw <- lapply(rep(1, n), unit, "null")
  
  fg <- lapply(seq_along(plots), function(i) gtable_frame(g = grobs[[i]], 
                                                          width = hw[[i]], height = hw[[i]]))
  
  spl <- split(fg, rep(1, n))
  
  rows <- lapply(spl, function(r) do.call(gridExtra::gtable_cbind, r))
  
  gt <- do.call(gridExtra::gtable_rbind, rows)
  
  invisible(gt)
}

#bal.plot
get.var.from.list.with.time <- function(var.name, covs.list) {
  var.name.in.covs <- vapply(covs.list, function(x) var.name %in% names(x), logical(1))
  var <- unlist(lapply(covs.list[var.name.in.covs], function(x) x[[var.name]]))
  times <- rep(var.name.in.covs, each = NCOL(covs.list[[1]]))
  return(list(var = var, times = times))
}

#print.bal.tab
print.data.frame_ <- function(x, ...) {
  if (is_not_null(x) && NROW(x) > 0 && NCOL(x) > 0) {
    print.data.frame(x, ...)
  }
}

#set.cobalt.options
acceptable.options <- function() {
  TF <- c(TRUE, FALSE)
  return(list(un = TF,
              continuous = c("raw", "std"),
              binary = c("raw", "std"),
              imbalanced.only = TF,
              disp.means = TF,
              disp.sds = TF,
              disp.v.ratio = TF,
              disp.ks = TF,
              disp.subclass = TF,
              disp.bal.tab = TF,
              cluster.summary = TF,
              cluster.fun = c("min", "mean", "max"),
              imp.summary = TF,
              imp.fun = c("min", "mean", "max"),
              multi.summary = TF,
              msm.summary = TF,
              target.summary = TF,
              int_sep = " * ",
              factor_sep = "_",
              center = TF))
}

#On attach
.onAttach <- function(...) {
  cobaltLib <- dirname(system.file(package = "cobalt"))
  version <- packageDescription("cobalt", lib.loc = cobaltLib)$Version
  BuildDate <- packageDescription("cobalt", lib.loc = cobaltLib)$Date
  
  foo <- paste0(" cobalt (Version ", version, ", Build Date: ", BuildDate, ")\n", 
                "   Please read the documentation at ?bal.tab to understand the default outputs.\n",
                "   Submit bug reports and feature requests to https://github.com/ngreifer/cobalt/issues\n",
                "   Install the development version (not guaranteed to be stable) with:\n",
                "     devtools::install_github(\"ngreifer/cobalt\")\n",
                "   Thank you for using cobalt!")
  packageStartupMessage(foo)
}

.onLoad <- function(libname, pkgname) {
  backports::import(pkgname)
}

#To pass CRAN checks:
utils::globalVariables(c("distance", "addl", "addl.list", "distance.list",
                         "quick", "treat", "Sample", "min.stat",
                         "max.stat", "mean.stat", "count"))

#Internal Utilities not in SHARED.R
`%+%` <- function(...) {
  a <- list(...)
  if (is.atomic(a[[1]])) do.call(crayon::`%+%`, a)
  else do.call(ggplot2::`%+%`, a)
}
strsplits <- function(x, splits, fixed = TRUE, ...) {
  #Link strsplit but takes multiple split values.
  #Only works for one string at a time (in x).
  for (split in splits) x <- unlist(strsplit(x, split, fixed = TRUE, ...))
  return(x[x != ""]) # Remove empty values
}
paste. <- function(..., collapse = NULL) {
  #Like paste0 but with sep = ".'
  paste(..., sep = ".", collapse = collapse)
}

#This document is shared across cobalt, WeightIt, and optweight as a symbolic link.
#Any edits will be automatically synced across all folders. Make sure functions work
#in all packages!
#The original file is in cobalt/R/.

#Strings
word_list <- function(word.list = NULL, and.or = c("and", "or"), is.are = FALSE, quotes = FALSE) {
  #When given a vector of strings, creates a string of the form "a and b"
  #or "a, b, and c"
  #If is.are, adds "is" or "are" appropriately
  L <- length(word.list)
  if (quotes) word.list <- vapply(word.list, function(x) paste0("\"", x, "\""), character(1L))
  if (L == 0) {
    out <- ""
    attr(out, "plural") = FALSE
  }
  else {
    word.list <- word.list[!word.list %in% c(NA_character_, "")]
    L <- length(word.list)
    if (L == 0) {
      out <- ""
      attr(out, "plural") = FALSE
    }
    else if (L == 1) {
      out <- word.list
      if (is.are) out <- paste(out, "is")
      attr(out, "plural") = FALSE
    }
    else {
      and.or <- match_arg(and.or)
      if (L == 2) {
        out <- paste(word.list, collapse = paste0(" ", and.or," "))
      }
      else {
        out <- paste(paste(word.list[seq_len(L-1)], collapse = ", "),
                     word.list[L], sep = paste0(", ", and.or," "))
        
      }
      if (is.are) out <- paste(out, "are")
      attr(out, "plural") = TRUE
    }
    
    
  }
  return(out)
}
firstup <- function(x) {
  #Capitalize first letter
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
expand.grid_string <- function(..., collapse = "") {
  return(apply(expand.grid(...), 1, paste, collapse = collapse))
}
num_to_superscript <- function(x) {
  nums <- setNames(c("\u2070",
                     "\u00B9",
                     "\u00B2",
                     "\u00B3",
                     "\u2074",
                     "\u2075",
                     "\u2076",
                     "\u2077",
                     "\u2078",
                     "\u2079"),
                   as.character(0:9))
  x <- as.character(x)
  splitx <- strsplit(x, "", fixed = TRUE)
  supx <- sapply(splitx, function(y) paste0(nums[y], collapse = ""))
  return(supx)
}
ordinal <- function(x) {
  if (!is.numeric(x) || !is.vector(x) || is_null(x)) stop("x must be a numeric vector.")
  if (length(x) > 1) return(vapply(x, ordinal, character(1L)))
  else {
    x0 <- abs(x)
    out <- paste0(x0, switch(substring(x0, nchar(x0), nchar(x0)),
                             "1" = "st",
                             "2" = "nd",
                             "3" = "rd",
                             "th"))
    if (sign(x) == -1) out <- paste0("-", out)
    
    return(out)
  }
}
round_df_char <- function(df, digits, pad = "0", na_vals = "") {
  if (NROW(df) == 0) return(df)
  nas <- is.na(df)
  if (!is.data.frame(df)) df <- as.data.frame.matrix(df, stringsAsFactors = FALSE)
  infs <- sapply(df, function(x) !is.na(x) & (x == Inf | x == -Inf), simplify = "array")
  rn <- rownames(df)
  cn <- colnames(df)
  df <- as.data.frame(lapply(df, function(col) {
    if (suppressWarnings(all(!is.na(as.numeric(as.character(col)))))) {
      as.numeric(as.character(col))
    } else {
      col
    }
  }), stringsAsFactors = FALSE)
  nums <- vapply(df, is.numeric, logical(1))
  o.negs <- sapply(1:NCOL(df), function(x) if (nums[x]) df[[x]] < 0 else rep(FALSE, length(df[[x]])))
  df[nums] <- round(df[nums], digits = digits)
  
  df[nas | infs] <- ""
  
  df <- as.data.frame(lapply(df, format, scientific = FALSE, justify = "none"), stringsAsFactors = FALSE)
  
  for (i in which(nums)) {
    if (any(grepl(".", df[[i]], fixed = TRUE))) {
      s <- strsplit(df[[i]], ".", fixed = TRUE)
      lengths <- lengths(s)
      digits.r.of.. <- sapply(seq_along(s), function(x) {
        if (lengths[x] > 1) nchar(s[[x]][lengths[x]])
        else 0 })
      df[[i]] <- sapply(seq_along(df[[i]]), function(x) {
        if (df[[i]][x] == "") ""
        else if (lengths[x] <= 1) {
          paste0(c(df[[i]][x], rep(".", pad == 0), rep(pad, max(digits.r.of..) - digits.r.of..[x] + as.numeric(pad != 0))),
                 collapse = "")
        }
        else paste0(c(df[[i]][x], rep(pad, max(digits.r.of..) - digits.r.of..[x])),
                    collapse = "")
      })
    }
  }
  
  df[o.negs & df == 0] <- paste0("-", df[o.negs & df == 0])
  
  # Insert NA placeholders
  df[nas] <- na_vals
  df[infs] <- "N/A"
  
  if (length(rn) > 0) rownames(df) <- rn
  if (length(cn) > 0) names(df) <- cn
  
  return(df)
}
text_box_plot <- function(range.list, width = 12) {
  full.range <- range(unlist(range.list))
  ratio = diff(full.range)/(width+1)
  rescaled.range.list <- lapply(range.list, function(x) round(x/ratio))
  rescaled.full.range <- round(full.range/ratio)
  d <- as.data.frame(matrix(NA_character_, ncol = 3, nrow = length(range.list),
                            dimnames = list(names(range.list), c("Min", paste(rep(" ", width + 1), collapse = ""), "Max"))),
                     stringsAsFactors = FALSE)
  d[,"Min"] <- vapply(range.list, function(x) x[1], numeric(1L))
  d[,"Max"] <- vapply(range.list, function(x) x[2], numeric(1L))
  for (i in seq_len(nrow(d))) {
    spaces1 <- rescaled.range.list[[i]][1] - rescaled.full.range[1]
    #|
    dashes <- max(0, diff(rescaled.range.list[[i]]) - 2)
    #|
    spaces2 <- max(0, diff(rescaled.full.range) - (spaces1 + 1 + dashes + 1))
    
    d[i, 2] <- paste0(paste(rep(" ", spaces1), collapse = ""), "|", paste(rep("-", dashes), collapse = ""), "|", paste(rep(" ", spaces2), collapse = ""))
  }
  return(d)
}
equivalent.factors <- function(f1, f2) {
  return(nunique(f1) == nunique(interaction(f1, f2, drop = TRUE)))
}
equivalent.factors2 <- function(f1, f2) {
  return(qr(matrix(c(rep(1, length(f1)), as.numeric(f1), as.numeric(f2)), ncol = 3))$rank == 2)
}
paste. <- function(..., collapse = NULL) {
  #Like paste0 but with sep = ".'
  paste(..., sep = ".", collapse = collapse)
}
wrap <- function(s, nchar, ...) {
  vapply(s, function(s_) {
    x <- strwrap(s_, width = nchar, ...)
    paste(x, collapse = "\n")
  }, character(1L))
}

#Numbers
check_if_zero <- function(x) {
  # this is the default tolerance used in all.equal
  tolerance <- .Machine$double.eps^0.5
  # If the absolute deviation between the number and zero is less than
  # the tolerance of the floating point arithmetic, then return TRUE.
  # This means, to me, that I can treat the number as 0 rather than
  # -3.20469e-16 or some such.
  abs(x - 0) < tolerance
}
between <- function(x, range, inclusive = TRUE, na.action = FALSE) {
  if (!all(is.numeric(x))) stop("x must be a numeric vector.", call. = FALSE)
  if (length(range) != 2) stop("range must be of length 2.", call. = FALSE)
  if (anyNA(range) || !is.numeric(range)) stop("range must contain numeric entries only.", call. = FALSE)
  range <- sort(range)
  
  if (anyNA(x)) {
    if (length(na.action) != 1 || !is.atomic(na.action)) stop("na.action must be an atomic vector of length 1.", call. = FALSE)
  }
  if (inclusive) out <- ifelse(is.na(x), na.action, x >= range[1] & x <= range[2])
  else out <- ifelse(is.na(x), na.action, x > range[1] & x < range[2])
  
  return(out)
}

#Statistics
binarize <- function(variable, zero = NULL, one = NULL) {
  nas <- is.na(variable)
  if (!is_binary(variable[!nas])) stop(paste0("Cannot binarize ", deparse(substitute(variable)), ": more than two levels."))
  if (is.character(variable)) variable <- factor(variable)
  variable.numeric <- as.numeric(variable)
  if (is_null(zero)) {
    if (is_null(one)) {
      if (0 %in% variable.numeric) zero <- 0
      else zero <- min(variable.numeric, na.rm = TRUE)
    }
    else {
      if (one %in% levels(variable)) zero <- levels(variable)[levels(variable) != one]
      else stop("The argument to \"one\" is not the name of a level of variable.", call. = FALSE)
    }
  }
  else {
    if (zero %in% levels(variable)) zero <- zero
    else stop("The argument to \"zero\" is not the name of a level of variable.", call. = FALSE)
  }
  
  newvar <- setNames(ifelse(!nas & variable.numeric == zero, 0L, 1L), names(variable))
  newvar[nas] <- NA_integer_
  return(newvar)
}
ESS <- function(w) {
  sum(w)^2/sum(w^2)
}
center <- function(x, at = NULL, na.rm = TRUE) {
  if (is.data.frame(x)) {
    x <- as.matrix.data.frame(x)
    type <- "df"
  }
  if (!is.numeric(x)) stop("x must be numeric.")
  else if (is.array(x) && length(dim(x)) > 2) stop("x must be a numeric or matrix-like (not array).")
  else if (!is.matrix(x)) {
    x <- matrix(x, ncol = 1)
    type <- "vec"
  }
  else type <- "matrix"
  if (is_null(at)) at <- colMeans(x, na.rm = na.rm)
  else if (length(at) %nin% c(1, ncol(x))) stop("at is not the right length.")
  out <- x - matrix(at, byrow = TRUE, ncol = ncol(x), nrow = nrow(x))
  if (type == "df") out <- as.data.frame.matrix(out)
  else if (type == "vec") out <- drop(out)
  return(out)
}
w.m <- function(x, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- rep(1, length(x))
  w[is.na(x)] <- NA_real_
  return(sum(x*w, na.rm=na.rm)/sum(w, na.rm=na.rm))
}
col.w.m <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) w <- 1
  w.sum <- colSums(w*!is.na(mat))
  return(colSums(mat*w, na.rm = na.rm)/w.sum)
}
col.w.v <- function(mat, w = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is.data.frame(mat)) {
      if (any(vapply(mat, is_, logical(1L), types = c("factor", "character")))) {
        stop("mat must be a numeric matrix.")
      }
      else mat <- as.matrix.data.frame(mat)
    }
    else if (is.numeric(mat)) {
      mat <- matrix(mat, ncol = 1)
    }
    else stop("mat must be a numeric matrix.")
  }
  
  if (is_null(bin.vars)) bin.vars <- rep(FALSE, ncol(mat))
  else if (length(bin.vars) != ncol(mat) || any(is.na(as.logical(bin.vars)))) {
    stop("bin.vars must be a logical vector with length equal to the number of columns of mat.", call. = FALSE)
  }
  bin.var.present <- any(bin.vars)
  non.bin.vars.present <- any(!bin.vars)
  
  var <- setNames(numeric(ncol(mat)), colnames(mat))
  if (is_null(w)) {
    if (non.bin.vars.present) {
      den <- colSums(!is.na(mat[, !bin.vars, drop = FALSE])) - 1
      var[!bin.vars] <- colSums(center(mat[, !bin.vars, drop = FALSE])^2, na.rm = na.rm)/den
    }
    if (bin.var.present) {
      means <- colMeans(mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else if (na.rm && any(is.na(mat))) {
    n <- nrow(mat)
    w <- array(w, dim = dim(mat))
    w[is.na(mat)] <- NA_real_
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    if (non.bin.vars.present) {
      x <- sqrt(w[, !bin.vars, drop = FALSE]) * center(mat[, !bin.vars, drop = FALSE], 
                                                       at = colSums(w[, !bin.vars, drop = FALSE] * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - colSums(w[, !bin.vars, drop = FALSE]^2, na.rm = na.rm))
    }
    if (bin.var.present) {
      means <- colSums(w[, bin.vars, drop = FALSE] * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  else {
    if (is_null(w)) w <- rep(1, nrow(mat))
    w <- w/sum(w)
    if (non.bin.vars.present) {
      x <- sqrt(w) * center(mat[, !bin.vars, drop = FALSE], 
                            at = colSums(w * mat[, !bin.vars, drop = FALSE], na.rm = na.rm))
      var[!bin.vars] <- colSums(x*x, na.rm = na.rm)/(1 - sum(w^2))
    }
    if (bin.var.present) {
      means <- colSums(w * mat[, bin.vars, drop = FALSE], na.rm = na.rm)
      var[bin.vars] <- means * (1 - means)
    }
  }
  return(var)
}
col.w.cov <- function(mat, y, w = NULL, na.rm = TRUE) {
  if (!is.matrix(mat)) {
    if (is_null(w)) return(cov(mat, y, use = if (na.rm) "pair" else "everything"))
    else mat <- matrix(mat, ncol = 1)
  }
  if (is_null(w)) {
    y <- array(y, dim = dim(mat))
    y[is.na(mat)] <- NA
    mat[is.na(y)] <- NA
    den <- colSums(!is.na(mat*y)) - 1
    cov <- colSums(center(mat, na.rm = na.rm)*center(y, na.rm = na.rm), na.rm = na.rm)/den
  }
  else if (na.rm && any(is.na(mat))) {
    n <- nrow(mat)
    w <- array(w, dim = dim(mat))
    w[is.na(mat)] <- NA_real_
    s <- colSums(w, na.rm = na.rm)
    w <- mat_div(w, s)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x*y, na.rm = na.rm)/(1 - colSums(w^2, na.rm = na.rm))
  }
  else {
    n <- nrow(mat)
    w <- w/sum(w)
    x <- w * center(mat, at = colSums(w * mat, na.rm = na.rm))
    cov <- colSums(x*y, na.rm = na.rm)/(1 - sum(w^2))
  }
  return(cov)
}
col.w.r <- function(mat, y, w = NULL, s.weights = NULL, bin.vars = NULL, na.rm = TRUE) {
  if (is_null(w) && is_null(s.weights)) return(cor(mat, y, w, use = if (na.rm) "pair" else "everything"))
  else {
    cov <- col.w.cov(mat, y = y, w = w, na.rm = na.rm)
    den <- sqrt(col.w.v(mat, w = s.weights, bin.vars = bin.vars, na.rm = na.rm)) * 
      sqrt(col.w.v(y, w = s.weights, na.rm = na.rm))
    return(cov/den)
  }
}
coef.of.var <- function(x, pop = TRUE, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  if (pop) sqrt(mean((x-mean(x))^2))/mean(x)
  else sd(x)/mean(x)
}
mean.abs.dev <- function(x, na.rm = TRUE) {
  if (na.rm) x <- x[!is.na(x)]
  mean(abs(x - mean(x)))
}
geom.mean <- function(y, na.rm = TRUE) {
  exp(mean(log(y[is.finite(log(y))]), na.rm = na.rm))
}
mat_div <- function(mat, vec) {
  mat/vec[col(mat)]
}
abs_ <- function(x, ratio = FALSE) {
  if (ratio) pmax(x, 1/x)
  else (abs(x))
}

#Formulas
is.formula <- function(f, sides = NULL) {
  res <- is.name(f[[1]])  && deparse(f[[1]]) %in% c( '~', '!') &&
    length(f) >= 2
  if (is_not_null(sides) && is.numeric(sides) && sides %in% c(1,2)) {
    res <- res && length(f) == sides + 1
  }
  return(res)
}

#treat/covs
get.covs.and.treat.from.formula <- function(f, data = NULL, terms = FALSE, sep = "", ...) {
  A <- list(...)
  
  #Check if data exists
  if (is_not_null(data)) {
    if (is.data.frame(data)) {
      data.specified <- TRUE
    }
    else {
      warning("The argument supplied to data is not a data.frame object. This may causes errors or unexpected results.", call. = FALSE)
      data <- environment(f)
      data.specified <- FALSE
    }
  }
  else {
    data <- environment(f)
    data.specified <- FALSE
  }
  
  env <- environment(f)
  
  tryCatch(tt <- terms(f, data = data),
           error = function(e) {
             if (conditionMessage(e) == "'.' in formula and no 'data' argument") {
               stop("'.' is not allowed in formulas.", call. = FALSE)
             }
             else stop(conditionMessage(e), call. = FALSE)
           })
  
  #Check if response exists
  if (is.formula(tt, 2)) {
    resp.vars.mentioned <- as.character(tt)[2]
    resp.vars.failed <- vapply(resp.vars.mentioned, function(v) {
      null_or_error(try(eval(parse(text=v)[[1]], data, env), silent = TRUE))
    }, logical(1L))
    
    if (any(resp.vars.failed)) {
      if (is_null(A[["treat"]])) stop(paste0("The given response variable, \"", as.character(tt)[2], "\", is not a variable in ", word_list(c("data", "the global environment")[c(data.specified, TRUE)], "or"), "."), call. = FALSE)
      tt <- delete.response(tt)
    }
  }
  else resp.vars.failed <- TRUE
  
  if (any(!resp.vars.failed)) {
    treat.name <- resp.vars.mentioned[!resp.vars.failed][1]
    treat <- eval(parse(text=treat.name)[[1]], data, env)
  }
  else {
    treat <- A[["treat"]]
    treat.name <- NULL
  }
  
  #Check if RHS variables exist
  tt.covs <- delete.response(tt)
  
  rhs.vars.mentioned.lang <- attr(tt.covs, "variables")[-1]
  rhs.vars.mentioned <- vapply(rhs.vars.mentioned.lang, deparse, character(1L))
  rhs.vars.failed <- vapply(rhs.vars.mentioned, function(v) {
    null_or_error(try(eval(parse(text=v)[[1]], data, env), silent = TRUE))
  }, logical(1L))
  
  if (any(rhs.vars.failed)) {
    stop(paste0(c("All variables in formula must be variables in data or objects in the global environment.\nMissing variables: ",
                  paste(rhs.vars.mentioned[rhs.vars.failed], collapse=", "))), call. = FALSE)
    
  }
  
  rhs.term.labels <- attr(tt.covs, "term.labels")
  rhs.term.orders <- attr(tt.covs, "order")
  
  rhs.df <- vapply(rhs.vars.mentioned, function(v) {
    d <- try(eval(parse(text=v)[[1]], data, env), silent = TRUE)
    is.data.frame(d) || is.matrix(d)
  }, logical(1L))
  
  if (any(rhs.df)) {
    if (any(rhs.vars.mentioned[rhs.df] %in% unlist(lapply(rhs.term.labels[rhs.term.orders > 1], function(x) strsplit(x, ":", fixed = TRUE))))) {
      stop("Interactions with data.frames are not allowed in the input formula.", call. = FALSE)
    }
    addl.dfs <- setNames(lapply(rhs.vars.mentioned[rhs.df], function(x) {as.data.frame(eval(parse(text=x)[[1]], data, env))}),
                         rhs.vars.mentioned[rhs.df])
    
    for (i in rhs.term.labels[rhs.term.labels %in% rhs.vars.mentioned[rhs.df]]) {
      ind <- which(rhs.term.labels == i)
      rhs.term.labels <- append(rhs.term.labels[-ind],
                                values = names(addl.dfs[[i]]),
                                after = ind - 1)
    }
    
    if (data.specified) data <- do.call("cbind", unname(c(addl.dfs, list(data))))
    else data <- do.call("cbind", unname(addl.dfs))
  }
  
  if (is_null(rhs.term.labels)) {
    new.form <- as.formula("~ 1")
    tt.covs <- terms(new.form)
    covs <- data.frame(Intercept = rep(1, if (is_null(treat)) 1 else length(treat)))
    if (is_not_null(treat.name) && treat.name == "Intercept") {
      names(covs) <- "Intercept_"
    }
  }
  else {
    new.form <- as.formula(paste("~", paste(rhs.term.labels, collapse = " + ")))
    tt.covs <- terms(new.form)
    attr(tt.covs, "intercept") <- 0
    
    #Get model.frame, report error
    mf.covs <- quote(stats::model.frame(tt.covs, data,
                                        drop.unused.levels = TRUE,
                                        na.action = "na.pass"))
    
    tryCatch({covs <- eval(mf.covs)},
             error = function(e) {stop(conditionMessage(e), call. = FALSE)})
    
    if (is_not_null(treat.name) && treat.name %in% names(covs)) stop("The variable on the left side of the formula appears on the right side too.", call. = FALSE)
  }
  
  if (s <- !identical(sep, "")) {
    if (!is.character(sep) || length(sep) > 1) stop("sep must be a string of length 1.", call. = FALSE)
    original.covs.levels <- setNames(vector("list", ncol(covs)), names(covs))
    for (i in names(covs)) {
      if (is.character(covs[[i]])) covs[[i]] <- factor(covs[[i]])
      if (is.factor(covs[[i]])) {
        original.covs.levels[[i]] <- levels(covs[[i]])
        levels(covs[[i]]) <- paste0(sep, original.covs.levels[[i]])
      }
    }
  }
  
  #Get full model matrix with interactions too
  covs.matrix <- model.matrix(tt.covs, data = covs,
                              contrasts.arg = lapply(Filter(is.factor, covs),
                                                     contrasts, contrasts=FALSE))
  
  if (s) {
    for (i in names(covs)) {
      if (is.factor(covs[[i]])) {
        levels(covs[[i]]) <- original.covs.levels[[i]]
      }
    }
  }
  
  if (!terms) attr(covs, "terms") <- NULL
  
  return(list(reported.covs = covs,
              model.covs = covs.matrix,
              treat = treat,
              treat.name = treat.name))
}
assign.treat.type <- function(treat) {
  #Returns treat with treat.type attribute
  nunique.treat <- nunique(treat)
  if (nunique.treat == 2) {
    treat.type <- "binary"
  }
  else if (nunique.treat < 2) {
    stop("The treatment must have at least two unique values.", call. = FALSE)
  }
  else if (is.factor(treat) || is.character(treat)) {
    treat.type <- "multinomial"
    treat <- factor(treat)
  }
  else {
    treat.type <- "continuous"
  }
  attr(treat, "treat.type") <- treat.type
  return(treat)
}
get.treat.type <- function(treat) {
  return(attr(treat, "treat.type"))
}
has.treat.type <- function(treat) {
  is_not_null(get.treat.type(treat))
}
process.s.weights <- function(s.weights, data = NULL) {
  #Process s.weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (is_null(data)) {
        stop("s.weights was specified as a string but there was no argument to data.", call. = FALSE)
      }
      else if (s.weights %in% names(data)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  return(s.weights)
}

#Uniqueness
nunique <- function(x, nmax = NA, na.rm = TRUE) {
  if (is_null(x)) return(0)
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (is.factor(x)) return(nlevels(x))
    else return(length(unique(x, nmax = nmax)))
  }
  
}
nunique.gt <- function(x, n, na.rm = TRUE) {
  if (missing(n)) stop("n must be supplied.")
  if (n < 0) stop("n must be non-negative.")
  if (is_null(x)) FALSE
  else {
    if (na.rm) x <- x[!is.na(x)]
    if (n == 1) !all_the_same(x)
    else if (length(x) < 2000) nunique(x) > n
    else tryCatch(nunique(x, nmax = n) > n, error = function(e) TRUE)
  }
}
all_the_same <- function(x) {
  if (is.double(x)) check_if_zero(abs(max(x) - min(x)))
  else !any(x != x[1])
}
is_binary <- function(x) !all_the_same(x) && all_the_same(x[x != x[1]])

#R Processing
is_ <- function(x, types, stop = FALSE, arg.to = FALSE) {
  s1 <- deparse(substitute(x))
  if (is_not_null(x)) {
    for (i in types) {
      if (i == "list") it.is <- is.vector(x, "list")
      else if (is_not_null(get0(paste.("is", i)))) {
        it.is <- get0(paste.("is", i))(x)
      }
      else it.is <- inherits(x, i)
      if (it.is) break
    }
  }
  else it.is <- FALSE
  
  if (stop) {
    if (!it.is) {
      s0 <- ifelse(arg.to, "The argument to ", "")
      s2 <- ifelse(any(types %in% c("factor", "character", "numeric", "logical")),
                   "vector", "")
      stop(paste0(s0, s1, " must be a ", word_list(types, and.or = "or"), " ", s2, "."), call. = FALSE)
    }
  }
  else {
    return(it.is)
  }
}
is_null <- function(x) length(x) == 0L
is_not_null <- function(x) !is_null(x)
clear_null <- function(x) {
  x[vapply(x, is_null, logical(1L))] <- NULL
  return(x)
}
probably.a.bug <- function() {
  fun <- paste(deparse(sys.call(-1)), collapse = "\n")
  stop(paste0("An error was produced and is likely a bug. Please let the maintainer know a bug was produced by the function\n",
              fun), call. = FALSE)
}
`%nin%` <- function(x, table) is.na(match(x, table, nomatch = NA_integer_))
null_or_error <- function(x) {is_null(x) || class(x) == "try-error"}
match_arg <- function(arg, choices, several.ok = FALSE) {
  #Replaces match.arg() but gives cleaner error message and processing
  #of arg.
  if (missing(arg))
    stop("No argument was supplied to match_arg.", call. = FALSE)
  arg.name <- deparse(substitute(arg))
  
  if (missing(choices)) {
    formal.args <- formals(sys.function(sysP <- sys.parent()))
    choices <- eval(formal.args[[as.character(substitute(arg))]],
                    envir = sys.frame(sysP))
  }
  
  if (is.null(arg))
    return(choices[1L])
  else if (!is.character(arg))
    stop(paste0("'", arg.name, "' must be NULL or a character vector"), call. = FALSE)
  if (!several.ok) {
    if (identical(arg, choices))
      return(arg[1L])
    if (length(arg) > 1L)
      stop(paste0("'", arg.name, "' must be of length 1"), call. = FALSE)
  }
  else if (is_null(arg))
    stop(paste0("'", arg.name, "' must be of length >= 1"), call. = FALSE)
  
  i <- pmatch(arg, choices, nomatch = 0L, duplicates.ok = TRUE)
  if (all(i == 0L))
    stop(paste0("'", arg.name, "' should be one of ", word_list(choices, and.or = "or", quotes = TRUE), "."),
         call. = FALSE)
  i <- i[i > 0L]
  if (!several.ok && length(i) > 1)
    stop("there is more than one match in 'match_arg'")
  choices[i]
}
last <- function(x) {
  x[[length(x)]]
}
len <- function(x, recursive = TRUE) {
  if (is.vector(x, "list")) sapply(x, len)
  else if (length(dim(x)) > 1) NROW(x)
  else length(x)
}

#Defunct; delete if everything works without them
.center <- function(x, na.rm = TRUE, at = NULL) {
  dimx <- dim(x)
  if (length(dimx) == 2L) x <- apply(x, 2, center, na.rm = na.rm, at = at)
  else if (length(dimx) > 2L) stop("x must be a numeric or matrix-like (not array).")
  else if (!is.numeric(x)) warning("x is not numeric and will not be centered.")
  else {
    if (is_null(at)) at <- mean(x, na.rm = na.rm)
    else if (!is.numeric(at)) stop("at must be numeric.")
    x <- x - at
  }
  return(x)
}
.w.v <- function(x, w = NULL, na.rm = TRUE) {
  .w.cov(x, x, w = w, na.rm = na.rm)
}
.w.cov <- function(x, y, w = NULL, na.rm = TRUE, type = 3) {
  
  if (length(x) != length(y)) stop("x and y must the same length")
  
  if (is_null(w)) w <- rep(1, length(x))
  else if (length(w) != length(x)) stop("weights must be same length as x and y")
  
  w[is.na(x) | is.na(y)] <- NA_real_
  
  wmx <- w.m(x, w, na.rm = na.rm)
  wmy <- w.m(y, w, na.rm = na.rm)
  
  wcov <- sum(w*(x - wmx)*(y - wmy), na.rm = na.rm) / .w.cov.scale(w, na.rm = na.rm, type = type)
  return(wcov)
}
.w.cov.scale <- function(w, type = 3, na.rm = TRUE) {
  
  sw <- sum(w, na.rm = na.rm)
  n <- sum(!is.na(w))
  vw1 <- sum((w - sw/n)^2, na.rm = na.rm)/n
  # vw2 <- sum((w - sw/n)^2, na.rm = na.rm)/(n-1)
  
  if (type == 1) sw
  else if (type == 2) sw - 1
  else if (type == 3) sw*(n-1)/n - vw1*n/sw
  # else if (type == 4) sw*(n-1)/n - vw2*n/sw
  
}
.w.r <- function(x, y, w = NULL, s.weights = NULL) {
  #Computes weighted correlation but using the unweighted (s.weighted) variances
  #in the denominator.
  if (is_null(s.weights)) s.weights <- rep(1, length(x))
  else if (length(s.weights) != length(x)) stop("s.weights must be same length as x and y")
  
  s.weights[is.na(x) | is.na(y)] <- NA_real_
  
  w_ <- w*s.weights
  
  r <- .w.cov(x, y, w_) / (sqrt(.w.v(x, s.weights) * .w.v(y, s.weights)))
  
  return(r)
}
.col.w.v <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    w <- rep(1, nrow(mat))
  }
  means <- col.w.m(mat, w, na.rm)
  w.scale <- apply(mat, 2, function(x) .w.cov.scale(w[!is.na(x)]))
  vars <- colSums(w*center(mat, at = means)^2, na.rm = na.rm)/w.scale
  
  return(vars)
}
.col.w.v.bin <- function(mat, w = NULL, na.rm = TRUE) {
  if (is_null(w)) {
    w <- rep(1, nrow(mat))
  }
  means <- col.w.m(mat, w, na.rm)
  vars <- means * (1 - means)
  return(vars)
}

#Functions to convert object to base.bal.tab input

x2base <- function(...) UseMethod("x2base")

x2base.matchit <- function(m, ...) {
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  
  if (any(class(m) == "matchit.subclass")) {
    subclass <- factor(m$subclass)
    method <- "subclassification"
  }
  else if (any(class(m) == "matchit.full")) {
    subclass <- NULL
    method <- "weighting"
  }
  else {
    subclass <- NULL
    method <- "matching"
  }
  
  weights <- data.frame(weights = m$weights)
  treat <- m$treat
  data <- A$data
  subset <- A$subset
  cluster <- A$cluster
  s.d.denom <- A$s.d.denom
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  
  if (is_not_null(m$model$model)) {
    if (nrow(m$model$model) == length(treat)) {
      covs <- data.frame(m$model$model[, names(m$model$model) %in% attributes(terms(m$model))$term.labels])
    }
    else {
      #Recreating covs from model object and m$X. Have to do this because when 
      #drop != NULL and reestimate = TRUE, cases are lost. This recovers them.
      
      order <- setNames(attr(m$model$terms, "order"),
                        attr(m$model$terms, "term.labels"))
      assign <- setNames(attr(m$X, "assign"), colnames(m$X))
      assign1 <- assign[assign %in% which(order == 1)] #Just main effects
      
      dataClasses <- attr(m$model$terms, "dataClasses")
      factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
      f0 <- setNames(lapply(factors.to.unsplit, 
                            function(x) {
                              if (dataClasses[x] == "factor")
                                list(levels = levels(m$model$model[[x]]),
                                     faclev = paste0(x, levels(m$model$model[[x]])))
                              else 
                                list(levels = unique(m$model$model[[x]]),
                                     faclev = paste0(x, unique(m$model$model[[x]])))
                            }),
                     factors.to.unsplit)
      covs <- as.data.frame(m$X[, names(assign1)])
      for (i in factors.to.unsplit) {
        covs <- unsplitfactor(covs, i, sep = "",
                              dropped.level = f0[[i]]$levels[f0[[i]]$faclev %nin% colnames(m$X)])
        if (dataClasses[i] == "logical") covs[[i]] <- as.logical(covs[[i]])
      }
    }
    
  }
  else if ("matchit.mahalanobis" %nin% class(m) && is_not_null(data)) {
    t.c <- get.covs.and.treat.from.formula(m$formula, data = data)
    covs <- t.c[["reported.covs"]]
  }
  else {
    covs <- data.frame(m$X)
  }
  m.data <- m$model$data
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(m.data) == cluster)) {
        cluster <- m.data[[cluster]]
      }
    }
    else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is_(subset, "logical")) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, m.data))
  }
  
  if (any(is.finite(m$distance))) {
    if (is_not_null(distance)) distance <- cbind(distance, distance = m$distance)
    else distance <- data.frame(distance = m$distance)
  }
  
  #Get s.d.denom
  estimand <- "ATT"
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method = method)
  
  ensure.equal.lengths <- TRUE
  vectors <- c("cluster", "treat", "subset", "subclass")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors, ifnotfound = list(NULL))), 
                        vapply(data.frames, 
                               function(x) NROW(get0(x)),
                               numeric(1L))), c(vectors, data.frames))
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to matchit()."), call. = FALSE)
  }
  
  if (any(c(is.na(covs), is.na(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$method <- method
  X$treat <- treat
  X$weights <- weights
  X$discarded <- m$discarded
  X$covs <- covs
  X$distance <- distance
  X$addl <- addl
  X$cluster <- factor(cluster)
  X$call <- m$call
  X$subclass <- factor(subclass)
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
}
x2base.ps <- function(ps, ...) {
  #stop.method
  #s.d.denom
  A <- list(...)
  
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
  if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
  
  if (is_not_null(A$stop.method)) {
    if (is.character(A$stop.method)) {
      rule1 <- names(ps$w)[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(names(ps$w)), x))), any)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- names(ps$w)
      }
    }
    else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(names(ps$w)))) {
      if (any(!A$stop.method %in% seq_along(names(ps$w)))) {
        message(paste0("Warning: There are ", length(names(ps$w)), " stop methods available, but you requested ", 
                       word_list(A$stop.method[!A$stop.method %in% seq_along(names(ps$w))], and.or = "and"),"."))
      }
      rule1 <- names(ps$w)[A$stop.method %in% seq_along(names(ps$w))]
    }
    else {
      warning("stop.method should be ", word_list(names(ps$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- names(ps$w)
    }
  }
  else {
    rule1 <- names(ps$w)
  }
  
  s <- names(ps$w)[match(tolower(rule1), tolower(names(ps$w)))]
  estimand <- ps$estimand
  
  weights <- data.frame(get.w(ps, s, estimand))
  treat <- ps$treat
  covs <- ps$data[, ps$gbm.obj$var.names, drop = FALSE]
  data <- A$data
  ps.data <- ps$data
  cluster <- A$cluster
  subset <- A$subset
  s.weights <- ps$sampw
  s.d.denom <- A$s.d.denom
  method <- rep("weighting", ncol(weights))
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(ps.data) == cluster)) {
        cluster <- ps.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, ps.data))
  }
  
  if (is_not_null(distance)) {
    if (length(s) == 1) {
      distance <- cbind(distance, prop.score = ps$ps[[s]])
    }
    else {
      distance <- cbind(distance, prop.score = ps$ps[s])
    }
  }
  else {
    if (length(s) == 1) {
      distance <- data.frame(prop.score = ps$ps[[s]])
    }
    else {
      distance <- data.frame(prop.score = ps$ps[s])
    }
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  
  X$weights <- weights
  X$treat <- treat
  X$distance <- distance
  X$addl <- addl
  X$covs <- covs
  X$call <- ps$parameters
  X$cluster <- factor(cluster)
  X$method <- method
  X$s.weights <- s.weights
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
}
x2base.mnps <- function(mnps, ...) {
  #stop.method
  #s.d.denom
  A <- list(...)
  
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  if (is_not_null(A)&& names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
  if (is_null(A$stop.method) == 0 && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
  
  if (is_not_null(A$stop.method)) {
    if (any(is.character(A$stop.method))) {
      rule1 <- mnps$stopMethods[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(mnps$stopMethods), x))), any)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- mnps$stopMethods
      }
    }
    else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(mnps$stopMethods))) {
      if (any(!A$stop.method %in% seq_along(mnps$stopMethods))) {
        message(paste0("Warning: There are ", length(mnps$stopMethods), " stop methods available, but you requested ", 
                       word_list(A$stop.method[!A$stop.method %in% seq_along(mnps$stopMethods)], and.or = "and"),"."))
      }
      rule1 <- mnps$stopMethods[A$stop.method %in% seq_along(mnps$stopMethods)]
    }
    else {
      warning("stop.method should be ", word_list(mnps$stopMethods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- mnps$stopMethods
    }
  }
  else {
    rule1 <- mnps$stopMethods
  }
  
  s <- mnps$stopMethods[match(tolower(rule1), tolower(mnps$stopMethods))]
  
  weights <- data.frame(get.w(mnps, s))
  treat <- mnps$treatVar
  covs <- mnps$data[mnps$psList[[1]]$gbm.obj$var.names]
  data <- A$data
  cluster <- A$cluster
  subset <- A$subset
  s.weights <- mnps$sampw
  s.d.denom <- A$s.d.denom
  focal <- mnps$treatATT
  method <- rep("weighting", ncol(weights))
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  mnps.data <- mnps$data
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(mnps.data) == cluster)) {
        cluster <- mnps.data[[cluster]]
      }
    }
    else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, mnps.data))
  }    
  # model.ps <- setNames(as.data.frame(lapply(mnps$psList, function(x) x$ps)), names(mnps$psList))
  # model.ps.combined <- numeric(length(treat))
  # for (i in levels(treat)) {
  #     model.ps.combined[treat == i] <- model.ps[treat == i, i]
  # }
  # if (length(distance) > 0) distance <- cbind(distance, prop.score = model.ps.combined)
  # else distance <- data.frame(prop.score = model.ps.combined)
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps()."), call. = FALSE)
  }
  
  #Get s.d.denom
  estimand <- mnps$estimand
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  
  X$weights <- weights
  X$treat <- treat
  X$distance <- distance
  X$addl <- addl
  X$covs <- covs
  X$call <- NULL
  X$cluster <- factor(cluster)
  X$s.weights <- mnps$sampw
  X$focal <- focal
  X$method <- method
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.ps.cont <- function(ps.cont, ...) {
  #stop.method
  #s.d.denom
  A <- list(...)
  
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]]
  if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
  
  if (is_not_null(A$stop.method)) {
    if (is.character(A$stop.method)) {
      rule1 <- names(ps.cont$w)[sapply(t(sapply(tolower(A$stop.method), function(x) startsWith(tolower(names(ps.cont$w)), x))), any)]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- names(ps.cont$w)
      }
    }
    else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(names(ps.cont$w)))) {
      if (any(!A$stop.method %in% seq_along(names(ps.cont$w)))) {
        message(paste0("Warning: There are ", length(names(ps.cont$w)), " stop methods available, but you requested ", 
                       word_list(A$stop.method[!A$stop.method %in% seq_along(names(ps.cont$w))], and.or = "and"),"."))
      }
      rule1 <- names(ps.cont$w)[A$stop.method %in% seq_along(names(ps.cont$w))]
    }
    else {
      warning("stop.method should be ", word_list(names(ps.cont$w), and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- names(ps.cont$w)
    }
  }
  else {
    rule1 <- names(ps.cont$w)
  }
  
  s <- names(ps.cont$w)[match(tolower(rule1), tolower(names(ps.cont$w)))]
  
  weights <- data.frame(get.w(ps.cont, s))
  treat <- ps.cont$treat
  covs <- ps.cont$data[, ps.cont$gbm.obj$var.names, drop = FALSE]
  data <- A$data
  ps.data <- ps.cont$data
  cluster <- A$cluster
  subset <- A$subset
  s.weights <- ps.cont$sampw
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(ps.data) == cluster)) {
        cluster <- ps.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, ps.data))
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to ps.cont()."), call. = FALSE)
  }
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  
  X$weights <- weights
  X$treat <- treat
  X$distance <- distance
  X$addl <- addl
  X$covs <- covs
  X$call <- ps.cont$parameters
  X$cluster <- factor(cluster)
  X$method <- rep("weighting", ncol(weights))
  X$s.weights <- s.weights
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "cont"
  
  return(X)
}
x2base.Match <- function(Match, ...) {
  
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Checks
  if (is_not_null(Match) && !is.list(Match)) {
    stop("'Match' object contains no valid matches")}
  
  #Get treat and covs
  data <- A$data
  t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
  
  #Initializing variables
  m <- Match
  estimand <- m$estimand
  method <- "matching"
  s.d.denom <- A$s.d.denom
  
  treat <- t.c[["treat"]]
  covs  <- t.c[["covs"]]
  
  weights <- data.frame(weights = get.w.Match(m))
  
  cluster <- A$cluster
  subset <- A$subset
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  dropped <- rep(FALSE, length(treat))
  if (is_not_null(m$index.dropped)) dropped[m$index.dropped] <- TRUE
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data))
  }
  
  ensure.equal.lengths <- TRUE
  covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
  vectors <- c("treat", "cluster", "subset")
  data.frames <- c(covs.data, "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in names(lengths)[names(lengths) != "weights"]) {
      if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to Match()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$treat <- treat
  X$weights <- weights
  X$discarded <- dropped
  X$distance <- NULL #NAs in distance bcause of incomplete list in Match object
  X$addl <- addl
  X$covs <- covs
  X$call <- NULL
  X$method <- method
  X$cluster <- factor(cluster)
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
}
x2base.formula <- function(formula, ...) {
  A <- list(...)
  
  if ("data" %in% names(A) && is_(A[["data"]], "mids")) {
    A[["data"]] <- imp.complete(A[["data"]])
    if (is_null(A[["imp"]])) A[["imp"]] <- A[["data"]][[".imp"]]
  }
  
  t.c <- get.covs.and.treat.from.formula(formula, A[["data"]], treat = A[["treat"]])
  covs <- t.c[["reported.covs"]]
  treat <- t.c[["treat"]]
  
  if (is_null(covs)) stop("The right hand side of the formula must contain covariates for which balance is to be assessed.", call. = FALSE)
  
  A[["covs"]] <- NULL
  A[["treat"]] <- NULL
  
  X <- do.call("x2base.data.frame", c(list(covs = covs, treat = treat), A))
  return(X)
}
x2base.data.frame <- function(covs, ...) {
  #treat
  #data
  #weights
  #distance
  #subclass
  #match.strata
  #addl
  #s.d.denom
  #method
  #cluster
  #estimand
  #imp
  
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  treat <- A$treat
  data <- A$data
  weights <- A$weights
  distance <- A$distance
  subclass <- A$subclass
  match.strata <- A$match.strata
  cluster <- A$cluster
  addl <- A$addl
  s.d.denom <- A$s.d.denom
  method <- A$method
  estimand <- A$estimand
  imp <- A$imp
  s.weights <- A$s.weights
  subset <- A$subset
  focal <- A$focal
  
  #Checks
  if (is_null(covs)) {
    stop("covs data.frame must be specified.", call. = FALSE)
  }
  is_(covs, "data.frame", stop = TRUE)
  
  #Process data
  if (is_not_null(data)) {
    if (is_(data, "mids")) {
      data <- imp.complete(data)
      if ("imp" %nin% names(A)) imp <- data[[".imp"]]
    }
    else if (!is_(data, "data.frame"))
    {
      warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
      data <- NULL
    }
  }
  
  specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
  if (is_not_null(weights)) {
    if (!is_(weights, c("character", "numeric", "data.frame", "list"))) {
      stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
    }
    specified["weights"] <- TRUE
  }
  if (is_not_null(subclass)){
    if (!is_(subclass, c("character", "numeric", "factor"))) {
      stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
    }
    specified["subclass"] <- TRUE
  }
  if (is_not_null(match.strata)) {
    if (!is_(match.strata, c("character", "numeric", "factor"))) {
      stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
    }
    specified["match.strata"] <- TRUE
  }
  
  #Getting method
  if (is_null(method)) {
    if (specified["match.strata"]) {
      if (sum(specified) > 1) {
        message(word_list(names(specified)[specified]), " are specified. Assuming \"matching\" and using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
        weights <- subclass <- NULL
      }
      method <- "matching"
    }
    else if (specified["subclass"]) {
      if (sum(specified) > 1) {
        message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
        weights <- match.strata <- NULL
      }
      method <- "subclassification"
      #weights <- rep(1, nrow(covs))
    }
    else if (specified["weights"]) {
      if (sum(specified) > 1) {
        message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
        match.strata <- subclass <- NULL
      }
      else {
        message("Assuming \"weighting\". If not, specify with an argument to method.")
      }
      method <- "weighting"
    }
    else {
      method <- "matching"
    }
  }
  else if (length(method) == 1) {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
    if (specified.method == "weighting") {
      if (specified["weights"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
          match.strata <- subclass <- NULL
        }
        method <- "weighting"
      }
      else if (specified["match.strata"]) {
        message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
        subclass <- NULL
        method <- "matching"
      }
      else if (specified["subclass"]) {
        message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
        method <- "subclassification"
        #weights <- rep(1, nrow(covs))
      }
      else {
        method <- "matching"
      }
    }
    else if (specified.method == "matching") {
      if (specified["match.strata"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
          weights <- subclass <- NULL
        }
        method <- "matching"
      }
      else if (specified["weights"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
          match.strata <- subclass <- NULL
        }
        method <- "matching"
      }
      else if (specified["subclass"]) {
        message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
        method <- "subclassification"
        #weights <- rep(1, nrow(covs))
      }
      else {
        method <- "matching"
      }
    }
    else if (specified.method == "subclassification") {
      if (specified["subclass"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
          weights <- match.strata <- NULL
        }
        method <- "subclassification"
        #weights <- rep(1, nrow(covs))
      }
      else if (specified["match.strata"]) {
        message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
        weights <- NULL
        method <- "matching"
      }
      else if (specified["weights"]) {
        message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
        method <- "weighting"
      }
    }
  }
  else {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    if (any(specified.method == "subclassification") || specified["subclass"]) {
      stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
    }
    else if (specified["match.strata"]) {
      stop("Only weights can be specified with mutiple methods.", call. = FALSE)
    }
    else if (!specified["weights"]) {
      warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
      method <- "matching"
    }
    else {
      #Matching and/or weighting with various weights
      method <- specified.method
      match.strata <- subclass <- NULL
    }
  }
  
  if (is_not_null(cluster) && !is_(cluster, c("character", "numeric", "factor"))) {
    stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
  }
  if (is_not_null(imp) && !is_(imp, c("character", "numeric", "factor"))) {
    stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
  }
  
  #Process treat
  if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
  
  
  if (is.character(treat) && length(treat)==1 && treat %in% names(data)) {
    treat <- data[[treat]]
  }
  else if (is_(treat, c("numeric", "logical", "factor")) || (is.character(treat) && length(treat) > 1)) {
    treat <- treat
  }
  else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
  
  if (sum(is.na(treat)) > 0)
    stop("Missing values exist in treat.", call. = FALSE)
  
  if (is_binary(treat)) {
    treat <- binarize(treat)
  }
  else if (is.character(treat)) {
    treat <- factor(treat)
  }
  
  #Process weights, addl, and distance
  for (i in c("weights", "addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data))
  }
  
  #Process subclass
  if (is_not_null(subclass)) {
    if (is_(subclass, c("numeric", "factor")) || (is.character(subclass) && length(subclass) > 1)) {
      subclass <- factor(subclass)
    }
    else if (is.character(subclass) && length(subclass)==1 && any(names(data) == subclass)) {
      subclass <- factor(data[[subclass]])
    }
    else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
  }
  
  #Process match.strata
  if (is_not_null(match.strata)) {
    if (is.character(match.strata) && length(match.strata)==1) {
      if (any(names(data) == match.strata)) {
        match.strata <- data[[match.strata]]
      }
      else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
    }
  }
  
  #Process sampling weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (any(names(data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
    if (anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  }
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is_(cluster, c("numeric", "factor")) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Process imp
  if (is_not_null(imp)) {
    if (is_(imp, c("numeric", "factor")) || (is.character(imp) && length(imp) > 1)) {
      imp <- imp
    }
    else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
      imp <- data[[imp]]
    }
    else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
    
    imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
    
    if (all_the_same(imp.lengths)) { #all the same
      for (i in vectors) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) { 
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            
            temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, temp.merge[[4]][order(temp.merge[[3]])])
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
      for (i in data.frames) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) {
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
    }
    else {
      problematic <- lengths > 0 & lengths != length(imp)
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
    }
    else ensure.equal.lengths <- FALSE
  }
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
  }
  
  #Turn match.strata into weights
  if (is_not_null(match.strata)) {
    weights <- data.frame(weights = match.strata2weights(match.strata = match.strata,
                                                         treat = treat,
                                                         covs = covs
    ))
  }
  
  if (is_not_null(weights)) {
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) !is.numeric(x), logical(1L)))) stop("All weights must be numeric.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    
    if (length(method) == 1) {
      method <- rep(method, ncol(weights))
    }
    else if (length(method) != ncol(weights)) {
      stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
  }
  
  #Check focal
  if (is_not_null(focal) && is.factor(treat)) {
    if (is.numeric(focal)) {
      if (focal <= nunique(treat)) focal <- levels(treat)[focal]
      else 
        stop(paste0("focal was specified as ", focal, 
                    ", but there are only ", levels(treat), " treatment groups."), call. = FALSE)
    }
    else {
      if (!any(levels(treat) == focal)) 
        stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
    }
  }
  
  #Get s.d.denom
  if (is_binary(treat) || !is.numeric(treat)) { #non-continuous
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method = method)
  }
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$method <- method
  X$covs <- covs
  X$weights <- weights
  X$treat <- treat
  X$distance <- distance
  X$subclass <- subclass
  X$cluster <- factor(cluster)
  X$call <- NULL
  X$addl <- addl
  X$imp <- factor(imp)
  X$s.weights <- s.weights
  X$focal <- focal
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.CBPS <- function(cbps.fit, ...) {
  #s.d.denom
  #cluster
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  #Checks
  
  treat <- model.response(model.frame(cbps.fit$terms, cbps.fit$data))
  covs <- cbps.fit$data[names(cbps.fit$data) %in% attributes(terms(cbps.fit))$term.labels]
  data <- A$data
  s.weights <- A$s.weights
  subset <- A$subset
  weights <- data.frame(weights = get.w(cbps.fit, use.weights = A$use.weights))
  cluster <- A$cluster
  s.d.denom <- A$s.d.denom
  estimand <- A$estimand
  method <- "weighting"
  c.data <- cbps.fit$data
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  #Process sampling weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (any(names(data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else if (any(names(c.data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
  }
  else s.weights <- rep(1, length(treat))
  
  weights <- weights/s.weights #Because CBPS weights contain s.weights in them
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(c.data) == cluster)) {
        cluster <- c.data[[cluster]]
      }
    }
    else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, c.data))
  }
  
  if (!any(class(cbps.fit) == "CBPSContinuous") && !is_binary(treat)) {
    #Multinomial
    # model.ps <- setNames(data.frame(cbps.fit$fitted.values), levels(treat))
    # model.ps.combined <- numeric(length(treat))
    # for (i in levels(treat)) {
    #     model.ps.combined[treat == i] <- model.ps[treat == i, i]
    # }
    # if (is_not_null(distance)) distance <- cbind(distance, prop.score = model.ps.combined)
    # else distance <- data.frame(prop.score = model.ps.combined)
  }
  else {
    if (all_the_same(cbps.fit$fitted.values)) {
      if (is_null(distance)) distance <- NULL
    }
    else {
      if (is_null(distance)) distance <- data.frame(prop.score = cbps.fit$fitted.values)
      else distance <- cbind(distance, prop.score = cbps.fit$fitted.values)
    }
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("cluster", "s.weights", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in names(lengths)[names(lengths) != "covs"]) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to CBPS()."), call. = FALSE)
  }
  
  #Get s.d.denom
  if (!any(class(cbps.fit) == "CBPSContinuous")) {
    if (is_binary(treat)) {
      X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
    }
    else {
      X$s.d.denom <- "pooled"
    }
  }
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$distance <- distance
  X$addl <- addl
  X$weights <- weights
  X$treat <- treat
  X$covs <- covs
  X$cluster <- factor(cluster)
  X$call <- cbps.fit$call
  X$s.weights <- s.weights
  X$method <- method
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.ebalance <- function(ebalance, ...) {
  #formula
  #data
  #treat
  #covs
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Get treat and covs
  data <- A$data
  t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
  
  #Initializing variables
  
  treat <- t.c[["treat"]]
  covs  <- t.c[["covs"]]
  cluster <- A$cluster
  subset <- A$subset
  weights <- data.frame(weights = get.w(ebalance, treat))
  method <- "weighting"
  s.d.denom <- A$s.d.denom
  estimand <- "ATT"
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data))
  }
  
  ensure.equal.lengths <- TRUE
  covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
  vectors <- c("treat", "cluster", "subset")
  data.frames <- c(covs.data, "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in names(lengths)[names(lengths) != "weights"]) {
      if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to ebalance()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$treat <- treat
  X$weights <- weights
  X$covs <- covs
  X$distance <- distance
  X$addl <- addl
  X$call <- NULL
  X$method <- method
  X$cluster <- factor(cluster)
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
}
x2base.ebalance.trim <- x2base.ebalance
x2base.optmatch <- function(optmatch, ...) {
  #formula
  #data
  #treat
  #covs
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Get treat and covs
  data <- A$data
  t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
  
  #Initializing variables
  treat <- binarize(t.c$treat)
  covs  <- t.c$covs
  distance <- A$distance
  subset <- A$subset
  cluster <- A$cluster
  s.d.denom <- A$s.d.denom
  estimand <- "ATT"
  method <- "matching"
  
  #Process match.strata (optmatch)
  if (length(optmatch) != length(treat) || length(optmatch) != nrow(covs)) {
    stop(paste0("The optmatch object must have the same length as ", ifelse(attr(t.c, "which")=="fd", "data", "covs"), "."), call. = FALSE)
  }
  a <- attributes(optmatch)
  
  d.reordered <- setNames(seq_len(nrow(covs)), rownames(covs))[a$names]
  
  weights <- data.frame(weights = get.w(optmatch))
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data))
  }    
  ensure.equal.lengths <- TRUE
  covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
  vectors <- c("treat", "cluster", "subset")
  data.frames <- c(covs.data, "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in names(lengths)[names(lengths) != "weights"]) {
      if (lengths[i] > 0 && lengths[i] != lengths["weights"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to optmatch()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$treat <- treat[d.reordered]
  X$distance <- distance[d.reordered, , drop = FALSE]
  X$covs <- covs[d.reordered, , drop = FALSE]
  X$weights <- weights
  X$addl <- addl[d.reordered, , drop = FALSE]
  X$call <- attr(optmatch, "call")
  X$method <- "matching"
  X$cluster <- factor(cluster[d.reordered])
  
  
  subset <- subset[d.reordered]
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
  
}
x2base.weightit <- function(weightit, ...) {
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  estimand <- weightit$estimand
  weights <- data.frame(weights = get.w(weightit))
  treat <- weightit$treat
  covs <- weightit$covs
  s.weights <- weightit$s.weights
  data <- A$data
  cluster <- A$cluster
  imp <- A$imp
  subset <- A$subset
  s.d.denom <- A$s.d.denom
  focal <- weightit$focal
  method <- rep("weighting", ncol(weights))
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  d.e.in.w <- vapply(c("covs", "exact", "by"), function(x) is_not_null(weightit[[x]]), logical(1L))
  if (any(d.e.in.w)) weightit.data <- do.call("cbind", unname(weightit[c("covs", "exact", "by")[d.e.in.w]]))
  else weightit.data <- NULL
  
  if (has.treat.type(treat)) {
    treat.type <- get.treat.type(treat)
  }
  else if (is_not_null(weightit$treat.type)) {
    treat.type <- weightit$treat.type
  }
  else {
    treat.type <- get.treat.type(assign.treat.type(treat))
  }
  
  if (is_null(covs)) stop("No covariates were specified in the weightit object.", call. = FALSE)
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
      stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(weightit.data) == cluster)) {
        cluster <- weightit.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process imp
  if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
    stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, weightit.data))
  }    
  
  if (is_not_null(distance)) distance <- cbind(distance, prop.score = weightit$ps)
  else if (is_not_null(weightit$ps)) distance <- data.frame(prop.score = weightit$ps)
  else distance <- NULL
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Process imp further
  if (is_not_null(imp)) {
    if (is.numeric(imp) || is.factor(imp) || (is.character(imp) && length(imp)>1)) {
      imp <- imp
    }
    else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
      imp <- data[[imp]]
    }
    else if (is.character(imp) && length(imp)==1 && any(names(weightit.data) == imp)) {
      imp <- weightit.data[[imp]]
    }
    else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
    
    imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
    
    if (all_the_same(imp.lengths)) { #all the same
      for (i in vectors) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) { 
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            
            temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, temp.merge[[4]][order(temp.merge[[3]])])
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
      for (i in data.frames) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) {
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
    }
    else {
      problematic <- lengths > 0 & lengths != length(imp)
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
    }
    else ensure.equal.lengths <- FALSE
  }
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
  }
  
  #Get s.d.denom
  if (treat.type != "continuous") {
    X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
  }
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$weights <- weights
  X$treat <- treat
  X$distance <- distance
  X$addl <- addl
  X$covs <- covs
  X$cluster <- factor(cluster)
  X$method <- method
  X$imp <- factor(imp)
  X$s.weights <- weightit$s.weights
  X$discarded <- weightit$discarded
  X$focal <- focal
  X$call <- weightit$call
  
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.designmatch <- function(dm, ...) {
  #formula
  #data
  #treat
  #covs
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Get treat and covs
  data <- A$data
  
  if (all(c("id_1", "id_2") %in% names(dm))) {
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs, needs.treat = FALSE)
    covs  <- t.c$covs
    
    if (anyDuplicated(c(dm[["id_1"]], dm[["id_2"]])) != 0) {
      stop("Some units are used more than once. Balance cannot be checked.", call. = FALSE)
    }
    treat <- rep(NA_real_, nrow(covs))
    treat[dm[["id_1"]]] <- 1
    treat[dm[["id_2"]]] <- 0
    
    in.matched <- !is.na(treat)
    weights <- NULL
  }
  else {
    t.c <- use.tc.fd(A$formula, data, A$treat, A$covs)
    covs  <- t.c$covs
    treat <- binarize(t.c$treat)
    in.matched <- rep(TRUE, length(treat))
    
    weights <- data.frame(weights = get.w.designmatch(dm, treat))
    
  }
  
  #Initializing variables
  distance <- A$distance
  subset <- A$subset
  cluster <- A$cluster
  estimand <- A$estimand
  s.d.denom <- A$s.d.denom
  method <- "matching"
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster) > 1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data))
  }    
  ensure.equal.lengths <- TRUE
  covs.data <- ifelse(attr(t.c, "which")=="fd", "data", "covs")
  vectors <- c("treat", "cluster", "subset")
  data.frames <- c(covs.data, "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in names(lengths)[names(lengths) != "treat"]) {
      if (lengths[i] > 0 && lengths[i] != lengths["treat"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original call to designmatch()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method)
  
  if (any(c(anyNA(covs), anyNA(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  
  X$treat <- treat
  X$distance <- distance
  X$covs <- covs
  X$weights <- weights
  X$addl <- addl
  X$call <- NULL
  X$method <- method
  X$cluster <- factor(cluster)
  
  X <- subset_X(X, in.matched)
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "binary"
  
  return(X)
  
}
.x2base.mimids <- function(mimids, ...) {
  
  nimp <- length(mimids[[2]]) - 1
  if (nimp == 1) {
    X <- x2base.matchit(mimids[[2]][-1][[1]])
  }
  else {
    Xs <- vector("list", nimp)
    for (i in 1:nimp) {
      Xs[[i]] <- x2base.matchit(mimids[[2]][-1][[i]])
      Xs[[i]][["imp"]] <- factor(rep(i, length(Xs[[i]][["treat"]])),
                                 levels = as.character(seq_along(Xs)))
    }
    n <- length(Xs[[1]][["treat"]])
    
    X <- setNames(vector("list", length(Xs[[1]])),
                  names(Xs[[1]]))
    for (x in names(X)) {
      if (is.data.frame(Xs[[1]][[x]]) || is.matrix(Xs[[1]][[x]])) {
        X[[x]] <- do.call("rbind", lapply(Xs, function(x_) x_[[x]]))
      }
      else if (length(Xs[[1]][[x]]) == n) {
        X[[x]] <- unlist(lapply(Xs, function(x_) x_[[x]]))
      }
      else {
        X[[x]] <- Xs[[1]][[x]]
      }
    }
    
    class(X) <- "imp"
  }
  
  if ("wimids" %in% class(mimids)) {
    X$weights[] <- unlist(lapply(mimids[[4]][-1], function(i) i[["inverse.weights"]]))
    X$s.d.denom <- "pooled"
    X$method <- "weighting"
  }
  
  return(X)
}
x2base.mimids <- function(mimids, ...) {
  
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  nimp <- length(mimids[[4]]) - 1
  
  if (any(class(mimids) == "wimids")) {
    subclass <- NULL
    method <- "weighting"
    estimand <- "ATE"
    weights <- do.call("rbind", lapply(mimids[[4]][-1], function(x) x["inverse.weights"]))
    m.distance <- data.frame(distance = unlist(lapply(mimids[[2]][-1], function(m) m[["distance"]])))
    if ("average.distance" %in% names(mimids[[4]][-1][[1]])) {
      m.distance$average.distance <- rep(rowMeans(do.call("cbind", lapply(mimids[[2]][-1], function(m) m[["distance"]]))), nimp)
    }
    discarded <- NULL
  } else {
    subclass <- NULL
    method <- "matching"
    estimand <- "ATT"
    weights <- do.call("rbind", lapply(mimids[[2]][-1], function(x) data.frame(weights = x[["weights"]])))
    m.distance <- data.frame(distance = unlist(lapply(mimids[[2]][-1], function(m) m[["distance"]])))
    if ("average.distance" %in% names(mimids[[4]][-1][[1]])) names(m.distance) <- "ave.distance"
    discarded <- if ("discarded" %in% names(mimids[[2]][-1][[1]])) unlist(lapply(mimids[[2]][-1], function(x) x[["discarded"]])) else NULL
  }
  
  #Process data
  data <- A$data
  
  if (is_not_null(data)) {
    if (is_(data, "mids")) {
      data <- imp.complete(data)
      if ("imp" %nin% names(A)) imp <- data[[".imp"]]
    }
    else if (!is_(data, "data.frame"))
    {
      warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
      data <- NULL
    }
  }
  
  m.data1 <- do.call("rbind", lapply(mimids[[2]][-1], function(m) m$model$data))
  m.data2 <- do.call("rbind", mimids[[4]][-1])
  
  subset <- A$subset
  cluster <- A$cluster
  s.d.denom <- A$s.d.denom
  imp <- m.data2[[".imp"]]
  weights[is.na(weights),] <- 0
  
  treat <- unlist(lapply(mimids[[2]][-1], function(m) m[["treat"]]))
  covs <- do.call("rbind", lapply(1:nimp, function(i) {
    m <- mimids[[2]][-1][[i]]
    if (is_not_null(m$model$model)) {
      if (nrow(m$model$model) == length(m$treat)) {
        covs <- data.frame(m$model$model[, names(m$model$model) %in% attributes(terms(m$model))$term.labels])
      }
      else {
        #Recreating covs from model object and m$X. Have to do this because when 
        #drop != NULL and reestimate = TRUE, cases are lost. This recovers them.
        
        order <- setNames(attr(m$model$terms, "order"),
                          attr(m$model$terms, "term.labels"))
        assign <- setNames(attr(m$X, "assign"), colnames(m$X))
        assign1 <- assign[assign %in% which(order == 1)] #Just main effects
        
        dataClasses <- attr(m$model$terms, "dataClasses")
        factors.to.unsplit <- names(dataClasses)[dataClasses %in% c("factor", "character", "logical")]
        f0 <- setNames(lapply(factors.to.unsplit, 
                              function(x) {
                                if (dataClasses[x] == "factor")
                                  list(levels = levels(m$model$model[[x]]),
                                       faclev = paste0(x, levels(m$model$model[[x]])))
                                else 
                                  list(levels = unique(m$model$model[[x]]),
                                       faclev = paste0(x, unique(m$model$model[[x]])))
                              }),
                       factors.to.unsplit)
        covs <- as.data.frame(m$X[, names(assign1)])
        for (i in factors.to.unsplit) {
          covs <- unsplitfactor(covs, i, sep = "",
                                dropped.level = f0[[i]]$levels[f0[[i]]$faclev %nin% colnames(m$X)])
          if (dataClasses[i] == "logical") covs[[i]] <- as.logical(covs[[i]])
        }
      }
      
    }
    else if (!anyNA(mimids[[4]][-1][[i]])) {
      t.c <- get.covs.and.treat.from.formula(m$formula, data = mimids[[4]][-1][[i]])
      covs <- t.c[["reported.covs"]]
    }
    else if (is_not_null(data)) {
      if (nrow(data) == length(mimids[[2]][-1][[i]][["treat"]])) {
        t.c <- get.covs.and.treat.from.formula(m$formula, data = data)
        covs <- t.c[["reported.covs"]]
      }
      else if (nrow(data) == length(treat)) {
        t.c <- get.covs.and.treat.from.formula(m$formula, 
                                               data = data[((i-1)*length(mimids[[2]][-1][[i]][["treat"]]) + 1):(i*length(mimids[[2]][-1][[i]][["treat"]])), , drop = FALSE])
        covs <- t.c[["reported.covs"]]
      }
      else {
        covs <- data.frame(m$X)
      }
    }
    else {
      covs <- data.frame(m$X)
    }
    return(covs)
  }))
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(m.data1) == cluster)) {
        cluster <- m.data1[[cluster]]
      }
    }
    else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is_(subset, "logical")) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl", "distance")) {
    assign(i, data.frame.process(i, A[[i]], treat, covs, data, m.data1, m.data2))
  }
  
  if (is_not_null(distance)) distance <- cbind(distance, m.distance)
  else distance <- m.distance
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal = NULL, method = method)
  
  ensure.equal.lengths <- TRUE
  vectors <- c("treat", "cluster", "subset")
  data.frames <- c("covs", "weights", "distance", "addl")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L))), c(vectors, data.frames))
  #Process imp further
  if (is_not_null(imp)) {
    
    imp.lengths <- vapply(unique(imp, nmax = nimp), function(i) sum(imp == i), numeric(1L))
    
    if (all_the_same(imp.lengths)) { #all the same
      for (i in vectors) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) { 
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            
            temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, temp.merge[[4]][order(temp.merge[[3]])])
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
      for (i in data.frames) {
        if (lengths[i] > 0 && lengths[i] != length(imp)) {
          if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
          if (lengths[i] == imp.lengths[1]) {
            temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                   order2 = seq_along(imp))
            temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                   get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
            )
            temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                by.y = 1:2, sort = FALSE)
            assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
          }
          else {
            problematic[i] <- TRUE
          }
        }
      }
    }
    else {
      problematic <- lengths > 0 & lengths != length(imp)
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
    }
    else ensure.equal.lengths <- FALSE
  }
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames[data.frames!="covs"])) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
  }
  
  if (any(c(is.na(covs), is.na(addl)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$method <- method
  X$treat <- treat
  X$weights <- weights
  X$discarded <- discarded
  X$covs <- covs
  X$distance <- distance
  X$addl <- addl
  X$cluster <- factor(cluster)
  X$imp <- factor(imp)
  X$call <- NULL
  X$subclass <- factor(subclass)
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- "imp"
  
  return(X)
}
x2base.wimids <- x2base.mimids

#MSMs wth multiple time points
x2base.iptw <- function(iptw, ...) {
  A <- list(...)
  
  X.names <- c("covs.list",
               "treat.list",
               "weights",
               "distance.list",
               "addl.list",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  if (is_not_null(A) && names(A)[1]=="" && is_null(A$stop.method)) A$stop.method <- A[[1]] #for bal.plot
  if (is_null(A$stop.method) && is_not_null(A$full.stop.method)) A$stop.method <- A$full.stop.method
  available.stop.methods <- names(iptw$psList[[1]]$ps)
  if (is_not_null(A$stop.method)) {
    if (any(is.character(A$stop.method))) {
      rule1 <- available.stop.methods[vapply(available.stop.methods, function(x) any(startsWith(tolower(x), tolower(A$stop.method))), logical(1L))]
      if (is_null(rule1)) {
        message(paste0("Warning: stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead."))
        rule1 <- available.stop.methods
      }
    }
    else if (is.numeric(A$stop.method) && any(A$stop.method %in% seq_along(available.stop.methods))) {
      if (any(!A$stop.method %in% seq_along(available.stop.methods))) {
        message(paste0("Warning: There are ", length(available.stop.methods), " stop methods available, but you requested ", 
                       word_list(A$stop.method[!A$stop.method %in% seq_along(available.stop.methods)], and.or = "and"),"."))
      }
      rule1 <- available.stop.methods[A$stop.method %in% seq_along(available.stop.methods)]
    }
    else {
      warning("stop.method should be ", word_list(available.stop.methods, and.or = "or", quotes = TRUE), ".\nUsing all available stop methods instead.", call. = FALSE)
      rule1 <- available.stop.methods
    }
  }
  else {
    rule1 <- available.stop.methods
  }
  
  s <- available.stop.methods[match(tolower(rule1), tolower(available.stop.methods))]
  estimand <- substr(toupper(s), nchar(s)-2, nchar(s))
  
  weights <- data.frame(get.w(iptw, s))
  treat.list <- lapply(iptw$psList, function(x) x$treat)
  covs.list <- lapply(iptw$psList, function(x) x$data[x$gbm.obj$var.names])
  subset <- A$subset
  data <- A$data
  cluster <- A$cluster
  ps.data <- iptw$psList[[1]]$data
  s.weights <- iptw$psList[[1]]$sampw
  ntimes <- iptw$nFits
  s.d.denom <- A$s.d.denom
  method <- rep("weighting", ncol(weights))
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  #Order covs.list
  all.covs <- unique(unlist(lapply(covs.list, names)))
  covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(ps.data) == cluster)) {
        cluster <- ps.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl.list", "distance.list")) {
    assign(i, list.process(i, A[[i]], ntimes, 
                           "the original call to iptw()",
                           treat.list,
                           covs.list,
                           data,
                           ps.data))
  }
  
  
  if (is_not_null(distance.list)) {
    for (ti in seq_along(distance.list)) {
      if (length(s) == 1) {
        distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw$psList[[ti]]$ps[[s]])
      }
      else {
        distance.list[[ti]] <- data.frame(distance[[ti]], prop.score = iptw$psList[[ti]]$ps[s])
      }
    }
    
  }
  else {
    distance.list <- vector("list", ntimes)
    for (ti in seq_along(distance.list)) {
      if (length(s) == 1) {
        distance.list[[ti]] <- data.frame(prop.score = iptw$psList[[ti]]$ps[[s]])
      }
      else {
        distance.list[[ti]] <- data.frame(prop.score = iptw$psList[[ti]]$ps[s])
      }
    }
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("weights")
  lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L)),
                        vapply(lists, function(x) {
                          if (is.null(get0(x))) 0 
                          else if (is.vector(get(x))) {
                            if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                            else max(lengths(get(x)))
                          }
                          else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                        }, numeric(1L))), c(vectors, data.frames, lists))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames, lists)) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to iptw()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat.list[[1]], focal = NULL, method)
  
  if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$weights <- weights
  X$treat.list <- lapply(treat.list, function(x) x)
  X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
  X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
  X$covs.list <- lapply(covs.list, function(x) x)
  X$call <- NULL
  X$cluster <- factor(cluster)
  X$method <- rep("weighting", ncol(weights))
  X$s.weights <- s.weights
  
  if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.data.frame.list <- function(covs.list, ...) {
  A <- list(...)
  X.names <- c("covs.list",
               "treat.list",
               "weights",
               "distance.list",
               "addl.list",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  covs.list <- covs.list
  treat.list <- A$treat.list
  data <- A$data
  weights <- A$weights
  distance.list <- A$distance.list
  cluster <- A$cluster
  addl.list <- A$addl.list
  s.d.denom <- A$s.d.denom
  method <- A$method
  imp <- A$imp
  s.weights <- A$s.weights
  subset <- A$subset
  focal <- A$focal
  ntimes <- length(covs.list)
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (is_not_null(s.weights) && any(vapply(s.weights, anyNA, logical(1L)))) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  #Checks
  if (is_null(covs.list)) {
    stop("covs.list must be specified.", call. = FALSE)
  }
  if (!is.list(covs.list)) {
    stop("covs.list must be a list of covariates for which balanced is to be assessed at each time point.", call. = FALSE)
  }
  if (any(!vapply(covs.list, is.data.frame, logical(1L)))) {
    stop("Each item in covs.list must be a data frame.", call. = FALSE)
  }
  
  if (is_not_null(data) && !is.data.frame(data)) {
    warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
    data <- NULL
  }
  
  specified <- setNames(rep(FALSE, 1), "weights")
  if (is_not_null(weights)) {
    if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
      stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
    }
    specified["weights"] <- TRUE
  }
  
  #Getting method
  if (is_null(method)) {
    if (specified["weights"]) {
      
      #message("Assuming \"weighting\". If not, specify with an argument to method.")
      
      method <- "weighting"
    }
    else {
      method <- "matching"
    }
  }
  else if (length(method) == 1) {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
    if (specified.method == "weighting") {
      if (specified["weights"]) {
        method <- "weighting"
      }
      else {
        method <- "matching"
      }
    }
    else {
      if (specified["weights"]) {
        warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
        method <- "matching"
      }
      else {
        method <- "matching"
      }
    }
  }
  else {
    specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
    if (any(specified.method == "subclassification") || specified["subclass"]) {
      warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
      method <- "matching"
    }
    else if (specified["match.strata"]) {
      warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
      method <- "matching"
    }
    else if (!specified["weights"]) {
      warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
      method <- "matching"
    }
    else {
      #Matching and/or weighting with various weights
      method <- specified.method
    }
  }
  
  if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
    stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
  }
  if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
    stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
  }
  
  #Order covs.list
  all.covs <- unique(unlist(lapply(covs.list, names)))
  covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
  
  #Process treat
  if (is_null(treat.list)) stop("treat.list must be specified.", call. = FALSE)
  if (!is.vector(treat.list)) {
    treat.list <- as.list(treat.list)
  }
  if (length(treat.list) != length(covs.list)) {
    stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
  }
  
  for (ti in seq_along(treat.list)) {
    if (is.numeric(treat.list[[ti]]) || is.factor(treat.list[[ti]]) || (is.character(treat.list[[ti]]) && length(treat.list[[ti]]) > 1)) {
      #treat.list[[ti]] <- treat.list[[ti]]
    }
    else if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1 && any(names(data) == treat.list[[ti]])) {
      names(treat.list)[ti] <- treat.list[[ti]]
      treat.list[[ti]] <- data[[treat.list[[ti]]]]
    }
    else stop("Each item in treat.list must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat.list[[ti]])) > 0)
      stop("Missing values exist in treat.list", call. = FALSE)
    
    if (is_binary(treat.list[[ti]])) {
      treat.list[[ti]] <- binarize(treat.list[[ti]])
    }
    else if (is.character(treat.list[[ti]])) {
      treat.list[[ti]] <- factor(treat.list[[ti]])
    }
  }
  
  #Process weights
  for (i in c("weights")) {
    assign(i, data.frame.process(i, A[[i]], do.call("cbind", treat.list), do.call("cbind", covs.list), data))
  }
  
  #Process addl and distance
  for (i in c("addl.list", "distance.list")) {
    assign(i, list.process(i, A[[i]], ntimes, 
                           "covs.list",
                           treat.list,
                           covs.list,
                           data
    ))
  }
  
  #Process sampling weights
  if (is_not_null(s.weights)) {
    if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
      stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
    }
    if (is.character(s.weights) && length(s.weights)==1) {
      if (any(names(data) == s.weights)) {
        s.weights <- data[[s.weights]]
      }
      else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
    }
    if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
  }
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
      cluster <- data[[cluster]]
    }
    else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("weights")
  lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is_null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L)),
                        vapply(lists, function(x) {
                          if (is_null(get0(x))) 0 
                          else if (is.vector(get(x))) {
                            if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                            else max(lengths(get(x)))
                          }
                          else max(vapply(get(x), function(y) if (is_not_null(y)) {if (is_not_null(nrow(y))) nrow(y) else length(y)} else 0, numeric(1L)))
                        }, numeric(1L))), c(vectors, data.frames, lists))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames, lists)) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs.list."), call. = FALSE)
  }
  
  if (is_not_null(weights)) {
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (any(vapply(weights, function(x) any(!is.finite(x)), logical(1L)))) stop("All weights must be numeric.", call. = FALSE)
    if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
    
    if (length(method) == 1) {
      method <- rep(method, ncol(weights))
    }
    else if (length(method) != ncol(weights)) {
      stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
    }
    
  }
  
  #Get s.d.denom
  X$s.d.denom <- rep("pooled", max(1, ncol(weights)))
  
  if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  if (is_null(s.weights)) s.weights <- rep(1, length(treat.list[[1]]))
  
  X$method <- method
  X$covs.list <- lapply(covs.list, function(x) x)
  X$weights <- weights
  X$treat.list <- lapply(treat.list, function(x) x)
  X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
  X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
  X$cluster <- factor(cluster)
  X$call <- NULL
  X$imp <- factor(imp)
  X$s.weights <- s.weights
  
  if (is_null(subset)) subset <- rep(TRUE, length(treat.list[[1]]))
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.formula.list <- function(formula.list, ...) {
  A <- list(...)
  A[["covs.list"]] <- NULL
  A[["treat.list"]] <- NULL
  
  treat.list <- covs.list <- vector("list", length(formula.list))
  for (i in seq_along(formula.list)) {
    t.c <- get.covs.and.treat.from.formula(formula.list[[i]], A[["data"]])
    covs.list[[i]] <- t.c[["reported.covs"]]
    treat.list[[i]] <- t.c[["treat"]]
    names(treat.list)[i] <- t.c[["treat.name"]]
  }
  
  X <- do.call("x2base.data.frame.list", c(list(covs.list, treat.list = treat.list), A))
  return(X)
}
x2base.CBMSM <- function(cbmsm, ...) {
  A <- list(...)
  
  X.names <- c("covs.list",
               "treat.list",
               "weights",
               "distance.list",
               "addl.list",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  ID <- sort(unique(cbmsm$id))
  times <- sort(unique(cbmsm$time))
  treat.list <- lapply(times, function(x) cbmsm$treat.hist[ID, x]) 
  covs <- cbmsm$data[names(cbmsm$data) %in% attributes(terms(cbmsm$model))$term.labels]
  weights <- data.frame(weights = get.w(cbmsm)[ID])
  ntimes <- length(times)
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  
  covs.list <- vector("list", ntimes)
  for (ti in times) {
    if (ti == 1) {
      covs.list[[ti]] <- setNames(data.frame(covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                                  paste0(names(covs), "0"))
    }
    else {
      covs.list[[ti]] <- cbind(covs.list[[ti - 1]], 
                               setNames(data.frame(cbmsm$y[cbmsm$time == ti - 1][ID], 
                                                   covs[cbmsm$time == ti, , drop = FALSE][ID, , drop = FALSE]),
                                        c(paste0("treat", ti - 1), paste0(names(covs), ti))))
    }
  }
  
  cluster <- A$cluster
  subset <- A$subset
  data <- A$data
  
  cbmsm.data <- cbmsm$data[cbmsm$time == 1, , drop = FALSE][ID, , drop = FALSE]
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
      stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(cbmsm.data) == cluster)) {
        cluster <- cbmsm.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl.list", "distance.list")) {
    assign(i, list.process(i, A[[i]], ntimes, 
                           "the original call to CBMSM()",
                           treat.list,
                           covs.list,
                           data,
                           cbmsm.data))
  }
  
  if (is_not_null(distance.list)) distance.list <- lapply(times, function(x) data.frame(distance.list[[x]], prop.score = cbmsm$fitted.values))
  else if (is_not_null(cbmsm$fitted.values)) distance.list <- lapply(times, function(x) data.frame(prop.score = cbmsm$fitted.values))
  else distance.list <- NULL
  
  ensure.equal.lengths <- TRUE
  vectors <- c("cluster", "subset")
  data.frames <- c("weights")
  lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L)),
                        vapply(lists, function(x) {
                          if (is.null(get0(x))) 0 
                          else if (is.vector(get(x))) {
                            if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                            else max(lengths(get(x)))
                          }
                          else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                        }, numeric(1L))), c(vectors, data.frames, lists))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames, lists)) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to CBMSM()."), call. = FALSE)
  }
  
  if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$weights <- weights
  X$treat.list <- lapply(treat.list, function(x) x)
  X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
  X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
  X$covs.list <- lapply(covs.list, function(x) x)
  X$call <- cbmsm$call
  X$cluster <- factor(cluster)
  X$method <- rep("weighting", ncol(weights))
  X$s.weights <- NULL
  X$s.d.denom <- "pooled"
  
  if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}
x2base.weightitMSM <- function(weightitMSM, ...) {
  A <- list(...)
  
  X.names <- c("covs.list",
               "treat.list",
               "weights",
               "distance.list",
               "addl.list",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass")
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  #Initializing variables
  estimand <- weightitMSM$estimand
  weights <- data.frame(weights = get.w(weightitMSM))
  treat.list <- weightitMSM$treat.list
  covs.list <- weightitMSM$covs.list
  if (is_null(covs.list)) stop("No covariates were specified in the weightit object.", call. = FALSE)
  covs <- do.call("cbind", covs.list)
  s.weights <- weightitMSM$s.weights
  data <- A$data
  cluster <- A$cluster
  imp <- A$imp
  subset <- A$subset
  ntimes <- length(treat.list)
  
  
  if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
  if (any(vapply(weights, function(x) any(x < 0), logical(1L)))) stop("Negative weights are not allowed.", call. = FALSE)
  if (is_not_null(s.weights) && anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
  
  weightitMSM.data <- weightitMSM$data
  d.e.in.w <- vapply(c("covs.list", "exact", "by"), function(x) is_not_null(weightitMSM[[x]]), logical(1L))
  if (any(d.e.in.w)) weightitMSM.data <- do.call("cbind", c(list(covs), weightitMSM[c("exact", "by")])[d.e.in.w])
  else weightitMSM.data <- NULL
  
  if (all(vapply(treat.list, has.treat.type, logical(1L)))) {
    treat.type <- vapply(treat.list, get.treat.type, character(1L))
  }
  else if (length(weightitMSM$treat.type) == length(treat.list)) {
    treat.type <- weightitMSM$treat.type
  }
  else {
    treat.type <- vapply(treat.list, function(treat) {
      get.treat.type(assign.treat.type(treat))
    }, character(1L))
  } 
  
  #Order covs.list
  all.covs <- unique(unlist(lapply(covs.list, names)))
  covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
  
  #Process cluster
  if (is_not_null(cluster)) {
    if (!is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
      stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
      cluster <- cluster
    }
    else if (is.character(cluster) && length(cluster)==1) {
      if (any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else if (any(names(weightitMSM.data) == cluster)) {
        cluster <- weightitMSM.data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in any given data set.", call. = FALSE)
    }
  }
  
  #Process subset
  if (is_not_null(subset)) {
    if (!is.logical(subset)) {
      stop("The argument to subset must be a logical vector.", call. = FALSE)
    }
    if (anyNA(subset)) {
      warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
      subset[is.na(subset)] <- FALSE
    }
  }
  
  #Process addl and distance
  for (i in c("addl.list", "distance.list")) {
    assign(i, list.process(i, A[[i]], ntimes, 
                           "the original call to weightitMSM()",
                           treat.list,
                           covs.list,
                           data,
                           weightitMSM.data))
  }
  
  if (is_not_null(distance.list)) distance.list <- lapply(seq_along(distance.list), function(x) data.frame(distance.list[[x]], prop.score = weightitMSM$ps.list[[x]]))
  else if (is_not_null(weightitMSM$ps.list)) distance.list <- lapply(seq_along(weightitMSM$ps.list), function(x) data.frame(prop.score = weightitMSM$ps.list[[x]]))
  else distance.list <- NULL
  
  ensure.equal.lengths <- TRUE
  vectors <- c("s.weights", "cluster", "subset")
  data.frames <- c("weights")
  lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
  problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
  lengths <- setNames(c(lengths(mget(vectors)), 
                        vapply(data.frames, 
                               function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                               }, numeric(1L)),
                        vapply(lists, function(x) {
                          if (is.null(get0(x))) 0 
                          else if (is.vector(get(x))) {
                            if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                            else max(lengths(get(x)))
                          }
                          else max(vapply(get(x), function(y) if (is_not_null(y)) nrow(y) else 0, numeric(1L)))
                        }, numeric(1L))), c(vectors, data.frames, lists))
  
  #Ensure all input lengths are the same.
  if (ensure.equal.lengths) {
    for (i in c(vectors, data.frames, lists)) {
      if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
        problematic[i] <- TRUE
      }
    }
  }
  if (any(problematic)) {
    stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as the original data set in the call to weightitMSM()."), call. = FALSE)
  }
  
  #Get s.d.denom
  X$s.d.denom <- rep("pooled", ncol(weights))
  
  if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
    warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
  }
  
  X$weights <- weights
  X$treat.list <- lapply(treat.list, function(x) x)
  X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
  X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
  X$covs.list <- lapply(covs.list, function(x) x)
  X$call <- NULL
  X$cluster <- factor(cluster)
  X$method <- rep("weighting", ncol(weights))
  X$s.weights <- s.weights
  
  if (is_null(subset)) subset <- rep(TRUE, lengths["treat.list"])
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}

x2base.default <- function(obj, ...) {
  
  A <- list(...)
  X.names <- c("covs",
               "treat",
               "weights",
               "distance",
               "addl",
               "s.d.denom",
               "call",
               "cluster",
               "imp",
               "s.weights",
               "focal",
               "discarded",
               "method",
               "subclass",
               "covs.list",
               "treat.list",
               "distance.list",
               "addl.list")
  
  X <- setNames(vector("list", length(X.names)),
                X.names)
  
  if (!is.list(obj)) stop("The input object must be an appropriate list, data.frame, formula, or the output of one of the supported packages.", call. = FALSE)
  
  Q <- list(treat = list(name = c("treat", "tr"), 
                         type = c("numeric", "character", "factor", "logical")),
            treat.list = list(name = c("treat.list", "treat", "tr"),
                              type = c("list", "data.frame")),
            covs = list(name = c("covs", "covariates"), 
                        type = c("data.frame")),
            covs.list = list(name = c("covs.list", "covs", "covariates"),
                             type = c("list")),
            formula = list(name = c("formula", "form"), 
                           type = c("formula")),
            formula.list = list(name = c("formula.list", "formula", "form"),
                                type = c("list")),
            data = list(name = c("data"),
                        type = c("data.frame", "mids")),
            weights = list(name = c("weights", "w", "wts"),
                           type = c("data.frame", "matrix", "numeric")),
            distance = list(name = c("distance", "ps", "pscore","p.score", "propensity.score"),
                            type = c("data.frame", "matrix", "numeric")),
            distance.list = list(name = c("distance.list", "ps.list", "distance", "ps"),
                                 type = c("list")),
            subclass = list(name = c("subclass", "strata"),
                            type = c("factor", "character", "numeric")),
            match.strata = list(name = c("match.strata"),
                                type = c("factor", "character", "numeric")),
            estimand = list(name = c("estimand", "target", "att", "ate"),
                            type = c("character", "numeric")),
            s.weights = list(name = c("s.weights", "sw", "sweights", "sampw"),
                             type = c("numeric")),
            focal = list(name = c("focal", "treatATT"), 
                         type = c("character", "numeric")),
            call = list(name = c("call"),
                        type = c("call")))
  
  P <- setNames(vector("list", length(Q)), names(Q))
  names(obj) <- tolower(names(obj))
  
  for (i in names(Q)) {
    if (is_null(A[[i]])) {
      for (j in Q[[i]][["name"]]) {
        if (is_null(P[[i]])) {
          if (is_not_null(obj[[j]])) {
            if (any(which.type <- vapply(Q[[i]][["type"]], function(x) is_(obj[[j]], x), logical(1L)))) {
              P[[i]] <- obj[[j]]
              attr(P[[i]], "name") <- j
              attr(P[[i]], "type") <- Q[[i]][["type"]][which.type]
            }
          }
        }
      }
      if (is_not_null(P[[i]])) {
        assign(i, P[[i]])
        A[[i]] <- P[[i]]
      }
    }
    assign(i, A[[i]])
  }
  
  msm <- FALSE
  
  #treat OK
  
  #treat.list
  if (is_not_null(treat.list)) {
    if (!all(sapply(treat.list, function(x) any(vapply(Q[["treat"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
      treat.list <- A[["treat.list"]]
    }
    msm <- TRUE
  }
  
  #covs 
  if (is_not_null(covs)) covs <- as.data.frame(covs)
  
  #covs.list
  if (is_not_null(covs.list)) {
    if (!all(sapply(covs.list, function(x) any(vapply(Q[["covs"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
      covs.list <- A[["covs.list"]]
    }
    msm <- TRUE
  }
  
  #formula OK
  
  #formula.list
  if (is_not_null(formula.list)) {
    if (!all(sapply(formula.list, function(x) any(vapply(Q[["formula"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
      formula.list <- A[["formula.list"]]
    }
    msm <- TRUE
  }
  
  #data
  if (is_not_null(data)) {
    if (is_(data, "mids")) {
      data <- imp.complete(data)
      if ("imp" %nin% names(A)) A[["imp"]] <- data[[".imp"]]
    }
    data <- as.data.frame(data)
  }
  
  #weights
  if (is_not_null(weights)) {
    if (is.vector(weights, "numeric")) weights <- data.frame(weights = weights)
    else weights <- as.data.frame(weights)
  }
  
  #distance
  if (is_not_null(distance)) {
    if (is.numeric(distance)) {
      if (is_not_null(attr(distance, "name"))) distance <- setNames(data.frame(distance),
                                                                    attr(distance, "name"))
      else distance <- data.frame(distance = distance)
    }
    else distance <- as.data.frame(distance)
  }
  
  #distance.list
  if (is_not_null(distance.list)) {
    if (!all(sapply(distance.list, function(x) any(vapply(Q[["distance"]][["type"]], function(c) is_(x, c), logical(1L)))))) {
      distance.list <- A[["distance.list"]]
    }
    #msm <- TRUE
  }
  
  #subclass
  if (is_not_null(subclass)) subclass <- factor(subclass)
  
  #match.strata
  if (is_not_null(match.strata)) match.strata <- factor(match.strata)
  
  #estimand
  if (is_not_null(estimand)) {
    estimand.name <- attr(estimand, "name")
    if (is_not_null(estimand.name) && toupper(estimand.name) == "ATT") {
      if (estimand == 0) estimand <- "ATE"
      else estimand <- "ATT"
    }
    else if (is_not_null(estimand.name) && toupper(estimand.name) == "ATE") {
      if (estimand == 0) estimand <- "ATT"
      else estimand <- "ATE"
    }
    else {
      if (tolower(estimand) %in% c("att", "treat", "treated", "tr", "t", "atet")) estimand <- "ATT"
      else if (tolower(estimand) %in% c("ate", "all")) estimand <- "ATE"
      else if (tolower(estimand) %in% c("atc", "control", "untreated", "u", "c", "ctrl", "atu", "atec", "ateu")) estimand <- "ATC"
      else estimand <- NULL
    }
  }
  
  #s.weights OK
  
  #focal OK
  
  #call OK
  
  #model (only to extract data)
  if (is_not_null(obj[["model"]])) {
    if (is_null(data) && "data" %in% names(obj[["model"]])) {
      data <- obj[["model"]][["data"]]
    }
  }
  
  if (!msm) {
    
    cluster <- A$cluster
    addl <- A[["addl"]]
    s.d.denom <- A$s.d.denom
    method <- A$method
    imp <- A$imp
    subset <- A$subset
    
    # if (length(distance) > 0 && !is.character(distance) && !is.numeric(distance) && !is.data.frame(distance)) {
    #     stop("The argument to distance must be a vector of distance scores or the (quoted) name of a variable in data that contains distance scores.", call. = FALSE)
    # }
    
    if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
      stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
      stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    t.c <- use.tc.fd(formula, data, treat, covs)
    
    treat <- t.c[["treat"]]
    covs  <- t.c[["covs"]]
    
    #Checks
    if (is_null(covs)) {
      stop("covs data.frame must be specified.", call. = FALSE)
    }
    if (!is.data.frame(covs)) {
      stop("covs must be a data.frame.", call. = FALSE)
    }
    if (is_not_null(data) && !is.data.frame(data)) {
      warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execution will halt.")
      data <- NULL
    }
    
    #Process treat
    if (is_null(treat)) stop("treat must be specified.", call. = FALSE)
    
    if (is.numeric(treat) || is.factor(treat) || (is.character(treat) && length(treat) > 1)) {
      treat <- treat
    }
    else if (is.character(treat) && length(treat)==1 && any(names(data) == treat)) {
      treat <- data[[treat]]
    }
    else stop("The argument to treat must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
    
    if (sum(is.na(treat)) > 0)
      stop("Missing values exist in treat.", call. = FALSE)
    
    if (nunique(treat) == 2) {
      treat <- binarize(treat)
    }
    else if (is.character(treat)) {
      treat <- factor(treat)
    }
    
    #Process weights, addl, and distance
    for (i in c("weights", "addl", "distance")) {
      assign(i, data.frame.process(i, A[[i]], treat, covs, data))
    }
    
    #Process subclass
    if (is_not_null(subclass)) {
      if (is.numeric(subclass) || is.factor(subclass) || (is.character(subclass) && length(subclass) > 1)) {
        subclass <- factor(subclass)
      }
      else if (is.character(subclass) && length(subclass)==1 && any(names(data) == subclass)) {
        subclass <- factor(data[[subclass]])
      }
      else stop("The name supplied to subclass is not the name of a variable in data.", call. = FALSE)
      
      subclass.weights <- match.strata2weights(match.strata = subclass,
                                               treat = treat,
                                               covs = covs)
    }
    
    #Process match.strata
    if (is_not_null(match.strata)) {
      if (is.numeric(match.strata) || is.factor(match.strata) || (is.character(match.strata) && length(match.strata) > 1)) {
        match.strata <- factor(match.strata)
      }
      else if (is.character(match.strata) && length(match.strata)==1 && any(names(data) == match.strata)) {
        match.strata <- factor(data[[match.strata]])
      }
      else stop("The name supplied to match.strata is not the name of a variable in data.", call. = FALSE)
      
      match.strata.weights <- match.strata2weights(match.strata = match.strata,
                                                   treat = treat,
                                                   covs = covs)
    }
    
    specified <- setNames(rep(FALSE, 3), c("match.strata", "subclass", "weights"))
    if (is_not_null(weights)) {
      if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
        stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
      }
      specified["weights"] <- TRUE
    }
    if (is_not_null(subclass)){
      if (!is.character(subclass) && !is.factor(subclass) && !is.numeric(subclass)) {
        stop("The argument to subclass must be a vector of subclass membership or the (quoted) name of a variable in data that contains subclass membership.", call. = FALSE)
      }
      specified["subclass"] <- TRUE
    }
    if (is_not_null(match.strata)) {
      if (!is.character(match.strata) && !is.numeric(match.strata) && !is.factor(match.strata)) {
        stop("The argument to match.strata must be a vector of match stratum membership or the (quoted) name of a variable in data that contains match stratum membership.", call. = FALSE)
      }
      specified["match.strata"] <- TRUE
    }
    
    #Getting method
    if (is_null(method)) {
      if (specified["match.strata"]) {
        if (specified["subclass"]) {
          if (isTRUE(all.equal(match.strata.weights, subclass.weights, 
                               check.attributes = FALSE))) {
            subclass.weights <- subclass <- NULL
            specified["subclass"] <- FALSE
          }
        }
        if (specified["weights"]) {
          if (ncol(weights) == 1 && isTRUE(all.equal(match.strata.weights, weights, 
                                                     check.attributes = FALSE))) {
            weights <- NULL
            specified["weights"] <- FALSE
          }
        }
        method <- "matching"
      }
      else if (specified["subclass"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Assuming \"subclassification\" and using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
          weights <- match.strata <- NULL
        }
        method <- "subclassification"
        #weights <- rep(1, nrow(covs))
      }
      else if (specified["weights"]) {
        if (sum(specified) > 1) {
          message(word_list(names(specified)[specified]), " are specified. Assuming \"weighting\" and using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
          match.strata <- subclass <- NULL
        }
        else {
          if (!any(c("optweight", "weightit") %in% class(obj))) {
            message("Assuming \"weighting\". If not, specify with an argument to method.")}
        }
        method <- "weighting"
      }
      else {
        method <- "matching"
      }
    }
    else if (length(method) == 1) {
      specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
      if (specified.method == "weighting") {
        if (specified["weights"]) {
          if (sum(specified) > 1) {
            message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
            match.strata <- subclass <- NULL
          }
          method <- "weighting"
        }
        else if (specified["match.strata"]) {
          message("method = \"weighting\" is specified, but no weights are present. Assuming \"matching\" and using match.strata instead.")
          subclass <- NULL
          method <- "matching"
        }
        else if (specified["subclass"]) {
          message("method = \"weighting\" is specified, but no weights are present. Assuming \"subclassification\" and using subclass instead.")
          method <- "subclassification"
          #weights <- rep(1, nrow(covs))
        }
        else {
          method <- "matching"
        }
      }
      else if (specified.method == "matching") {
        if (specified["match.strata"]) {
          if (sum(specified) > 1) {
            message(word_list(names(specified)[specified]), " are specified. Using match.strata and ignoring ", word_list(names(specified)[specified & names(specified)!="match.strata"]), ".")
            weights <- subclass <- NULL
          }
          method <- "matching"
        }
        else if (specified["weights"]) {
          if (sum(specified) > 1) {
            message(word_list(names(specified)[specified]), " are specified. Using weights and ignoring ", word_list(names(specified)[specified & names(specified)!="weights"]), ".")
            match.strata <- subclass <- NULL
          }
          method <- "matching"
        }
        else if (specified["subclass"]) {
          message("method = \"matching\" is specified, but no weights or match.strata are present. Assuming \"subclassification\" and using subclass instead.")
          method <- "subclassification"
          #weights <- rep(1, nrow(covs))
        }
        else {
          method <- "matching"
        }
      }
      else if (specified.method == "subclassification") {
        if (specified["subclass"]) {
          if (sum(specified) > 1) {
            message(word_list(names(specified)[specified]), " are specified. Using subclass and ignoring ", word_list(names(specified)[specified & names(specified)!="subclass"]), ".")
            weights <- match.strata <- NULL
          }
          method <- "subclassification"
          #weights <- rep(1, nrow(covs))
        }
        else if (specified["match.strata"]) {
          message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"matching\" and using match.strata instead.")
          weights <- NULL
          method <- "matching"
        }
        else if (specified["weights"]) {
          message("method = \"subclassification\" is specified, but no subclass is present. Assuming \"weighting\" and using weights instead.")
          method <- "weighting"
        }
      }
    }
    else {
      specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
      if (any(specified.method == "subclassification") || specified["subclass"]) {
        stop("Subclassification cannot be specified along with other methods.", call. = FALSE)
      }
      else if (specified["match.strata"]) {
        stop("Only weights can be specified with mutiple methods.", call. = FALSE)
      }
      else if (!specified["weights"]) {
        warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
        method <- "matching"
      }
      else {
        #Matching and/or weighting with various weights
        method <- specified.method
        match.strata <- subclass <- NULL
      }
    }
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
      if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
        stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
      }
      if (is.character(s.weights) && length(s.weights)==1) {
        if (any(names(data) == s.weights)) {
          s.weights <- data[[s.weights]]
        }
        else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
      }
      if (anyNA(s.weights)) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    }
    
    #Process cluster
    if (is_not_null(cluster)) {
      if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
        cluster <- cluster
      }
      else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
      if (!is.logical(subset)) {
        stop("The argument to subset must be a logical vector.", call. = FALSE)
      }
      if (anyNA(subset)) {
        warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
        subset[is.na(subset)] <- FALSE
      }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("treat", "subclass", "match.strata", "cluster", "s.weights", "subset")
    data.frames <- c("covs", "weights", "distance", "addl")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames))), c(vectors, data.frames))
    lengths <- setNames(c(lengths(mget(vectors, ifnotfound = list(NULL))), 
                          vapply(data.frames, 
                                 function(x) NROW(get0(x)), 
                                 numeric(1L))), c(vectors, data.frames))
    #Process imp
    if (is_not_null(imp)) {
      if (is_(imp, c("numeric", "factor")) || (is.character(imp) && length(imp) > 1)) {
        imp <- imp
      }
      else if (is.character(imp) && length(imp)==1 && any(names(data) == imp)) {
        imp <- data[[imp]]
      }
      else stop("The name supplied to imp is not the name of a variable in data.", call. = FALSE)
      
      imp.lengths <- vapply(unique(imp), function(i) sum(imp == i), numeric(1L))
      
      if (all_the_same(imp.lengths)) { #all the same
        for (i in vectors) {
          if (lengths[i] > 0 && lengths[i] != length(imp)) { 
            if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
            if (lengths[i] == imp.lengths[1]) {
              temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                     order2 = seq_along(imp))
              
              temp.var <- data.frame(sort(imp), rep(seq_len(lengths[i]), length(imp.lengths)),
                                     get(i)[rep(seq_len(lengths[i]), length(imp.lengths))]
              )
              temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                  by.y = 1:2, sort = FALSE)
              assign(i, temp.merge[[4]][order(temp.merge[[3]])])
            }
            else {
              problematic[i] <- TRUE
            }
          }
        }
        for (i in data.frames) {
          if (lengths[i] > 0 && lengths[i] != length(imp)) {
            if (nunique.gt(imp.lengths, 1)) stop("The number of units in each imputation must be the same unless other inputs provide an observation for each unit in each imputation.", call. = FALSE)
            if (lengths[i] == imp.lengths[1]) {
              temp.imp <- data.frame(imp = imp, order = rep(seq_len(lengths[i]), length(imp.lengths)),
                                     order2 = seq_along(imp))
              temp.var <- data.frame(sort(imp),rep(seq_len(lengths[i]), length(imp.lengths)),
                                     get(i)[rep(seq_len(lengths[i]), length(imp.lengths)), , drop = FALSE]
              )
              temp.merge <- merge(temp.imp, temp.var, by.x = c("imp", "order"), 
                                  by.y = 1:2, sort = FALSE)
              assign(i, setNames(temp.merge[order(temp.merge[[3]]), -c(1:3), drop = FALSE], names(get(i))))
            }
            else {
              problematic[i] <- TRUE
            }
          }
        }
      }
      else {
        problematic <- lengths > 0 & lengths != length(imp)
      }
      if (any(problematic)) {
        stop(paste0(word_list(names(problematic)[problematic]), " must have the same number of observations as imp."), call. = FALSE)
      }
      else ensure.equal.lengths <- FALSE
    }
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
      for (i in c(vectors, data.frames[data.frames!="covs"])) {
        if (lengths[i] > 0 && lengths[i] != lengths["covs"]) {
          problematic[i] <- TRUE
        }
      }
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs."), call. = FALSE)
    }
    
    #Turn match.strata into weights
    if (is_not_null(get0("match.strata.weights"))) {
      weights <- data.frame(weights = match.strata.weights)
    }
    
    if (is_not_null(weights)) {
      if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
      if (any(vapply(weights, function(x) !is.numeric(x), logical(1L)))) {
        stop("All weights must be numeric.", call. = FALSE)
      }
      if (length(method) == 1) {
        method <- rep(method, ncol(weights))
      }
      else if (length(method) != ncol(weights)) {
        stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
      }
      
    }
    
    #Check focal
    if (is_not_null(focal) && is.factor(treat)) {
      if (is.numeric(focal)) {
        if (focal <= nunique(treat)) focal <- levels(treat)[focal]
        else 
          stop(paste0("focal was specified as ", focal, 
                      ", but there are only ", levels(treat), " treatment groups."), call. = FALSE)
      }
      else {
        if (!any(levels(treat) == focal)) 
          stop(paste0("The name specified to focal is not the name of any treatment group."), call. = FALSE)
      }
    }
    
    #Get s.d.denom
    if (is_binary(treat) || !is.numeric(treat)) { #non-continuous
      X$s.d.denom <- get.s.d.denom(s.d.denom, estimand, weights, treat, focal, method)
    }
    
    if (any(c(is.na(covs), is.na(addl)))) {
      warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    if (is_null(subset)) subset <- rep(TRUE, length(treat))
    
    X$method <- method
    X$covs <- covs
    X$weights <- weights
    X$treat <- treat
    X$distance <- distance
    X$subclass <- subclass
    X$cluster <- factor(cluster)
    X$call <- call
    X$addl <- addl
    X$imp <- factor(imp)
    X$s.weights <- s.weights
    X$focal <- focal
    
  }
  else {
    cluster <- A$cluster
    addl.list <- A$addl.list
    s.d.denom <- A$s.d.denom
    method <- A$method
    imp <- A$imp
    subset <- A$subset
    
    if (any(vapply(weights, anyNA, logical(1L)))) stop("NAs are not allowed in the weights.", call. = FALSE)
    if (is_not_null(s.weights) && any(vapply(s.weights, anyNA, logical(1L)))) stop("NAs are not allowed in the sampling weights.", call. = FALSE)
    
    initial.list.lengths <- c(length(formula.list), length(covs.list), length(treat.list))
    if (!all_the_same(initial.list.lengths[initial.list.lengths != 0])) stop("The lists in the object were not the same length.", call. = FALSE)
    ntimes.guess <- max(initial.list.lengths)
    
    if (is_null(treat.list)) treat.list <- vector("list", length(ntimes.guess)) 
    if (is_null(covs.list)) covs.list <- vector("list", length(ntimes.guess)) 
    for (i in seq_len(ntimes.guess)) {
      t.c <- use.tc.fd(formula.list[[i]], data, treat.list[[i]], covs.list[[i]])
      
      treat.list[[i]] <- t.c[["treat"]]
      covs.list[[i]]  <- t.c[["covs"]]
      if (is_not_null(t.c[["treat.name"]])) names(treat.list)[i] <- t.c[["treat.name"]]
    }
    
    ntimes <- length(covs.list)
    
    #Checks
    if (is_null(covs.list)) {
      stop("A covariate list must be specified.", call. = FALSE)
    }
    if (!is.list(covs.list)) {
      stop("covs.list must be a list of covariates for which balanced is to be assessed at each time point.", call. = FALSE)
    }
    if (any(!vapply(covs.list, function(x) is.data.frame(x), logical(1L)))) {
      stop("Each item in covs.list must be a data frame.", call. = FALSE)
    }
    
    if (is_not_null(data) && !is.data.frame(data)) {
      warning("The argument to data is not a data.frame and will be ignored. If the argument to treat is not a vector, the execuction will halt.")
      data <- NULL
    }
    
    specified <- setNames(rep(FALSE, 1), "weights")
    if (is_not_null(weights)) {
      if (!is.character(weights) && !is.numeric(weights) && !is.data.frame(weights) && !is.list(weights)) {
        stop("The argument to weights must be a vector, list, or data frame of weights or the (quoted) names of variables in data that contain weights.", call. = FALSE)
      }
      specified["weights"] <- TRUE
    }
    
    #Getting method
    if (is_null(method)) {
      if (specified["weights"]) {
        
        #message("Assuming \"weighting\". If not, specify with an argument to method.")
        
        method <- "weighting"
      }
      else {
        method <- "matching"
      }
    }
    else if (length(method) == 1) {
      specified.method <- match_arg(method, c("weighting", "matching", "subclassification"))
      if (specified.method == "weighting") {
        if (specified["weights"]) {
          method <- "weighting"
        }
        else {
          method <- "matching"
        }
      }
      else {
        if (specified["weights"]) {
          warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
          method <- "matching"
        }
        else {
          method <- "matching"
        }
      }
    }
    else {
      specified.method <- match_arg(method, c("weighting", "matching", "subclassification"), several.ok = TRUE)
      if (any(specified.method == "subclassification") || specified["subclass"]) {
        warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
        method <- "matching"
      }
      else if (specified["match.strata"]) {
        warning("Only weighting is allowed with multiple treatment time points. Assuming weighting instead.", call. = FALSE)
        method <- "matching"
      }
      else if (!specified["weights"]) {
        warning("Multiple methods were specified, but no weights were provided. Providing unadjusted data only.", call. = FALSE)
        method <- "matching"
      }
      else {
        #Matching and/or weighting with various weights
        method <- specified.method
      }
    }
    
    if (is_not_null(cluster) && !is.character(cluster) && !is.numeric(cluster) && !is.factor(cluster)) {
      stop("The argument to cluster must be a vector of cluster membership or the (quoted) name of a variable in data that contains cluster membership.", call. = FALSE)
    }
    if (is_not_null(imp) && !is.character(imp) && !is.numeric(imp) && !is.factor(imp)) {
      stop("The argument to imp must be a vector of imputation IDs or the (quoted) name of a variable in data that contains imputation IDs.", call. = FALSE)
    }
    
    #Order covs.list
    all.covs <- unique(unlist(lapply(covs.list, names)))
    covs.list <- lapply(covs.list, function(x) x[all.covs[all.covs %in% names(x)]])
    
    #Process treat
    if (is_null(treat.list)) stop("A treatment list must be specified.", call. = FALSE)
    if (!is.vector(treat.list)) {
      treat.list <- as.list(treat.list)
    }
    if (length(treat.list) != length(covs.list)) {
      stop("treat.list must be a list of treatment statuses at each time point.", call. = FALSE)
    }
    
    for (ti in seq_along(treat.list)) {
      if (is.numeric(treat.list[[ti]]) || is.factor(treat.list[[ti]]) || (is.character(treat.list[[ti]]) && length(treat.list[[ti]]) > 1)) {
        #treat.list[[ti]] <- treat.list[[ti]]
      }
      else if (is.character(treat.list[[ti]]) && length(treat.list[[ti]])==1 && any(names(data) == treat.list[[ti]])) {
        names(treat.list)[ti] <- treat.list[[ti]]
        treat.list[[ti]] <- data[[treat.list[[ti]]]]
      }
      else stop("Each item in treat.list must be a vector of treatment statuses or the (quoted) name of a variable in data that contains treatment status.", call. = FALSE)
      
      if (sum(is.na(treat.list[[ti]])) > 0)
        stop("Missing values exist in treat.list", call. = FALSE)
      
      if (is_binary(treat.list[[ti]])) {
        treat.list[[ti]] <- binarize(treat.list[[ti]])
      }
      else if (is.character(treat.list[[ti]])) {
        treat.list[[ti]] <- factor(treat.list[[ti]])
      }
    }
    
    #Process weights
    for (i in c("weights")) {
      assign(i, data.frame.process(i, A[[i]], do.call("cbind", treat.list), do.call("cbind", covs.list), data))
    }
    
    #Process addl and distance
    for (i in c("addl.list", "distance.list")) {
      assign(i, list.process(i, A[[i]], ntimes, 
                             "covs.list",
                             treat.list,
                             covs.list,
                             data
      ))
    }
    
    #Process sampling weights
    if (is_not_null(s.weights)) {
      if (!(is.character(s.weights) && length(s.weights) == 1) && !is.numeric(s.weights)) {
        stop("The argument to s.weights must be a vector or data frame of sampling weights or the (quoted) names of variables in data that contain sampling weights.", call. = FALSE)
      }
      if (is.character(s.weights) && length(s.weights)==1) {
        if (any(names(data) == s.weights)) {
          s.weights <- data[[s.weights]]
        }
        else stop("The name supplied to s.weights is not the name of a variable in data.", call. = FALSE)
      }
    }
    
    #Process cluster
    if (is_not_null(cluster)) {
      if (is.numeric(cluster) || is.factor(cluster) || (is.character(cluster) && length(cluster)>1)) {
        cluster <- cluster
      }
      else if (is.character(cluster) && length(cluster)==1 && any(names(data) == cluster)) {
        cluster <- data[[cluster]]
      }
      else stop("The name supplied to cluster is not the name of a variable in data.", call. = FALSE)
    }
    
    #Process subset
    if (is_not_null(subset)) {
      if (!is.logical(subset)) {
        stop("The argument to subset must be a logical vector.", call. = FALSE)
      }
      if (anyNA(subset)) {
        warning("NAs were present in subset. Treating them like FALSE.", call. = FALSE)
        subset[is.na(subset)] <- FALSE
      }
    }
    
    ensure.equal.lengths <- TRUE
    vectors <- c("s.weights", "cluster", "subset")
    data.frames <- c("weights")
    lists <- c("treat.list", "distance.list", "addl.list", "covs.list")
    problematic <- setNames(rep(FALSE, length(c(vectors, data.frames, lists))), c(vectors, data.frames, lists))
    lengths <- setNames(c(lengths(mget(vectors)), 
                          vapply(data.frames, 
                                 function(x) {if (is.null(get0(x))) 0 else nrow(get(x))
                                 }, numeric(1L)),
                          vapply(lists, function(x) {
                            if (is.null(get0(x))) 0 
                            else if (is.vector(get(x))) {
                              if (is.data.frame(get(x)[[1]]) || is.matrix(get(x)[[1]])) max(vapply(get(x), nrow, numeric(1L)))
                              else max(lengths(get(x)))
                            }
                            else max(vapply(get(x), function(y) if (is_not_null(y)) {if (is_not_null(nrow(y))) nrow(y) else length(y)} else 0, numeric(1L)))
                          }, numeric(1L))), c(vectors, data.frames, lists))
    
    #Ensure all input lengths are the same.
    if (ensure.equal.lengths) {
      for (i in c(vectors, data.frames, lists)) {
        if (lengths[i] > 0 && lengths[i] != lengths["covs.list"]) {
          problematic[i] <- TRUE
        }
      }
    }
    if (any(problematic)) {
      stop(paste0(word_list(names(problematic[problematic])), " must have the same number of observations as covs.list."), call. = FALSE)
    }
    
    if (is_not_null(weights)) {
      if (any(vapply(weights, function(x) any(!is.finite(x)), logical(1L)))) {
        stop("All weights must be numeric.", call. = FALSE)
      }
      if (length(method) == 1) {
        method <- rep(method, ncol(weights))
      }
      else if (length(method) != ncol(weights)) {
        stop("Valid inputs to method must have length 1 or equal to the number of valid sets of weights.", call. = FALSE)
      }
      
    }
    
    #Get s.d.denom
    X$s.d.denom <- rep("pooled", max(1, ncol(weights)))
    
    if (any(vapply(c(covs.list, addl.list), anyNA, logical(1L)))) {
      warning("Missing values exist in the covariates. Displayed values omit these observations.", call. = FALSE)
    }
    
    if (is_null(s.weights)) s.weights <- rep(1, length(treat.list[[1]]))
    if (is_null(subset)) subset <- rep(TRUE, length(treat.list[[1]]))
    
    X$method <- method
    X$covs.list <- lapply(covs.list, function(x) x)
    X$weights <- weights
    X$treat.list <- lapply(treat.list, function(x) x)
    X$distance.list <- if (is_not_null(distance.list)) lapply(distance.list, function(x) x) else NULL
    X$addl.list <- if (is_not_null(addl.list)) lapply(addl.list, function(x) x) else NULL
    X$cluster <- factor(cluster)
    X$call <- call
    X$imp <- factor(imp)
    X$s.weights <- s.weights
  }
  
  X <- subset_X(X, subset)
  X <- setNames(X[X.names], X.names)
  
  class(X) <- get.X.class(X)
  
  return(X)
}

#Utility functions
f.build <- function(y, rhs) {
  if (missing(rhs)) stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
  else if ((is.data.frame(rhs) || is.matrix(rhs)) && is_not_null(colnames(rhs))) {
    vars <- paste0("`", gsub("`", "", colnames(rhs)), "`")
  }
  else if (is.character(rhs) && is_not_null(rhs)) {
    vars <- paste0("`", gsub("`", "", rhs), "`")
  }
  else stop("Right hand side argument to f.build() must be a vector of variable names or a data set with named variables.", call. = FALSE)
  
  if (missing(y) || identical(y, "")) y <- NULL
  else if (!is.character(y) || length(y) > 1) stop ("Response argument to f.build() must be the quoted name of the response variable.", call. = FALSE)
  
  f <- reformulate(vars, y)
  return(f)
}
splitfactor <- function(data, var.name, replace = TRUE, sep = "_", drop.level = NULL, drop.first = TRUE, drop.singleton = FALSE, drop.na = TRUE, check = TRUE) {
  #Splits factor into multiple (0, 1) indicators, replacing original factor in dataset. 
  #Retains all categories unless only 2 levels, in which case only the second level is retained.
  #If variable only has one level, will delete.
  #var.name= the name of the variable to split when data is specified
  #data=data set to be changed
  
  if (is.data.frame(data)) {
    data <- as.data.frame(data)
    if (check) {
      factor.names <- names(data)[vapply(data, function(x) is.factor(x) || is.character(x), logical(1L))]
      if (missing(var.name)) {
        var.name <- factor.names
      }
      else if (is.character(var.name)) {
        if (any(var.name %in% factor.names)) {
          if (any(!var.name %in% factor.names)) {
            not.in.factor.names <- var.name[!var.name %in% factor.names]
            warning(paste(word_list(not.in.factor.names, "and", is.are = TRUE), 
                          "not the name(s) of factor variable(s) in data and will not be split."), 
                    call. = FALSE)
          }
          var.name <- var.name[var.name %in% factor.names]
        }
        else {
          stop("No names in var.name are names of factor variables in data.", call. = FALSE)
        }
      }
      else {
        stop("var.name must be a character vector of the name(s) of factor variable(s) in data.", call. = FALSE)
      }
      if (is_null(factor.names)) {
        stop("There are no factor variables to split in data.", call. = FALSE)
      }
    }
    else {
      if (missing(var.name) || !is.character(var.name)) {
        stop("var.name must be a character vector of the names of variables in data.", call. = FALSE)
      }
      else {
        if (any(var.name %in% names(data))) {
          if (any(var.name %nin% names(data))) {
            not.in.data.names <- var.name[!var.name %in% names(data)]
            warning(paste(word_list(not.in.data.names, "and", is.are = TRUE), 
                          "not the name(s) of variable(s) in data and will not be split."), 
                    call. = FALSE)
          }
          var.name <- var.name[var.name %in% names(data)]
        }
        else {
          stop("No names in var.name are names of variables in data.", call. = FALSE)
        }
      }
    }
    
  }
  else if (is.atomic(data)) {
    dep <- deparse(substitute(data))
    data <- data.frame(data)
    if (missing(var.name)) {
      names(data) <- dep
    }
    else if (is.vector(var.name) && (is.atomic(var.name) || is.factor(var.name))) {
      if (is_null(var.name)) {
        names(data) <- dep
      }
      else if (length(var.name) == 1) {
        names(data) <- var.name
      }
      else {
        warning("Only using the first item of var.name.", call. = FALSE)
        names(data) <- var.name[1]
      }
    }
    else {
      stop("var.name must be an atomic or factor vector of length 1 with the stem of the new variable.", call. = FALSE)
    }
    var.name <- names(data)
  }
  else {
    stop("data must a be a data.frame or an atomic vector.", call. = FALSE)
  }
  
  if (is_not_null(drop.level) && length(var.name) > 1) {
    warning("drop.level cannot be used with multiple entries to var.name. Ignoring drop.level.", call. = FALSE)
    drop.level <- NULL
  }
  drop.na <- setNames(rep(drop.na, length(var.name)), var.name)
  for (v in var.name) {
    drop <- character(0)
    x <- factor(data[names(data) == v][[1]], exclude = NULL)
    na.level <- is.na(levels(x))
    levels(x) <- paste0(sep, levels(x))
    data[names(data) == v][[1]] <- x
    
    skip <- FALSE
    if (nlevels(x) > 1) {
      k <- model.matrix(as.formula(paste0("~`", v, "`- 1")), data = data)
      colnames(k) <- gsub("`", colnames(k), replacement = "")
      
      if (any(na.level)) {
        if (drop.na[v]) {
          k[k[,na.level] == 1,] <- NA_real_
        }
      }
      else drop.na[v] <- FALSE
      
    }
    else {
      if (drop.singleton) {
        data <- data[names(data)!=v]
        skip <- TRUE
      }
      else {
        k <- matrix(1, ncol = 1, nrow = length(x))
        colnames(k) <- paste0(v, levels(x)[1])
      }
    }
    
    if (!skip) {
      if (is_not_null(drop.level)) {
        if (is.character(drop.level) && length(drop.level) == 1 && drop.level %in% levels(x)) {
          drop <- drop.level
        }
        else {
          stop(paste("drop must be the name of a level of", v, "which is to be dropped."), call. = FALSE)
        }
      }
      else {
        if ((ncol(k) == 2 && (drop.first == "if2" || drop.first == TRUE)) ||
            (ncol(k) > 2 && drop.first == TRUE)) {
          drop <- levels(x)[1]
        }
      }
      
      dropl <- rep(FALSE, ncol(k))
      if (is_not_null(drop)) {
        dropl[!na.level & levels(x) %in% drop] <- TRUE
      }
      if (drop.na[v]) dropl[na.level] <- TRUE
      
      k <- k[,!dropl, drop = FALSE]
      
      if (ncol(data) == 1) {
        data <- data.frame(k, row.names = rownames(data))
      }
      else if (replace) {
        if (match(v, names(data)) == 1){
          data <- cbind(k, data[names(data)!=v], row.names = rownames(data))
        }
        else if (match(v, names(data)) == ncol(data)) {
          data <- cbind(data[names(data)!=v], k, row.names = rownames(data))
        }
        else {
          where <- match(v, names(data))
          data <- cbind(data[1:(where-1)], k, data[(where+1):ncol(data)], row.names = rownames(data))
        }
      }
      else {
        data <- data.frame(data, k, row.names = rownames(data))
      }
      
    }
    
  }
  
  return(data)
}
unsplitfactor <- function(data, var.name, replace = TRUE, sep = "_", dropped.level = NULL, dropped.na = TRUE) {
  
  if (!is.data.frame(data)) stop("data must be a data.frame containing the variables to unsplit.", call = FALSE)
  if (!is.character(var.name)) stop("var.name must be a string containing the name of the variables to unsplit.", call. = FALSE)
  if (is_not_null(dropped.level) && length(var.name) > 1) {
    warning("dropped.level cannot be used with multiple var.names and will be ignored.", call. = FALSE, immediate. = TRUE)
    dropped.level <- NULL
  }
  
  if (!is.character(var.name)) stop("var.name must be a character vector containing the name of the variable to unsplit.", call. = FALSE)
  if (length(sep) > 1 || !is.character(sep)) stop("sep must be a character vector of length 1 containing the seperating character in the names of the split variables.", call. = FALSE)
  if (length(dropped.level) > 1 && !is.atomic(dropped.level)) {
    warning("dopped.level must be an atomic vector of length 1 containing the value of the dropped category of the split variable. It will be ignored.", call. = FALSE, immediate. = TRUE)
    dropped.level <- NULL
  }
  not.the.stem <- character(0)
  
  for (v in var.name) {
    dropped.level0 <- dropped.level
    var.to.combine <- data[startsWith(names(data), paste0(v, sep))]
    if (is_null(var.to.combine)) {
      not.the.stem <- c(not.the.stem, paste0(v, sep))
      next
    }
    
    if (any(rowSums(apply(var.to.combine, 2, is.na)) %nin% c(0, ncol(var.to.combine)))) {
      stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the <NA> pattern.", call. = FALSE)
    }
    NA.column <- character(0)
    
    if (!isTRUE(dropped.na)) {
      NA.column <- paste0(v, sep, ifelse(dropped.na == FALSE, "NA", dropped.na))
      if (NA.column %in% names(var.to.combine)) {
        var.to.combine[var.to.combine[[NA.column]] == 1,] <- NA_real_
        var.to.combine <- var.to.combine[names(var.to.combine) != NA.column]
      }
      else {
        stop(paste("There is no variable called", word_list(NA.column, quotes = TRUE), "to generate the NA values."), call. = FALSE)
      }
    }
    var.sum <- rowSums(var.to.combine)
    if (isTRUE(all.equal(unique(var.sum), 1))) {
      #Already unsplit
    }
    else if (isTRUE(all.equal(sort(unique(var.sum)), c(0, 1)))) {
      #Missing category
      
      if (is_null(dropped.level)) {
        k.levels0 <- sapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep), fixed = TRUE)[[1]][2])
        
        if (suppressWarnings(all(!is.na(as.numeric(k.levels0))))) {
          dropped.level0 <- as.character(min(as.numeric(k.levels0)) - 1)
          dropped.name <- paste0(v, sep, dropped.level0)
        }
        else {
          message("The dropped category will be set to NA.")
          dropped.name <- dropped.level0 <- NA_character_
        }
        
      }
      else dropped.name <- paste0(v, sep, dropped.level)
      var.to.combine <- setNames(data.frame(1-var.sum, var.to.combine),
                                 c(dropped.name, names(var.to.combine)))
      
    }
    else {
      stop("The variables in data selected based on var.name and sep do not seem to form a split variable based on the row sums.", call. = FALSE)
    }
    
    k.levels <- vapply(names(var.to.combine), function(x) strsplit(x, paste0(v, sep), fixed = TRUE)[[1]][2], character(1L))
    
    k <- rep(NA_character_, nrow(data))
    for (i in seq_along(k.levels)) {
      k <- ifelse(var.to.combine[[i]] == 1, k.levels[i], k)
    }
    
    k <- factor(k, levels = k.levels)
    
    
    if (replace) {
      where <- which(names(data) %in% c(names(var.to.combine), NA.column))
      
      data[[min(where)]] <- k
      remove.cols <- where[where!=min(where)]
      if (is_not_null(remove.cols)) data <- data[-remove.cols]
      names(data)[min(where)] <- v
    }
    else {
      data <- cbind(data, setNames(data.frame(k), v))
    }
  }
  
  if (is_not_null(not.the.stem)) warning(paste0(word_list(not.the.stem, is.are = TRUE, quotes = TRUE), " not the stem of any variables in data and will be ignored. Ensure var.name and sep are correct."), call. = FALSE)
  
  return(data)
}

var.names <- function(b, type, file = NULL, minimal = FALSE) {
  if (is_not_null(b[["print.options"]][["co.names"]])) {
    if (minimal) vars <- unique(unlist(lapply(b[["print.options"]][["co.names"]], function(x) x[["component"]][x[["type"]] == "base"])))
    else vars <- vapply(b[["print.options"]][["co.names"]], function(x) paste(x[["component"]], collapse = ""), character(1))
  }
  else {
    vars <- NULL
    var.containers <- c(quote(b[["Balance"]]),
                        quote(b[["Cluster.Balance"]][[1]][["Balance"]]),
                        quote(b[["Subclass.Balance"]][[1]]),
                        quote(b[["Imputation.Balance"]][[1]][["Balance"]]),
                        quote(b[["Imputation.Balance"]][[1]][["Cluster.Balance"]][[1]][["Balance"]]),
                        quote(b[["Pair.Balance"]][[1]]),
                        quote(b[["Time.Balance"]][[1]][["Balance"]]))
    for (i in var.containers) {
      obj <- eval(i)
      if (is_not_null(obj)) {
        vars <- rownames(obj)
        break
      }
      else obj <- NULL
    }
    if (is_null(vars)) stop("No variable names were found in the object. It is probably not a bal.tab object.", call. = FALSE)
    if (minimal) warning("minimal is being set to FALSE because the part of the object required for it to be TRUE is missing.", call. = FALSE)
  }
  
  if (is_not_null(file)) {
    if (!endsWith(file, ".csv")) stop("The filename in file must end in \".csv\".", call. = FALSE)
  }
  
  if (missing(type)) {
    if (is_not_null(file)) type <- "df"
    else type <- "vec"
  }
  else {
    possible.types <- c("df", "vec")
    type <- possible.types[pmatch(type, possible.types)]
  }
  
  if (is.na(type)) stop("type must be \"df\" or \"vec\"")
  else if (type == "df") {
    out <- data.frame(old = vars, new = vars, stringsAsFactors = FALSE, row.names = NULL)
  }
  else {
    out <- setNames(vars, vars)
  }
  
  if (is_not_null(file)) {
    if (type == "df") {
      write.csv(out, file = file, row.names = FALSE)
      invisible(out)
    }
    else {
      warning("Only type = \"df\" is compatible with a file name.", call. = FALSE)
      out
    }
  }
  else out
}

set.cobalt.options <- function(..., default = FALSE) {
  opts <- list(...)
  if (is_not_null(opts) && (is_null(names(opts)) ||  "" %in% names(opts))) {
    stop("All arguments must be named.", call. = FALSE)
  }
  # if ("continuous" %in% names(opts)) names(opts)[names(opts) == "continuous"] <- "cont"
  # if ("binary" %in% names(opts)) names(opts)[names(opts) == "binary"] <- "bin"
  
  multiple.allowed <- c("cluster.fun", "imp.fun")
  any.string.allowed <- c("int_sep", "factor_sep")
  
  if (any(duplicates <- table(names(opts)) > 1)) {
    stop(paste0(word_list(names(duplicates)[duplicates], is.are = TRUE), " present more than once in the input to set.cobalt.options."), call. = FALSE)
  }
  
  if (any(names(opts) %nin% names(acceptable.options()))) {
    warning(paste("The following are not acceptable options and will be ignored:", word_list(unique(names(opts)[names(opts) %nin% names(acceptable.options())]))), call. = FALSE, immediate. = TRUE)
    opts <- opts[names(opts) %in% names(acceptable.options())]
  }
  
  if (default) {
    return.to.default <- setdiff(names(acceptable.options()), names(opts))
  }
  else return.to.default <- NULL
  
  multiple.opts <- NULL
  bad.opts <- NULL
  for (i in names(opts)) {
    if (is_null(opts[[i]])) {
      return.to.default <- c(return.to.default, i)
      opts[[i]] <- NULL
    }
    else {
      if (length(opts[[i]]) > 1 && i %nin% multiple.allowed) multiple.opts <- c(multiple.opts, i)
      if (mode(opts[[i]]) != mode(acceptable.options()[[i]]) || 
          (!(is.character(opts[[i]]) && is.character(acceptable.options()[[i]]) && (i %in% any.string.allowed || !is.na(pmatch(opts[[i]], acceptable.options()[[i]])))) &&
           !all(opts[[i]] %in% acceptable.options()[[i]]))) bad.opts <- c(bad.opts, i)
    }
  }
  
  if (is_not_null(opts)) {
    both.opts <- intersect(multiple.opts, bad.opts)
    multiple.opts <- multiple.opts[multiple.opts %nin% both.opts]
    bad.opts <- bad.opts[bad.opts %nin% both.opts]
    problematic.opts <- setNames(vector("list", 3), c("multiple", "bad", "both"))
    problematic.opts[["multiple"]] <- setNames(lapply(multiple.opts, function(i) {
      paste(i, "must be of length 1.")
    }), multiple.opts)
    problematic.opts[["bad"]] <- setNames(lapply(bad.opts, function(i) {
      if (i %in% any.string.allowed) paste0(i, " must be a character string.")
      else paste0(i, " must be ", word_list(acceptable.options()[[i]], quotes = is.character(acceptable.options()[[i]]), and.or = "or"), ".")
    }), bad.opts)
    problematic.opts[["both"]] <- setNames(lapply(both.opts, function(i) {
      if (i %in% any.string.allowed) paste0(i, " must be a character string of length 1.")
      else paste0(i, " must be one of ", word_list(acceptable.options()[[i]], quotes = is.character(acceptable.options()[[i]]), and.or = "or"), ".")
    }), both.opts)
    
    problems <- do.call("c", unname(problematic.opts))
    problems <- problems[names(opts)[names(opts) %in% names(problems)]]
    if (is_not_null(problems)) {
      stop(do.call("paste", c(list(""), problems, list("\nNo options will be set.", sep = "\n"))), call. = FALSE)
    }
    
    names(opts) <- paste0("cobalt_", names(opts))
    options(opts)
  }
  
  if (is_not_null(return.to.default)) {
    options(setNames(replicate(length(return.to.default), NULL), paste0("cobalt_", return.to.default)))
  }
  # if ("continuous" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "continuous"] <- "cont"
  # if ("binary" %in% names(opts)) names(acceptable.options)[names(acceptable.options) == "binary"] <- "bin"
}
get.cobalt.options <- function(...) {
  opts <- list(...)
  
  opts <- clear_null(opts)
  if (is_null(opts)) opts <- names(acceptable.options())
  else {
    if (!all(vapply(opts, is.character, logical(1L)))) {
      stop("All arguments must be strings containing the name of an option to return.", call. = FALSE)
    }
    opts <- do.call("c", opts)
    if (any(not.in.accept <- opts %nin% names(acceptable.options()))) {
      plural <- sum(not.in.accept) > 1
      stop(paste0(word_list(opts[not.in.accept], is.are = TRUE, quotes = TRUE),
                  " not", ifelse(plural, "", " an"), " acceptable option", 
                  ifelse(plural, "s", ""), "."), call. = FALSE)
    }
  }
  
  out <- setNames(lapply(paste0("cobalt_", opts), getOption), opts)
  return(out)
  
}

love.plot <- function(x, stats = "mean.diffs", abs = TRUE, agg.fun = NULL, 
                      var.order = NULL, drop.missing = TRUE, drop.distance = FALSE, 
                      threshold = NULL, line = FALSE, stars = "none", grid = TRUE, 
                      limits = NULL, colors = NULL, shapes = NULL, alpha = 1, size = 3, 
                      wrap = 30, var.names = NULL, title, sample.names, labels = FALSE,
                      position = "right", themes = NULL, ...) {
  
  #Replace .all and .none with NULL and NA respectively
  .call <- match.call(expand.dots = TRUE)
  .alls <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.all)), logical(1L))
  .nones <- vapply(seq_along(.call), function(x) identical(.call[[x]], quote(.none)), logical(1L))
  if (any(c(.alls, .nones))) {
    .call[.alls] <- expression(NULL)
    .call[.nones] <- expression(NA)
    return(eval(.call))
  }
  
  #Re-call bal.tab with disp.v.ratio or disp.ks if stats = "v" or "k".
  if (!exists(deparse(substitute(x)))) { #if x is not an object (i.e., is a function call)
    mc <- match.call()
    replace.args <- function(m) {
      #m_ is bal.tab call or list (for do.call)
      if (any(sapply(stats, function(x) pmatch(x, "variance.ratios", 0L) != 0L))) {
        if (!isTRUE(eval(m[["disp.v.ratio"]]))) {
          m[["un"]] <- TRUE
          m[["disp.v.ratio"]] <- TRUE
        }
      }
      if (any(sapply(stats, function(x) pmatch(x, "ks.statistics", 0L) != 0L))) {
        if (!isTRUE(eval(m[["disp.ks"]]))) {
          m[["un"]] <- TRUE
          m[["disp.ks"]] <- TRUE
        }
      }
      
      if (any(names(m) == "cluster")) {
        m[["cluster.summary"]] <- TRUE
        if (any(names(m) == "cluster.fun")) m[["cluster.fun"]] <- NULL
      }
      if (any(names(m) == "imp")) {
        m[["imp.summary"]] <- TRUE
        if (any(names(m) == "imp.fun")) m[["imp.fun"]] <- NULL
      }
      
      m[["abs"]] <- abs
      
      return(m)
    }
    
    if (deparse(mc[["x"]][[1]]) %in% c("bal.tab", methods("bal.tab"))) { #if x i bal.tab call
      mc[["x"]] <- replace.args(mc[["x"]])
      x <- eval(mc[["x"]])
      
    }
    else if (deparse(mc[["x"]][[1]]) == "do.call") { #if x is do.call
      d <- match.call(eval(mc[["x"]][[1]]), mc[["x"]])
      if (deparse(d[["what"]]) %in% c("bal.tab", methods("bal.tab"))) {
        d[["args"]] <- replace.args(d[["args"]])
        x <- eval(d)
      }
    }
  }
  
  if ("bal.tab" %nin% class(x)) {
    #Use bal.tab on inputs first, then love.plot on that
    m <- match.call()
    m.b <- m; m.b[[1]] <- quote(bal.tab); names(m.b)[2] <- ''
    
    m.l <- m; 
    m.l[["x"]] <- m.b
    
    return(eval(m.l))
  }
  
  args <- list(...)
  
  if (any(class(x) == "bal.tab.cont")) {
    stats <- "correlations"
  }
  else stats <- match_arg(stats, c("mean.diffs", "variance.ratios", "ks.statistics"), several.ok = TRUE)
  
  which.stat <- c(mean.diffs = "Diff", variance.ratios = "V.Ratio", ks.statistics = "KS", correlations = "Corr")[stats]
  which.stat2 <- c(Diff = "Mean Difference", V.Ratio = "Variance Ratio", KS = "Kolmogorov-Smirnov Statistic", Corr = "Correlation")[which.stat]
  
  #shape (deprecated)
  #un.color (deprecated)
  #adj.color (deprecated)
  #cluster.fun (deprecated)
  #star_char
  
  p.ops <- c("which.cluster", "which.imp", "which.treat", "which.time", "disp.subclass")
  for (i in p.ops) {
    if (i %in% names(args)) attr(x, "print.options")[[i]] <- args[[i]]
  }
  
  #Using old argument names
  if (is_not_null(args$cluster.fun) && is_null(agg.fun)) agg.fun <- args$cluster.fun
  if (is_not_null(args$no.missing)) drop.missing <- args$no.missing
  
  Agg.Fun <- NULL
  subtitle <- NULL
  facet <- NULL
  
  cluster.names.good <- NULL
  imp.numbers.good <- NULL
  treat.names.good <- NULL; disp.treat.pairs <- NULL
  time.names.good <- NULL
  subclass.names <- NULL
  
  #NA = aggregate, NULL = individual
  null.cluster <- is_null(attr(x, "print.options")$which.cluster)
  na.cluster <- !null.cluster && is.na(attr(x, "print.options")$which.cluster)
  null.imp <- is_null(attr(x, "print.options")$which.imp)
  na.imp <- !null.imp && is.na(attr(x, "print.options")$which.imp)
  null.treat <- is_null(attr(x, "print.options")$which.treat)
  na.treat <- !null.treat && is.na(attr(x, "print.options")$which.treat)
  null.time <- is_null(attr(x, "print.options")$which.time)
  na.time <- !null.time && is.na(attr(x, "print.options")$which.time)
  
  config <- "agg.none"
  
  #Get B and config
  if (any(class(x) == "bal.tab.msm")) {
    #Get time.names.good
    time.names <- names(x[["Time.Balance"]])
    if (na.time) {
      config <- "agg.time"
      time.names.good <- NULL
    }
    else if (null.time) {
      config <- "agg.none"
      time.names.good <- setNames(rep(TRUE, length(time.names)), time.names)
    }
    else if (is.character(attr(x, "print.options")$which.time)) {
      time.names.good <- sapply(time.names, function(x) any(sapply(attr(x, "print.options")$which.time, function(y) identical(x, y))))
      if (any(time.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
      }
    }
    else if (is.numeric(attr(x, "print.options")$which.time)) {
      time.names.good <- setNames(seq_along(x$Time.Balance) %in% attr(x, "print.options")$which.time, time.names)
      if (any(time.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.time are valid treatment names or time point indices.", call. = FALSE)
      }
    }
    else stop("The argument to which.time must be .none, .all, or the names of treatments or indices of time points.", call. = FALSE)
    
    
    #Get B from x
    if (config == "agg.none") {
      B <- do.call("rbind", lapply(names(x[["Time.Balance"]])[time.names.good], 
                                   function(t) cbind(x[["Time.Balance"]][[t]][["Balance"]],
                                                     time = t,
                                                     variable.names = rownames(x[["Time.Balance"]][[t]][["Balance"]]))))
      facet <- "time"
    }
    else if (config == "agg.time") {
      if (is_null(x[["Balance.Across.Times"]])) {
        stop("Cannot aggregate across time periods without a balance summary across time periods.\nThis may be because multinomial treatments were used, multiple treatment types were used,\nor quick was set to TRUE and msm.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      #Agg.Fun <- switch(match_arg(agg.fun), mean = "Mean", max = "Max", range = "Range")
      Agg.Fun <- "Max"
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Time Points"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Time Points")
      }
      B <- cbind(x[["Balance.Across.Times"]],
                 variable.names = rownames(x[["Balance.Across.Times"]]))
    }
  }
  else if (any(class(x) == "bal.tab.imp.cluster")) {
    #Get imp.numbers.good
    imp.numbers <- seq_along(x[["Imputation.Balance"]])
    if (na.imp) {
      imp.numbers.good <- NULL
    }
    else if (null.imp) {
      imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
    }
    else if (is.numeric(attr(x, "print.options")$which.imp)) {
      imp.numbers.good <- setNames(imp.numbers %in% attr(x, "print.options")$which.imp, imp.numbers)
    }
    else stop("The argument to which.imp must be .none, .all, or the indices of imputations.", call. = FALSE)
    
    #Get cluster.names.good
    cluster.names <- names(x[["Cluster.Balance.Across.Imputations"]])
    if (na.cluster) {
      cluster.names.good <- NULL
    }
    else if (null.cluster) {
      cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
    }
    else if (is.character(attr(x, "print.options")$which.cluster)) {
      cluster.names.good <- sapply(cluster.names, function(n) any(sapply(attr(x, "print.options")$which.cluster, function(y) identical(n, y))))
    }
    else if (is.numeric(attr(x, "print.options")$which.cluster)) {
      cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% attr(x, "print.options")$which.cluster, cluster.names)
    }
    else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
    
    #Set configuration type of B using which.imp and which.cluster
    if (na.imp) { #aggregate over all imps
      if (na.cluster) {
        config <- "agg.all"
      }
      else { #1, #6
        if (any(cluster.names.good)) {
          config <- "agg.imp"
        }
        else {
          stop("Make sure the arguments to which.cluster are valid names or indices of clusters.", call. = FALSE)
        }
        
      }
    }
    else if (sum(imp.numbers.good) == 1) {
      if (na.cluster) {
        config <- "agg.cluster"
      }
      else { 
        if (any(cluster.names.good)) {
          config <- "agg.none"
        }
        else {
          stop("Make sure the arguments to which.cluster are valid names or indices of clusters.", call. = FALSE)
        }
      }
    }
    else if (sum(imp.numbers.good) > 1) {
      if (na.cluster) {
        config <- "agg.cluster"
      }
      else if (sum(cluster.names.good) == 1) {
        config <- "agg.none"
      }
      else {
        stop("At least one of which.cluster or which.imp must be .none or of length 1.", call. = FALSE)
      }
    }
    else {
      stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
    }
    
    #Get B from x based on configuration
    if (config == "agg.none") {
      B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good],
                                   function(i) do.call("rbind", lapply(names(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]])[cluster.names.good],
                                                                       function(y) cbind(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]][[y]][["Balance"]],
                                                                                         cluster = y,
                                                                                         imp = paste("Imputation:", i),
                                                                                         variable.names = rownames(x[["Imputation.Balance"]][[i]][["Cluster.Balance"]][[y]][["Balance"]]))))))
      if (sum(imp.numbers.good) == 1) {
        facet <- "cluster"
        subtitle <- paste("Imputation:", imp.numbers[imp.numbers.good])
      }
      else {
        facet <- "imp"
        subtitle <- paste("Cluster:", cluster.names[cluster.names.good])
      }
    }
    else if (config == "agg.imp") {
      if (is_null(x[["Cluster.Balance.Across.Imputations"]])) {
        stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                        mean = "Mean", max = "Max", range = "Range")
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Imputations"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Imputations", sep = " ")
      }
      B <- do.call("rbind", lapply(names(x[["Cluster.Balance.Across.Imputations"]])[cluster.names.good], 
                                   function(c) cbind(x[["Cluster.Balance.Across.Imputations"]][[c]][["Cluster.Balance"]], 
                                                     cluster = c, 
                                                     variable.names = rownames(x[["Cluster.Balance.Across.Imputations"]][[c]][["Cluster.Balance"]]))))
      facet <- "cluster"
    }
    else if (config == "agg.cluster") {
      if (is_null(x[["Imputation.Balance"]][[1]][["Cluster.Summary"]])) {
        stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                        mean = "Mean", max = "Max", range = "Range")
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Clusters"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters", sep = " ")
      }
      B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good], 
                                   function(i) cbind(x[["Imputation.Balance"]][[i]][["Cluster.Summary"]], 
                                                     imp = paste("Imputation:", i), 
                                                     variable.names = rownames(x[["Imputation.Balance"]][[i]][["Cluster.Summary"]]))))
      facet <- "imp"
    }
    else if (config == "agg.all") {
      if (is_null(x[["Balance.Across.Imputations"]])) {
        stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and cluster.summary or imp.summary were set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      #Cluster.Fun <- switch(match_arg(cluster.fun), mean = "Mean", max = "Max", range = "Range")
      Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                        mean = "Mean", max = "Max", range = "Range")
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Clusters and Imputations"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters and Imputations", sep = " ")
      }
      B <- cbind(x[["Balance.Across.Imputations"]],
                 variable.names = row.names(x[["Balance.Across.Imputations"]]))
      facet <- NULL
    }
  }
  else if (any(class(x) == "bal.tab.imp")) {
    #Get imp.numbers.good
    imp.numbers <- seq_along(x[["Imputation.Balance"]])
    if (na.imp) {
      config <- "agg.imp"
      imp.numbers.good <- NULL
    }
    else if (null.imp) {
      config <- "agg.none"
      imp.numbers.good <- setNames(rep(TRUE, length(imp.numbers)), imp.numbers)
    }
    else if (is.numeric(attr(x, "print.options")$which.imp)) {
      imp.numbers.good <- setNames(imp.numbers %in% attr(x, "print.options")$which.imp, imp.numbers)
      if (any(imp.numbers.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.imp are valid imputation indices.", call. = FALSE)
      }
    }
    else stop("The argument to which.imp must be .none, .all, or the indices of imputations.", call. = FALSE)
    
    #Get B from x
    if (config == "agg.none") {
      B <- do.call("rbind", lapply(names(x[["Imputation.Balance"]])[imp.numbers.good], 
                                   function(i) cbind(x[["Imputation.Balance"]][[i]][["Balance"]],
                                                     imp = paste("Imputation:", i),
                                                     variable.names = rownames(x[["Imputation.Balance"]][[i]][["Balance"]]))))
      facet <- "imp"
    }
    else if (config == "agg.imp") {
      if (is_null(x[["Balance.Across.Imputations"]])) {
        stop("Cannot aggregate across imputations without a balance summary across imputations.\nThis may be because quick was set to TRUE and imp.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                        mean = "Mean", max = "Max", range = "Range")
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Imputations"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Imputations")
      }
      B <- cbind(x[["Balance.Across.Imputations"]],
                 variable.names = rownames(x[["Balance.Across.Imputations"]]))
    }
  }
  else if (any(class(x) == "bal.tab.cluster")) {
    #Get cluster.names.good
    cluster.names <- names(x[["Cluster.Balance"]])
    if (na.cluster) {
      config <- "agg.cluster"
      cluster.names.good <- NULL
    }
    else if (null.cluster) {
      config <- "agg.none"
      cluster.names.good <- setNames(rep(TRUE, length(cluster.names)), cluster.names)
    }
    else if (is.character(attr(x, "print.options")$which.cluster)) {
      cluster.names.good <- sapply(cluster.names, function(c) any(sapply(attr(x, "print.options")$which.cluster, function(y) identical(c, y))))
      if (any(cluster.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
      }
    }
    else if (is.numeric(attr(x, "print.options")$which.cluster)) {
      cluster.names.good <- setNames(seq_along(x$Cluster.Balance) %in% attr(x, "print.options")$which.cluster, cluster.names)
      if (any(cluster.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
      }
    }
    else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
    
    #Get B from x
    if (config == "agg.none") {
      B <- do.call("rbind", lapply(names(x[["Cluster.Balance"]])[cluster.names.good], 
                                   function(c) cbind(x[["Cluster.Balance"]][[c]][["Balance"]],
                                                     cluster = c,
                                                     variable.names = rownames(x[["Cluster.Balance"]][[c]][["Balance"]]))))
      facet <- "cluster"
    }
    else if (config == "agg.cluster") {
      if (is_null(x[["Cluster.Summary"]])) {
        stop("Cannot aggregate across clusters without a balance summary across clusters.\nThis may be because quick was set to TRUE and cluster.summary set to FALSE in the original bal.tab() call.", call. = FALSE)
      }
      Agg.Fun <- switch(tolower(match_arg(agg.fun, c("range", "max", "mean"))), 
                        mean = "Mean", max = "Max", range = "Range")
      if (Agg.Fun == "Range") {
        subtitle <- "Range Across Clusters"
      }
      else {
        subtitle <- paste(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), "Across Clusters")
      }
      B <- cbind(x[["Cluster.Summary"]],
                 variable.names = rownames(x[["Cluster.Summary"]]))
    }
  }
  else if (any(class(x) == "bal.tab.multi")) {
    #Get cluster.names.good
    treat.names <- attr(x, "print.options")$treat.names
    if (na.treat) {
      config <- "agg.pair"
      treat.names.good <- NULL
    }
    else if (null.treat) {
      config <- "agg.none"
      treat.names.good <- rep(TRUE, length(treat.names))
    }
    else if (is.character(attr(x, "print.options")$which.treat)) {
      treat.names.good <- treat.names %in% attr(x, "print.options")$which.treat
      if (any(treat.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.treat are valid treatment names or indices.", call. = FALSE)
      }
    }
    else if (is.numeric(attr(x, "print.options")$which.treat)) {
      treat.names.good <- seq_along(treat.names) %in% attr(x, "print.options")$which.treat
      if (any(treat.names.good)) {
        config <- "agg.none"
      }
      else {
        stop("Make sure the arguments to which.cluster are valid cluster names or indices.", call. = FALSE)
      }
    }
    else stop("The argument to which.cluster must be .none, .all, or the names or indices of clusters.", call. = FALSE)
    
    
    if (na.treat) {
      disp.treat.pairs <- character(0)
    }
    else {
      if (attr(x, "print.options")$pairwise) {
        if (sum(treat.names.good) == 1) {
          disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) any(attr(x[["Pair.Balance"]][[p]], "print.options")$treat.names == attr(x, "print.options")$treat.names[treat.names.good]))]
        }
        else {
          disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) all(attr(x[["Pair.Balance"]][[p]], "print.options")$treat.names %in% attr(x, "print.options")$treat.names[treat.names.good]))]
        }
      }
      else {
        if (sum(treat.names.good) == 1) {
          disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
            treat.names <- attr(x[["Pair.Balance"]][[p]], "print.options")$treat.namestreat.names
            any(treat.names[treat.names != "Others"] == attr(x, "print.options")$treat.names[treat.names.good])})]
        }
        else {
          disp.treat.pairs <- names(x[["Pair.Balance"]])[sapply(names(x[["Pair.Balance"]]), function(p) {
            treat.names <- attr(x[["Pair.Balance"]][[p]], "print.options")$treat.namestreat.names
            all(treat.names[treat.names != "Others"] %in% attr(x, "print.options")$treat.names[treat.names.good])})]
        }
      }
      
    }
    
    #Get B from x
    if (config == "agg.none") {
      B <- do.call("rbind", lapply(disp.treat.pairs,
                                   function(p) cbind(x[["Pair.Balance"]][[p]][["Balance"]],
                                                     treat.pair = p,
                                                     variable.names = rownames(x[["Pair.Balance"]][[p]][["Balance"]]))))
      facet <- "treat.pair"
    }
    else if (config == "agg.pair") {
      #Agg.Fun <- switch(match_arg(agg.fun), mean = "Mean", max = "Max", range = "Range")
      Agg.Fun <- "Max"
      if (Agg.Fun == "Range") {
        subtitle <- paste0("Range Across Treatment", ifelse(attr(x, "print.options")$pairwise, " Pairs", "s"))
      }
      else {
        subtitle <- paste0(ifelse(Agg.Fun == "Mean", "Average", Agg.Fun), " Across Treatment", ifelse(attr(x, "print.options")$pairwise, " Pairs", "s"))
      }
      B <- cbind(x[["Balance.Across.Pairs"]],
                 variable.names = rownames(x[["Balance.Across.Pairs"]]))
    }
  }
  else if (any(class(x) == "bal.tab.subclass")) {
    if (any(stats == "variance.ratios")) stop("Variance ratios not currently supported for subclassification.", call. = FALSE)
    if (any(stats == "ks.statistics")) stop("KS statistics not currently supported for subclassification.", call. = FALSE)
    if (any(class(x) == "bal.tab.cont")) stop("Continuous treatments not currently supported for subclassification.", call. = FALSE)
    subclass.names <- names(x[["Subclass.Balance"]])
    sub.B <- do.call("cbind", lapply(subclass.names, function(s) {
      sub <- x[["Subclass.Balance"]][[s]]
      sub.B0 <- setNames(sub[endsWith(names(sub), ".Adj")],
                         gsub(".Adj", paste0(".Subclass ", s), names(sub)[endsWith(names(sub), ".Adj")]))
      return(sub.B0) }))
    B <- cbind(x[["Balance.Across.Subclass"]], sub.B, variable.names = row.names(x[["Balance.Across.Subclass"]]))
    if (attr(x, "print.options")$disp.subclass) attr(x, "print.options")$weight.names <- c("Adj", paste("Subclass", subclass.names))
    else attr(x, "print.options")$weight.names <- "Adj"
    subtitle <- "Across Subclasses"
  }
  else {
    B <- cbind(x[["Balance"]], variable.names = row.names(x[["Balance"]]))
  }
  
  if (is_not_null(facet) && length(stats) > 1) {
    stop("stats can only have a length of 1 when faceting by other dimension (e.g., cluster, treatment).", call. = FALSE)
  }
  
  if (is_not_null(agg.fun) && config == "agg.none") {
    warning("No aggregation will take place, so agg.fun will be ignored. Remember to set 'which.<ARG> = .none' to aggregate across <ARG>.", call. = FALSE)
  }
  
  #Process abs
  if (config == "agg.none") {
    if (!abs && !attr(x, "print.options")[["abs"]]) {
      abs <- FALSE
    }
    else {
      if (!abs && attr(x, "print.options")[["abs"]])
        warning("abs is TRUE in the bal.tab object; ignoring abs in the call to love.plot().", call. = FALSE)
      abs <- TRUE
    }
  }
  else if (config %in% c("agg.time", "agg.pair")) {
    abs <- TRUE
  }
  else if (config %in% c("agg.imp", "agg.cluster", "agg.all")) {
    abs <- attr(x, "print.options")[["abs"]]
  }
  
  #Process variable names
  if (is_not_null(var.names)) {
    if (is.data.frame(var.names)) {
      if (ncol(var.names)==1) {
        if (is_not_null(row.names(var.names))) {
          new.labels <- setNames(unlist(as.character(var.names[,1])), rownames(var.names))
        }
        else warning("var.names is a data.frame, but its rows are unnamed.", call. = FALSE)
      }
      else {
        if (all(c("old", "new") %in% names(var.names))) {
          new.labels <- setNames(unlist(as.character(var.names[,"new"])), var.names[,"old"])
        }
        else {
          if (ncol(var.names)>2) warning("Only using first 2 columns of var.names", call. = FALSE)
          new.labels <- setNames(unlist(as.character(var.names[,2])), var.names[,1])
        }
      } 
    }
    else if (is.atomic(var.names) || is.factor(var.names)) {
      if (is_not_null(names(var.names))) {
        new.labels <- setNames(as.character(var.names), names(var.names))
      }
      else warning("var.names is a vector, but its values are unnamed.", call. = FALSE)
    }
    else if (is.list(var.names)) {
      if (all(sapply(var.names, function(x) is.character(x) || is.factor(x)))) {
        if (is_not_null(names(var.names))) {
          new.labels <- unlist(var.names) #already a list
        }
        else warning("var.names is a list, but its values are unnamed.", call. = FALSE)
      }
      else warning("var.names is a list, but its values are not the new names of the variables.", call. = FALSE)
    }
    else warning("Argument to var.names is not one of the accepted structures and will be ignored.\n  See help(love.plot) for details.", immediate.=TRUE, call. = FALSE)
    
    co.names <- attr(x, "print.options")[["co.names"]]
    seps <- attr(co.names, "seps")
    for (i in names(co.names)) {
      comp <- co.names[[i]][["component"]]
      type <- co.names[[i]][["type"]]
      
      if (i %in% names(new.labels) && !is.na(new.labels[i])) {
        co.names[[i]][["component"]] <- new.labels[i]
        co.names[[i]][["type"]] <- "base"
      }
      else {
        if ("isep" %in% type) {
          named.vars <- character(sum(type == "isep") + 1)
          sep.inds <- c(which(type == "isep"), length(comp) + 1)
          named.vars <- lapply(seq_along(sep.inds), function(k) {
            inds <- (if (k == 1) seq(1, sep.inds[k] - 1) 
                     else seq(sep.inds[k-1] + 1, sep.inds[k] - 1))
            var <- comp[inds]
            var.is.base <- type[inds] == "base"
            pasted.var <- paste(var, collapse = "")
            if (pasted.var %in% names(new.labels)) return(new.labels[pasted.var])
            else return(paste(ifelse(var.is.base & var %in% names(new.labels) & !is.na(new.labels[var]), new.labels[var], var), collapse = ""))
          })
          co.names[[i]][["component"]] <- do.call("paste", c(unname(named.vars), list(sep = seps["int"])))
        }
        else co.names[[i]][["component"]] <- ifelse(type == "base" & comp %in% names(new.labels) & !is.na(new.labels[comp]), new.labels[comp], comp)
      }
    }
    
    recode.labels <- setNames(names(co.names), 
                              vapply(co.names, function(x) paste0(x[["component"]], collapse = ""), character(1L)))
    
    B[["variable.names"]] <- do.call(f.recode, c(list(B[["variable.names"]]), as.list(recode.labels)))
  }
  
  distance.names <- as.character(unique(B[["variable.names"]][B[["Type"]] == "Distance"], nmax = sum(B[["Type"]] == "Distance")))
  if (drop.distance) {
    B <- B[B[["variable.names"]] %nin% distance.names, , drop = FALSE]
  }
  
  #Process variable order
  if (is_not_null(var.order) && "love.plot" %nin% class(var.order)) {
    if (attr(x, "print.options")$nweights == 0) {
      ua <- c("Unadjusted", "Alphabetical")
      names(ua) <- c("unadjusted", "alphabetical")
    }
    else if (attr(x, "print.options")$nweights == 1) {
      ua <- c("Adjusted", "Unadjusted", "Alphabetical")
      names(ua) <- c("adjusted", "unadjusted", "alphabetical")
    }
    else {
      ua <- c("Unadjusted", attr(x, "print.options")$weight.names, "Alphabetical")
      names(ua) <- c("unadjusted", attr(x, "print.options")$weight.names, "alphabetical")
    }
    var.order <- ua[match_arg(var.order, tolower(ua))]
  }
  
  #Process sample names
  ntypes <- if (is_null(subclass.names)) length(attr(x, "print.options")$weight.names) + 1 else 2
  if (!missing(sample.names)) {
    if (!is.vector(sample.names, "character")) {
      warning("The argument to sample.names must be a character vector. Ignoring sample.names.", call. = FALSE)
      sample.names <- NULL
    }
    else if (length(sample.names) %nin% c(ntypes, ntypes - 1)) {
      warning("The argument to sample.names must contain as many names as there are sample types, or one fewer. Ignoring sample.names.", call. = FALSE)
      sample.names <- NULL
    }
  }
  else sample.names <- NULL
  
  #Process limits
  if (is_not_null(limits)) {
    if (!is.vector(limits, "list")) {
      limits <- list(limits)
    }
    if (any(vapply(limits, 
                   function(l) !is.vector(l, "numeric") || length(l) %nin% c(0L, 2L), 
                   logical(1L)))) {
      warning("limits must be a list of numeric vectors of legnth 2. Ignoring limits.", call. = FALSE)
      limits <- NULL
    }
    
    if (is_not_null(names(limits))) {
      names(limits) <- stats[pmatch(names(limits), stats, duplicates.ok = TRUE)]
      limits <- limits[!is.na(names(limits))]
    }
    else {
      names(limits) <- stats[1:length(limits)]
    }
  }
  
  #Setting up appearance
  
  #Alpha (transparency)
  if (is.numeric(alpha[1]) && 
      !is.na(alpha[1]) && 
      between(alpha[1], c(0,1))) alpha <- alpha[1]
  else {
    warning("The argument to alpha must be a number between 0 and 1. Using 1 instead.", call. = FALSE)
    alpha <- 1
  }
  
  #Color
  if (is_not_null(args[["colours"]])) colors <- args[["colours"]]
  
  if (is_null(colors)) {
    if (shapes.ok(shapes, ntypes) && length(shapes) > 1 && length(shapes) == ntypes) {
      colors <- rep("black", ntypes)
    }
    else colors <- gg_color_hue(ntypes)
  }
  else {
    if (length(colors) == 1) colors <- rep(colors, ntypes)
    else if (length(colors) > ntypes) {
      colors <- colors[seq_len(ntypes)]
      warning(paste("Only using first", ntypes, "value", if (ntypes > 1) "s " else " ", "in colors."), call. = FALSE)
    }
    else if (length(colors) < ntypes) {
      warning("Not enough colors were specified. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(ntypes)
    }
    
    if (!all(sapply(colors, isColor))) {
      warning("The argument to colors contains at least one value that is not a recognized color. Using default colors instead.", call. = FALSE)
      colors <- gg_color_hue(ntypes)
    }
    
  }
  # colors[] <- vapply(colors, col_plus_alpha, character(1L), alpha = alpha)
  fill <- colors
  
  #Shapes
  if (is_null(shapes)) {
    shapes <- assign.shapes(colors)
  }
  else {
    #check shapes
    if (!shapes.ok(shapes, ntypes)) {
      warning(paste("The argument to shape must be", ntypes, "valid shape", if (ntypes > 1) "s." else ".", " See ?love.plot for more information.\nUsing default shapes instead."), call. = FALSE)
      shapes <- assign.shapes(colors)
    }
    else if (length(shapes) == 1) shapes <- rep(shapes, ntypes)
  }
  
  #Size
  if (is.numeric(size[1])) size <- size[1]
  else {
    warning("The argument to size must be a number. Using 3 instead.", call. = FALSE)
    size <- 3
  }
  stroke0 <- stroke <- rep(0, ntypes)
  size0 <- size <- rep(size, ntypes)
  
  shapes.with.fill <- grepl("filled", shapes, fixed = TRUE)
  stroke[shapes.with.fill] <- size[shapes.with.fill]/3
  size[shapes.with.fill] <- size[shapes.with.fill]* .58
  
  # stroke <- .8*size
  
  if (is_not_null(facet)) {
    if (is_not_null(var.order) && "love.plot" %nin% class(var.order) && tolower(var.order) != "alphabetical" && (sum(cluster.names.good) > 1 || sum(imp.numbers.good) > 1 || length(disp.treat.pairs) > 1 || sum(time.names.good) > 1)) {
      warning("var.order cannot be set with faceted plots (unless \"alphabetical\"). Ignoring var.order.", call. = FALSE)
      var.order <- NULL
    }
  }
  
  agg.range <- is_not_null(Agg.Fun) && Agg.Fun == "Range"
  
  #Process thresholds
  stat2threshold <- c(mean.diffs = "m.threshold",
                      variance.ratios = "v.threshold",
                      ks.statistics = "ks.threshold",
                      correlations = "r.threshold")
  thresholds <- setNames(sapply(stats, function(i) {
    if (is_not_null(attr(x, "print.options")[[stat2threshold[i]]])) {
      return(attr(x, "print.options")[[stat2threshold[i]]])
    }
    else {
      return(NA_real_)
    }
  }), stats)
  
  if (is_not_null(threshold)) {
    if (!all(is.na(threshold)) && !is.numeric(threshold)) {
      stop("threshold must be numeric.", call. = FALSE)
    }
    
    if (is_not_null(names(threshold))) {
      names(threshold) <- stats[pmatch(names(threshold), stats, duplicates.ok = TRUE)]
      threshold <- threshold[!is.na(names(threshold))]
    }
    else {
      names(threshold) <- stats[1:length(threshold)]
    }
    
    thresholds[names(threshold)] <- as.numeric(threshold)
  }
  thresholds <- thresholds[!is.na(thresholds)]
  
  #Title
  if (missing(title)) title <- "Covariate Balance"
  else title <- as.character(title)
  # if (missing(subtitle)) subtitle <- as.character(subtitle)
  
  #Process themes
  if (is_not_null(themes)) {
    if (!is.vector(themes, "list")) {
      themes <- list(themes)
    }
    if (any(vapply(themes, 
                   function(t) !all(c("theme", "gg") %in% class(t)), 
                   logical(1L)))) {
      warning("themes must be a list of \"theme\" objects. Ignoring themes.", call. = FALSE)
      themes <- NULL
    }
    
    if (is_not_null(names(themes))) {
      names(themes) <- stats[pmatch(names(themes), stats, duplicates.ok = TRUE)]
      themes <- themes[!is.na(names(themes))]
    }
    else {
      names(themes) <- stats[1:length(themes)]
    }
  }
  
  plot.list <- setNames(vector("list", length(stats)), stats)
  for (s in stats) {
    variable.names <- as.character(B[["variable.names"]])
    if (s == "mean.diffs") {
      binary <- attr(x, "print.options")$binary
      continuous <- attr(x, "print.options")$continuous
      #All std, no std, some std
      if ((binary == "std" || sum(B[["Type"]] == "Binary") == 0) && 
          (continuous == "std" || sum(B[["Type"]] != "Binary") == 0)) {
        xlab.diff <- "Standardized Mean Differences"
      } 
      else if ((binary == "raw" || sum(B[["Type"]] == "Binary") == 0) && 
               (continuous == "raw" || sum(B[["Type"]] != "Binary") == 0)) {
        xlab.diff <- "Mean Differences"
      }
      else {
        stars <- match_arg(stars, c("none", "std", "raw"))
        if (stars == "none") {
          warning("Standardized mean differences and raw mean differences are present in the same plot. \nUse the 'stars' argument to distinguish between them and appropriately label the x-axis.", call. = FALSE)
          xlab.diff <- "Mean Differences"
        }
        else {
          if (length(args$star_char) == 1 && is.character(args$star_char)) star_char <- args$star_char
          else star_char <- "*"
          
          vars_to_star <- setNames(rep(FALSE, nrow(B)), variable.names)
          if (stars == "std") {
            if (binary == "std") vars_to_star[variable.names[B[["Type"]] == "Binary"]] <- TRUE
            if (continuous == "std") vars_to_star[variable.names[B[["Type"]] != "Binary"]] <- TRUE
            xlab.diff <- "Mean Differences"
          }
          else if (stars == "raw") {
            if (binary == "raw") vars_to_star[variable.names[B[["Type"]] == "Binary"]] <- TRUE
            if (continuous == "raw") vars_to_star[variable.names[B[["Type"]] != "Binary"]] <- TRUE
            xlab.diff <- "Standardized Mean Differences"
          }
          variable.names[vars_to_star[variable.names]] <- paste0(variable.names[vars_to_star[variable.names]], star_char)
        }
      }
    }
    
    #Get SS
    if (agg.range) {
      SS <- do.call("rbind", 
                    lapply(c("Un", attr(x, "print.options")$weight.names),
                           function(w) data.frame(var = variable.names,
                                                  type = B[["Type"]],
                                                  min.stat = B[[paste.("Min", which.stat[s], w)]],
                                                  max.stat = B[[paste.("Max", which.stat[s], w)]],
                                                  mean.stat = B[[paste.("Mean", which.stat[s], w)]],
                                                  Sample = ifelse(w == "Un", "Unadjusted", 
                                                                  ifelse(w == "Adj", "Adjusted", w)),
                                                  row.names = NULL)))
      
      if (is_not_null(facet)) {
        if ("cluster" %in% facet) {
          SS$cluster <- rep(B[["cluster"]], 1 + attr(x, "print.options")$nweights)
        }
        if ("imp" %in% facet) {
          SS$imp <- rep(B[["imp"]], 1 + attr(x, "print.options")$nweights)
        }
        if ("treat.pair" %in% facet) {
          SS$treat.pair <- rep(B[["treat.pair"]], 1 + attr(x, "print.options")$nweights)
        }
      }
      
      
      if (all(sapply(SS[c("min.stat", "max.stat", "mean.stat")], is.na))) 
        stop(paste("No balance statistics to display. This can occur when", 
                   switch(s, variance.ratios = "disp.v.ratio", ks.statistics = "disp.ks"), 
                   "= FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
      
      missing.stat <- all(is.na(SS[["mean.stat"]]))
      if (missing.stat) stop(paste0(word_list(firstup(tolower(which.stat2))), 
                                    " cannot be displayed. This can occur when ", 
                                    word_list(paste.("disp", tolower(which.stat[s])), and.or = "and", is.are = TRUE), 
                                    " FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
      
      gone <- character(0)
      for (i in levels(SS$Sample)) {
        if (all(sapply(SS[SS$Sample == i, c("min.stat", "max.stat", "mean.stat")], is.na))) {
          gone <- c(gone, i)
          if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
          SS <- SS[SS[["Sample"]] != i,]
        }
      }
      
      dec <- FALSE
      
      if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
      
      if (is_not_null(var.order)) {
        if ("love.plot" %in% class(var.order)) {
          old.vars <- levels(var.order$data$var)
          old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
          if (any(SS[["var"]] %nin% old.vars)) {
            warning("The love.plot object in var.order doesn't have the same variables as the current input. Ignoring var.order.", call. = FALSE)
            var.order <- NULL
          }
          else {
            SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
          }
        }
        else if (tolower(var.order) == "alphabetical") {
          if ("time" %in% facet) {
            covnames0 <- vector("list", length(unique(SS[["time"]])))
            for (i in seq_along(covnames0)) {
              if (i == 1) {
                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
              }
              else {
                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
              }
            }
            covnames <- unlist(covnames0)
          }
          else covnames <- sort(levels(SS[["var"]]))
          SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
          
        }
        else if (var.order %in% ua) {
          if (var.order %in% gone) {
            warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
            var.order <- NULL
          }
          else {
            v <- as.character(SS[["var"]][order(SS[["mean.stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
            
            SS[["var"]] <- factor(SS[["var"]], 
                                  levels=c(v[v %nin% distance.names], 
                                           sort(distance.names, decreasing = TRUE)))
          }
        }
        
      }
      if (is_null(var.order)) {
        covnames <- as.character(unique(SS[["var"]]))
        SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
      }
      # SS[, "Sample"] <- factor(SS[, "Sample"], levels = c("Adjusted", "Unadjusted"))
      SS[["Sample"]] <- factor(SS[["Sample"]])
      if (s == "mean.diffs" && any(base::abs(SS[["max.stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
      if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["min.stat"]]),]
      SS[["stat"]] <- SS[["mean.stat"]]
    }
    else {
      SS <- do.call("rbind", 
                    lapply(c("Un", attr(x, "print.options")$weight.names),
                           function(w) data.frame(var = variable.names,
                                                  type = B[["Type"]],
                                                  stat = B[[ifelse(is_null(Agg.Fun), paste.(which.stat[s], w),
                                                                   paste.(Agg.Fun, which.stat[s], w))]],
                                                  Sample = ifelse(w == "Un", "Unadjusted", 
                                                                  ifelse(w == "Adj", "Adjusted", w)),
                                                  row.names = NULL)))
      
      if (is_not_null(facet)) {
        for (i in c("cluster", "imp", "treat.pair", "time")) {
          if (i %in% facet) {
            SS[[i]] <- rep(B[[i]], 1 + attr(x, "print.options")$nweights)
          }
        }
      }
      
      missing.stat <- all(is.na(SS[["stat"]]))
      if (missing.stat) stop(paste0(word_list(firstup(tolower(which.stat2))), 
                                    " cannot be displayed. This can occur when ", 
                                    word_list(paste.("disp", tolower(which.stat[s])), and.or = "and"), 
                                    " are FALSE and quick = TRUE in the original call to bal.tab()."), call. = FALSE)
      
      gone <- character(0)
      for (i in levels(SS$Sample)) {
        if (all(is.na(SS[["stat"]][SS$Sample==i]))) {
          gone <- c(gone, i)
          if (i == "Unadjusted") warning("Unadjusted values are missing. This can occur when un = FALSE and quick = TRUE in the original call to bal.tab().", call. = FALSE, immediate. = TRUE)
          SS <- SS[SS[["Sample"]]!=i,]
        }
      }
      
      if (abs) {
        if (s == "variance.ratios") SS[["stat"]] <- pmax(SS[["stat"]], 1/SS[["stat"]])
        else if (s == "mean.diffs") SS[["stat"]] <- base::abs(SS[["stat"]])
      }
      dec <- FALSE
      
      if (is_not_null(plot.list[[1]])) var.order <- plot.list[[1]]
      
      if (is_not_null(var.order)) {
        if ("love.plot" %in% class(var.order)) {
          old.vars <- levels(var.order$data$var)
          old.vars[endsWith(old.vars, "*")] <- substr(old.vars[endsWith(old.vars, "*")], 1, nchar(old.vars[endsWith(old.vars, "*")])-1)
          if (any(SS[["var"]] %nin% old.vars)) {
            warning("The love.plot object in var.order doesn't have the same variables as the current input. Ignoring var.order.", call. = FALSE)
            var.order <- NULL
          }
          else {
            SS[["var"]] <- factor(SS[["var"]], levels = old.vars[old.vars %in% SS[["var"]]])
          }
        }
        else if (tolower(var.order) == "alphabetical") {
          if ("time" %in% facet) {
            covnames0 <- vector("list", length(unique(SS[["time"]])))
            for (i in seq_along(covnames0)) {
              if (i == 1) {
                covnames0[[i]] <- sort(levels(SS[["var"]][SS[["time"]] == i]))
              }
              else {
                covnames0[[i]] <- sort(setdiff(levels(SS[["var"]][SS[["time"]] == i]), unlist(covnames0[seq_along(covnames0) < i])))
              }
            }
            covnames <- unlist(covnames0)
          }
          else covnames <- sort(levels(SS[["var"]]))
          SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[!covnames %in% distance.names]), sort(distance.names, decreasing = TRUE)))
          
        }
        else if (var.order %in% ua) {
          if (var.order %in% gone) {
            warning(paste0("var.order was set to \"", tolower(var.order), "\", but no ", tolower(var.order), " ", tolower(which.stat2), "s were calculated. Ignoring var.order."), call. = FALSE, immediate. = TRUE)
            var.order <- NULL
          }
          else {
            v <- as.character(SS[["var"]][order(SS[["stat"]][SS[["Sample"]]==var.order], decreasing = dec, na.last = FALSE)])
            
            SS[["var"]] <- factor(SS[["var"]], 
                                  levels=c(v[v %nin% distance.names], 
                                           sort(distance.names, decreasing = TRUE)))
          }
        }
        
      }
      if (is_null(var.order)) {
        covnames <- as.character(unique(SS[["var"]])) #Don't use levels here to preserve original order
        SS[["var"]] <- factor(SS[["var"]], levels = c(rev(covnames[covnames %nin% distance.names]), sort(distance.names, decreasing = TRUE)))
      }
      SS[["Sample"]] <- factor(SS[["Sample"]])
      if (s == "mean.diffs" && any(base::abs(SS[["stat"]]) > 5, na.rm = TRUE)) warning("Large mean differences detected; you may not be using standardized mean differences for continuous variables.", call.=FALSE)
      if (length(stats) == 1 && drop.missing) SS <- SS[!is.na(SS[["stat"]]),]
    }
    
    SS <- SS[order(SS[["var"]], na.last = FALSE),]
    SS[["var"]] <- factor(SS[["var"]])
    
    #Make the plot
    #library(ggplot2)
    
    baseline.xintercept <- switch(s, 
                                  "mean.diffs" = 0, 
                                  "variance.ratios" = 1, 
                                  "ks.statistics" = 0, 
                                  "correlations" = 0)
    threshold.xintercepts <- NULL
    if (s == "correlations") {
      if (abs) {
        xlab <- "Absolute Treatment-Covariate Correlations"
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
      }
      else {
        xlab <- "Treatment-Covariate Correlations"
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = -thresholds[s], upper = thresholds[s])
      }
      scale_Statistics <- scale_x_continuous
    }
    else if (s == "mean.diffs") {
      if (abs) {
        xlab <- paste("Absolute", xlab.diff)
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
      }
      else {
        xlab <- xlab.diff
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = -thresholds[s], upper = thresholds[s])
      }
      scale_Statistics <- scale_x_continuous
    }
    else if (s == "variance.ratios") {
      xlab <- "Variance Ratios"
      if (abs) {
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = max(thresholds[s], 1/thresholds[s]))
      }
      else {
        if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = min(c(thresholds[s], 1/thresholds[s])), upper = max(c(thresholds[s], 1/thresholds[s])))
      }
      scale_Statistics <- scale_x_log10
    }
    else if (s == "ks.statistics") {
      xlab <- "Kolmogorov-Smirnov Statistics"
      if (s %in% names(thresholds)) threshold.xintercepts <- c(lower = base::abs(thresholds[s]))
      scale_Statistics <- scale_x_continuous
    }
    
    apply.limits <- FALSE
    SS[["on.border"]] <- FALSE
    if (is_not_null(limits[[s]])) {
      if (limits[[s]][2] < limits[[s]][1]) {
        limits[[s]] <- c(limits[[s]][2], limits[[s]][1])
      }
      
      if (limits[[s]][1] >= baseline.xintercept) limits[[s]][1] <- baseline.xintercept - .05*limits[[s]][2]
      if (limits[[s]][2] <= baseline.xintercept) limits[[s]][2] <- baseline.xintercept - .05*limits[[s]][1]
      
      if (identical(scale_Statistics, scale_x_log10)) limits[[s]][limits[[s]] <= 1e-2] <- 1e-2
      
      if (agg.range) {
        # for (i in c("min.stat", "mean.stat", "max.stat")) {
        #     if (any(SS[[i]] < limits[[s]][1], na.rm = TRUE)) {
        #         if (i == "mean.stat") SS[["on.border"]][SS[[i]] < limits[[s]][1]] <- TRUE
        #         if (i == "mean.stat") SS[[i]][SS[[i]] < limits[[s]][1]] <- limits[[s]][1]
        #         else SS[[i]][SS[[i]] < limits[[s]][1]] <- limits[[s]][1]
        #     }
        #     if (any(SS[[i]] > limits[[s]][2], na.rm = TRUE)) {
        #         if (i == "mean.stat") SS[["on.border"]][SS[[i]] > limits[[s]][2]] <- TRUE
        #         if (i == "mean.stat") SS[[i]][SS[[i]] > limits[[s]][2]] <- limits[[s]][2]
        #         else SS[[i]][SS[[i]] > limits[[s]][2]] <- limits[[s]][2]
        #         
        #         # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
        #     }
        # }
        
        if (any(SS[["mean.stat"]] < limits[[s]][1], na.rm = TRUE)) {
          SS[["on.border"]][SS[["mean.stat"]] < limits[[s]][1]] <- TRUE
          SS[["mean.stat"]][SS[["mean.stat"]] < limits[[s]][1]] <- limits[[s]][1]
          SS[["max.stat"]][SS[["max.stat"]] < limits[[s]][1]] <- limits[[s]][1]
          SS[["min.stat"]][SS[["min.stat"]] < limits[[s]][1]] <- limits[[s]][1]
        }
        if (any(SS[["mean.stat"]] > limits[[s]][2], na.rm = TRUE)) {
          SS[["on.border"]][SS[["mean.stat"]] > limits[[s]][2]] <- TRUE
          SS[["mean.stat"]][SS[["mean.stat"]] > limits[[s]][2]] <- limits[[s]][2]
          SS[["max.stat"]][SS[["max.stat"]] > limits[[s]][2]] <- limits[[s]][2]
          SS[["min.stat"]][SS[["min.stat"]] > limits[[s]][2]] <- limits[[s]][2]
          # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
        }
        # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
      }
      else {
        if (any(SS[["stat"]] < limits[[s]][1], na.rm = TRUE)) {
          SS[["on.border"]][SS[["stat"]] < limits[[s]][1]] <- TRUE
          SS[["stat"]][SS[["stat"]] < limits[[s]][1]] <- limits[[s]][1]
        }
        if (any(SS[["stat"]] > limits[[s]][2], na.rm = TRUE)) {
          SS[["on.border"]][SS[["stat"]] > limits[[s]][2]] <- TRUE
          SS[["stat"]][SS[["stat"]] > limits[[s]][2]] <- limits[[s]][2]
          # warning("Some points will be removed from the plot by the limits.", call. = FALSE)
        }
      }
      
      apply.limits <- TRUE
    }
    
    if (is_not_null(sample.names)) {
      if (length(sample.names) == ntypes - 1) {
        levels(SS$Sample)[-1] <- sample.names
      }
      else if (length(sample.names) == ntypes) {
        levels(SS$Sample) <- sample.names
      }
    }
    
    lp <- ggplot(aes(y = var, x = stat, group = Sample), data = SS) +
      theme(panel.background = element_rect(fill = "white"),
            axis.text.x = element_text(color = "black"),
            axis.text.y = element_text(color = "black"),
            panel.border = element_rect(fill = NA, color = "black"),
            plot.background = element_blank(),
            legend.background = element_blank(),
            legend.key = element_blank()
      ) +
      scale_shape_manual(values = shapes) +
      scale_size_manual(values = size) +
      scale_discrete_manual(aesthetics = "stroke", values = stroke) +
      scale_fill_manual(values = fill) +
      scale_color_manual(values = colors) +
      labs(y = NULL, x = wrap(xlab, wrap))
    
    lp <- lp + geom_vline(xintercept = baseline.xintercept,
                          linetype = 1, color = "gray5")
    
    if (is_not_null(threshold.xintercepts)) {
      lp <- lp + geom_vline(xintercept = threshold.xintercepts,
                            linetype = 2, color = "gray8")
    }
    
    if (agg.range) {
      position.dodge <- ggstance::position_dodgev(.5*(size0[1]/3))
      if (line == TRUE) { #Add line except to distance
        f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
        lp <- lp + ggplot2::layer(geom = "path", data = f, 
                                  position = position.dodge, 
                                  stat = "identity",
                                  mapping = aes(x = mean.stat, color = Sample), 
                                  params = list(size = size0[1]*.8/3, na.rm = TRUE,
                                                alpha = alpha))
      }
      
      lp <- lp +
        geom_errorbarh(aes(y = var, xmin = min.stat, xmax = max.stat,
                           color = Sample), position = position.dodge,
                       size = size0[1]*.8/3,
                       alpha = alpha, 
                       height = 0,
                       na.rm = TRUE) +
        geom_point(aes(y = var, 
                       x = mean.stat, 
                       shape = Sample,
                       size = Sample,
                       stroke = Sample,
                       color = Sample),
                   fill = "white", na.rm = TRUE,
                   alpha = alpha,
                   position = position.dodge)
      
    }
    else {
      if (is_null(subclass.names) || !attr(x, "print.options")$disp.subclass) {
        if (line == TRUE) { #Add line except to distance
          f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
          lp <- lp + ggplot2::layer(geom = "path", data = f(SS),
                                    position = "identity", stat = "identity",
                                    mapping = aes(color = Sample),
                                    params = list(size = size0[1]*.8/3,
                                                  na.rm = TRUE,
                                                  alpha = alpha))
        }
        lp <- lp + geom_point(data = SS, aes(shape = Sample,
                                             size = Sample,
                                             stroke = Sample,
                                             color = Sample),
                              fill = "white", 
                              na.rm = TRUE,
                              alpha = alpha)
        
      }
      else {
        SS.u.a <- SS[SS$Sample %in% c("Unadjusted", "Adjusted"),]
        SS.u.a$Sample <- factor(SS.u.a$Sample)
        if (line == TRUE) { #Add line except to distance
          f <- function(q) {q[["stat"]][q$type == "Distance"] <- NA; q}
          lp <- lp + ggplot2::layer(geom = "path", data = f(SS.u.a),
                                    position = "identity", stat = "identity",
                                    mapping = aes(color = Sample),
                                    params = list(size = size*.8,
                                                  na.rm = TRUE,
                                                  alpha = alpha))
        }
        lp <- lp + geom_point(data = SS.u.a,
                              aes(shape = Sample,
                                  size = Sample,
                                  stroke = Sample,
                                  color = Sample),
                              fill = "white",
                              na.rm = TRUE)
        lp <- lp + geom_text(data = SS[SS$Sample %nin% c("Unadjusted", "Adjusted"),],
                             mapping = aes(label = gsub("Subclass ", "", Sample)),
                             size = 2.5*size0[1]/3, na.rm = TRUE)
      }
      
      
    }
    
    if (!drop.distance && is_not_null(distance.names)) {
      lp <- lp + geom_hline(linetype = 1, color = "black",
                            yintercept = nunique(SS[["var"]]) - length(distance.names) + .5)
    }
    if (apply.limits) {
      lp <- lp + scale_Statistics(limits = limits[[s]], expand = c(0, 0))
    }
    else {
      lp <- lp + scale_Statistics()
    }
    
    if (isFALSE(grid)) {
      lp <- lp + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank())
    }
    else {
      lp <- lp + theme(panel.grid.major = element_line(color = "gray87"),
                       panel.grid.minor = element_line(color = "gray90"))
    }
    
    if (is_not_null(facet)) {
      lp <- lp + facet_grid(f.build(".", facet), drop = FALSE) + labs(x = xlab)
    }
    
    class(lp) <- c(class(lp), "love.plot")
    plot.list[[s]] <- lp
  }
  
  if (length(stats) > 1 || isTRUE(args$use.grid)) {
    
    if (!is.character(position) || length(position) > 1) {
      position <- NA_character_
    }
    position <- match_arg(as.character(position), 
                          c("right", "left", "top", "bottom", "none"))
    
    #Process labels
    if (isTRUE(labels)) labels <- LETTERS[seq_along(plot.list)]
    else if (is_null(labels) || isFALSE(labels)) labels <- NULL
    else if (!is.atomic(labels) || length(labels) != length(plot.list)) {
      warning("labels must be TRUE or a string with the same length as stats. Ignoring labels.", call. = FALSE)
      labels <- NULL
    }
    else labels <- as.character(labels)
    
    # p <- ggpubr::ggarrange(plotlist = plot.list, common.legend = TRUE, legend = position, 
    #                align = "hv", nrow = 1)
    # if (is_not_null(subtitle)) {
    #     p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(subtitle, size = 11))
    # }
    # p <- ggpubr::annotate_figure(p, top = ggpubr::text_grob(title, size = 13.2))
    # 
    # P <- attr(P, "plots")
    
    plots.to.combine <- plot.list
    for (i in seq_along(plots.to.combine)) {
      if (i > 1) {
        plots.to.combine[[i]] <- plots.to.combine[[i]] + 
          theme(axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                legend.position = "none")
      }
      else {
        plots.to.combine[[i]] <- plots.to.combine[[i]] + theme(legend.position = "none")
      }
      
      if (is_not_null(labels)) {
        plots.to.combine[[i]] <- plots.to.combine[[i]] + labs(title = labels[i])
      }
      
      if (is_not_null(themes[[stats[i]]])) {
        plots.to.combine[[i]] <- plots.to.combine[[i]] + themes[[stats[i]]]
      }
    }
    
    g <- ggarrange_simple(plots = plots.to.combine, nrow = 1)
    title.grob <- grid::textGrob(title, gp=grid::gpar(fontsize=13.2))
    subtitle.grob <- grid::textGrob(subtitle, gp=grid::gpar(fontsize=13.2))
    
    if (position == "none") {
      p <- gridExtra::arrangeGrob(grobs = list(g), nrow = 1)
    }
    else {
      legg <- ggplot2::ggplotGrob(plots.to.combine[[1]] + theme(legend.position = position))
      leg <- legg$grobs[[which(legg$layout$name == "guide-box")]]
      
      if (position == "left") {
        p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 1, 
                                    widths = grid::unit.c(sum(leg$widths), grid::unit(1, "npc") - sum(leg$widths)))
      }
      else if (position == "right") {
        p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 1, 
                                    widths = grid::unit.c(grid::unit(1, "npc") - sum(leg$widths), sum(leg$widths)))
      }
      else if (position == "top") {
        p <- gridExtra::arrangeGrob(grobs = list(leg, g), nrow = 2,
                                    heights = grid::unit.c(sum(leg$heights), grid::unit(1, "npc") - sum(leg$heights)))
      }
      else if (position == "bottom") {
        p <- gridExtra::arrangeGrob(grobs = list(g, leg), nrow = 2,
                                    heights = grid::unit.c(grid::unit(1, "npc") - sum(leg$heights), sum(leg$heights)))
      }
    }
    
    if (is_not_null(subtitle)) {
      p <- gridExtra::arrangeGrob(p, top = subtitle.grob)
    }
    p <- gridExtra::arrangeGrob(p, top = title.grob)
    
    grid::grid.newpage()
    grid::grid.draw(p)
    
    attr(p, "plots") <- plot.list
    class(p) <- c(class(p), "love.plot")
    
    return(invisible(p))
  }
  else {
    
    p <- plot.list[[1]] + 
      labs(title = title, subtitle = subtitle) +
      theme(plot.title = element_text(hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5),
            legend.position = position)
    
    if (is_not_null(themes[[1]])) {
      p <- p + themes[[1]]
    }
    
    return(p)
    
  }
  
}
plot.bal.tab <- love.plot


scale_discrete_manual <- function(aesthetics, ..., values, breaks = waiver()) {
  manual_scale(aesthetics, values, breaks, ...)
}


manual_scale <- function(aesthetic, values = NULL, breaks = waiver(), ..., limits = NULL) {
  # check for missing `values` parameter, in lieu of providing
  # a default to all the different scale_*_manual() functions
  if (is_missing(values)) {
    values <- NULL
  } else {
    force(values)
  }
  
  if (is.null(limits)) {
    limits <- names(values)
  }
  
  # order values according to breaks
  if (is.vector(values) && is.null(names(values)) && !is.waive(breaks) &&
      !is.null(breaks) && !is.function(breaks)) {
    if (length(breaks) <= length(values)) {
      names(values) <- breaks
    } else {
      names(values) <- breaks[1:length(values)]
    }
  }
  
  pal <- function(n) {
    if (n > length(values)) {
      abort(glue("Insufficient values in manual scale. {n} needed but only {length(values)} provided."))
    }
    values
  }
  discrete_scale(aesthetic, "manual", pal, breaks = breaks, limits = limits, ...)
}
