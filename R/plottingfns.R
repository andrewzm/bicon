#' @title Plot meshes and observations
#' @description This function is the general one used for plotting/visualising objects mesh objects.
#' @param x a mesh object
#' @param ... other parameters used to configure the plot. These include \code{max} (upperbound on colour scale) and \code{min} (lowerbound on colour scale)
#' @return a ggplot2 object
#' @export
#' @examples 
#' ## See vignette min_max_T
setGeneric("plot")


#' @title Plots an interpolated field from a mesh
#' @docType methods
#' @description This function takes a mesh and a column name in the mesh, to generate an interpolated
#' field of the indicated column
#' @param x a mesh object
#' @param y a character indicated which column to plot
#' @param ds the resolution of the plot (defaulted to 40 x 40)
#' @param max upperbound on colour scale
#' @param min lowerbound on colour scale
#' @return a ggplot2 object
#' @export
#' @examples 
#' ## See vignette min_max_T
setGeneric("plot_interp", function(x,y,ds,...) standardGeneric("plot_interp"))

.check_plot_args <- function(args,z=0) {
  
  if ("max" %in% names(args)) {
    stopifnot(is.numeric(args$max))
  } else {
    args$max <- max(z,na.rm=T) 
  }
  
  if ("min" %in% names(args)) {
    stopifnot(is.numeric(args$min))
  } else {
    args$min <- min(z,na.rm=T)
  }
  
  if("leg_title" %in% names(args)) {
    stopifnot(is.character(args$leg_title))
  } else {
    args$leg_title <- ""
  }
  
  if ("plot_dots" %in% names(args)) {
    stopifnot(is.logical(args$plot_dots))
  } else {
    args$plot_dots <- T
  }
  
  if ("g" %in% names(args)) {
    stopifnot(is.ggplot(args$g))    
  } else {
    args$g <- NULL
  }
  
  if("pt_size" %in% names(args)) {
    stopifnot(is.numeric(args$pt_size))
  } else {
    args$pt_size <- 1
  }
  
  if("size" %in% names(args)) {
    args$pt_size <- args$size
  } 
  
  if(!("palette" %in% names(args))) {
    args$palette <- "RdYlBu"
  } 
  
  if(!("reverse" %in% names(args))) {
    args$reverse <- FALSE
  } 
  
  return(args)
}


#' @rdname plot_interp
#' @aliases plot_interp,FEBasis,character-method
setMethod("plot_interp",signature(x = "FEBasis",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,ds=40,...) {
            args <- list(...)
            Mesh <- x
            df <- getDf(Mesh)
            Mesh['to_plot'] <- df[y]
            args <- .check_plot_args(args,df[y])
            zhi <- args$max
            zlo <- args$min 
            leg_title <- args$leg_title
            palette <- args$palette
            reverse <- args$reverse
            
            mesh_grid <- interp(Mesh["x"],Mesh["y"],Mesh["to_plot"],
                                xo=seq(min(Mesh["x"]),max(Mesh["x"]),length=ds), 
                                yo=seq(min(Mesh["y"]),max(Mesh["y"]),length=ds))
            expanded_grid <- cbind(expand.grid(mesh_grid$x,mesh_grid$y),as.vector(mesh_grid$z))
            names(expanded_grid) <- c("x","y","z")
            g <- OverlayPlot(expanded_grid,GGstd = NULL,leg_title=leg_title,zlo=zlo,zhi=zhi,palette=palette,reverse=reverse)
            return(g)
          })

#' @rdname plot_interp
#' @aliases plot_interp,FEBasis-method
setMethod("plot",signature(x = "FEBasis"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)
            args <- .check_plot_args(args)
            plot_dots <- args$plot_dots
            g <- args$g            
            g <- ggplotGraph(as.matrix(x@pars$K - diag(diag(x@pars$K))),fixed=x@pars$p,point_size=0.5,g=g,plot_dots=plot_dots) 
            return(g)
          })

#' @rdname plot
#' @aliases plot,FEBasis,character-method
setMethod("plot",signature(x = "FEBasis",y="character"),  # GRBF basis with mean offset as last weight
          function(x,y,...) {
            args <- list(...)
            
            xlocs <- x@pars$p[,1]
            ylocs <- x@pars$p[,2]
            if(!(y %in% names(x@pars$vars))) stop("Did not find label in mesh data frame")
            z <- x@pars$vars[[y]]
            df <- data.frame(x = xlocs, y=ylocs, z=z)
            args <- .check_plot_args(args,z)
            pt_size <- args$pt_size
            df$p1 <- args$max 
            df$p2 <- args$min 
            leg_title <- args$leg_title
            
            g <- LinePlotTheme() + geom_point(data=df,aes(x,y,colour=pmax(pmin(z,p1),p2)),size=pt_size) + 
              scale_colour_gradient2(low="red",mid="light yellow",high="blue",
                                     guide=guide_legend(title=leg_title,reverse=T))
            return(g)
          })


ggplotGraph <- function(Q,fixed=NA,point_size=0.5,g = NULL,plot_dots=T) {
  
  net <- network(Q,directed = F,density=0.3)
  m <- as.matrix.network.adjacency(net) # get sociomatrix
  plotcord <- data.frame(X1=fixed[,1],X2=fixed[,2])
  edglist <- as.matrix.network.edgelist(net)
  edges <- data.frame(plotcord[edglist[,1],], plotcord[edglist[,2],])
  colnames(edges) <-  c("X1","Y1","X2","Y2")
  edges$midX  <- (edges$X1 + edges$X2) / 2
  edges$midY  <- (edges$Y1 + edges$Y2) / 2
  if (class(g)[1] == "NULL") {
    pnet <- ggplot()  +
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),
                   data=edges, size = 0.3, colour="dark grey") +
      scale_colour_brewer(palette="Set1") +
      scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
      # discard default grid + titles in ggplot2
      theme(panel.background = element_blank()) + theme(legend.position="none")+
      theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
      theme( legend.background = element_rect(colour = NA)) +
      theme(panel.background = element_rect(fill = "white", colour = NA)) +
      theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
  } else {
    pnet <- g  +
      geom_segment(aes(x=X1, y=Y1, xend = X2, yend = Y2),
                   data=edges, size = 0.3, colour="dark grey") +
      scale_colour_brewer(palette="Set1")
    #scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) 
    # discard default grid + titles in ggplot2
    #theme(panel.background = element_blank(), legend.position="none",
    #      legend.background = element_rect(colour = NA), panel.background = element_rect(fill = "white", colour = NA),
    #      panel.grid.minor = element_blank(), panel.grid.major = element_blank()) 
  }
  if(plot_dots) pnet <- pnet + geom_point(aes(X1, X2),colour="red", data=plotcord,size=point_size)
  return(pnet)
  
}

## Plot a spatial field from vertices
PlotSpatField <- function(p,x,my_title = '',ylab="",zlim=c(0,0),palette=brewer.pal(11,"RdBu"),res=200)  {
  surf <- interp(x=p[,1],y=p[,2],z=x,xo=seq(min(p[,1]), max(p[,1]), length = res),yo=seq(min(p[,2]), max(p[,2]), length = res))
  if ((zlim[1] == 0) && (zlim[2] == 0)) zlim <- range(surf$z,na.rm=T)
  zlen <- zlim[2] - zlim[1]
  
  surf$z[which(surf$z < zlim[1])] = zlim[1]
  surf$z[which(surf$z > zlim[2])] = zlim[2]
  
  
  image(x=surf$x, y=surf$y, z=surf$z,col=palette, axes=F, zlim=zlim, xlab="", ylab=ylab)
  image.plot( zlim=zlim, col=palette, legend.only=TRUE, horizontal =TRUE)
  grid()
  box()
  # filled.contour(surf$x, surf$y, surf$z,color.palette=terrain.colors,xlab='x',ylab='y')
  title(my_title)
  return(surf)
}


OverlayPlot <- function(GG,GGstd = NULL,leg_title="",zlo = -0.6,zhi = 0.6,alphalo = 0.2, alphahi = 1,palette="RdYlBu",reverse=F) {
  GG$zlo = zlo
  GG$zhi = zhi
  if (!is.null(GGstd)) {
    GGstd$alphalo = alphalo
    GGstd$alphahi = alphahi
  }
  
  g <- LinePlotTheme() + geom_tile(data=subset(GG,!(is.na(z))),aes(x=x,y=y,fill=pmin(pmax(z,zlo),zhi))) +
    #scale_fill_gradient2(low=(muted("red")),mid="light yellow",high=(muted("blue")),limits=c(zlo,zhi),
    #                     guide=guide_colourbar(title=leg_title)) + 
    scale_fill_distiller(palette=palette,reverse=reverse,limits=c(zlo,zhi),
                         guide=guide_colourbar(title=leg_title,position=c(-1500,500)))
  
  if (!is.null(GGstd)) {
    g <- g + geom_tile(data=GGstd,aes(x=x,y=y, alpha=-pmax(pmin(z,alphahi),alphalo)),fill="white") +scale_alpha(guide = 'none')
  }
  
  
  g <- g + theme(legend.position="right") + 
    theme(legend.key.width=unit(4,"lines")) + theme(legend.key.height=unit(4,"lines")) +
    theme(text = element_text(size=40))
  return(g)
}

#' @title Line-plot theme
#' 
#' @description Formats a ggplot object for neat plotting.
#'
#' @return Object of class \code{ggplot}
#' @keywords ggplot
#' @export
#' @examples
#' \dontrun{
#' X <- data.frame(x=runif(100),y = runif(100), z = runif(100))
#' LinePlotTheme() + geom_point(data=X,aes(x,y,colour=z))
#'}
LinePlotTheme <- function () 
{
    ggplot() + theme(panel.background = element_rect(fill = "white", colour = "black"), 
                     text = element_text(size = 20), 
                     panel.grid.major = element_line(colour = "light gray", size = 0.05), 
                     panel.border = element_rect(fill = NA, colour = "black"))
}

scale_colour_brewer <- function(..., type = "seq", palette = 1) {
    discrete_scale("colour", "brewer", brewer_pal(type, palette), ...)
}

scale_fill_brewer <- function(..., type = "seq", palette = 1) {
    discrete_scale("fill", "brewer", brewer_pal(type, palette), ...)
}

scale_colour_distiller <- function(..., type = "seq", palette = 1, values = NULL, space = "Lab", na.value = "grey50") {
    # warn about using a qualitative brewer palette to generate the gradient
    type <- match.arg(type, c("seq", "div", "qual"))
    if (type == "qual") {
        warning("Using a discrete colour palette in a continuous scale.\n  Consider using type=\"seq\" or type=\"div\" instead")
    }
    continuous_scale("colour", "distiller",
                     gradient_n_pal(brewer_pal(type, palette)(6), values, space), na.value = na.value, ...)
    # NB: 6 colours per palette gives nice gradients; more results in more saturated colours which do not look as good
}

scale_fill_distiller <- function(..., type = "seq", palette = 1, values = NULL, space = "Lab", na.value = "grey50",reverse=F) {
    type <- match.arg(type, c("seq", "div", "qual"))
    if (type == "qual") {
        warning("Using a discrete colour palette in a continuous scale.\n  Consider using type=\"seq\" or type=\"div\" instead")
    }
    pal <- brewer_pal(type, palette)(6)
    if(reverse) pal <- rev(pal)
    
    continuous_scale("fill", "distiller",
                     gradient_n_pal(pal, values, space), na.value = na.value, ...)
}

