weight.edges <-
function(nb.object, coords, distmat=NULL, alpha=2, beta=NULL, max.d=NULL,
	unit.angle="degrees", rot.angle=0, rm.same.y=TRUE, plot.connexions=TRUE)
{
### Olivier Gauthier, Pierre Legendre and F. Guillaume Blanchet - August 2013
##########################################################################################
	### General check-up
	if(ncol(coords)!=3){
		stop("'coords' needs to have three columns")
	}
	if(is.character(coords[,1])){
		stop("The first column of 'coords' needs to be numeric")
	}else{
		coords<-as.matrix(coords)
	}
	
	if(rot.angle == 0){
		if (unit.angle == "degrees" | unit.angle == "radians") {
			rot.angle<-0
		}else{
			stop("Units for angles must be either 'degrees' or 'radians'")
		}
	}
	
	if(is.null(distmat)){
		links.length.mat <- as.matrix(dist(coords[,2:3]))
	}else{
		links.length.mat <- as.matrix(distmat)
	}
	if(is.null(max.d)){	
		max.d <- max(links.length.mat)/2
		print("The default value for 'max.d' was used. The weights may not be optimum")
	}
		
	### Construct the links matrix
	links <- listw2sn(nb2listw(nb.object, zero.policy = TRUE))[, 1:2]
	links <- rm.double.link(links)
	
	### Rotate the coordinates
	if(missing(unit.angle)){
		unit.angle<-"degrees"
	}
	
	if(rot.angle != 0){
		if (unit.angle == "degrees") {
			rot.angle <- pi/180 * rot.angle
		}else{
			if (unit.angle == "radians") {
				rot.angle <- rot.angle
			}else{
				stop("Units for angles must be either 'degrees' or 'radians'")
			}
		}
	}
	
	coords[, 2:3] <- round(rotation(as.matrix(coords[, 2:3]), rot.angle), digits = 8)
	
	### Remove edges perpendicular to the directional gradient
	if(rm.same.y == TRUE) {
		links <- remove.same.y(coords = coords, link = links)
	}
	
	low.y <- which(coords[, 3] == min(coords[, 3]))
	
	### Plot connexion diagram
	if(plot.connexions) {
		xy.range<-apply(coords[,2:3],2,range)
		
		xy.range.min<-xy.range[1,2]-((xy.range[2,2]-xy.range[1,2])/5)
		xy.range.max<-xy.range[2,2]

		par(mar=c(1,1,1,1))
		plot(coords[, 2:3], pch = 20, asp = 1, cex = 0.5, axes = FALSE, xlab = "", ylab = "", ylim=c(xy.range.min,xy.range.max))
		segments(x0 = coords[links[, 1], 2], y0 = coords[links[, 1], 3], x1 = coords[links[, 2], 2], y1 = coords[links[, 2], 3], col = "red")
		text(coords[,2:3],labels=as.character(coords[,1]),pos=2,col="red")

		site0<-c(mean(coords[,2]),xy.range.min)
		
		points(site0[1],site0[2],pch=19,col="blue")
		segments(x0=site0[1],y0=site0[2],x1=coords[low.y,2],y1=coords[low.y,3],col="blue")
		par(mar=c(5,4,4,2))
	}
	
	### Reorganize the edges 
	n.low.y <- length(low.y)
	links <- rbind(cbind(rep(0, n.low.y), low.y), as.matrix(links))
	nrow.links <- nrow(links)
	links <- cbind(1:nrow.links, links)
	links <- r.order.link(nrow.links, links, coords)
	links <- links[, 2:3]
	
	### Remove edges directly linked with site 0
	linkstorm<-unique(which(links==0,arr.ind=TRUE)[,1])
	links<-links[-linkstorm,]
	nrow.links <- nrow(links)
	
	### Find the length of each edge
	links.length <- numeric(length=nrow.links)
	for(i in 1:nrow.links){
		links.length[i] <- links.length.mat[links[i,1],links[i,2]]
	}
	
	### Calculate the weight associated to each edge
	if (is.null(beta)) {
		# First weighting function, Legendre & Legendre (2012, eq. 14.3)
		w <- 1 - (links.length/max.d)^alpha
		} else {
		# Second weighting function, Legendre & Legendre (2012, eq. 14.4)
		w <- 1 / links.length^beta
		}
	w
}
