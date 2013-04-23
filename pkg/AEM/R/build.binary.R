`build.binary` <-
function (nb.object, coords, unit.angle="degrees", rot.angle = 0, rm.same.y = TRUE, plot.connexions = TRUE) {
	require(spdep)
	if(is.character(coords[,1])){
		stop("The first column of 'coords' needs to be numeric")
	}else{
		coords<-as.matrix(coords)
	}
	link <- listw2sn(nb2listw(nb.object, zero.policy = TRUE))[, 1:2]
	link <- rm.double.link(link)
	n <- nrow(coords)
	if(missing(unit.angle)){
		unit.angle<-"degrees"
	}
	if (rot.angle == 0){
		if (unit.angle == "degrees" | unit.angle == "radians") {
			rot.angle<-0
		}else{
			stop("Unit for angles be either 'degrees' or 'radians'")
		}
	}
	if (rot.angle != 0){
		if (unit.angle == "degrees") {
			rot.angle <- pi/180 * rot.angle
		}else{
			if (unit.angle == "radians") {
				rot.angle <- rot.angle
			}else{
				stop("Unit for angles be either 'degrees' or 'radians'")
			}
		}
	}
	coords[, 2:3] <- round(rotation(as.matrix(coords[, 2:3]), rot.angle), digits = 8)
	if (rm.same.y == TRUE) {
		link <- remove.same.y(coords = coords, link = link)
	}
	
	low.coord <- which(coords[, 3] == min(coords[, 3]))

	if (plot.connexions) {
		xy.range<-apply(coords[,2:3],2,range)
		
		xy.range.min<-xy.range[1,2]-((xy.range[2,2]-xy.range[1,2])/5)
		xy.range.max<-xy.range[2,2]

		par(mar=c(1,1,1,1))
		plot(coords[, 2:3], pch = 20, asp = 1, cex = 0.5, axes = FALSE, 
			xlab = "", ylab = "", ylim=c(xy.range.min,xy.range.max))
		segments(x0 = coords[link[, 1], 2], y0 = coords[link[, 
			1], 3], x1 = coords[link[, 2], 2], y1 = coords[link[, 
			2], 3], col = "red")
		text(coords[,2:3],labels=as.character(coords[,1]),pos=2,col="red")

		site0<-c(mean(coords[,2]),xy.range.min)
		
		points(site0[1],site0[2],pch=19,col="blue")
		segments(x0=site0[1],y0=site0[2],x1=coords[low.coord,2],y1=coords[low.coord,3],col="blue")
		par(mar=c(5,4,4,2))
	}
	n.low.coord <- length(low.coord)
	link <- rbind(cbind(rep(0, n.low.coord), low.coord), as.matrix(link))
	nrow.link <- nrow(link)
	link <- cbind(1:nrow.link, link)
	points.order <- sort(coords[, 3], decreasing = FALSE, index.return = TRUE)$ix
	points.order <- c(0, points.order)
	link <- r.order.link(nrow.link, link, coords)
	link <- cbind(link, link[, 2])
	link.tmp <- link
	link2fac <- as.factor(link[, 2])
	matR <- mat.or.vec(1, n * nrow.link)
	mat <- .C("buildbinary", as.integer(nrow(link)), as.integer(link), 
		as.integer(points.order), as.integer(length(points.order)), 
		as.integer(n), matres = as.integer(matR), PACKAGE = "AEM")$matres
	res.mat <- matrix(mat, nrow = n, ncol = nrow.link, byrow = TRUE)
	res.link <- link[, 2:3]
	colnames(res.link) <- c("from", "to")
	
	res<-list(res.mat, res.link)
	names(res)<-c("se.mat", "edges")
	return(res)
}

