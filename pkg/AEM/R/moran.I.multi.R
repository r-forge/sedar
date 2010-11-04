`moran.I.multi` <-
function(eigenvector.mat,link,weight,scaled=FALSE,normalize=FALSE,na.rm = FALSE,test.type="permutation",nperm=999,plot.res=TRUE){
	
	if(!is.matrix(link)){
		link<-as.matrix(link)
	}
	if(mode(link)!="numeric"){
		stop("'link' should be numeric")
	}
	
	if(missing(weight)){
		weight<-rep(1,nrow(link))
	}
	
	#CC# Creat the distance matrix for the first class of distance (first neighbour)
	nsite<-nrow(eigenvector.mat)
	mat.W<-matrix(0,ncol=nsite,nrow=nsite)
	
	for(i in 1:nrow(link)){
		mat.W[link[i,1],link[i,2]]<-weight[i]
	}
	
	#CC# Build the result object
	n.eigen<-ncol(eigenvector.mat)
	res.mat<-matrix(nrow=n.eigen,ncol=2)
	
	#CC# Calculate a Moran's for each eigenvectors
	#CC# Extract the p-value and the observed Moran's I
	for(i in 1:n.eigen){
		#CC# Define alternative hypotheses
		moran.val<-moran.I.basic(eigenvector.mat[,i],mat.W,scaled=scaled)
		alter <- ifelse(moran.val$observed > moran.val$expected, "greater", "less")

		res<-moran.I.uni(eigenvector.mat[,i],mat.W,scaled=scaled,normalize=normalize,na.rm = na.rm,test.type=test.type,nperm=nperm,alternative=alter)
		res.mat[i,2]<-res$p.value
		res.mat[i,1]<-res$observed
	}
	
	colnames(res.mat)<-c("observed","p.value")
	
	#CC# Calculate the expected of the Moran's I under the null hypothesis
	res.expe<-res$expected
	
	#CC# draw a plot of the results
	if(plot.res){
		plot(res.mat[,1],xlab="Eigenvector",ylab="Moran's I",type="n",las=1)
		abline(h=res.expe,col="red",lwd=2)
		points(res.mat[,1])
		points(which(res.mat[,2]<=0.05),res.mat[which(res.mat[,2]<=0.05),1],pch=19)
		legend("topright",c("Not significant Moran's I","Significant Moran's I",expression(paste("Expected value under ", H[0]))),col=c("black","black","red"),lty=c(0,0,1),lwd=c(0,0,2),pch=c(1,19,-1))
	}
	
	#CC# Return results
	result<-list(res.mat,res.expe)
	names(result)<-c("res.mat","expected")
	return(result)
}

