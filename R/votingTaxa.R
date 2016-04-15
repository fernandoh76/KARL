#' votingTaxa
#'
#' This function helps to re-assign misclassified instances in default classification models by calculating distances to neighboring points in the PCA space.
#' @param ids Is a vector containing identifiers for misclassified genomes.
#' @param nproc Is the number of processors to perform calculations.
#' @param plim Is the maximum number of Principal Components taken for distance calculation.
#' @keywords PCA enzymes classification KARL
#' @export
#' @examples
#' load('evalActinobacteria.df.i3.Rdata')
#' actinoTprF<-evalActinobacteria.df.i3@perfectErrors$TprF
#' actino.vt7<-votingTaxa(ids=actinoTprF,nproc=10,plim=7)

votingTaxa<-function(ids,nproc,plim=3){

	require(doParallel,quietly=T)

	#############################
	# Define internal functions #
	#############################

	chunker<-function(m,n){
		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}

	getdist<-function(y){

		q<-query
		s<-subjects[schunks[[y]],]

		r<-apply(s,1,function(j){
			
			rb<-rbind(q,j)
			dt<-sqrt(sum((rb[1,]-rb[2,])^2))
		})

		return(r)
	}

	#################################
	# Parallel distance calculation #
	#################################

	distdata<-prcompEC$x[,1:plim]
	rnames<-as.vector(alldata$id)
	rownames(distdata)<-rnames

	errids<-which(rnames%in%ids)
	
	result<-list()

	for (e in 1:length(errids)){
		
		query<-distdata[errids[e],]
		subjects<-distdata[-errids[e],]

		sdim<-dim(subjects)
		sseq<-seq(1:sdim[1])

		schunks<-chunker(sseq,nproc)

		registerDoParallel(cl=nproc,cores=nproc)

		distances<-foreach(y=1:nproc) %dopar% {
			
			getdist(y)	
		}

		finaldist<-unlist(distances)

		result[[e]]<-finaldist
	}
	
	names(result)<-rnames[errids]

	return(result)
}
