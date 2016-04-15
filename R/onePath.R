#' onePath
#'
#' Retrieves the frequency of each enzyme to a certain pathway in a set of taxa.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param taxa A vector of names for desired taxa.
#' @keywords pid Identification code for desired pathway.
#' @export
#' @examples
#' tone<-c('Gammaproteobacteria','Mollicutes','Spirochaetia','Actinobacteria')
#' pone<-onePath(rank='class',taxa=tone,pid='00550')

onePath<-function(rank,taxa,pid){

	lev<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,lev])

	efreq<-NULL

	for (x in 1:length(taxa)){

		t1<-alldata[which(l1==taxa[x]),9:dim(alldata)[2]]
		a1<-apply(t1,2,sum)/dim(t1)[1]	
		
		efreq<-rbind(efreq,a1)

	}

	rownames(efreq)<-taxa
	colenz<-colnames(efreq)

	l<-which(names(li.pathways)==paste('ec',pid,sep=''))
	e1<-li.pathways[[l]]
	result<-efreq[,which(colenz%in%e1)]	

	rownames(result)<-taxa

	return(result)
}
