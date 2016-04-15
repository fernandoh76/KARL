#' consPath
#'
#' This function calculates the completeness of each metabolic pathway for a certain taxon.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param taxa The taxon name, or a vector of names.
#' @keywords pathways enzymes KARL
#' @export
#' @examples
#' taxco<-c('Streptomyces','Acinetobacter','Mycoplasma','Escherichia')
#' consPath(rank='genus',taxa=taxco)

consPath<-function(rank,taxa,pid){

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

	result<-NULL

	for (l in 1:length(li.pathways)){

		e1<-li.pathways[[l]]
		w1<-efreq[,which(colenz%in%e1)]	

		if (length(dim(w1)==2)){	
		
			a1<-apply(w1,1,sum)/dim(w1)[2]
			result<-cbind(result,a1)

		} else {
	
			result<-cbind(result,w1)
		}

	}

	rownames(result)<-taxa
	colnames(result)<-gsub('ec','path:map',names(li.pathways))

	return(result)
}
