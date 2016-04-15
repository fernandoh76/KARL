#' cumsumPath
#'
#' This function plots a cumulative curve for pathway completeness.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param taxa A vector of taxa names.
#' @param pid The pathway id.
#' @param ecol A vector of colors of same length as taxa.
#' @keywords pathways enzymes KARL
#' @export
#' @examples
#' taxcum<-c('Actinobacteria','Tenericutes','Firmicutes','Proteobacteria','Bacteroidetes')
#' taxcol<-c('red','green','blue','violet','orange')
#' cumsumPath(rank='phylum',taxa=taxcum,pid='00550',ecol=taxcol)

cumsumPath<-function(rank,taxa,pid,ecol){

	lev<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,lev])
	
	einp<-li.pathways[[grep(pid,names(li.pathways))]]
	allcol<-colnames(alldata)

	etaxa<-NULL

	for (x in 1:length(taxa)){

		t1<-alldata[which(l1==taxa[x]),which(allcol%in%einp)]
		a1<-apply(t1,1,sum)/dim(t1)[2]	
		
		etaxa[[x]]<-a1

	}

	names(etaxa)<-taxa

	intervals<-seq(0,1,.05)

	cumintervals<-function(j){
		
		aux<-cumsum(rev(table(cut(j,breaks=intervals))))/length(j)
	
		return(aux)
	}

	ecum<-lapply(etaxa,cumintervals)

	for (e in 1:length(ecum)){

		if (e==1){

			plot(ecum[[e]],type='l',col=ecol[e],lwd=2,xaxt='n',
				xlab='Percentage of genomes',
				ylab='Fraction of enzymes present',
				main=paste('Pathway',pid)
				)

			axis(1,at=c(1,5,10,15,20),las=2,labels=c(0,25,50,75,100))

		} else if (e>1){
	
			lines(ecum[[e]],col=ecol[e],lwd=2)
		}

	}
	
	legend('topleft',legend=taxa,fill=ecol,border='white',cex=.7)

}
