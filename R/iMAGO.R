#' iMAGO
#'
#' Calculates MAGO index for re-classification of wrongly classified genomes and outputs membership probabilities.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param vclist Is the output of votingTaxa() function.
#' @keywords MAGO classification PCA reassignment KARL
#' @export
#' @examples
#' actino.imago<-iMAGO(vclist=actino.v7,rank='phylum')

iMAGO<-function(vclist,rank){

	tabrank<-sort(table(alldata[,rank]))
	namrank<-names(tabrank)

	ids<-names(vclist)

	probabilities<-list()

	for (v in 1:length(vclist)){
	
		vcdist<-vclist[[v]]
		maxdis<-max(vcdist)
		maxrnd<-round(maxdis)

		if (maxrnd<maxdis){

			maxrnd<-maxrnd+1
		}

		sdis<-seq(0,maxrnd)
		bands<-length(sdis)-1
		radius<-seq(0.5,bands,1)

		matdist<-matrix(nrow=length(tabrank),ncol=bands,0)
		rownames(matdist)<-namrank

		l1<-sdis[1]
		u1<-sdis[2]

		i1<-names(which(vcdist>l1 & vcdist<=u1))

		v1<-alldata[which(alldata$id%in%i1),rank]
		t1<-table(v1)
		n1<-names(t1)

		for (b in 1:bands){

			lower<-sdis[b]
			upper<-sdis[b+1]
			
			inband<-names(which(vcdist>lower & vcdist<=upper))

			voters<-alldata[which(alldata$id%in%inband),rank]
			tabvot<-table(voters)
			namvot<-names(tabvot)

			namlen<-length(which(namvot%in%n1))

			if (namlen==length(namvot)){
			
				for (w in 1:length(namvot)){

					xt<-which(namrank==namvot[w])
					yt<-b

					mago<-(tabvot[w]/tabrank[xt])/radius[b]
					matdist[xt,yt]<-mago
				}

				n1<-namvot

			} else {

				break
			}	
		}

		sumdist<-apply(matdist,1,sum)
		pdist<-NULL
		
		for (s in 1:length(sumdist)){

			pr<-sumdist[s]/sum(sumdist)
			pdist<-c(pdist,pr)
		}

		names(pdist)<-names(sumdist)

		probabilities[[v]]<-sort(pdist)

	}

	return(probabilities)
}
