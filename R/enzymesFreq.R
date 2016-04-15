#' enzymesFreq
#'
#' This function calculates and plots the frequency of each enzyme in a pair of given taxa.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax1 The name of first taxon.
#' @param tax2 The name of second taxon.
#' @keywords frequency enzymes KARL
#' @export
#' @examples
#' enzymesFreq(rank='class',tax='Actinobacteria',tax2='Gammaproteobacteria')

enzymesFreq<-function(rank,tax1,tax2){

	x<-which(colnames(alldata)==rank)

	l1<-as.vector(alldata[,x])

	s1<-alldata[which(l1==tax1),9:dim(alldata)[2]]
	s2<-alldata[which(l1==tax2),9:dim(alldata)[2]]

	a1<-apply(s1,2,sum)/dim(s1)[1]
	a2<-apply(s2,2,sum)/dim(s2)[1]
	a3<-rbind(a1,a2)

	colors<-rep('grey',dim(a3)[2])

	for (x in 1:dim(a3)[2]){

		if (a3[1,x]>=.9 & a3[2,x]>=.9){

			colors[x]<-'orange'

		} else if (a3[1,x]<=.1 & a3[2,x]<=.1){

			colors[x]<-'black'

		} else if (a3[2,x]>=.9 & a3[1,x]<=.1){

			colors[x]<-'seagreen3'

		} else if (a3[2,x]<=.1 & a3[1,x]>=.9){

			colors[x]<-'orchid2'

		}
	}			

	plot(a1,a2,pch=19,main=rank,xlab=tax1,ylab=tax2,col=colors)

	abline(h=.1,col='red',lty=2)
	abline(v=.1,col='red',lty=2)
	abline(h=.9,col='red',lty=2)
	abline(v=.9,col='red',lty=2)
}
