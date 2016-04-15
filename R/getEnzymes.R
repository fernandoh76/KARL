#' getEnzymes
#'
#' Outputs a list with taxon-specific, shared and absent enzymes between taxa.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax1 The first taxon name.
#' @param tax2 The second taxon name.
#' @param dwtax1 Is the lower frequency cut-off for tax1 (default 0.1).
#' @param dwtax2 Is the lower frequency cut-off for tax2 (default 0.1).
#' @param uptax1 Is the lower frequency cut-off for tax1 (default 0.9).
#' @param uptax2 Is the lower frequency cut-off for tax2 (default 0.9).
#' @keywords specific compare enzymes KARL
#' @export
#' @examples
#' enzymes<-getEnzymes(rank='class',tax1='Actinobacteria',tax2='Gammaproteobacteria')

getEnzymes<-function(rank,tax1,tax2,uptax1=.9,dwtax1=.1,uptax2=.9,dwtax2=.1){

	x<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,x])

	s1<-alldata[which(l1==tax1),9:dim(alldata)[2]]
	s2<-alldata[which(l1==tax2),9:dim(alldata)[2]]

	a1<-apply(s1,2,sum)/dim(s1)[1]
	a2<-apply(s2,2,sum)/dim(s2)[1]
	a3<-rbind(a1,a2)

	intax1<-which(a3[1,]>=uptax1 & a3[2,]<=dwtax2)
	intax2<-which(a3[2,]>=uptax2 & a3[1,]<=dwtax1)

	inboth<-which(a3[1,]>=uptax1 & a3[2,]>=uptax2)
	noboth<-which(a3[1,]<=dwtax1 & a3[2,]<=dwtax2)

	listEnzymes<-list(inboth,noboth,intax1,intax2)

	n1<-paste('in',tax1,sep='')
	n2<-paste('in',tax2,sep='')

	names(listEnzymes)<-c('inBoth','noBoth',n1,n2)

	return(listEnzymes)
}
