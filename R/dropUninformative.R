#' dropUninformative
#'
#' Feature selection by removing uninformative enzymes based on correlation and co-occurrence.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name.
#' @param dwtax The lower frequency cut-off.
#' @param uptax The upper frequency cut-off.
#' @param cf The correlation cut-off.
#' @keywords feature selection enzymes correlation frequency KARL
#' @export
#' @examples
#' streptococcus.droped<-dropUninformative(rank='genus',tax='Streptococcus',dwtax=0.1,uptax=0.9,cf=0.9)

dropUninformative<-function(rank,tax,dwtax=.1,uptax=.9,cf=.9){
	
	require(caret,quietly=T)

	x<-which(colnames(alldata)==rank)

	l1<-as.vector(alldata[,x])

	s1<-alldata[which(l1==tax),9:dim(alldata)[2]]
	s2<-alldata[which(l1!=tax),9:dim(alldata)[2]]

	a1<-apply(s1,2,sum)/dim(s1)[1]
	a2<-apply(s2,2,sum)/dim(s2)[1]
	a3<-rbind(a1,a2)

	inboth<-which(a3[1,]>=uptax & a3[2,]>=uptax)
	noboth<-which(a3[1,]<=dwtax & a3[2,]<=dwtax)

	outside<-names(c(inboth,noboth))

	filtdata<-alldata[,-which(colnames(alldata)%in%outside)]

	cordata<-filtdata[,9:dim(filtdata)[2]]
	corMat<-cor(cordata)
	
	highCor<-findCorrelation(corMat,cutoff=cf)
	
	finaldata<-cbind(filtdata[,1:8],cordata[,-highCor])

	return(finaldata)

}
