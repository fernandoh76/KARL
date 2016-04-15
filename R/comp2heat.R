# Uses alldata.Rdata object
# requires package GPLOTS

#' comp2heat
#'
#' This function produces a heatmap using enzymes presence/absence comparing two selected taxa.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax1 The first taxon name.
#' @param tax2 The second taxon name.
#' @param reduce Takes logical (True or False) to reduce rows when presence/absence vectors are identical.
#' @keywords heatmap presence/absence enzymes KARL
#' @export
#' @examples
#' comp2heat(rank='genus',tax1='Actinobacillus',tax2='Bartonella',reduce=FALSE)

comp2heat<-function(rank,tax1,tax2,reduce=F){
	
	require(gplots)
	
	options(expressions=10000)
	
	x<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,x])

	if (reduce==F){

		s1<-alldata[which(l1==tax1),9:dim(alldata)[2]]
		s2<-alldata[which(l1==tax2),9:dim(alldata)[2]]
		he<-as.matrix(rbind(s1,s2))

		d1<-dim(s1)[1]
		d2<-dim(s2)[1]
		
		rowcol<-c(rep('forestgreen',d1),rep('firebrick',d2))
	
	} else if (reduce==T){
		
		s1<-unique(alldata[which(l1==tax1),9:dim(alldata)[2]])
		s2<-unique(alldata[which(l1==tax2),9:dim(alldata)[2]])
		he<-as.matrix(rbind(s1,s2))

		d1<-dim(s1)[1]
		d2<-dim(s2)[1]
		
		rowcol<-c(rep('forestgreen',d1),rep('firebrick',d2))		
	}
	
	hecol<-colorRampPalette(c('grey','black'))(2)
	
	heatmap.2(he,
		labRow=F,
		labCol=F,
		RowSideColors=rowcol,
		col=hecol,
		trace='none',
		key=F)
				
	legend(x=0,y=1,
		legend=c('Present','Absent'),
		fill=c('black','grey'),
		border=F,
		cex=.7
		)

	legend(x=0,y=0.90,
		legend=c(tax1,tax2),
		fill=c('forestgreen','firebrick'),
		border=F,
		cex=.7
		)
}
