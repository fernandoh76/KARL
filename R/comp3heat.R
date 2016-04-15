# Uses alldata.Rdata object
# Produces a heatmap considering enzymes P/A for pairs of taxa at the same rank
# RANK takes 'phylum', 'class', 'order' or 'genus'
# TAX1 is one selected taxon and TAX2 the second selected taxon
# requires package GPLOTS

#' comp3heat
#'
#' This function produces a heatmap using enzymes presence/absence comparing two selected taxa and one selected outgroup.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax1 The first taxon name.
#' @param tax2 The second taxon name.
#' @param outgroup The third taxon as outgroup.
#' @param reduce Takes logical (True or False) to reduce rows when presence/absence vectors are identical.
#' @keywords heatmap outgroup presence/absence enzymes KARL
#' @export
#' @examples
#' comp3heat(rank='genus',tax1='Xylella',tax2='Dyella',outgroup='Nevskia',reduce=FALSE)

comp3heat<-function(rank,tax1,tax2,outgroup,reduce=F){
	
	require(gplots)
	
	options(expressions=10000)
	
	x<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,x])

	if (reduce==F){

		s1<-alldata[which(l1==tax1),9:dim(alldata)[2]]
		s2<-alldata[which(l1==tax2),9:dim(alldata)[2]]
		s3<-alldata[which(l1==outgroup),9:dim(alldata)[2]]
		he<-as.matrix(rbind(s1,s2,s3))

		d1<-dim(s1)[1]
		d2<-dim(s2)[1]
		d3<-dim(s3)[1]
		
		rowcol<-c(rep('forestgreen',d1),rep('firebrick',d2),rep('blue',d3))
	
	} else if (reduce==T){
		
		s1<-unique(alldata[which(l1==tax1),9:dim(alldata)[2]])
		s2<-unique(alldata[which(l1==tax2),9:dim(alldata)[2]])
		s3<-alldata[which(l1==outgroup),9:dim(alldata)[2]]
		he<-as.matrix(rbind(s1,s2,s3))

		d1<-dim(s1)[1]
		d2<-dim(s2)[1]
		d3<-dim(s3)[1]
				
		rowcol<-c(rep('forestgreen',d1),rep('firebrick',d2),rep('blue',d3))		
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
		legend=c(tax1,tax2,outgroup),
		fill=c('forestgreen','firebrick','blue'),
		border=F,
		cex=.7
		)
}
