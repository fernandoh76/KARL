# Produces a heatmap considering enzymes P/A in one taxon
# RANK takes 'phylum', 'class', 'order', 'family' or 'genus'
# TAX is the selected taxon name 
# LROWS takes T or F to decide plotting rownames or not (genome names)
# REDUCE takes T or F and reduce heatmap rows when presence/absence vectors are identical.

#' comp1heat
#'
#' This function produces a heatmap using enzymes presence/absence in a single taxon.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The selected taxon name.
#' @param lrows Takes logical (True or False) to decide plotting rownames (genome names) or not.
#' @param reduce Takes logical (True or False) to reduce rows when presence/absence vectors are identical.
#' @keywords heatmap presence/absence enzymes KARL
#' @export
#' @examples
#' comp1heat(rank='genus',tax='Salmonella',lrows=FALSE,reduce=TRUE)
#' comp1heat(rank='genus',tax='Actinobacillus',lrows=T,reduce=F)

comp1heat<-function(rank,tax,lrows=F,reduce=F){
	
	require(gplots)
	options(expressions=10000)

	maintext<-paste(rank,'::',tax,sep='')

	x<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,x])
	he<-as.matrix(alldata[which(l1==tax),9:dim(alldata)[2]])
	
	hecol<-colorRampPalette(c('grey','black'))(2)
	
	if (reduce==F){
	
		if (lrows==T){	
			
			ini<-paste(strsplit(tax,'')[[1]][1],'.',sep='')
			rs<-alldata[which(l1==tax),'genome']
			rows<-gsub(tax,ini,rs)

			heatmap.2(he,
				main=maintext,
				labRow=rows,
				cexRow=.4,
				labCol=F,
				col=hecol,
				trace='none',
				margins=c(10,10),
				key=F)
			
		} else if (lrows==F){
						
			rows<-F

			heatmap.2(he,
				main=maintext,
				labRow=rows,
				labCol=F,
				col=hecol,
				trace='none',
				key=F)
		}
		
	} else if (reduce==T){
		
		he<-unique(he)
	
		if (lrows==T){	
			
			ini<-paste(strsplit(tax,'')[[1]][1],'.',sep='')
			rs<-alldata[which(l1==tax),'genome']
			rows<-gsub(tax,ini,rs)

			heatmap.2(he,
				main=maintext,
				labRow=rows,
				cexRow=.4,
				labCol=F,
				col=hecol,
				trace='none',
				margins=c(10,10),
				key=F)
			
		} else if (lrows==F){
						
			rows<-F

			heatmap.2(he,
				main=maintext,
				labRow=rows,
				labCol=F,
				col=hecol,
				trace='none',
				key=F)
		}
	
	}	
	
	legend('topleft',
		legend=c('Present','Absent'),
		fill=c('black','grey'),
		border=F,
		cex=.7
		)	
}
