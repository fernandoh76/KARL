#' taxChild
#'
#' Outputs the children taxa immediately below a certain rank and taxon.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name.
#' @param counting Is a logical (T or F) for counting the number of children or not.
#' @keywords taxa children KARL
#' @export
#' @examples
#' taxChild(rank='class',tax='Actinobacteria',counting=FALSE)

taxChild<-function(rank,tax,counting=F){
	
	x<-which(colnames(alldata)==rank)
	a1<-which(as.vector(alldata[,x])==tax)
	
	if (counting==F){
	
		children<-unique(as.vector(alldata[a1,x+1]))
		
		return(children)
		
	} else if (counting==T){
		
		children<-sort(table(as.vector(alldata[a1,x+1])))
		
		return(children)
	}
}
