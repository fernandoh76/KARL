#' enzymesVenn
#'
#' This function takes the output from getEnzymes() and draws a Venn diagram.
#' @param elist Is the resulting object from getEnzymes().
#' @keywords frequency enzymes Venn diagram KARL
#' @export
#' @examples
#' enzymesVenn(elist=enzymes)

enzymesVenn<-function(elist){

	require(VennDiagram,quietly=T)
	require(seqinr,quietly=T)

	counts<-unlist(lapply(elist,length))
	ncounts<-names(counts)

	nBoth<-counts[1]
	nNo<-counts[2]

	n3<-counts[3]
	n4<-counts[4]

	nams3<-c2s(s2c(ncounts[3])[-c(1,2)])
	nams4<-c2s(s2c(ncounts[4])[-c(1,2)])

	colors<-c('forestgreen','firebrick','grey')

	x11()

	draw.triple.venn(

		area1=n3+nBoth,
		area2=n4+nBoth,
		area3=nNo,

		n12=nBoth,
		n23=0,
		n13=0,
		n123=0,

		category=c(nams3,nams4,'Absent'),
		col=colors,
		cat.col=colors,
		lwd=5)
}
	
