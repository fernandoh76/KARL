#' findPath
#'
#' This function identifies the 5-number code for KEGG Pathways given a keyword.
#' @param keyword A keyword to link to KEGG.
#' @keywords KEGG pathway KARL
#' @export
#' @examples
#' findPath(keyword='peptidoglycan')

findPath<-function(keyword){

	cmd<-paste('wget -O tmppath http://rest.kegg.jp/find/pathway/',keyword,sep='')

	system(cmd,ignore.stderr=T)

	pinfo<-read.table('tmppath',sep='\t',header=F)
	colnames(pinfo)<-c('pathway','description')

	return(pinfo)

	system('rm -rf tmppath')
}
