#' predictOne
#'
#' Predicts the membership of a genome or set to a particular taxon.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name.
#' @param input A presence/absence vector or data frame.
#' @keywords predict taxon enzymes presence/absence KARL.
#' @export
#' @examples
#' input.pred<-predictOne(rank='phylum',tax='Proteobacteria',input=inputdata)

predictOne<-function(rank,tax,input){

	model.file<-paste(rank,tax,'model','Rdata',sep='.')
	model.path<-paste(system.file('models',rank,package='KARL'),model.file,sep='/')

	load(model.path)

	model.pred<-predict(get(gsub('.Rdata','',model.file)),input)
	
	model.pred[which(model.pred==0)]<-'Yes'
	model.pred[which(model.pred==1)]<-'No'

	model.pred<-as.data.frame(model.pred)
	colnames(model.pred)<-tax
	rownames(model.pred)<-rownames(input)

	return(model.pred)
}
