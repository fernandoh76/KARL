# DATA is alldata by default. Can be a subset of selected enzymes.
# rank takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
# TAX is the selected taxon to classify at the desired taxonomic rank.
# PREFIX is for the name of model object.

#' buildSVM
#'
#' This function allows to build the classification model for a given taxon using all enzymes or a subset.
#' @param data The alldata object by default, can be a subset of selected enzymes (tipically by feature selection).
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax Is the selected taxon to classify at the desired taxonomic rank.
#' @keywords SVM classification SMO RWeka KARL
#' @export
#' @examples
#' buildSVM(rank='phylum',tax='Cyanobacteria',prefix='Cyanobacteria')

buildSVM<-function(data=alldata,rank,tax,prefix){

	options(java.parameters='-Xmx10g')

	library('RWeka')
	library('rJava')

	svmdata<-data
	lev<-which(colnames(svmdata)==rank)
	
	nas<-which(is.na(svmdata[,lev])==T)

	if (length(nas)>0){
		
		svmdata<-svmdata[-nas,]
	}

	svmdata$category[svmdata[,lev]==tax]<-T
	svmdata$category[svmdata[,lev]!=tax]<-F
	
	model.svm<-SMO(
		category~.,
		data=svmdata[,9:dim(svmdata)[2]],
		control=Weka_control(K=list('PolyKernel',E=1))
		)

	.jcache(model.svm$classifier)

	nmodel<-paste(prefix,'model',sep='.')

	assign(nmodel,model.svm)
	
	save(list=nmodel,file=paste(nmodel,'.Rdata',sep=''))

	print(paste(tax,'model done!'))
	print('#-#-#-#-#')

	return(model.svm)
}
