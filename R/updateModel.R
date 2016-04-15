#' updateModel
#'
#' Updates a selected classification model without the need of re-installing the package.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name.
#' @param mode Takes 'de' (default) for building the model with the whole data, or 'fs' for building the model after feature selection procedure.
#' @param replace Takes a logical (defualt is TRUE) indicating if the original model is replaced by the updated one or not and keeped in the workind directory.
#' @keywords update model KARL
#' @export
#' @examples
#' updateModel(rank='genus',tax='Escherichia',replace=T)

updateModel<-function(rank,tax,mode='de',replace=T){

	if (mode!='de' | mode!='fs'){
	
		stop('Parameter "mode" was wrongly set!')
	
	} else{

		print(paste('Updating model for ',rank,' ',tax,'!',sep=''))

		if (mode=='de'){

			print('--> Proceed without feature selection...')
			print('--> Building model with default data...')

			bluildSVM(rank=rank,tax=tax,prefix=paste(rank,tax,sep='.'))
			
			if (replace==T){

				oname<-paste(prefix,'model.Rdata',sep='.')
				opath<-gsub(oname,'',system.file(rank,oname,package='KARL'))

				system(paste('mv',oname,opath))

				print('--> Updating successful!')
			
			} else if (replace==F){

				print('--> Updating successful!')
			}

		} else if (mode=='fs'){

			print('--> Starting feature selection...')
			print('--> Droping uninformative enzymes...')

			new.droped<-dropUninformative(rank=rank,tax=tax)

			print('--> Ranking most informative enzymes...')

			new.ranked<-rankSlope(idata=new.droped,rank=rank,tax=tax,plotting=F)

			print('--> Building model with selected data...')

			buildSVM(data=new.ranked,tax=tax,rank=rank,prefix=paste(rank,tax,sep='.'))

			if (replace==T){				

				oname<-paste(prefix,'model.Rdata',sep='.')
				opath<-gsub(oname,'',system.file(rank,oname,package='KARL'))

				system(paste('mv',oname,opath))

				print('--> Updating successful!')
			
			} else if (replace==F){

				print('--> Updating successful!')
			}
		}

		print('Updating finished!')
	}
}
