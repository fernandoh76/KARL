#' updateAll
#'
#' Updates all classification models without the need of re-installing the package.
#' @param mode Takes 'de' (default) for building the model with the whole data, or 'fs' for building the model after feature selection procedure.
#' @param replace Takes a logical (defualt is TRUE) indicating if the original model is replaced by the updated one or not and keeped in the workind directory.
#' @keywords update model KARL
#' @export
#' @examples
#' updateAll(mode='de',replace=T)

updateAll<-function(mode='de',replace=T){

	if (mode!='de' | mode!='fs'){
	
		stop('Parameter "mode" was wrongly set!')
	
	} else{

		print('Updating all models!')

		levels<-c('domain','phylum','class','order','family','genus')

		if (mode=='de'){

			print('--> Starting feature selection...')
			print('--> Building models with default data...')

			for (l in levels){

				print(paste('--> Building',l,'models...'))

				taxa<-as.vector(unique(alldata[,l]))

				if (replace==T){

					for (x in taxa){

						buildSVM(level=l,tax=x,prefix=paste(level,tax,sep='.'))

						oname<-paste(level,tax,'model.Rdata',sep='.')
						opath<-gsub(oname,'',system.file(level,oname,package='KARL'))
						system(paste('mv',oname,opath))

						print(paste('		-->',x))
					}
				
					print(paste('--> Updating',l,'successful!'))
				
				} else if (replace==F){

					system(paste('mkdir',l))

					for (x in taxa){

						buildSVM(level=l,tax=x,prefix=paste(level,tax,sep='.'))

						print(paste('		-->',x))
					}

					system(paste('mv *model.Rdata ./',l,sep=''))

					print(paste('--> Updating',l,'successful!'))
				}
			}
	
			print('--> Updating finished!')
		
		} else if (mode=='fs'){

			print('--> Proceed to feature selection...')
			print('--> Building models with default data...')

			for (l in levels){

				print(paste('--> Building',l,'models...'))

				taxa<-as.vector(unique(alldata[,l]))

				if (replace==T){

					for (x in taxa){

						new.droped<-dropUninformative(level=level,tax=tax)
						new.ranked<-rankSlope(idata=new.droped,level=level,tax=tax,plotting=F)
						
						buildSVM(data=new.ranked,tax=tax,level=level,prefix=paste(level,tax,sep='.'))

						oname<-paste(level,tax,'model.Rdata',sep='.')
						opath<-gsub(oname,'',system.file(level,oname,package='KARL'))
						system(paste('mv',oname,opath))

						print(paste('		-->',x))
					}
				
					print(paste('--> Updating',l,'successful!'))
				
				} else if (replace==F){

					system(paste('mkdir',l))

					for (x in taxa){

						new.droped<-dropUninformative(level=level,tax=tax)
						new.ranked<-rankSlope(idata=new.droped,level=level,tax=tax,plotting=F)
						
						buildSVM(data=new.ranked,tax=tax,level=level,prefix=paste(level,tax,sep='.'))

						print(paste('		-->',x))
					}

					system(paste('mv *model.Rdata ./',l,sep=''))

					print(paste('--> Updating',l,'successful!'))
				}
			}
	
			print('--> Updating finished!')
		}
	}
}
