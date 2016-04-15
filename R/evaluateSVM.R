#' evaluateSVM
#'
#' This function performs a m-repeated n-fold cross-validation and stores detailed results for each genome in a class called evalModel.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name.
#' @param folds Is the number of folds in cross-validation.
#' @param iters Is the number of iterations in the repeated cross-validation.
#' @param parallel Is a logical indication if doing parallelization or sequential computing. Number of cores will be equal to number of folds. calculation.
#' @keywords evaluate model SVM cross-validation KARL
#' @export
#' @examples
#' evalCyanobacteria.df.i3<-evaluateSVM(rank='phylum',tax='Cyanobacteria',folds=10,parallel=T,iters=3)

evaluateSVM<-function(idata=alldata,rank,tax,folds=10,iters=3,parallel=T){
	
	require(doParallel,quietly=T)

	#############################
	# Define internal functions #
	#############################

	# Generates data chunks for cross-validation

	chunker<-function(m,n){
		s<-split(m,cut(seq_along(m),n,labels=F))
		return(s)
	}

	# Merges a list of dataframes to parse results

	merge.dflist<-function(x,y){
		merge(x,y,all=TRUE,by='id')
	}

	# Implements n-fold cross-validation

	crossval<-function(y){

		options(java.parameters='-Xmx10g')
		require(RWeka,quietly=T)

		allchunks<-mapply(c,truechunks,falsechunks,SIMPLIFY=F)
			
		testy<-allchunks[[y]]
		trainy<-unlist(allchunks[-y])

		test<-randomdata[which(randomdata$id%in%testy),]
		train<-randomdata[which(randomdata$id%in%trainy),]

		print(paste('Training model with training set',y))
	
		train.model<-SMO(category~.,
		             data=train[,9:dim(train)[2]],
			     control=Weka_control(K=list('PolyKernel',E=1))
			     )

		print(paste('Evaluating test set',y))
		evaluation<-predict(train.model,test[,(9:dim(test)[2]-1)])
		valids<-as.data.frame(cbind(as.vector(test$id),evaluation!=1))
		colnames(valids)<-c('id','prediction')
			
		return(valids)
	}

	######################
	# Define new classes #
	######################

	# Define class 'evalModel' to store results

	setClass('evalModel',slots=c(
		 perfectMatches='list',
		 variousMatches='list',
		 perfectErrors='list',
	         variousErrors='list',
	         confusionMatrix='matrix'),
		 where=.GlobalEnv
		)

	##################
	# Start analysis #
	##################
	
	lev<-which(colnames(idata)==rank)
	ntax<-length(which(idata[,lev]==tax))

	if (ntax<=1){
		notenough<-paste('Less than 2 organisms belonging to',tax,
				 ': unsufficient information to proceed!'
				)

		stop(notenough)
	
	} else {
		
		whichna<-which(is.na(idata[,lev])==T)

		if (length(whichna)>0){
			
			evaldata<-idata[-whichna,]

		} else if (length(whichna)==0){

			evaldata<-idata

		}

		evaldata$category[evaldata[,lev]==tax]<-T
		evaldata$category[evaldata[,lev]!=tax]<-F

		truedata<-evaldata[which(evaldata$category==T),]
		falsedata<-evaldata[which(evaldata$category==F),]

		dimtrue<-dim(truedata)
		dimfalse<-dim(falsedata)

		if (folds>=ntax){
			
			foldms<-paste('Not enough ',tax,' members (',dimtrue[1],
				       ') to perform ',folds,'-fold cross-validation!',sep='')
			print(foldms)

			folds<-dimtrue[1]

			start<-paste(iters,'-repeated, ',folds,
                             '-fold cross-validation started!',
                             sep='')
		} else {

			start<-paste(iters,'-repeated, ',folds,
                             '-fold cross-validation started!',
                             sep='')
		}

		print(start)

		# Starts repeated cross-validation

		t1<-round(as.vector(proc.time()[3])/60,digits=1)

		result<-NULL

		for (x in 1:iters){

			print(paste('Repetition',x,'started!'))

			# Chunks false data

			randomfalse<-falsedata[sample(dimfalse[1]),]

			catfalse<-cbind(
      			         as.vector(randomfalse$id),as.vector(randomfalse$category)
				)

			colnames(catfalse)<-c('id','category')
			catfalse<-as.data.frame(catfalse)
	
			fseq<-seq(1,dimfalse[1])
			fchunks<-chunker(fseq,folds)

			falsechunks<-lapply(fchunks,function(x){as.vector(falsedata[x,'id'])})

			# Chunks true data
			
			randomtrue<-truedata[sample(dimtrue[1]),]
			
			cattrue<-cbind(
      			         as.vector(randomtrue$id),as.vector(randomtrue$category)
				)

			colnames(cattrue)<-c('id','category')
			cattrue<-as.data.frame(cattrue)
	
			tseq<-seq(1,dimtrue[1])
			tchunks<-chunker(tseq,folds)

			truechunks<-lapply(tchunks,function(x){as.vector(truedata[x,'id'])})
		
			# Starts parallel computing for cross-validation

			randomdata<-rbind(randomfalse,randomtrue)
			randomcat<-rbind(catfalse,cattrue)

			if (parallel==T){

				registerDoParallel(makeCluster(folds))

				print(paste('Starting parallel computing for training',rank,tax))

				partial<-foreach(y=1:folds) %dopar% {
			
					crossval(y)
	
				}
			
			} else if (parallel==F){

				registerDoParallel(makeCluster(folds))

				print(paste('Starting sequential computing for training',rank,tax))

				partial<-foreach(y=1:folds) %do% {
			
					crossval(y)
	
				}
			}
			
			stopCluster(folds)

			print('Training and evaluation finished')

			# Store results

			result[[x]]<-do.call(rbind,partial)
			
			print(paste('Repetition',x,'finished!'))

		}

		t2<-round(as.vector(proc.time()[3])/60,digits=1)
		tf<-t2-t1

		print(paste('Model evaluation for',rank,tax,'elapsed',tf,'minutes'))
		print('Parsing results')

		# Parse results

		output<-Reduce(merge.dflist,result)
		output<-merge(output,randomcat,by='id')
		
		ites<-paste('ite',seq(1,iters),sep='')
		colnames(output)<-c('id',ites,'category')

		dio<-dim(output)[2]
			
		tabeval<-apply(output[,2:(dio-1)],1,table)		
		tableng<-unlist(lapply(tabeval,length))
	
		perfect<-output[which(tableng==1),]
		various<-output[which(tableng!=1),]

		pdim<-dim(perfect)
		vdim<-dim(various)
	
		if (pdim[1]>0 & vdim[1]>0) {

			# Perfect

			pFT<-which(perfect$category==F & perfect$ite1==T)
			pTF<-which(perfect$category==T & perfect$ite1==F)
			pFF<-which(perfect$category==F & perfect$ite1==F)
			pTT<-which(perfect$category==T & perfect$ite1==T)

			idpFT<-as.vector(perfect[pFT,'id'])
			idpTF<-as.vector(perfect[pTF,'id'])
			idpFF<-as.vector(perfect[pFF,'id'])
			idpTT<-as.vector(perfect[pTT,'id'])
	
			perfectErrors<-list(idpFT,idpTF)
			names(perfectErrors)<-c('FprT','TprF')

			perfectMatches<-list(idpFF,idpTT)
			names(perfectMatches)<-c('FprF','TprT')
	
			# Various

			tabvar<-apply(various[,2:(vdim[2]-1)],1,table)
			rowvar<-which(rownames(tabvar)=='TRUE')
	
			coff<-iters/2

			variousDecision<-rep('FALSE',vdim[1])
			variousDecision[which(tabvar[rowvar,]>coff)]<-'TRUE'

			dfDecision<-cbind(various[,c('id','category')],variousDecision)
	
			vTF<-which(dfDecision$variousDecision==F & dfDecision$category==T)
			vFT<-which(dfDecision$variousDecision==T & dfDecision$category==F)
			vTT<-which(dfDecision$variousDecision==T & dfDecision$category==T)
			vFF<-which(dfDecision$variousDecision==F & dfDecision$category==F)

			idvTF<-as.vector(various[vTF,'id'])
			idvFT<-as.vector(various[vFT,'id'])
			idvTT<-as.vector(various[vTT,'id'])
			idvFF<-as.vector(various[vFF,'id'])

			variousErrors<-list(idvTF,idvFT)
			names(variousErrors)<-c('TprF','FprT')

			variousMatches<-list(idvFF,idvTT)
			names(variousMatches)<-c('FprF','TprT')

			# Confusion matrix
	
			pmSum<-lapply(perfectMatches,length)
			vmSum<-lapply(variousMatches,length)
			peSum<-lapply(perfectErrors,length)
			veSum<-lapply(variousErrors,length)

			cmFF<-sum(pmSum$FprF,vmSum$FprF)
			cmTT<-sum(pmSum$TprT,vmSum$TprT)
			cmFT<-sum(peSum$FprT,veSum$FprT)
			cmTF<-sum(peSum$TprF,veSum$TprF)

			confusionMatrix<-rbind(c(cmFF,cmFT),c(cmTF,cmTT))
			colnames(confusionMatrix)<-c('prF','prT')
			rownames(confusionMatrix)<-c('F','T')

			# Generates object with final results

			resultEval<-new('evalModel',
					perfectMatches=perfectMatches,
					variousMatches=variousMatches,
				        perfectErrors=perfectErrors,
					variousErrors=variousErrors,
				        confusionMatrix=confusionMatrix
					)
	
			return(resultEval)

		} else if (pdim[1]>0 & vdim[1]==0) {
		
			# Perfect

			pFT<-which(perfect$category==F & perfect$ite1==T)
			pTF<-which(perfect$category==T & perfect$ite1==F)
			pFF<-which(perfect$category==F & perfect$ite1==F)
			pTT<-which(perfect$category==T & perfect$ite1==T)

			idpFT<-as.vector(perfect[pFT,'id'])
			idpTF<-as.vector(perfect[pTF,'id'])
			idpFF<-as.vector(perfect[pFF,'id'])
			idpTT<-as.vector(perfect[pTT,'id'])
	
			perfectErrors<-list(idpFT,idpTF)
			names(perfectErrors)<-c('FprT','TprF')

			perfectMatches<-list(idpFF,idpTT)
			names(perfectMatches)<-c('FprF','TprT')
	
			# Various
	
			variousErrors<-list(NULL,NULL)
			names(variousErrors)<-c('TprF','FprT')

			variousMatches<-list(NULL,NULL)
			names(variousMatches)<-c('FprF','TprT')
	
			# Confusion matrix
	
			pmSum<-lapply(perfectMatches,length)
			vmSum<-lapply(variousMatches,length)
			peSum<-lapply(perfectErrors,length)
			veSum<-lapply(variousErrors,length)

			cmFF<-sum(pmSum$FprF,vmSum$FprF)
			cmTT<-sum(pmSum$TprT,vmSum$TprT)
			cmFT<-sum(peSum$FprT,veSum$FprT)
			cmTF<-sum(peSum$TprF,veSum$TprF)

			confusionMatrix<-rbind(c(cmFF,cmFT),c(cmTF,cmTT))
			colnames(confusionMatrix)<-c('prF','prT')
			rownames(confusionMatrix)<-c('F','T')

			# Generates object with final results

			finalresult<-new('evalModel',perfectMatches=perfectMatches,
				         perfectErrors=perfectErrors,
		                         variousMatches=variousMatches,
		                         variousErrors=variousErrors,
		                         confusionMatrix=confusionMatrix)
		
			return(finalresult)	

		} else if (pdim[1]==0 & vdim[1]>0) {
	
			# Perfect

			perfectErrors<-list(NULL,NULL)
			names(perfectErrors)<-c('FprT','TprF')
	
			perfectMatches<-list(NULL,NULL)
			names(perfectMatches)<-c('FprF','TprT')
			
			# Various
	
			tabvar<-apply(various[,2:(vdim[2]-1)],1,table)
			rowvar<-which(rownames(tabvar)=='TRUE')
	
			coff<-iters/2
	
			variousDecision<-rep('FALSE',vdim[1])
			variousDecision[which(tabvar[rowvar,]>coff)]<-'TRUE'
	
			dfDecision<-cbind(various[,c('id','category')],variousDecision)
		
			vTF<-which(dfDecision$variousDecision==F & dfDecision$category==T)
			vFT<-which(dfDecision$variousDecision==T & dfDecision$category==F)
			vTT<-which(dfDecision$variousDecision==T & dfDecision$category==T)
			vFF<-which(dfDecision$variousDecision==F & dfDecision$category==F)
	
			idvTF<-as.vector(various[vTF,'id'])
			idvFT<-as.vector(various[vFT,'id'])
			idvTT<-as.vector(various[vTT,'id'])
			idvFF<-as.vector(various[vFF,'id'])
	
			variousErrors<-list(idvTF,idvFT)
			names(variousErrors)<-c('TprF','FprT')
	
			variousMatches<-list(idvFF,idvTT)
			names(variousMatches)<-c('FprF','TprT')
	
			# Confusion matrix
		
			pmSum<-lapply(perfectMatches,length)
			vmSum<-lapply(variousMatches,length)
			peSum<-lapply(perfectErrors,length)
			veSum<-lapply(variousErrors,length)
	
			cmFF<-sum(pmSum$FprF,vmSum$FprF)
			cmTT<-sum(pmSum$TprT,vmSum$TprT)
			cmFT<-sum(peSum$FprT,veSum$FprT)
			cmTF<-sum(peSum$TprF,veSum$TprF)
	
			confusionMatrix<-rbind(c(cmFF,cmFT),c(cmTF,cmTT))
			colnames(confusionMatrix)<-c('prF','prT')
			rownames(confusionMatrix)<-c('F','T')
	
			# Generates object with final results
	
			finalresult<-new('evalModel',perfectMatches=perfectMatches,
				         perfectErrors=perfectErrors,
		                         variousMatches=variousMatches,
		                         variousErrors=variousErrors,
		                         confusionMatrix=confusionMatrix)
		
			return(finalresult)
	
		} else {
			stop('No input to calculate evaluation: model error!')
		}
	}

	#####################
	# Finishes analysis #
	#####################		
}
