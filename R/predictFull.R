#' predictFull
#'
#' Predicts the whole lineage of a genome or set from domain to genus.
#' @param input A presence/absence vector or data frame.
#' @keywords predict lineage enzymes presence/absence KARL.
#' @export
#' @examples
#' full.pred<-predictFull(input=inputdata)

predictFull<-function(input){

	##### Define internal function ####

	parsePrediction<-function(x){

		colna<-colnames(x)
		rowna<-rownames(x)

		parsed<-list()	
	
		positives<-NA
		negatives<-NA
		ambiguous<-NA

		for (n in 1:dim(x)[1]){
	
			w<-which(x[n,]=='Yes')
			le<-length(w)

			if (le==1){
			
				positives<-c(positives,paste(rowna[n],colna[w],sep='-'))
		
			} else if (le==0){

				negatives<-c(negatives,rowna[n])
		
			} else if (le>1){

				ambiguous<-c(ambiguous,paste(rowna[n],paste(colna[w],collapse='/'),sep='-'))
			}
		}

		parsed[[1]]<-positives
		parsed[[2]]<-negatives
		parsed[[3]]<-ambiguous
		
		names(parsed)<-c('positives','negatives','ambiguous')

		return(parsed)
	}		

	##### Output object #####

	din<-dim(input)[1]

	rown<-rownames(input)

	full.pred<-matrix(nrow=din,ncol=6,'Unknown')

	rownames(full.pred)<-rown

	colnames(full.pred)<-c('domain','phylum','class','order','family','genus')

	##### Genus prediction #####	

	rank<-'genus'

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL

	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input)))
	}

	result.rank<-as.data.frame(result.rank)

	parsed.result<-parsePrediction(result.rank)

	gposi<-parsed.result$positives
	gposi<-gposi[-is.na(gposi)]

	if (length(gposi)>0){

		for (p in gposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			gen<-spl[[1]][2]
		
			full.pred[pos,1:6]<-as.vector(as.matrix(alldata[which(alldata$genus==gen)[1],3:8]))
		}
	}

	gambi<-parsed.result$ambiguous
	gambi<-gambi[-is.na(gambi)]

	if (length(gambi)>0){

		for (p in gambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			gen<-spl[[1]][2]
		
			full.pred[pos,6]<-gen
		}

		
	}

	gnext<-unique(c(parsed.result$negatives,parsed.result$ambiguous))
	gnext<-gnext[-is.na(gnext)]
	gnext<-rown[which(rown%in%gnext)]

	#### Family prediction ####

	rank<-'family'

	input2<-input[gnext,]

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL

	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input2)))
	}

	parsed.result<-parsePrediction(result.rank)
		
	fposi<-parsed.result$positives
	fposi<-fposi[-is.na(fposi)]

	if (length(fposi)>0){

		for (p in fposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			fam<-spl[[1]][2]
		
			full.pred[pos,1:5]<-as.vector(as.matrix(alldata[which(alldata$family==fam)[1],3:7]))
		}
	}

	fambi<-parsed.result$ambiguous
	fambi<-fambi[-is.na(fambi)]

	if (length(fambi)>0){

		for (p in fambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			fam<-spl[[1]][2]
		
			full.pred[pos,5]<-fam
		}
	}

	fnext<-unique(c(parsed.result$negatives,parsed.result$ambiguous))
	fnext<-fnext[-is.na(fnext)]
	fnext<-rown[which(rown%in%fnext)]
	
	#### Order prediction ####
	
	rank<-'order'

	input3<-input[fnext,]

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL
	
	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input3)))
	}

	parsed.result<-parsePrediction(result.rank)
		
	oposi<-parsed.result$positives
	oposi<-oposi[-is.na(oposi)]

	if (length(oposi)>0){

		for (p in oposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			ord<-spl[[1]][2]
		
			full.pred[pos,1:4]<-as.vector(as.matrix(alldata[which(alldata$order==ord)[1],3:6]))
		}
	}

	oambi<-parsed.result$ambiguous
	oambi<-oambi[-is.na(oambi)]

	if (length(oambi)>0){

		for (p in oambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			ord<-spl[[1]][2]
		
			full.pred[pos,4]<-ord
		}
	}

	onext<-unique(c(parsed.result$negatives,parsed.result$ambiguous))
	onext<-onext[-is.na(onext)]
	onext<-rown[which(rown%in%onext)]

	#### Class prediction ####

	rank<-'class'

	input4<-input[onext,]

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL
	
	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input4)))
	}

	parsed.result<-parsePrediction(result.rank)
		
	cposi<-parsed.result$positives
	cposi<-cposi[-is.na(cposi)]

	if (length(cposi)>0){

		for (p in cposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			cla<-spl[[1]][2]
		
			full.pred[pos,1:3]<-as.vector(as.matrix(alldata[which(alldata$class==cla)[1],3:5]))
		}
	}

	cambi<-parsed.result$ambiguous
	cambi<-cambi[-is.na(cambi)]

	if (length(cambi)>0){

		for (p in cambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			cla<-spl[[1]][2]
		
			full.pred[pos,3]<-cla
		}
	}

	cnext<-unique(c(parsed.result$negatives,parsed.result$ambiguous))
	cnext<-cnext[-is.na(cnext)]
	cnext<-rown[which(rown%in%cnext)]
	
	#### Phylum prediction ####

	rank<-'phylum'

	input5<-input[cnext,]

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL
	
	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input5)))
	}

	parsed.result<-parsePrediction(result.rank)
		
	pposi<-parsed.result$positives
	pposi<-pposi[-is.na(pposi)]

	if (length(pposi)>0){

		for (p in pposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			phy<-spl[[1]][2]
		
			full.pred[pos,1:2]<-as.vector(as.matrix(alldata[which(alldata$phylum==phy)[1],3:4]))
		}
	}

	pambi<-parsed.result$ambiguous
	pambi<-pambi[-is.na(pambi)]

	if (length(pambi)>0){

		for (p in pambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			phy<-spl[[1]][2]
		
			full.pred[pos,2]<-phy
		}
	}

	pnext<-unique(c(parsed.result$negatives,parsed.result$ambiguous))
	pnext<-pnext[-is.na(pnext)]
	pnext<-rown[which(rown%in%pnext)]
	
	#### Domain prediction ####

	rank<-'domain'

	input6<-input[pnext,]

	inrank<-sort(names(which(table(alldata[,rank])>1)))

	result.rank<-NULL
	
	for (l in 1:length(inrank)){

		result.rank<-cbind(result.rank,as.matrix(predictOne(rank=rank,tax=inrank[l],input=input6)))
	}

	parsed.result<-parsePrediction(result.rank)
		
	dposi<-parsed.result$positives
	dposi<-dposi[-is.na(dposi)]

	if (length(dposi)>0){

		for (p in dposi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			dom<-spl[[1]][2]
		
			full.pred[pos,1]<-as.vector(as.matrix(alldata[which(alldata$domain==dom)[1],3]))
		}
	}

	dambi<-parsed.result$ambiguous
	dambi<-dambi[-is.na(dambi)]

	if (length(dambi)>0){

		for (p in dambi){

			spl<-strsplit(p,split='-')
			pos<-which(rown==spl[[1]][1])
			dom<-spl[[1]][2]
		
			full.pred[pos,1]<-dom
		}
	}

	full.pred<-as.data.frame(full.pred)

	return(full.pred)

}
