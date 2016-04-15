#' rankSlope
#'
#' Calculates Information Gain for each enzyme and rank the values using the output from dropUninformative().
#' @param idata Takes the output from dropUninformative() function.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax The taxon name. The values of rank and tax must be the same than in the respective dropUninformative() run.
#' @param plotting Is a logical indicating if producing plots or not.
#' @keywords Information Gain enzymes presence/absence KARL
#' @export
#' @examples
#' streptococcus.ranked<-rankSlope(idata=streptococcus.droped,rank='genus',tax='Streptococcus',plotting=T)

rankSlope<-function(idata,rank,tax,plotting=T){
	
	print(Sys.time())

	require(RWeka,quietly=T)
	require(segmented,quietly=T)

	d<-which(colnames(idata)==rank)

	resp<-as.vector(idata[,d])
	resp[which(resp!=tax)]<-'No'

	gdata<-cbind(idata[,9:dim(idata)[2]],resp)
	
	rankEC<-GainRatioAttributeEval(resp~.,data=gdata)
	rankEC<-rankEC[-which(rankEC<=0.01)]

	if (length(rankEC)==0){
		
		print('No informative enzymes!')
		print(paste('  --> Retaining the original set of',dim(idata)[2]-8,'enzymes'))

		print(Sys.time())
		print('---------------------------------')

		return(idata)

	} else {
	
		sorted.rank<-as.vector(rev(sort(rankEC)))
	
		named.rank<-names(rev(sort(rankEC)))
		inx<-seq(1:length(sorted.rank))
		glm.rank<-glm(inx~sorted.rank)

		lx<-length(inx)

		if (lx>10){

			print('More than 10 informative enzymes!')
			print(paste('  --> Doing Davies Test over',lx,'enzymes'))

			ss<-c(10:lx)

			best<-c()

			for (s in ss){
				davies<-davies.test(obj=glm.rank,seg.Z=~sorted.rank,k=s)
				best<-c(best,as.vector(davies$statistic))
			}

			sel<-summary(best)[4]

			if (plotting==T){

				x11()

		 		plot(density(best),
			     		col='grey',
		            		lwd=3,
		             		xlab='Gain Ratio',
		             		main=paste("Davies's slope change test for",rank,tax)
		            		)

				legend('topright',
			       		legend=c("Davies' test cutoff"),
		               		lty=2,
		               		lwd=4,
		               		col='red',
			       		cex=.7
			       		)

				abline(v=sel,col='red',lty=2,lwd=2)
		
				x11()
	
				selcol<-rep('grey',length(sorted.rank))
				selcol[which(sorted.rank>=sel)]<-'black'

				plot(sorted.rank,
			     	type='l',
		             	col='black',
			     	main=paste('Most informative attributes for',rank,tax),
		             	ylab='Ranked Information Gain Ratio',
			     	xlab='Ranking'
			    	)

				legend('topright',
			       		legend=c("Davies' test cutoff"),
		               		lty=2,
		               		lwd=4,
		               		col='red',
			       		cex=.7
			       		)

				points(sorted.rank,pch=20,col=selcol)
				abline(h=sel,lty=2,col='red')

				Sys.sleep(12)
				graphics.off()
			}
	
			selected<-named.rank[which(sorted.rank>=sel)]
	
			sdata<-cbind(idata[,1:8],idata[,which(colnames(idata)%in%selected)])
	
			print(paste('  --> Retaining',length(selected),'enzymes'))
			print(Sys.time())
			print('---------------------------------')

			return(sdata)
	
		} else if (lx<=10){

			print('Less than 10 informative enzymes!')
			print(paste('  --> Skipping Davies Test over',lx,'enzymes'))
			print(paste('  --> Retaining',lx,'enzymes'))
	
			sdata<-cbind(idata[,1:8],idata[,which(colnames(idata)%in%named.rank)])
	
			print(Sys.time())
			print('---------------------------------')		
	
			return(sdata)
	
		}
	}
}
