#' prcompTax
#'
#' Plots Principal Components in 2D and 3D and allows taxa coloring.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param taxa A vector of names for desired taxa.
#' @param colors A vector of colors of same length as taxa.
#' @param axes A vector of length 2 or 3 indicating which axes to be plotted.
#' @param legpos The position of the legend.
#' @keywords PCA taxa enzymes presence/absence KARL.
#' @export
#' @examples
#' prcompTax(rank='genus',taxa=c('Streptomyces','Lactobacillus'),axes=c(1,2),colors=c('red','darkgreen'),legpos='topright')

prcompTax<-function(rank,taxa,colors,axes=c(1,2),legpos='topright'){

	require(rgl,quietly=T)

	'%notin%'<-Negate('%in%')
	taxa.ranks<-c('domain','phylum','class','order','family','genus')	

	if (rank%notin%taxa.ranks){
		stop('Argument "rank" was wrongly specified.')

	} else if (missing(taxa)){
		stop('Argument "taxa" was not specified.')

	} else if (missing(taxa) & missing(rank)){
		stop('Arguments "taxa" and "rank" were not specified.')

	} else {

		prx<-prcompEC$x
	
		lev<-which(colnames(alldata)==rank)
		lea<-length(axes)

		if (missing(colors)){

			if (lea==2){
				
				plot(prx[,axes[1]],
		     		     prx[,axes[2]],
		     		     pch=19,
  		     		     xlab=paste('PC',axes[1],sep=''),
		     	             ylab=paste('PC',axes[2],sep='')
				)

			} else if (lea==3){

				plot3d(prx[,axes[1]],
		       		       prx[,axes[2]],
		                       prx[,axes[3]],
				       col='grey',
		                       xlab=paste('PC',axes[1],sep=''),
		                       ylab=paste('PC',axes[2],sep=''),
		                       zlab=paste('PC',axes[3],sep='')
				)

			} else if (lea>3 | lea<2){
				stop('Argument "axes" was wrongly specified.')
			}

		} else {

			if (lea==2){
				plot(prx[,axes[1]],
		     		     prx[,axes[2]],
		     		     pch=19,
	             		     col='lightgrey',
  		     		     xlab=paste('PC',axes[1],sep=''),
		     	             ylab=paste('PC',axes[2],sep='')
				)
				
				for (x in 1:length(taxa)){
					orgs<-which(alldata[,lev]==taxa[x])
				
					if (length(orgs)==0){
						warning(paste('Taxon',taxa[x],'not present.'))
					} else {
						points(prx[orgs,axes[1]],
						     prx[orgs,axes[2]],
						     pch=19,
						     col=colors[x],
					             xlab=paste('PC',axes[1],sep=''),
						     ylab=paste('PC',axes[2],sep='')
						)
					}
				}

				legend(legpos,
					legend=taxa,
					col=colors,
                                	cex=.7,
					pch=19
				)

			} else if (lea==3){

				prcolors<-rep('lightgrey',dim(alldata)[1])
			
				for (x in 1:length(taxa)){
					orgs<-which(alldata[,lev]==taxa[x])
	
					if (length(orgs)==0){
						w<-paste('Taxon',
						         taxa[x],'not present.')
	
						warning(w,call.=F)
			
					} else {
						prcolors[orgs]<-colors[x]
					}
	
				}

				plot3d(prx[,axes[1]],
		       		       prx[,axes[2]],
		                       prx[,axes[3]],
		                       col=prcolors,
				       type='s',
				       radius=.2,
				       lwd=0,
		                       xlab=paste('PC',axes[1],sep=''),
		                       ylab=paste('PC',axes[2],sep=''),
		                       zlab=paste('PC',axes[3],sep='')
				)
				
				legend3d(legpos,
					legend=taxa,
					col=colors,
					pch=19,
					cex=.7
				)

			} else if (lea>3 | lea<2){
				stop('Argument "axes" was wrongly specified.')
			}
		}			
	}
}

