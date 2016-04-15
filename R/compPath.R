# Uses 'alldata' and 'koenzymes' objects

#' compPath
#'
#' This function calculates the frequency for each enzyme in a pair of given taxa and displays these frequencies over a selected pathway. Outputs a PNG file.
#' @param rank Takes 'domain', 'phylum', 'class', 'order', 'family' or 'genus'.
#' @param tax1 The first taxon name.
#' @param tax2 The second taxon name.
#' @param pid Takes a 5-number code for any KEGG Pathway identifier. Can use findPath() function to find the code given a keyword.
#' @keywords pathways KEGG enzymes KARL
#' @export
#' @examples
#' compPath(rank='class',tax1='Betaproteobacteria',tax2='Mollicutes',pid='00550')

compPath<-function(rank,tax1,tax2,pid){
	
	require(pathview)
	require(png)

	x<-which(colnames(alldata)==rank)
	l1<-as.vector(alldata[,x])

	s1<-alldata[which(l1==tax1),9:dim(alldata)[2]]
	s2<-alldata[which(l1==tax2),9:dim(alldata)[2]]

	a1<-apply(s1,2,sum)/dim(s1)[1]
	a2<-apply(s2,2,sum)/dim(s2)[1]
	a3<-rbind(a1,a2)
	
	kofreq<-NULL
	a3col<-colnames(a3)
	konames<-gsub('ec:','',names(koenzymes))
	
	for (a in 1:length(a3col)){
		grp<-which(konames==a3col[a])
		le<-length(grp)
		if (le>0){
			aux<-NULL
			for (l in 1:le){
				aux<-cbind(aux,a3[,a])
			}
			colnames(aux)<-gsub('ko:','',koenzymes[grp])
			kofreq<-cbind(kofreq,aux)
		}
	}
	
	pathview(gene.data=t(kofreq),
		pathway.id=pid,
		species='ko',
		kegg.native=T,
		limit=list(gene=c(0,1)),
		low=list(gene='grey'),
		mid=list(gene='orange'),
		high=list(gene='red'),
		bins=list(gene=10),
		out.suffix='intermed',
		)
	
	img<-readPNG(list.files(pattern='intermed.multi.png'))
	system('rm -rf *.intermed.multi.png')

	h<-dim(img)[1]
	w<-dim(img)[2]
	
	outsuf<-paste('path',pid,sep=':')
	outname<-paste(outsuf,tax1,'vs',tax2,'png',sep='.')

	png(outname,width=w,height=h)
	par(mar=c(0,0,0,0),xpd=NA,mgp=c(0,0,0),oma=c(0,0,0,0),ann=F)
	plot.new()
	plot.window(0:1,0:1)

	usr<-par('usr')
	rasterImage(as.raster(img),xleft=0,xright=1,ybottom=0,ytop=1)

	text(.98,.89,paste('Left::',tax1,sep=''),cex=1.5,adj=c(1,NA))
	text(.98,.86,paste('Right::',tax2,sep=''),cex=1.5,adj=c(1,NA))
		
	dev.off()

}
