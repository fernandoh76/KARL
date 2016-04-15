#' link2path
#'
#' Identifies and plots significant pathways between pairs of taxa.
#' @param elist is a list of enzymes resulting from getEnzymes().
#' @param conf is a p-value for confidence cut-off.
#' @param col1 color for first taxon.
#' @param col2 color for second taxon.
#' @keywords pathway significance KARL
#' @export
#' @examples
#' enzymes<-getEnzymes(tax1='Helicobacteraceae',tax2='Enterococcaceae',level='family')
#' sig.paths<-link2path(elist=enzymes,conf=0.01,col1='blue',col2='maroon')

link2path<-function(elist,conf=.01,col1='darkgreen',col2='maroon'){
	
	require(seqinr)

	options(warn=-1)

	system('touch link2path_1.out')
	system('touch link2path_2.out')

	elist1<-names(elist[[3]])
	elist2<-names(elist[[4]])

	tax1<-c2s(s2c(names(elist[3]))[-c(1,2)])
	tax2<-c2s(s2c(names(elist[4]))[-c(1,2)])
	
	for (e in elist1){

		cmd<-paste('wget -O tmppath http://rest.kegg.jp/link/pathway/',e,sep='')

		system(cmd,ignore.stderr=T)
		system('cat tmppath >> link2path_1.out')
	}

	for (e in elist2){

		cmd<-paste('wget -O tmppath http://rest.kegg.jp/link/pathway/',e,sep='')

		system(cmd,ignore.stderr=T)
		system('cat tmppath >> link2path_2.out')
	}

	g1<-read.table('link2path_1.out',sep='\t',header=F)
	g2<-read.table('link2path_2.out',sep='\t',header=F)

	system('rm -rf tmppath link2path_*')

	a1<-g1[grep('path:map',as.vector(g1$V2)),]
	a2<-g2[grep('path:map',as.vector(g2$V2)),]
	
	ec1<-gsub('ec:','',as.vector(a1$V1))
	pi1<-gsub('path:map','',as.vector(a1$V2))
	co1<-rep('1',length(ec1))

	ec2<-gsub('ec:','',as.vector(a2$V1))
	pi2<-gsub('path:map','',as.vector(a2$V2))
	co2<-rep('2',length(ec2))

	df<-as.data.frame(rbind(cbind(ec1,pi1,co1),cbind(ec2,pi2,co2)))
	
	colnames(df)<-c('enzyme','pathway','color')

	unique.enzymes<-unique(df$enzyme)
	unique.pathway<-unique(df$pathway)

	mat<-matrix(nrow=length(unique.enzymes),ncol=length(unique.pathway),0)
	
	colnames(mat)<-unique.pathway
	rownames(mat)<-unique.enzymes

	for (x in 1:dim(df)[1]){
		
		wx<-which(unique.enzymes==as.vector(df[x,1]))
		wy<-which(unique.pathway==as.vector(df[x,2]))
		wz<-as.vector(df[x,3])

		mat[wx,wy]<-wz
	}

	total<-apply(mat,2,function(x){length(which(x!=0))})
	ex.g1<-apply(mat,2,function(x){length(which(x==1))})
	ex.g2<-apply(mat,2,function(x){length(which(x==2))})

	counts<-rbind(total,ex.g1,ex.g2)

	pvals<-c()

	for (y in 1:dim(counts)[2]){
		
		grs<-counts[2:3,y]
		obs<-counts[1,y]				

		pvals<-c(pvals,prop.test(grs,c(obs,obs),correct=T)$p.value)
	}

	names(pvals)<-colnames(counts)

	sorted.pvals<-sort(pvals[which(pvals<conf)])

	signi<-format(sorted.pvals,scientific=T,digits=2)

	signi.paths<-names(signi)
	
	system('touch signipaths')

	for (s in signi.paths){

		cmd<-paste('wget -O tmppath http://rest.kegg.jp/find/pathway/',s,sep='')

		system(cmd,ignore.stderr=T)
		system('cat tmppath >> signipaths')
	}

	signi.tab<-read.table('signipaths',sep='\t',header=F)
	
	system('rm -rf signipaths tmppath')


	c1<-gsub('path:map','',as.vector(signi.tab$V1))
	c2<-as.vector(signi.tab$V2)
	c3<-as.data.frame(cbind(c1,c2))

	c4<-as.data.frame(cbind(names(signi),signi))
	
	c5<-merge(c3,c4,by.x='c1',by.y='V1')
	colnames(c5)<-c('pid','description','p-value')
	
	

	aux1<-which(unlist(lapply(apply(mat,1,unique),function(x){'1'%in%x}))==T)
	aux2<-which(unlist(lapply(apply(mat,1,unique),function(x){'2'%in%x}))==T)

	colos<-rep('pink',dim(mat)[2])
	colos[which(colnames(mat)%in%names(signi))]<-'red'

	colcol<-c(rep(col1,length(aux1)),rep(col2,length(aux2)))

	auz<-colnames(mat)
	'%notin%'<-Negate('%in%')
	auz[which(auz%notin%names(signi))]<-''

	mycolors<-colorRampPalette(c('white',col1,col2))(3)
	
	heatmap.2(apply
		(mat,1,as.numeric),
		dendrogram='none',
		Rowv=F,Colv=F,
		tracecol=NULL,
		key=F,labCol='',
		labRow=auz,
		col=mycolors,
		colRow=as.vector(df[,3]),
		RowSideColors=colos,
		ColSideColors=colcol,
		offsetRow=-53.5,
		cexRow=.8)

	legend(x=0.22,y=0.9,
	       legend=c(tax1,tax2),
	       fill=c(col1,col2),
	       cex=.7,border=F)
	
	return(c5)
}
