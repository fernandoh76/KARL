#' updateData
#'
#' Updates all data used in Explorer and Predictor modules without the need of re-installing the package.
#' @param route is the full ftp path where genomes are stored at PATRIC database.
#' @param replace takes a logical indicating if replacing old objects by new ones. If FALSE new objects are stored in a 'tmpdir' in the working directory.
#' @keywords update enzymes presence/absence KARL
#' @export
#' @examples
#' patric<-'ftp://ftp.patricbrc.org/patric2/patric3/genomes'
#' updateData(route=patric,replace=F)

updateData<-function(route){
	
	print('Starting data updating...')

	if (missing(route)){

		stop('"route" argument not specified with no default.')

	} else {	

		print('[1/4] Getting new genomes.')

		system('mkdir tmpdir')
		setwd('./tmpdir')

		require('KEGGREST',quietly=T)
		require('RCurl',quietly=T)

		current<-rownames(alldata)

		dirs<-getURL(route,verbose=F,dirlistonly=T)
		allgenomes<-unlist(strsplit(dirs,'\n'))
		
		news<-setdiff(allgenomes,current)
		
		if (length(news)==0){

			stop('	--> The are not new genomes to add.')

		} else {

			with.pathway<-c()
			without.pathway<-c()
	
			for (n in news){

				pfile<-paste(n,'.PATRIC.pathway.tab',sep='')

				cmd<-paste('wget ',route,n,'/',pfile,sep='')
				system(cmd,ignore.stderr=T,ignore.stdout=T)

				lfile<-list.files('.',pattern='pathway.tab')

				if (pfile%in%lfile){

					tab<-read.csv(pfile,sep='\t',header=T)
					gname<-as.vector(tab[1,2])

					with.pathway<-rbind(with.pathway,c(p,gname))

					print(paste(p,'With pathway annotation',sep=' --> '))

				} else {

					without.pathway<-c(without.pathway,p)

					print(paste(p,'Without pathway annotation',sep=' --> '))
				}
			}
		
			colnames(with.pathway)<-c('id','genome')

			with.pathway<-as.data.frame(with.pathway)

			print(paste('New genomes with pathway annotation -->',dim(with.pathway)[1]))
			print(paste('New genomes without pathway annotation -->',dim(without.pathway)[1]))
			
			print('[2/4] Parsing pathway files.')
	
			new.listEC<-list()
			new.listPW<-list()
	
			for (w in 1:dim(with.pathway)[1]){

				p1<-paste(with.pathway[w,1],'.PATRIC.pathway.tab',sep='')

				print(p1)

				t1<-read.csv(p1,sep='\t',header=T,dec=',')

				new.listEC[[g]]<-unique(as.vector(t1$ec_number))
				new.listPW[[g]]<-unique(as.vector(t1$pathway_id))
			}

			names(new.listEC)<-as.vector(with.pathway[,1])
			names(new.listPW)<-as.vector(with.pathway[,1])

			numEC<-unlist(lapply(new.listEC,length))
			numPW<-unlist(lapply(new.listPW,length))

			allenzymes<-sort(unique(unlist(new.listEC)))

			new.dataEC<-NULL
	
			for (l in 1:length(new.listEC)){

				genzymes<-rep(0,length(allenzymes))	
				genzymes[which(allenzymes%in%new.listEC[[l]])]<-1

				new.dataEC<-rbind(new.dataEC,genzymes)
			}

			colnames(new.dataEC)<-allenzymes
			rownames(new.dataEC)<-names(new.listEC)

			u<-as.vector(with.pathway[,2])
			gs<-unlist(lapply(u,function(x){strsplit(x,' ')[[1]][1]}))
			ugs<-unique(gs)

			print('[3/4] Seeking for new taxa')

			taxa.genus<-as.vector(taxa$genus)

			'%notin%'<-Negate('%in%')
			notin.taxa<-ugs[which(ugs%notin%taxa.genus)]

			library(taxize)

			taxa.levels<-c('phylum','class','order','family','genus')

			plus.taxa<-NULL

			for (n in notin.taxa){

				plus.taxa<-rbind(plus.taxa,tax_name(n,db='ncbi',get=taxa.levels))
			}

			alltaxa<-rbind(taxa[,3:8],plus.taxa[,3:8])
			ualltaxa<-unique(alltaxa)

			taxa.withpathway<-NULL
	
			for (g in 1:length(gs)){

				this<-which(ualltaxa$genus==gs[g])

				if (length(this)>0){

					taxa.withpathway<-rbind(taxa.withpathway,ualltaxas[this,])

				} else {

					taxa.withpathway<-rbind(taxa.withpathway,rep('NA',5))
				}
			}

			new.genomes<-cbind(with.pathway,taxa.withpathway)

			nalist<-apply(taxa.withpathway,1,unique)
			nare<-which(lapply(nalist,unique)=='NA')

			for (e in nare){

				ck<-as.vector(new.genomes[e,2])
				tx<-tax_name(ck,db='ncbi',get=taxa.levels)

				new.genomes[e,3:8]<-tx[3:8]
			}

			newdata<-cbind(newgenomes,new.dataEC)
			rownames(newdata)<-NULL

			print('[4/4] Updating enzymes matrix')
			
			###############
			# New alldata #
			###############

			enzymes<-colnames(alldata)
			newenzymes<-colnames(newdata)

			addenzymes<-setdiff(newenzymes,enzymes)

			orderdata<-newdata[,1:8]
						
			for (z in enzymes){

				orderdata<-cbind(orderdata,newdata[,which(newenzymes==z)])
			}

			new.alldata<-rbind(alldata,orderdata)

			for (a in addenzymes){

				new<-c(rep(0,dim(alldata)[2]),newdata[,which(newenzymes==a)])
				new.alldata<-cbind(new.alldata,new)
			}

			colnames(new.alldata)<-c('id','genome',taxa.levels,enzymes,newenzymes)

			############
			# New taxa #
			############

			new.taxa<-new.alldata[,1:7]

			################
			# New prcompEC #
			################

			new.prcompEC<-prcomp(new.alldata[,9:dim(new.alldata)[2]])

			#################
			# New koenzymes #
			#################

			enzymes<-colnames(new.alldata)[-c(1:9)]

			new.koenzymes<-c()

			for (e in enzymes){

				new.koenzymes<-c(new.koenzymes,keggLink('ko',e))
			}

			###################
			# New li.pathways #
			###################

			new.pathways<-list.files('.',pattern='PATRIC.pathway.tab')

			new.path.enzymes<-NULL

			for (p in new.pathways){
				
				pa<-read.csv(p,sep='\t',header=T)
				new.path.enzymes<-rbind(new.path.enzymes,pa)
			}

			unipath<-unique(new.path.enzymes$pathway_id)

			new.li.pathways<-NULL

			for (x in 1:length(unipath)){
	
				ecun<-unique(as.vector(new.path.enzymes[which(new.path.enzymes$pathway_id==unipath[x]),'ec_number']))
				new.li.pathways[[x]]<-ecun
			}

			names(new.li.pathways)<-unipath

			keys<-unique(c(names(li.pathways),names(new.li.pathways)))
			this.li.pathways<-lapply(setNames(mapply(c,li.pathways[keys],new.li.pathways[keys]),keys),unique)

			##################
			# Save new files #
			##################
	
			if (replace==T){

				alldata.path<-system.file('data','alldata.Rdata',package='KARL')
				alldata<-new.alldata
				save(alldata,file=alldata.path)

				koenzymes.path<-system.file('data','koenzymes.Rdata',package='KARL')
				koenzymes<-new.koenzymes
				save(koenzymes,file=koenzymes.path)

				li.pathways.path<-system.file('data','li.pathways.Rdata',package='KARL')
				li.pathways<-this.li.pathways
				save(li.pathways,file=li.pathways.path)

				prcompEC.path<-system.file('data','prcompEC.Rdata',package='KARL')
				prcompEC<-new.prcompEC
				save(prcompEC,file=prcompEC.path)

				taxa.path<-system.file('data','taxa.Rdata',package='KARL')
				taxa<-new.taxa
				save(taxa,file=taxa.path)

				new.genomes.path<-gsub('taxa.Rdata','new.genomes.Rdata',taxa.path)
				save(new.genomes,file=new.genomes.path)

				system('rm -rf tmpdir')
			}
		}
	}
}	
