#' assignR
#'
#' This function creates a presence/absence vector of enzymes for an input genome.
#' @param genes Is the result from findR function, either in nucleotides or amino acids.
#' @param type Takes 'prot' or 'nucl' depending on genes parameter datatype.
#' @param qc Is the query coverage treshold selected to consider a good enough hit.
#' @param id Is the identity treshold selected to consider a good enough hit.
#' @param thrds Is the number of threads to be used in Blast searches.
#' @keywords blastp enzymes KARL presence/absence
#' @export
#' @examples
#' p.a.vector<-assignR(genes='pegs.fasta',type='prot',qc=90,id=50,thrds=8)

assignR<-function(genes,type,qc=90,id=50,thrds=1){
	
	require(seqinr,quietly=T)

	###################################
        # internal function: fastasplit() #
        ###################################

	fastasplit<-function(f,n){

		fasta<-read.fasta(f)
		seqs<-lapply(getSequence(fasta),toupper)
		vecu<-seq(1:length(seqs))
		
		s<-split(vecu,cut(seq_along(vecu),n,labels=F))

		for (x in 1:length(s)){
		
			write.fasta(seqs[s[[x]]],names=paste('name',s[[x]],sep=''),file.out=paste('splitted_',x,'.faa',sep=''))
		}
	}

	##################################

	enzymes<-colnames(alldata)[9:dim(alldata)[2]]

	dbspath<-system.file('databases',package='KARL')
	blastp<-paste(system.file('blast',package='KARL'),'blastp',sep='/')

	result<-c()
	
	if (type=='prot'){

		fastasplit(f=genes,n=thrds)

		for (y in enzymes){

			cmd1<-paste("seq 1 ",thrds,
				    " | xargs -P ",thrds," -I {} ",blastp,"-query splitted_{}.faa -db ",dbspath,"/",y,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

			cmd2<-'cat partial_* > blout'
	
			system(cmd1,wait=T)
			system(cmd2,wait=T)

			print(paste('Searching for enzyme ',y,sep=''))
		
			if (file.info('blout')$size==0){

				result<-c(result,0)

				system('rm -rf blout partial_*')

			} else {

				blout<-read.table('blout',sep='\t',header=F)

				system('rm -rf blout partial_*')

				colnames(blout)<-c('qid','sid','pid','gaps','len','qlen')

				blout$qcov<-(blout$len-blout$gaps)/blout$qlen*100

				hits<-which(blout$qcov>=qc & blout$pid>=id)

				if (length(hits)>0){

					result<-c(result,1)

				} else {

					result<-c(result,0)
				}
			}
		}

	} else {
		
		ffn<-read.fasta(genes)

		sec<-lapply(getSequence(ffn),toupper)

		faa<-lapply(sec,function(j){translate(j,numcode=11)})

		write.fasta(faa,names=names(ffn),file.out='aux.faa')

		fastasplit(f='aux.faa',n=thrds)

		for (y in enzymes){

			cmd1<-paste("seq 1 ",thrds,
				    " | xargs -P ",thrds," -I {} ",blastp,"-query splitted_{}.faa -db ",dbspath,"/",y,
				    " -evalue 0.001 -outfmt '6 qseqid sseqid pident gaps length qlen' -out partial_{}",sep='')

			cmd2<-'cat partial_* > blout'
	
			system(cmd1,wait=T)
			system(cmd2,wait=T)

			print(paste('Searching for enzyme ',y,sep=''))
		
			if (file.info('blout')$size==0){

				result<-c(result,0)

				system('rm -rf blout partial_* aux.faa')

			} else {

				blout<-read.table('blout',sep='\t',header=F)

				system('rm -rf blout partial_* splitted_* aux.faa')

				colnames(blout)<-c('qid','sid','pid','gaps','len','qlen')

				blout$qcov<-(blout$len-blout$gaps)/blout$qlen*100

				hits<-which(blout$qcov>=qc & blout$pid>=id)

				if (length(hits)>0){

					result<-c(result,1)

				} else {

					result<-c(result,0)
				}
			}
		}			
	}

	system('rm -rf splitted_*')

	result2<-data.frame(t(result))
	colnames(result2)<-colnames(alldata)[9:dim(alldata)[2]]
	rownames(result2)<-strsplit(genes,'.')[[1]][1]

	return(result2)
}
