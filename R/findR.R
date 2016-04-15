#' findR
#'
#' This function performs genes prediction using Prodigal and/or Glimmer.
#' @param genome takes the file name for contigs or completed genome.
#' @param complete takes a logical indicating if the genome is closed or not.
#' @param program takes 'prodigal', 'glimmer' or 'both' for choosing the prediction algorithm. When 'both' is selected the intersection of individual searches is outputed using blastp.
#' @param outpre is the prefix for output files.
#' @keywords prodigal glimmer genes KARL
#' @export
#' @examples
#' pred.genes<-findR(genome='contigs.fasta',complete=F,program='both',outpre='pegs')



findR<-function(genome,complete=F,program='both',outpre){
	
	print('# Gene finding started #')
	print('------------------------')

	makeblastdb<-paste(system.file('blast',package='KARL'),'makeblastdb',sep='/')
	blast<-paste(system.file('blast',package='KARL'),'blastp',sep='/')
	blastdbcmd<-paste(system.file('blast',package='KARL'),'blastdbcmd',sep='/')
	prodigal<-system.file('prodigal','prodigal',package='KARL')
	glimmer<-system.file('glimmer/scripts','g3-iterated.csh',package='KARL')
	extract<-system.file('glimmer/bin','extract',package='KARL')

	glpath<-system.file('glimmer/scripts',package='KARL')
	glbin<-system.file('glimmer/bin',package='KARL')
	elph<-system.file('glimmer',package='KARL')
	
	Sys.setenv("awkpath"=glpath)
	Sys.setenv("glimmerpath"=glbin)
	Sys.setenv("elphbin"=elph)

	require(seqinr)
	
	######################################
	# Define internal function jointPred #
	######################################

	# Outputs the intersection between Glimmer and Prodigal predictions, it is used
	# when PROGRAM parameter in FINDR function is set to 'both'

	jointPred<-function(outfile){
	
		print('[3/3] Finding joint prediction')
		pdbp<-paste(makeblastdb,'-in prodigal.faa -dbtype prot -hash_index -parse_seqids -title prdb -out prdb')
		ndbp<-paste(makeblastdb,'-in prodigal.ffn -dbtype nucl -hash_index -parse_seqids -title prdbn -out prdbn')
		pdbg<-paste(makeblastdb,'-in glimmer.faa -dbtype prot -hash_index -parse_seqids -title gldb -out gldb')
	
		system(pdbp,ignore.stdout=T)
		system(ndbp,ignore.stdout=T)
	
		system(pdbg,ignore.stdout=T)	

		blast1<-paste(blast," -query prodigal.faa -db gldb -evalue 0.001", 
			" -outfmt '6 qseqid sseqid pident length gaps qlen' -out blast1",
			sep='')

		blast2<-paste(blast," -query glimmer.faa -db prdb -evalue 0.001", 
			" -outfmt '6 qseqid sseqid pident length gaps qlen' -out blast2",
			sep='')
	
		system(blast1)
		system(blast2)

		blout1<-read.table('blast1',sep='\t',header=F)
		blout2<-read.table('blast2',sep='\t',header=F)

		colnames(blout1)<-c('qid','sid','pid','len','gaps','qlen')
		colnames(blout2)<-c('qid','sid','pid','len','gaps','qlen')

		blout1$qcov<-(blout1$len-blout1$gaps)/blout1$qlen
		blout2$qcov<-(blout2$len-blout2$gaps)/blout2$qlen


		b1<-blout1[which(blout1$pid>=90 & blout1$qcov>=.90),c(1,2)]
		b2<-blout2[which(blout2$pid>=90 & blout2$qcov>=.90),c(2,1)]
	
		p1<-apply(b1,1,function(x){paste(x,collapse='_')})
		p2<-apply(b2,1,function(x){paste(x,collapse='_')})	
	
		ptab<-table(c(p1,p2))
		bids<-names(which(ptab>1))
		ids<-unlist(lapply(bids,function(z){strsplit(z,'_')[[1]][1]}))

		cat(ids,sep='\n',file='ids.txt')

		system(paste(blastdbcmd,'-entry_batch ids.txt -db prdb > bothx.faa'))
		system(paste(blastdbcmd,'-entry_batch ids.txt -db prdbn > bothx.ffn'))
		system('rm -rf ids.txt blast1 blast2 prdb* gldb* glimmer* prodigal*')

		bothn<-read.fasta('bothx.ffn')
		bothp<-read.fasta('bothx.faa')
	
		seqn<-lapply(getSequence(bothn),toupper)
		seqp<-lapply(getSequence(bothp),toupper)
		namn<-paste('both',seq(1:length(seqn)),sep='')
		namp<-paste('both',seq(1:length(seqp)),sep='')
	
		write.fasta(seqn,names=namn,file.out=paste(outfile,'ffn',sep='.'))
		write.fasta(seqp,names=namp,file.out=paste(outfile,'faa',sep='.'))
		system('rm -rf bothx.*')

	}

	if (program=='both'){

		if (complete==F){
			
			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){

				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')

			pcmd<-paste(prodigal,
				' -i genome.fasta ', 
				'-d prodigal.ffn -o prodigal.out &> plog',
				sep=''
				)
			
			gcmd<-paste(glimmer,
				' genome.fasta glimmer &> glog',
				sep=''
				)

			print('[1/3] Running Prodigal')
			system(pcmd)

			print('[2/3] Running Glimmer')
			system(gcmd)
						
			system(paste(extract,'genome.fasta glimmer.run1.predict &> elog > glimm.ffn'))
			system('rm -rf prodigal.out glimmer.* genome.fasta plog glog elog')

			prod<-read.fasta('prodigal.ffn')
			pffn<-lapply(getSequence(prod),toupper)
			pfaa<-lapply(pffn,function(k){translate(k,numcode=11)})			
			pnam<-paste('prodigal',seq(1:length(prod)),sep='')
			write.fasta(pffn,names=pnam,file.out='prodigal.ffn')
			write.fasta(pfaa,names=pnam,file.out='prodigal.faa')
			
			glim<-read.fasta('glimm.ffn')
			gffn<-lapply(getSequence(glim),toupper)		
	
			ns<-grep(c2s(nstr),lapply(gffn,c2s))

			if (length(ns)>0){

				gffn<-gffn[-ns]
				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out='glimmer.ffn')
				write.fasta(gfaa,names=gnam,file.out='glimmer.faa')

			} else {

				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out='glimmer.ffn')
				write.fasta(gfaa,names=gnam,file.out='glimmer.faa')
			}
			
			system('rm -rf glimm.ffn')
			jointPred(outfile=outpre)
		
			print('-------------------------')
			print('# Gene finding finished #')
			
		} else if (complete==T){
			
			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){

				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')

			pcmd<-paste(prodigal,
				' -i genome.fasta -c ', 
				'-d prodigal.ffn -o prodigal.out &> glog',
				sep=''
				)
			
			gcmd<-paste(glimmer,
				' genome.fasta glimmer &> plog',
				sep=''
				)

			print('[1/3] Running Prodigal')			
			system(pcmd)

			print('[2/3] Running Glimmer')	
			system(gcmd)
			
			system(paste(extract,'genome.fasta glimmer.run1.predict &> elog > glimm.ffn'))
			system('rm -rf prodigal.out glimmer.* genome.fasta glog plog elog')

			prod<-read.fasta('prodigal.ffn')
			pffn<-lapply(getSequence(prod),toupper)
			pfaa<-lapply(pffn,function(k){translate(k,numcode=11)})			
			pnam<-paste('prodigal',seq(1:length(prod)),sep='')
			write.fasta(pffn,names=pnam,file.out='prodigal.ffn')
			write.fasta(pfaa,names=pnam,file.out='prodigal.faa')			
	
			glim<-read.fasta('glimm.ffn')
			gffn<-lapply(getSequence(glim),toupper)
			
			ns<-grep(c2s(nstr),lapply(gffn,c2s))

			if (length(ns)>0){

				gffn<-gffn[-ns]
				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out='glimmer.ffn')
				write.fasta(gfaa,names=gnam,file.out='glimmer.faa')

			} else {

				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out='glimmer.ffn')
				write.fasta(gfaa,names=gnam,file.out='glimmer.faa')
			}
			
			system('rm -rf glimm.ffn')
			jointPred(outfile=outpre)

			print('-------------------------')
			print('# Gene finding finished #')	
		}	
		
	} else if (program=='glimmer'){
	
		if (complete==F){
		
			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){
	
				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')
			
			gcmd<-paste(glimmer,
				' genome.fasta glimmer &> glog',
				sep=''
				)

			print('[1/1] Running Glimmer')

			system(gcmd)			
	
			system(paste(extract,'genome.fasta glimmer.run1.predict &> elog > glimm.ffn'))
			system('rm -rf glimmer.* genome.fasta elog')
	
			glim<-read.fasta('glimm.ffn')
			system('rm -rf glimm.ffn glog')
			gffn<-lapply(getSequence(glim),toupper)

			ns<-grep(c2s(nstr),lapply(gffn,c2s))

			if (length(ns)>0){

				gffn<-gffn[-ns]
				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out=paste(outpre,'ffn',sep='.'))
				write.fasta(gfaa,names=gnam,file.out=paste(outpre,'faa',sep='.'))

			} else {

				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out=paste(outpre,'ffn',sep='.'))
				write.fasta(gfaa,names=gnam,file.out=paste(outpre,'faa',sep='.'))
			}

			print('-------------------------')
			print('# Gene finding finished #')

		} else if (complete==T){

			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){
				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')

			gcmd<-paste(glimmer,
				' genome.fasta glimmer &> glog',
				sep=''
				)

			print('[1/1] Running Glimmer')

			system(gcmd)			
	
			system(paste(extract,'genome.fasta glimmer.run1.predict &> elog > glimm.ffn'))
			system('rm -rf glimmer.* genome.fasta elog')
	
			glim<-read.fasta('glimm.ffn')
			system('rm -rf glimm.ffn glog')
			gffn<-lapply(getSequence(glim),toupper)

			ns<-grep(c2s(nstr),lapply(gffn,c2s))

			if (length(ns)>0){

				gffn<-gffn[-ns]
				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out=paste(outpre,'ffn',sep='.'))
				write.fasta(gfaa,names=gnam,file.out=paste(outpre,'faa',sep='.'))

			} else {

				gfaa<-lapply(gffn,function(i){translate(i,numcode=11)})
				gnam<-paste('glimmer',seq(1:length(glim)),sep='')
				write.fasta(gffn,names=gnam,file.out=paste(outpre,'ffn',sep='.'))
				write.fasta(gfaa,names=gnam,file.out=paste(outpre,'faa',sep='.'))
			}

			print('-------------------------')
			print('# Gene finding finished #')		
		}
	
	} else if (program=='prodigal'){

		if (complete==F){

			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){

				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')

			pcmd<-paste(prodigal,
				' -i genome.fasta ', 
				'-d prodigal.ffn -o prodigal.out &> plog',
				sep=''
				)

			print('[1/1] Running Prodigal')
			system(pcmd)

			system('rm -rf prodigal.out genome.fasta plog')

			prod<-read.fasta('prodigal.ffn')
			system('rm -rf prodigal.ffn')
			pffn<-lapply(getSequence(prod),toupper)
			pfaa<-lapply(pffn,function(k){translate(k,numcode=11)})			
			pnam<-paste('prodigal',seq(1:length(prod)),sep='')
			write.fasta(pffn,names=pnam,file.out=paste(outpre,'ffn',sep='.'))
			write.fasta(pfaa,names=pnam,file.out=paste(outpre,'faa',sep='.'))

			print('-------------------------')
			print('# Gene finding finished #')			
				
		} else if (complete==T){
			
			nstr<-rep('N',50)
			contigs<-read.fasta(genome)
			fasta<-lapply(getSequence(contigs),toupper)
			geno<-c()
			
			for (f in 1:length(fasta)){
				geno<-c(geno,nstr,fasta[[f]])
			}

			write.fasta(geno,names='genome',file.out='genome.fasta')

			pcmd<-paste(prodigal,
				' -i genome.fasta -c ', 
				'-d prodigal.ffn -o prodigal.out &> plog',
				sep=''
				)

			print('[1/1] Running Prodigal')
			system(pcmd)

			system('rm -rf prodigal.out genome.fasta plog')

			prod<-read.fasta('prodigal.ffn')
			system('rm -rf prodigal.ffn')
			pffn<-lapply(getSequence(prod),toupper)
			pfaa<-lapply(pffn,function(k){translate(k,numcode=11)})			
			pnam<-paste('prodigal',seq(1:length(prod)),sep='')
			write.fasta(pffn,names=pnam,file.out=paste(outpre,'ffn',sep='.'))
			write.fasta(pfaa,names=pnam,file.out=paste(outpre,'faa',sep='.'))

			print('-------------------------')
			print('# Gene finding finished #')			
		}

	} else {
		stop('Argument "program" was wrongly set, please set it to "both", "prodigal" or "glimmer".')
	}
}

