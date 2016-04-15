#' identifiR
#'
#' Identifies a genome based on 16S similarity with SILVA database.
#' @param inseq Input set of genomic contigs or closed genome.
#' @param taxout Filename for taxonomic report.
#' @param geneout Filename for the resulting 16S sequence found in genome.
#' @keywords 16S Rnammer SILVA KARL
#' @export
#' @examples
#' identifiR(inseq='contigs.fasta',taxout='taxreport.txt',geneout='contigs16S.fasta')

identifiR<-function(inseq,taxout,geneout){

	rnammer<-system.file('rnammer','rnammer',package='KARL')
	blastn<-system.file('blast','blastn',package='KARL')
	blastdbcmd<-system.file('blast','blastdbcmd',package='KARL')
	silva<-gsub('.nsq','',system.file('silva','silva.nsq',package='KARL'))
	
	require(seqinr)

	rnam<-paste(rnammer,' -S bact -m ssu -f ',geneout,' < ',inseq,sep='')
	print('Searching 16S gene...')

	system(rnam,ignore.stdout=TRUE)

	if (file.info(geneout)$size>0){

		print('16S gene FOUND in contigs!')
		print('Proceeding to taxonomic identification...')

		steen<-read.fasta(geneout)
		s<-lapply(getSequence(steen),toupper)

		write.fasta(s[[1]],names='query1',file.out='tmp.fasta')

		system(paste(blastn,"-query tmp.fasta -db",silva,"-evalue 0.001 -outfmt '6 qseqid sseqid pident length gaps slen evalue' -out tmp.blout"))
		system('rm -rf tmp.fasta')

		if (file.info('tmp.blout')$size>0){

			blout<-read.table('tmp.blout',sep='\t',header=FALSE)

			system('rm -rf tmp.blout')

			colnames(blout)<-c('qid','sid','pid','lgth','gaps','slen','evalue')
			blout$qcov<-(blout$lgth-blout$gaps)/blout$slen
			bl<-blout[which(blout$qcov>=.98),]
			bl<-bl[order(blout$pid,decreasing=TRUE),]
			hit<-bl[1,]

			if (hit$pid>=99){

				id<-as.vector(bl[1,'sid'])
				dbm<-paste(blastdbcmd,'-entry',id,'-db',silva,'> closest16S.fasta')

				system(dbm)

				hdr<-system("grep '>' closest16S.fasta",intern=TRUE)
				uhdr<-unlist(strsplit(hdr,';'))
				uhdr[1]<-strsplit(uhdr[1],' ')[[1]][2]

				cat('Your strain was identified at species level with identity > 99%',file=taxout,sep='\n\n',append=TRUE)
				cat('###### Taxonomy for the closest 16S sequence ######',file=taxout,sep='\n',append=TRUE)			
				cat(uhdr,file=taxout,sep='\n',append=TRUE)

				print(paste('Identification finished, results in',taxout))

			} else if (hit$pid>=97 & hit$pid<99){

				id<-as.vector(bl[1,'sid'])
				dbm<-paste(blastdbcmd,'-entry',id,'-db',silva,'> closest16S.fasta')

				system(dbm,ignore.stdout=TRUE)

				hdr<-system("grep '>' closest16S.fasta",intern=TRUE)
				uhdr<-unlist(strsplit(hdr,';'))
				uhdr[1]<-strsplit(uhdr[1],' ')[[1]][2]
				cat('Your strain was identified at genus with identity < 99%',file=taxout,sep='\n',append=TRUE)
				cat('It probably represents a (new) species not present in the current 16S database',file=taxout,sep='\n\n',append=TRUE)
				cat('###### Taxonomy for the closest 16S sequence ######',file=taxout,sep='\n',append=TRUE)			
				cat(uhdr,file=taxout,sep='\n',append=TRUE)

				print(paste('Identification finished, results in',taxout))

			} else if (hit$pid<97){

				id<-as.vector(bl[1,'sid'])
				dbm<-paste(blastdbcmd,'-entry',id,'-db',silva,'> closest16S.fasta')
				system(dbm,ignore.stdout=TRUE)

				hdr<-system("grep '>' closest16S.fasta",intern=TRUE)
				uhdr<-unlist(strsplit(hdr,';'))
				uhdr[1]<-strsplit(uhdr[1],' ')[[1]][2]

				cat('Your strain presented < 97% identity with all 16S sequences in the database',file=taxout,sep='\n\n',append=TRUE)
				cat('It probably represents a new species',file=taxout,sep='\n',append=TRUE)
				cat('###### Taxonomy for the closest 16S sequence ######',file=taxout,sep='\n',append=TRUE)			
				cat(uhdr,file=taxout,sep='\n',append=TRUE)

				print(paste('Identification finished, results in', taxout))

			}		
		} else {

			print('No significant hits in 16S database!')
		}

	} else {

		print('No 16S genes found in contigs!')	
	}
}
