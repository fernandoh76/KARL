### De novo reads assembly using SPAdes
### r1 is the name of fastq file containing reads 1
### r2 is the name of fastq file containing reads 2 (if left empty single read sequencing is assumed)
### Mode 'fast' runs SPAdes by default (faster). Mode 'accurate' runs SPAdes with wider K-mer range and using '--careful' flag (slower).
### out is the name of the output contigs file
### report is the name for the assembly report
### path is the full path were SPAdes is installed

#' assemblR
#'
#' This function allows to automatically assemble a genome from single or pair-end reads using SPAdes.
#' @param r1 The filename for fastq containing reads 1.
#' @param r2 The filename for fastq containing reads 2 (left empty if single reads are assumed).
#' @param mode Takes 'fast' for SPAdes running by default (faster), or 'accurate' for SPAdes running with a wider K-mer range and using the '--careful' flag (slower).
#' @param out The name of the output file with contigs.
#' @param report The name of the assembly report file.
#' @keywords SPAdes assembly KARL
#' @export
#' @examples
#' assemblR(r1='reads1.fastq',r2='reads2.fastq',mode='accurate',out='contigs.fasta',report='assembly_report.txt')

assemblR<-function(r1,r2,mode,out,report){

	require(seqinr)

	spades<-system.file('spades/bin', 'spades.py', package = 'KARL')

	'%notin%'<-Negate('%in%')

	e1<-exists('r1')
	e2<-exists('r2')

	if (mode=='fast'){

		if (e1==TRUE & e2==TRUE){

			print('Assembling PAIRED-END reads in FAST mode...')

			spades<-paste(spades,' -1 ',r1,' -2 ',r2,' -o tmp',sep='')

			system(spades,ignore.stdout=TRUE)

			system(paste('mkdir ',out,sep=''))
			system(paste('mv tmp/*.log ',out,sep=''),ignore.stdout=TRUE)
			system(paste('mv tmp/contigs.fasta ',out,sep=''),ignore.stdout=TRUE)

			system('rm -rf tmp',ignore.stdout=TRUE)

		} else if (e1==TRUE & e2==FALSE){

			print('Assembling SINGLE-END reads in FAST mode...')

			spades<-paste(spades,' -1 ',r1,' -o ',out,sep='')

			system(spades,ignore.stdout=TRUE)

			system(paste('mkdir ',out,sep=''))
			system(paste('mv tmp/*.log ',out,sep=''),ignore.stdout=TRUE)
			system(paste('mv tmp/contigs.fasta ',out,sep=''),ignore.stdout=TRUE)

			system('rm -rf tmp',ignore.stdout=TRUE)

		} else if (e1==FALSE & e2==FALSE){

			stop('Error in reads file names. Please set at least the r1 parameter.')

		}

	} else if (mode=='accurate'){

		if (e1==TRUE & e2==TRUE){

			print('Assembling PAIRED-END reads in ACCURATE mode...')

			spades<-paste(spades,' --careful -k 21,27,33,39,45,51,57,63,69,75,81,87,93,99 -1 ',r1,' -2 ',r2,' -o tmp',sep='')

			system(spades,ignore.stdout=TRUE)

			system(paste('mkdir ',out,sep=''))
			system(paste('mv tmp/*.log ',out,sep=''),ignore.stdout=TRUE)
			system(paste('mv tmp/contigs.fasta ',out,sep=''),ignore.stdout=TRUE)

			system('rm -rf tmp',ignore.stdout=TRUE)

		} else if (e1==TRUE & e2==FALSE){

			print('Assembling SINGLE-END reads in FAST mode...')

			spades<-paste(spades,' --careful -k 21,27,33,39,45,51,57,63,69,75,81,87,93,99 -1 ',r1,' -o tmp',sep='')

			system(spades,ignore.stdout=TRUE)
			system(spades,ignore.stdout=TRUE)

			system(paste('mkdir ',out,sep=''))
			system(paste('mv tmp/*.log ',out,sep=''),ignore.stdout=TRUE)
			system(paste('mv tmp/contigs.fasta ',out,sep=''),ignore.stdout=TRUE)

			system('rm -rf tmp',ignore.stdout=TRUE)

		} else if (e1==FALSE & e2==FALSE){

			stop('Error in reads file names. Please set at least the r1 parameter.')
			
		}

	} else if (mode%notin%c('fast','accurate')){

		stop('Parameter mode is wrong or not specified: please provide fast or accurate')
	}

	contigs<-read.fasta(paste(out,'/contigs.fasta',sep=''))

	seqs<-getSequence(contigs)

	contig_num<-length(seqs)

	lens<-unlist(lapply(seqs,length))

	contig_min<-min(lens)
	contig_max<-max(lens)
	contig_avg<-round(mean(lens))

	lens_max<-round(length(which(lens>=1000))/contig_num*100)
	lens_min<-round(length(which(lens<1000))/contig_num*100)

	contig_sum<-sum(lens)

	cons<-rev(sort(lens))

	n50<-cons[which(cumsum(cons)>=sum(cons)/2)[1]]

	rpt<-paste(out,'/',report,sep='')

	cat('#### REPORT ####',file=rpt,sep='\n')

	cat('\n',file=rpt,sep='\n',append=TRUE)

	cat(paste('Genome length: ',contig_sum,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Number of contigs: ',contig_num,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Percentage of contigs >= 1000 bp: ',lens_max,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Percentage of contigs < 1000 bp: ',lens_min,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Min contig length: ',contig_min,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Max contig length: ',contig_max,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('Average contig length: ',contig_avg,sep=''),file=rpt,sep='\n',append=TRUE)

	cat(paste('N50: ',n50,sep=''),file=rpt,sep='\n',append=TRUE)

	print('ASSEMBY HAS FINISHED')
}
