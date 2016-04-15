#' initializeKARL
#'
#' This function generates and stores models and databases for prediction.
#' @param data default is 'alldata' object.
#' @param rank takes 'domain', 'phylum', 'class', 'order', 'family', 'genus' or 'all' (default: 'all').
#' @param min takes the minimum number of genomes per taxon to build a model (default: 2).
#' @keywords models prediction
#' @export
#' @examples
#' initializeKARL(rank='all',min=2)

initializeKARL<-function(

	             rank='all',
		     silva.url='http://www.arb-silva.de/fileadmin/silva_databases/current/Exports/SILVA_123_SSURef_Nr99_tax_silva.fasta.gz',
		     min=2

                       ){

	cdir<-getwd()

	print('### Initialization started... ###')
	print('### Initializing models...    ###')

	'%notin%'<-Negate('%in%')

	ranks<-c('domain','phylum','class','order','family','genus','all')

	if (rank%notin%ranks){

		stop("*** ERROR *** Parameter 'rank' has been wrongly set")

	} else if (rank=='all'){

		for (r in ranks){

			print(paste('--> Initializing ',r,' models!',sep=''))
			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))

			tdir<-paste(mdir,'/',r,sep='')
			system(paste('mkdir',tdir))
			setwd(tdir)
	
			taxa<-names(which(table(alldata[,r])>=min))

			for (x in taxa){

				bld<-buildSVM(rank=r,tax=x,prefix=x)
			}
		}
	
	} else if (rank=='domain'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))
			
			tdir<-paste(mdir,'/domain',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing domain models!')
	
			taxa<-names(which(table(alldata[,'domain'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='domain',tax=x,prefix=x)
			}

			setwd(cdir)
	
	} else if (rank=='phylum'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))

			tdir<-paste(mdir,'/phylum',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing phylum models!')
	
			taxa<-names(which(table(alldata[,'phylum'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='phylum',tax=x,prefix=x)
			}			
	
			setwd(cdir)

	} else if (rank=='class'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))	
		
			tdir<-paste('mkdir ',mdir,'/class',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing class models!')

			taxa<-names(which(table(alldata[,'class'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='class',tax=x,prefix=x)
			}

			setwd(cdir)

	} else if (rank=='order'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))

			tdir<-paste('mkdir ',mdir,'/order',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing order models!')
	
			taxa<-names(which(table(alldata[,'order'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='order',tax=x,prefix=x)
			}

			setwd(cdir)

	} else if (rank=='family'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))		
	
			tdir<-paste('mkdir ',mdir,'/family',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing family models!')
	
			taxa<-names(which(table(alldata[,'family'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='family',tax=x,prefix=x)
			}

			setwd(cdir)

	} else if (rank=='genus'){

			mdir<-gsub('blast','models',system.file('blast',package='KARL'))
			system(paste('mkdir',mdir))
			
			tdir<-paste('mkdir ',mdir,'/genus',sep='')
			system(paste('mkdir',tdir))

			setwd(tdir)

			print('--> Initializing genus models!')
	
			taxa<-names(which(table(alldata[,'genus'])>=min))

			for (x in taxa){

				bld<-buildSVM(rank='genus',tax=x,prefix=x)
			}

			setwd(cdir)
	}
	print('### Initializing databases... ###')

	db1<-system.file('databases1',package='KARL')
	db2<-system.file('databases2',package='KARL')

	cmddb<-paste('mv ',db2,'/* ',db1,sep='')
	system(cmddb)
	
	system(paste('rm -rf',db2))

	dbfinal<-gsub('databases1','databases',db1)

	cmdfinal<-paste('mv ',db1,' ',dbfinal,sep='')
	system(cmdfinal)

	makeblastdb<-system.file('blast/makeblastdb',package='KARL')

	setwd(dbfinal)

	faas<-list.files(pattern='.faa')

	for (f in faas){

		tit<-gsub('.faa','',f)
		cmd<-paste(makeblastdb,' -in ',f,' -dbtype prot -title ',tit,' -out ',tit,sep='')
		system(cmd,ignore.stderr=T,ignore.stdout=T)
	}

	silva<-gsub('blast','silva',system.file('blast',package='KARL'))
	system(paste('mkdir',silva))
	setwd(silva)

	cmdsilva<-paste('curl -O',silva)
	system(cmdsilva,ignore.stderr=T,ignore.stdout=T)
	
	system('gunzip -d *.gz',ignore.stderr=T,ignore.stdout=T)
	system('grep ">" *.fasta > headers.txt',ignore.stderr=T,ignore.stdout=T)
	
	cmdsilva2<-paste(makeblastdb,' -in *.fasta -dbtype nucl -hash_index -parse_seqids -title silvadb -out silvadb',sep='')
	system(cmdsilva2,ignore.stderr=T,ignore.stdout=T)

	
	ddir<-system.file('data',package='KARL')
	setwd(ddir)

	prcompEC<-prcomp(alldata[,9:dim(alldata)[2]])
	save(prcompEC,file='prcompEC.RData')

	setwd(cdir)
	
	print('### Initialization finished!  ###')
}
