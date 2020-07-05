dir.create('nulldistr')

run_distr<-function(){
for(i in 1:1000){
	
	source('largeModel_code_par.r')

	iter = gsub(' |:','',timestamp(prefix='',suffix=''))

	dir.create(paste0('nulldistr/iter_',iter))

	system(paste0('mv data_out/Gene_Choice_aggstats.rda nulldistr/iter_',iter,'/'))
	system(paste0('mv data_out/Gene_Choice.rda nulldistr/iter_',iter,'/'))

	# clean up
	#system('rm *_secretor.*')

	# #### save Gene Linkage Scores
	# system(paste0('mv gene_linkage_score.* nulldistr/iter_',iter,'/'))

	# #### save Best Gene / Order
	# system(paste0('mv best_gene.* nulldistr/iter_',iter,'/'))
	# system(paste0('mv best_gene_order.* nulldistr/iter_',iter,'/'))

	# #### save Model Score
	# system(paste0('mv Model_Scores.desc nulldistr/iter_',iter,'/'))
	# system(paste0('mv Model_Scores.bin nulldistr/iter_',iter,'/'))

	# #### save Model Score Groups
	# system(paste0('mv Model_Score_Groups.desc nulldistr/iter_',iter,'/'))
	# system(paste0('mv Model_Score_Groups.bin nulldistr/iter_',iter,'/'))

	# #### save Gene Choice
	# system(paste0('mv data_out/Gene_Choice.* nulldistr/iter_',iter,'/'))
	# system(paste0('mv data_out/Structure_Choice.* nulldistr/iter_',iter,'/'))

}
}

#load_null_distr=function(){
	library(abind)
	iters = system('ls nulldistr/',intern=T)
	distrL=lapply(iters,function(i){
		load(paste0('nulldistr/',i,'/Gene_Choice_aggstats.rda'))
		out=data_allr[,c('PROP','GLS','MSC')]
		nam=data_allr[,c("dataset","Linkage","Link","Gene")]
		rownames(out) = apply(nam,1,function(x) paste(x,collapse='_'))
		out
	})
	distr=abind(distrL,along=3)
	save(distr,file='nulldistr.rda')
#	return(distr)
#}

# par(mfrows=c(3,1))
# hist(tmp[,1,],main='Prop')
# hist(tmp[,2,],main='Gene-Linkage Score')
# hist(tmp[,3,],main='Model Score Contribution')
