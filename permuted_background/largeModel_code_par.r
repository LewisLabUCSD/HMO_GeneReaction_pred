
if(F){
	#### run code
	cli=15
	setwd('../HMO_biosynthesis.RAND')
	source('calc_null_distr.r')
	load_null_distr()
}


source('../HMO_biosynthesis/functions.R')
source('../HMO_biosynthesis/unittests.r')

## parallelize


dbg = F
run_GLS=T
run_GLbin=T
run_LS=T
run_MS=T
run_groups=T
run_inclusion=F
run_check=F
vis=F

cooc = F
colink=F

run_SC=T
run_Struc_CBX=F
anyl_CBX=F
run_GC=T

anyl_GLSvMS=F
write_scores=F

cluster=F
network=F
CBX_extend=F

new_sels = F
gene_v_str=F

RAND=F

runPCA=F

inclusion_rank=F

#gene_spaggetti=F
quantile_norm=F

cluster_vis=F

skip_secretor=F

#cli = 15

rm_genes = c('POFUT1b','POFUT2a','POFUT2b','FUT8')



#if(gene_spaggetti){
	library(openxlsx)
	# gene expression spaggetti plots
	
	# split data1,data2,secretor,nonsecretor
	d1_orig = read.csv('../HMO_biosynthesis/data/raw/data_HMO.data_used_gene1.csv') #sheetName='Data_used_gene1')
	d1_orig = d1_orig[,!colnames(d1_orig) %in% c("X","Mean","organism","syn","ilmn2")]

	d2_orig = read.csv('../HMO_biosynthesis/data/raw/data_HMO.data_used_gene2.csv') #sheetName='Data_used_gene1')

	if(quantile_norm){ # additional normalziation (skipped)
		library(preprocessCore)
		d1 = cbind(d1_orig[,1:2],normalize.quantiles(as.matrix(d1_orig[,-(1:2)]))); colnames(d1) = colnames(d1_orig)
		d2 = cbind(d2_orig[,1:2],normalize.quantiles(as.matrix(d2_orig[,-(1:2)]))); colnames(d2) = colnames(d2_orig)
	}else{
		d1=d1_orig; d2=d2_orig
	}
	gc = read.csv('../HMO_biosynthesis/data/raw/data_HMO.gene_correspondence.csv')[,1:4] #sheetName='Data_used_gene1')
	gc[gc=='']=NA
	gc=na.omit(gc)
	tmp = rbind( cbind( melt(d1,id.vars=c('ilmn','gene')) , ds='d1') , cbind( melt(d2,id.vars=c('ilmn','gene')) , ds='d2') )
	tmp$expression = as.numeric(as.character(tmp$value))
	count = rowSums(table(tmp$gene,tmp$ds))
	tmp = tmp[as.character(tmp$gene) %in% names(count)[count>107],]
	g=ggplot(data=tmp,aes(x=gene,y=expression+1,color=ds))+geom_boxplot()+scale_y_log10()+
		theme(axis.text.x = element_text(angle = 45, hjust = 1,size=)) 
	#ggsave(plot=g,filename='figures/gene_expression.pdf')
	tmp2= as.data.frame( do.call(rbind,strsplit(as.character(tmp$variable),'\\.'))[,1:2] )
	colnames(tmp2) = c('subject','time_days')
	tmp2[,2] = ifelse( grepl('h',tmp2[,2]) , as.numeric(gsub('h','',tmp2[,2]))/24 , as.numeric(gsub('D','',tmp2[,2])) )
	tmp = cbind( tmp, tmp2)
	g2=grid.arrange(grobs=lapply(c('d1','d2'),function(ds){ 
		ggplot(data=tmp[tmp$ds==ds,],aes(x=time_days,y=expression+1,color=subject))+geom_point()+geom_line()+scale_y_log10()+
			theme(axis.text.x = element_text(angle = 45, hjust = 1,size=)) + facet_wrap(~gene)
	}))
	#ggsave(plot=g2,filename='figures/gene_expression.bytime.pdf',height=10,width=15)
	g3=grid.arrange(grobs=lapply(c('d1','d2'),function(ds){ 
		ggplot(data=tmp[tmp$ds==ds,],aes(x=gene,y=expression+1,color=subject))+geom_boxplot()+scale_y_log10()+
			theme(axis.text.x = element_text(angle = 45, hjust = 1,size=)) 
	}))
	#ggsave(plot=g3,filename='figures/gene_expression.bysubject.pdf',height=10,width=15)
	rm_genes = unique( c(rm_genes, 
		as.character(d1[apply( d1[,-(1:2)] , 1 , quantile)[4,]<=0,'gene']) ,
		as.character(d2[apply( d2[,-(1:2)] , 1 , quantile)[4,]<=0,'gene']) ) )
#}


stat = c('R','P')[1]
DS = c('data1','data2','data1_secretor','data2_secretor')

# load scores codes
combinations = do.call(rbind,strsplit(gsub('_sec','-sec',system('ls ../HMO_biosynthesis/data/scores_split/*.mat',intern=T)),'_|\\.'))
unq_tmp=unique(paste(combinations[,7],ifelse(combinations[,8]=='mat',NA,combinations[,8]),sep='_'))
tmp=do.call(rbind,strsplit(unq_tmp,'_'))

# load model codes
combinations2 = do.call(rbind,strsplit(gsub('_sec','-sec',system('ls ../HMO_biosynthesis/data/models_split/*.mat',intern=T)),'_|\\.'))
unq_tmp2=unique(paste(combinations2[,6],ifelse(combinations2[,7]=='mat',NA,combinations2[,7]),sep='_'))
tmp2=do.call(rbind,strsplit(unq_tmp,'_'))

# check common and duplicated
rm_fail = c( unq_tmp[duplicated(tmp[,1]) & tmp[,2]=='NA'] , union(unq_tmp,unq_tmp2)[!union(unq_tmp,unq_tmp2) %in% intersect(unq_tmp,unq_tmp2)])
print( paste('removed:' , paste( rm_fail , collapse=',')))
combinations_unq <<- gsub('_NA','',unq_tmp[!unq_tmp %in% rm_fail])

if(dbg){ combinations_unq <<- combinations_unq[1:10] }

### Load Model-Correlation Score & Compute Gene-Linkage-Model Score
linkage_orig = load_linkages()

### 8/1/19, added ST3GAL1 bc it is gal-galnac specific, and B3GALNT2 bc it is a glcnac transferase (not relevant)
linkage = linkage_update(linkage_orig)

#linkage = lapply(names(linkage_orig),function(xn){ 
#	x = linkage_orig[[xn]]
#	#if(grepl('L3_|L9',xn)){
#	if(grepl('L3_',xn)){
#	 	rm_genes_i = c(rm_genes,"FUT2")
#	}else if(grepl('L4_',xn)){
#		rm_genes_i = c(rm_genes,"ST3GAL1a")
#	}else if(grepl('L6_',xn)){
#		rm_genes_i = c(rm_genes,"B3GALNT2")
#	}else if(grepl('L9_',xn)){
#		rm_genes_i = c(rm_genes,"FUT2",'FUT4','FUT11','FUT6b')
#	}else{
#		rm_genes_i = rm_genes
#	}
#	x[!x%in%rm_genes_i]
#} )
names(linkage) = names(linkage_orig)
gl_names_orig = gsub('_gene[1-9]|_gene','',paste(names(unlist(linkage_orig)),unlist(linkage_orig),sep='_'))
gl_names = gsub('_gene[1-9]|_gene','',paste(names(unlist(linkage)),unlist(linkage),sep='_'))

# functions
which.max.na = function(x) ifelse( all(is.na(x)),NA,which.max(x) )
max.na = function(x) ifelse( all(is.na(x)),NA,max(x, na.rm=T) )
mean.na = function(x) ifelse( all(is.na(x)),NA,mean(x, na.rm=T) )

save(combinations_unq,stat,gl_names,linkage,which.max.na,max.na,mean.na,DS,linkage_orig,gl_names_orig,file='combinations_unq.rda')

##################
# Scoring
ptm <- proc.time()


if(run_GLS){
print('calculate gene-linkage score: run_GLS')

cl<-makeCluster(spec = cli)
registerDoParallel(cl = cl)

## calculate mean
Gene_Linkage_Scores = lapply(DS,function(dataset){
	# calc N
	#breaks = do.call(cbind,lapply( combinations_unq , function(cbx){
	#breaks = unlist(mclapply( combinations_unq , function(cbx){
	breaks<-unlist(foreach(cbx = combinations_unq,.errorhandling='pass',.packages='R.matlab') %dopar% {
		load('combinations_unq.rda')
		dim( readMat(paste0('../HMO_biosynthesis/data/scores_split/',stat,'_',dataset,'_',cbx,'.mat'))[[1]] )[2]
#		as.double(unlist(  lapply(r,function(x) dim(x)[2])) )
	})#))
	#model_count =  unique( apply(breaks,1,sum) )
	model_count =  sum(breaks)
	break_idx = c(0,cumsum(breaks))
	# initialize gene-linkage correlation matrix
	gene_linkage_raw = big.matrix(nrow=model_count,ncol= length(unlist(linkage)),type='double',
		separated = FALSE, backingfile = paste0("gene_linkage_raw.",dataset,".bin"), descriptorfile = paste0("gene_linkage_raw.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(gene_linkage_raw) = gl_names
	options(bigmemory.allow.dimnames=FALSE)

	gene_linkage_shuffle = big.matrix(nrow=model_count,ncol= length(unlist(linkage)),type='double',
		separated = FALSE, backingfile = paste0("gene_linkage_shuffle.",dataset,".bin"), descriptorfile = paste0("gene_linkage_shuffle.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(gene_linkage_shuffle) = gl_names
	options(bigmemory.allow.dimnames=FALSE)

	gene_linkage_score = big.matrix(nrow=model_count,ncol= length(unlist(linkage)),type='double',
		separated = FALSE, backingfile = paste0("gene_linkage_score.",dataset,".bin"), descriptorfile = paste0("gene_linkage_score.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(gene_linkage_raw) = gl_names
	options(bigmemory.allow.dimnames=FALSE)
#	desc = describe(gene_linkage_raw)

#	if(dataset==DS[1]){
#		inclusion = big.matrix(ncol=model_count,nrow= 221,type='integer',
#			separated = FALSE, backingfile = "inclusion.bin", descriptorfile = "inclusion.desc")
#	}

	#i=1
	# populate gene-linkage correlation matrix
#res<-foreach(i = 1:n,.combine = rbind) %dopar% {
#require(bigmemory)
#data<-attach.big.matrix(datadesc)
#	for(ci in 1:length(combinations_unq)){
	if(grepl('secretor',dataset)&skip_secretor){
		return(describe(gene_linkage_score))
	}

	# perfect minimal hash Fredman, Komlos & Szemeredi 1948; allows a perfect mapping to a limited domain
	hash = function(x,n,p,k){ if(n<p){ (k*x %% p) %% n } else { stop('n must be less than p') } }
	p=2038074743
	seed = sample.int(2^10,1)

	res<-foreach(ci = 1:length(combinations_unq),.packages=c('R.matlab','bigmemory')) %dopar% {
		require(bigmemory)
		load('combinations_unq.rda')
		cbx = combinations_unq[ci]
		print(paste('populating',cbx))
		r=readMat(paste0('../HMO_biosynthesis/data/scores_split/',stat,'_',dataset,'_',cbx,'.mat'))
		r=r[c(2:10,1)]
		if( any( sapply(r,nrow) != sapply(linkage_orig,length) )){stop('linkage names dont make raw linkage scores')}
		imp=readMat(paste0('../HMO_biosynthesis/data/scores_split/model_impossible_',cbx,'.mat'))[[1]]
#		j=(i+dim(r[[1]])[2]-1)
		i = break_idx[ci]+1
		j = break_idx[ci+1]
		tmp = do.call(cbind,lapply(r,function(x) t(x) ) )
		desc = paste0("gene_linkage_raw.",dataset,".desc")
		data<-attach.big.matrix(desc)
		# desc2 = paste0("gene_linkage_shuffle.",dataset,".desc")
		# data2<-attach.big.matrix(desc2)
		#data<-gene_linkage_raw
		for(k in 1:ncol(data)){
			k_raw = which( gl_names[k] == gl_names_orig )
			tmp[imp,k_raw] = NA
			idx = (i:j)+nrow(data)*(k-1)
			idx_hash = hash(idx,n=prod(dim(data)),p=p,k=seed)+1
			data[idx_hash] = tmp[,k_raw]
			#data[i:j,k] = tmp[,k_raw]
			# data2[i:j,k] = tmp[,k_raw]
		}

#		if(dataset==DS[1]){
#			solCombi=readMat(paste0('data/models_split/solCombi_',cbx,'.mat'))$solCombi
#			incl = attach.big.matrix('inclusion.desc')
#			for(j in 1:nrow(incl)){
#				incl[j,i:j]=incl[j,]
#			}
#		}
		NULL
#		i=j+1
	}

	########################## RANDOMIZE: gene_linkage_raw matrix
	# 	desc = paste0("gene_linkage_raw.",dataset,".desc")
	# 	gene_linkage_raw<-attach.big.matrix(desc)
	# 	desc2 = paste0("gene_linkage_shuffle.",dataset,".desc")
	# 	gene_linkage_shuffle<-attach.big.matrix(desc2)
	# n=prod(dim(gene_linkage_shuffle))
	# idx_rand = sample.int(n,size=n/2,replace=F,useHash=T)
	# step=1e4
	# idx_max = prod(dim(gene_linkage_raw))
	# idx_rand = 1:idx_max
	# idx_seq = seq(1,idx_max,step)
	# low = 1
	# for(high in idx_seq[-1]){
	# 	idx_rand_i = sample(idx_rand,step)
	# 	gene_linkage_shuffle[low:(high-1)] = gene_linkage_raw[idx_rand_i]
	# 	low=high
	# 	idx_rand = idx_rand[-which(idx_rand%in%idx_rand_i)]
	# 	gc(reset=T)
	# }
	# gc(reset=T)
	####################################
	###########################################

	# calc sd
	linkage_sd = sapply( names(linkage) , function(l) sd(na.omit(as.vector( gene_linkage_raw[,grepl(gsub('_gene','',l),gl_names)] ))) )
	# calc mean
	linkage_mean = sapply( names(linkage) , function(l) mean(na.omit(as.vector( gene_linkage_raw[,grepl(gsub('_gene','',l),gl_names)] ))) )
	save(linkage_sd,linkage_mean,file='stats.rda')

	# normalized gene linkage score
	jnk<-foreach(i = 1:ncol(gene_linkage_raw),.packages=c('bigmemory')) %dopar% {
#	for(i in 1:ncol(gene_linkage_raw)){
		require(bigmemory)
		load('combinations_unq.rda')
		load('stats.rda')
		data_raw   = attach.big.matrix(paste0("gene_linkage_raw.",dataset,".desc"))
		data_score = attach.big.matrix(paste0("gene_linkage_score.",dataset,".desc"))
		l = paste( paste(strsplit( gl_names[i] , '_')[[1]][1:2],collapse='_'),'gene',sep='_') 
		data_score[,i] = (data_raw[,i]-linkage_mean[[l]]) / linkage_sd[[l]]
		NULL
	}

	return(describe(gene_linkage_score))
})
names(Gene_Linkage_Scores) = DS
Gene_Linkage_Scores = lapply( Gene_Linkage_Scores, function(d){	
	if(is.big.matrix(d)){d=describe(d)}
	x = attach.big.matrix(d)
	print(dim(x))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(x) = gl_names
	options(bigmemory.allow.dimnames=FALSE)
	return(describe(x))
})
stopCluster(cl)
}else{
	Gene_Linkage_Scores = lapply( DS,function(x) dget( paste0("gene_linkage_score.",x,".desc") ) ) 
}

#stop()

##############

if(run_GLbin){
print('calculate linkage genes: run_GLbin')

cl<-makeCluster(spec = 5)
registerDoParallel(cl = cl)

### get inclusion for best genes
Best_Genes = lapply(Gene_Linkage_Scores,function(glsM){ 
	desc = glsM

	dataset=strsplit(desc@description$filename,split='\\.')[[1]][2]
	print(dataset)
	gc(reset=T)

	# initi linkage score mx
	best_genes = big.matrix(nrow=desc@description$nrow,ncol= length(gl_names),type='integer',
		separated = FALSE, backingfile = paste0("best_gene.",dataset,".bin"), descriptorfile = paste0("best_gene.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(best_genes) = gl_names # gsub('_gene','',names(linkage))
	options(bigmemory.allow.dimnames=FALSE)
	desc_genes = describe(best_genes)

	best_genes_order = big.matrix(nrow=desc@description$nrow,ncol= length(gl_names),type='integer',
		separated = FALSE, backingfile = paste0("best_gene_order.",dataset,".bin"), descriptorfile = paste0("best_gene_order.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(best_genes) = gl_names # gsub('_gene','',names(linkage))
	options(bigmemory.allow.dimnames=FALSE)
	desc_genes_order = describe(best_genes_order)

	if(grepl('secretor',dataset)&skip_secretor){
		return(desc_genes)
	}

	# get max gene
	tmp<- foreach(l = names(linkage),.packages=c('R.matlab','bigmemory')) %dopar% {
		require(bigmemory)
		gc(reset=T)
		load('combinations_unq.rda')
		gls = attach.big.matrix(desc)
		l_idx = grepl(gsub('_gene','',l),gl_names)
		mx_gene= attach.big.matrix(desc_genes )
		mx_gene_order= attach.big.matrix(desc_genes_order )
		print(l)
		print(gl_names[which(l_idx)])
		na_mods = c()
		if(sum(l_idx)==1){ 
			mx_gene[,l_idx]= rep(1,nrow(mx_gene)) 
			mx_gene_order[,l_idx]= rep(1,nrow(mx_gene_order)) 
		}else{
			tmpi = which(l_idx)[ unlist( apply( gls[,l_idx] , 1 , which.max.na ) ) ]
			tmpj =  apply( gls[,l_idx] , 1 , order , decreasing=T ) 
			#print(1)
			#mx_gene[cbind(1:length(tmpi),tmpi)] = integer(1)
			#print(1)
			na_mods = is.na(tmpi)
			#print(1)
			count=1
			for(i in which(l_idx)){
				print(i)
				mx_gene[na_mods,i] = NA
				mx_gene[!na_mods & tmpi==i,i]=1
				mx_gene_order[na_mods,i] = NA
				mx_gene_order[,i] = tmpj[count,]
				count=count+1
			}
			
		}
		NULL
	}#) 
	#colnames(mx_gene) = gsub('_gene','',names(linkage))
#	keep=sapply( mx_gene,function(x) length(levels(x))>1)
	#mx_gene[is.na(mx_gene)] = 'NA'
	
	return(desc_genes)
#	tmp=Matrix( data.matrix( dummy.data.frame( mx_gene,drop=F,omit.constants=F) ) , sparse=T)
	#tmp = do.call(cbind,apply( mx_gene , 2 , function(mx) { dummy( mx,drop=F) } ))
	#options(na.action="na.pass")
	#tmp=model.matrix(as.formula(paste('~',paste(colnames(mx_gene),collapse='+'),'+0')),data=mx_gene)
})

print('Gene_Best computed')
print(proc.time()-ptm)
stopCluster(cl)
}else{
	Best_Genes = lapply( DS,function(x) dget( paste0("best_gene.",x,".desc") ) ) 
}

##############


if(run_LS){
print('calculate linkage score: run_LS')
# calc gene-linkage score matrix
ptm = proc.time()
cl<-makeCluster(spec = cli)
registerDoParallel(cl = cl)
Linkage_Scores =  lapply(Gene_Linkage_Scores,function(glsM){ 
	desc = glsM
	dataset=strsplit(desc@description$filename,split='\\.')[[1]][2]
	print(dataset)
	gc(reset=T)
	
	# initi linkage score mx
	link_score = big.matrix(nrow=desc@description$nrow,ncol= length(linkage),type='double',
		separated = FALSE, backingfile = paste0("linkage_score.",dataset,".bin"), descriptorfile = paste0("linkage_score.",dataset,".desc"))
	options(bigmemory.allow.dimnames=TRUE)
	colnames(link_score) = names(linkage)
	options(bigmemory.allow.dimnames=FALSE)
	desc_link = describe(link_score)


	if(grepl('secretor',dataset)&skip_secretor){
		return(describe(desc_genes))
	}

	#score = do.call( cbind , mclapply( names(linkage) , function(l) {
	tmp<-foreach(l = names(linkage),.packages=c('bigmemory')) %dopar% {
		require(bigmemory)
		load('combinations_unq.rda')
		gls = attach.big.matrix(desc)
		link = attach.big.matrix(desc_link)
		l_idx = grepl(gsub('_gene','',l),gl_names)
		#if( !sum(l_idx) %in% c('1','3','4','5','7','9')){stop(paste('incorrect index in',l))}

		if(sum(l_idx)==1){
			link[,names(linkage)==l] = gls[,l_idx]
		}else{
#n=1e4
#ptm=proc.time(); tmp=apply( attach.big.matrix(Gene_Linkage_Scores[[1]])[,l_idx] , 1 , max.na  ) ; print(proc.time()-ptm)
#ptm=proc.time(); tmp=biganalytics::apply( as.big.matrix( attach.big.matrix(Gene_Linkage_Scores[[1]])[1:n,l_idx] ) , 1 , max.na  ) ; print(proc.time()-ptm)
#ptm=proc.time(); tmp=biganalytics::apply( attach.big.matrix(Gene_Linkage_Scores[[1]])  , 1 , function(x) max.na(x[l_idx])  ) ; print(proc.time()-ptm)
#			link[,names(linkage)==l] = biganalytics::apply( as.big.matrix(gls[,l_idx]) , 1 , max.na  )
			link[,names(linkage)==l] = apply( gls[,l_idx] , 1 , max.na  ) 
		}
		NULL
	} #))
	#colnames(score) = names(linkage)
	gc(reset=T)
	describe( link_score )
}) 
print('Linkage_Scores computed')
print(proc.time()-ptm)

stopCluster(cl)
}else{
	Linkage_Scores = lapply( DS,function(x) dget( paste0("linkage_score.",x,".desc") ) ) 
}

############################3

if(run_MS){
print('calculate model score: run_MS')
#### compute model scores
cl<-makeCluster(spec = 5)
registerDoParallel(cl = cl)

Model_Scores = big.matrix(nrow=Linkage_Scores[[1]]@description$nrow,ncol= length(DS),type='double',
	separated = FALSE, backingfile = 'Model_Scores.bin', descriptorfile = 'Model_Scores.desc')
options(bigmemory.allow.dimnames=TRUE)
colnames(Model_Scores) = DS
options(bigmemory.allow.dimnames=FALSE)
tmp<-foreach(ds = DS,.packages=c('bigmemory')) %dopar% {
	load('combinations_unq.rda')
	link = attach.big.matrix( paste0("linkage_score.",ds,".desc") ) 
	#dataset=ds #strsplit(ls@description$filename,split='\\.')[[1]][2]
	model_scores = attach.big.matrix('Model_Scores.desc')
	model_scores[,DS==ds] = biganalytics::apply( link , 1 , mean.na  ) 
	NULL
}
print('Model_Scores computed')
print( proc.time()-ptm)
stopCluster(cl)
}else{
	Model_Scores = attach.big.matrix('Model_Scores.desc')
}

##############

if(run_groups){
print('calculate model groups: run_groups')

cl<-makeCluster(spec = 5)
registerDoParallel(cl = cl)

Model_Score_Groups = big.matrix(nrow=Linkage_Scores[[1]]@description$nrow,ncol= length(DS),type='integer',
	separated = FALSE, backingfile = 'Model_Score_Groups.bin', descriptorfile = 'Model_Score_Groups.desc')
options(bigmemory.allow.dimnames=TRUE)
colnames(Model_Scores) = DS
options(bigmemory.allow.dimnames=FALSE)

tmp<-foreach(ds = DS,.packages=c('bigmemory')) %dopar% {
	load('combinations_unq.rda')
	score = attach.big.matrix('Model_Scores.desc')[,DS==ds]
	MSG = attach.big.matrix('Model_Score_Groups.desc')
	z = scale(score)
	p = pnorm(z,lower.tail=F); 
	q = p.adjust(p,'fdr'); 
	lab = ifelse(p<.05,1,0); 
#	lab = ifelse(q<.1,1,0); 
	MSG[,DS==ds] = lab
	NULL
}
##
stopCluster(cl)
}else{
	Model_Score_Groups = attach.big.matrix('Model_Score_Groups.desc')
}

# if(vis){
# n=1e6
# if(dbg){
# 	jpeg('figures/sigmoid_models.sample.png',height=1000,width=1000)
# }else{
# 	png('figures/sigmoid_models.png',height=1000,width=1000)
# }
# par(mfrow=c(4,5))
# print('visualize')

# for(ds in DS){
# 	score = Model_Scores[,DS==ds]
# 	score2= Model_Scores[,DS== (ds2<-ifelse(grepl('1',ds),gsub('1','2',ds),gsub('2','1',ds))) ]
# 	lab = Model_Score_Groups[,DS==ds]
# 	lab2 = Model_Score_Groups[,DS== ds2 ]

# 	tmp = data.frame(score=score[idx<-sample(1:length(score),n)],lab=lab[idx])
# 	tmp = tmp[order(tmp$score,na.last=T),]
# 	tmp2 = data.frame(score2=score2[idx],lab=lab2[idx])
# 	tmp2 = tmp2[order(tmp2$score,na.last=T),]
# 	plot( 1:nrow(tmp) ,tmp$score,col=c('black','red')[factor(tmp$lab,levels=c('0','1'))] ,main=ds )
# 	plot( 1:nrow(tmp) ,tmp$score,col=c('black','red')[factor(tmp2$lab,levels=c('0','1'))] ,main='significant in validation set' )
# 	plot( 1:nrow(tmp) ,tmp$score,col=c('black','red')[factor(tmp$lab*tmp2$lab,levels=c('0','1'))] ,main='intersection with validation set' )
# 	#gg=qplot( x=tmp$score , y=tmp2$score , col=factor(tmp$lab+tmp2$lab) ) + geom_point() + stat_smooth(method='lm') + xlab(ds) + ylab(ds2)
# 	#ggsave(plot=gg,filename=paste0('figures/compare.',ds,'.png'))
# 	qqnorm(sample(score,n))
# 	hist(score)
# }

# dev.off()
# try(dev.off())
# }

# check scores
#ni=1

# if(run_check){
# 	di=sample(1:length(DS),1)
# 	tmp00=attach.big.matrix(paste0("gene_linkage_raw.",DS[di],".desc"))
# 	# get mean of sd
# 	mn<-lapply(names(linkage),function(l){ mean(as.vector(unlist(tmp00[,grepl(gsub('_gene','',l),gl_names)])),na.rm=T)  })
# 	sd<-lapply(names(linkage),function(l){ sd(as.vector(unlist(tmp00[,grepl(gsub('_gene','',l),gl_names)])),na.rm=T)  })
# 	names(mn)=names(sd)=names(linkage)

# check=do.call(rbind,lapply(sample(1:48e6,20),function(ni){
# 	print(ni)
# 	tmp0=tmp00[ni,]

# 	tmp1a=attach.big.matrix(Gene_Linkage_Scores[[di]])[ni,]
# 	print( all( tmp1a == (tmp1b<-unlist(lapply(names(linkage),function(l) ((tmp0[grepl(gsub('_gene','',l),gl_names)]-mn[[l]])/sd[[l]])) ) ) ))
# 	tmp2a=attach.big.matrix(Linkage_Scores[[di]])[ni,]
# 	print( tmp2a == (tmp2b<-unlist(lapply(names(linkage),function(l){ max.na(tmp1a[grepl(gsub('_gene','',l),gl_names)]) }))) )
# 	Model_Scores[ni,di] == mean(tmp2a)
# 	Model_Scores[ni,di] == mean(tmp2b)
# 	cbind(ni,tmp2a,tmp2b)
# }))
# if(!all(check[,2]==check[,3],na.rm=T)){stop('Model_Scores don"t validate')}
# }

####################
### structure choice
if(run_inclusion){
	print('load inclusion: run_inclusion')
	cl<-makeCluster(spec = cli)
	registerDoParallel(cl = cl)

	#inclusion =  do.call(cbind,mclapply( combinations_unq , function(cbx){
	#inclusion<-foreach(cbx = combinations_unq,.combine=cbind,.packages=c('Matrix','R.matlab')) %dopar% {
	inclusion_desc<-foreach(cbx = combinations_unq,.packages=c('R.matlab','bigmemory')) %dopar% {
		require(bigmemory)
		describe( as.big.matrix( readMat(paste0('../HMO_biosynthesis/data/models_split/solCombi_',cbx,'.mat'))$solCombi ,type='integer',
			separated = FALSE, backingfile = paste0('.inclusion_tmp',cbx,'.bin'), descriptorfile = paste0('.inclusion_tmp',cbx,'.desc')) )
	#		r = Matrix( readMat(paste0('data/models_split/solCombi_',cbx,'.mat'))$solCombi , sparse=T)
	}#))

	inclusion = big.matrix(ncol=sum(breaks<-unlist(lapply(inclusion_desc,function(x) x@description$ncol))),nrow=inclusion_desc[[1]]@description$nrow ,type='integer',
		separated = FALSE, backingfile = 'inclusion.rxn.bin', descriptorfile = 'inclusion.rxn.desc')

	options(bigmemory.allow.dimnames=TRUE)
	#colnames(inclusion) = annotate_inclusion()
	options(bigmemory.allow.dimnames=FALSE)

	break_idx = c(0,cumsum(breaks))
	save(break_idx,file='breaks_cum.rda')
	tmp<-foreach( ci = 1:length(combinations_unq) , .packages=c('bigmemory')) %dopar% {
	#for( ci in 1:length(combinations_unq)){ # , .packages=c('bigmemory')) %do% {
		require(bigmemory)
		load('combinations_unq.rda')
		cbx = combinations_unq[ci]
		print(paste('populating',cbx))
		mat = attach.big.matrix( paste0('.inclusion_tmp',cbx,'.desc') )
		inc_out = attach.big.matrix( paste0('inclusion.rxn.desc') )
		load('breaks_cum.rda')
		i = break_idx[ci]+1
		j = break_idx[ci+1]
		if( (j-i+1) != dim(mat)[2]){stop('wrong dimensions')}
		for(rxn in 1:nrow(mat)){
			inc_out[rxn,i:j] = mat[rxn,]
		}
		NULL
	}

	stopCluster(cl)
	#inclusion = attach.big.matrix('inclusion.desc')
}else{
	inclusion = attach.big.matrix( paste0('../HMO_biosynthesis/inclusion.rxn.desc') )
}

if(run_SC){
#rxn = grepl('-HMO',rownames(inclusion))

cl<-makeCluster(spec = cli)
registerDoParallel(cl = cl)

# background co-occurance
rxn_full =  annotate_inclusion()
rxn = !grepl('HMO',rxn_full)
lst=list()


## calculate relevant proportions
Structure_Choice = lapply(DS,function(dataset){
	# get selected models and scores
	# get selected models and scores
	model_score_i = Model_Scores[,DS==dataset]
	infeasible = is.na(Model_Scores[,DS==dataset])

	groups_i = Model_Score_Groups[,DS==dataset]
	selected_models = which( groups_i==1 & !infeasible)
	selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible)

	save(DS,dataset,MAD,invert_ds,file='tmp.rda')
	print(dataset)

	if(grepl('secretor',dataset)&skip_secretor){
		return(list(0))
	}

	print('calc proportions')
	out<<-data.frame(
		rxn_name = annotate_inclusion(),
		background_size_complete= dim(inclusion)[2], 
		background_size= dim(inclusion)[2]-sum(infeasible),
		selection_size = length( selected_models ),
		selection_inter_size = length( selected_models_inter ),
		background_success=foreach(idx=1:nrow(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,!infeasible] ,na.rm=T) },
		selection_success=foreach(idx=1:nrow(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,selected_models ] ,na.rm=T) },
		selection_inter_success=foreach(idx=1:nrow(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,selected_models_inter] ,na.rm=T) }
	)	
	out$background_failure=out$background_size-out$background_success
	out$selection_failure=out$selection_size-out$selection_success
	out$selection_inter_failure=out$selection_inter_size-out$selection_inter_success
	print(tail(out))

	print('enrichment')
	out$enrichment = 		phyper( out$selection_success-1 , 		m = out$background_success , n = out$background_failure , k = out$selection_size ,lower.tail=F)
	out$enrichment_inter = phyper( out$selection_inter_success-1 , m = out$background_success , n = out$background_failure , k = out$selection_inter_size ,lower.tail=F)
	out$p_background = out$background_success / out$background_size
	out$p_selection = out$selection_success / out$selection_size
	out$p_selection_inter = out$selection_inter_success / out$selection_inter_size
	print(tail(out))
	gc(reset=T)
	print('mean - background')
	#out$t_test_background = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_background = foreach(idx = 1:nrow(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory); load('tmp.rda'); i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,] ; model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset];
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
#		tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
	print(tail(out))
	gc(reset=T)
	print('mean - selected')
#	out$t_test_selected = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_selected = foreach(idx = 1:nrow(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory); load('tmp.rda'); i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,] ; model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1
		i = i[selected_models]; model_score_i = model_score_i[selected_models]
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
#		tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
	print(tail(out))
	gc(reset=T)
	print('mean - selected_inter')
#	out$t_test_selected = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_selected_inter = foreach(idx = 1:nrow(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory); load('tmp.rda'); i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx,] ; model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1 & (attach.big.matrix('Model_Score_Groups.desc')[,DS==invert_ds(dataset)])==1
		i = i[selected_models]; model_score_i = model_score_i[selected_models]
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
#		tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
#	print('wilcox test')
###	out$t_wilcox_background = do.call(rbind, apply(inclusion,1,function(i){
#	out$t_wilcox_background = foreach(idx = 1:nrow(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
#		require(bigmemory); load('tmp.rda'); i = attach.big.matrix('inclusion.rxn.desc')[idx,] ; model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
#		if( length(table( sample( i , 1e6)))<2 ){ if( length(table(i ))<2 ) {return(NA)}}
#		tmp=NA; try( tmp<-tidy( wilcox.test( model_score_i[i==1] , model_score_i[i==0]  ,na.action="na.omit") ) ); tmp
#	} #) )
###	out$t_wilcox_selected = do.call(rbind, apply(inclusion,1,function(i){
#	out$t_wilcox_selected = foreach(idx = 1:nrow(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
#		require(bigmemory); load('tmp.rda'); i = attach.big.matrix('inclusion.rxn.desc')[idx,] ; model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset];  
#		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1
#		if( length(table( sample( i[selected_models] , 1e6)))<2 ){ if( length(table(i[selected_models] ))<2 ) {return(NA)}}
#		tmp=NA; try( tmp<-tidy( wilcox.test( model_score_i[selected_models & i==1] , model_score_i[selected_models & i==0]  ,na.action="na.omit") ) ); tmp
#	} #) )
	print(tail(out))
#	out$glm_background =    do.call(rbind,apply(inclusion,1,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i ~ i  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))
#	out$glm_selected =      do.call(rbind,apply(inclusion,1,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i[selected_models] ~ i[selected_models]  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))

	# select co-occurance

	save(out,file=paste0('data_out/Structure_Choice.',dataset,'.rda'))
	out
})
stopCluster(cl)

if(skip_secretor){
	Structure_Choice[[3]] = Structure_Choice[[1]]
	Structure_Choice[[4]] = Structure_Choice[[2]]
}

names(Structure_Choice) = DS

Structure_Choice=lapply(Structure_Choice, function(sc){
	idx=which( unlist(lapply(sc,is.data.frame) ) | unlist(lapply(sc,is.matrix) ) )
	cbind( sc[,-idx], do.call(cbind, lapply(idx,function(i){ 
		tmp=sc[,i]
		colnames(tmp) = paste(colnames(sc)[i],colnames(tmp),sep='.')
		tmp
	} ) ) )
})

#lapply(DS,function(ds) write.csv(Structure_Choice[[ds]],file=paste0('data_out/Structure_Choice.',ds,'.csv')))
#library(xlsx)
try(write.xlsx(Structure_Choice,file=paste0('data_out/Structure_Choice.xlsx')))
save(Structure_Choice,file=paste0('data_out/Structure_Choice.rda'))

if(! 'rxn_name'%in%colnames(Structure_Choice[[1]])){
	Structure_Choice=lapply(Structure_Choice,function(sc){
	cbind(rxn_name=annotate_inclusion(),sc)
	})
}

if(run_Struc_CBX){
	cl<-makeCluster(spec = 5)
	registerDoParallel(cl = cl)

	inclusion=attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
	background_tab = table(apply(inclusion[!grepl('HMO',annotate_inclusion()),!is.na(Model_Scores[,1])],2,paste,collapse='_') )

	Structure_CBX = mclapply(DS,function(dataset){
		print(dataset)
		# get selected models and scores
		model_score_i = Model_Scores[,DS==dataset]
		infeasible = is.na(model_score_i)

		groups_i = Model_Score_Groups[,DS==dataset]
		selected_models = which( groups_i==1 & !infeasible)
		selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible) 
		save(DS,dataset,MAD,file='tmp.rda')
		print('calc cbx')
			rxn_names = annotate_inclusion()[!grepl('HMO',ann<-annotate_inclusion())]
			infeasible_cbx_freq=table(apply(inclusion[!grepl('HMO',ann),infeasible],2,paste,collapse='_') )
			selection_cbx_freq=table(apply(inclusion[!grepl('HMO',ann),selected_models ],2,paste,collapse='_'))
			selection_inter_cbx_freq=table(apply(inclusion[!grepl('HMO',ann),selected_models_inter ],2,paste,collapse='_'))
		print('merge cbx')
		out = merge(as.matrix(selection_cbx_freq),as.matrix(selection_inter_cbx_freq),by='row.names',all=T)
		out = merge(out,as.matrix(infeasible_cbx_freq),by.x='Row.names',by.y='row.names',all=T)
		colnames(out) =c( 'CBX' , paste(dataset,c('selection','intersection','infeasible'),sep='_'))

		save(out,file=paste0('data_out/Structure_Choice_cbx.',dataset,'.rda'))
		out
	})
	out = merge( Structure_CBX[[1]] , Structure_CBX[[2]] , by='CBX',all=T)
	out = merge( out , Structure_CBX[[3]] , by='CBX',all=T)
	out = merge( out , Structure_CBX[[4]] , by='CBX',all=T)
	bt = melt( background_tab)
	colnames(bt) = c('CBX','background')
	out = merge(out,bt,by='CBX',all=T)

	save(out,Structure_CBX,background_tab,file='data_out/Structure_Choice_cbx.rda')
	write.csv(out,file='data_out/Structure_Choice_cbx.csv')

	stopCluster(cl)
	}
}

if(CBX_extend){
	load('data_out/Structure_Choice_cbx.rda')
	rxn_names = annotate_inclusion()[!grepl('HMO',ann<-annotate_inclusion())]
	cbx_include = do.call(rbind,lapply( out$CBX , function(x) strsplit(x,'_')[[1]] ) )
	colnames(cbx_include) = gsub('_rxn','',rxn_names)

	tmp=list()
	for(i in c(1,100,1000,2000,4000)){
		tmp[[as.character(i)]]= melt( sapply( grep('_selection|data2_intersection|data2_secretor_intersection',colnames(out),value=T) , function(n) {
			apply( cbx_include[!is.na(out[,n]) & out[,n]>i,] , 2 , function(x) sum(as.numeric(x))/length(x) )
		}))
	}
	m=melt(tmp)
	colnames(m) =c('HMOi','DS','jnk','proportion','Min_Model_Count')
	m$HMO = gsub('[1-9]','',m$HMOi)
	m$subtype = ifelse( grepl('secretor',m$DS) , 'secretor' , 'all')
	m$ds = ifelse( grepl('1',m$DS) , 'data1' , 'data2')
	m$sel = ifelse( grepl('intersection',m$DS) , 'intersection' , 'selection')
	m$ds = ifelse( m$sel=='intersection', 'data1/2',m$ds)

	g=grid.arrange( grobs=lapply( unique(m$HMO) , function(h){ 
			ggplot( data=m[m$HMO==h & m$Min_Model_Count==1,],aes(x=HMOi,y=proportion,fill=paste(ds,sel))) + geom_bar(stat='identity',position='dodge') +
				facet_grid( ~ subtype , scale='free') + 
				theme(axis.text.x = element_text(angle = 90, hjust = 1))
		}) )
	ggsave(plot=g,filename='figures/proportions_by_combination.pdf')
}

#}else{
#}
#if(F){
	load('data_out/Structure_Choice.rda')
	load('data_out/Structure_Choice_cbx.rda')

	pdf('figures/Structure_Choice.pdf',height=15)
	par(mfrow=c(3,1),mar=c(10,8,2,2))
	lapply(Structure_Choice,function(sc){
		barplot( t( data.matrix( sc[idx<-!grepl('HMO',sc$rxn_name),c('p_background','p_selection','p_selection_inter')] ) ),names.arg=sc$rxn_name[idx],
			beside=T,las=2,ylab='proportion',legend.text=c('Background','Selection','Intersection'))

		barplot( t( data.matrix( sc[idx,grepl('.mean1',colnames(sc))]  - sc[idx,grepl('.mean0',colnames(sc))] )) ,
			names.arg=sc$rxn_name[idx],beside=T,las=2,ylab='Mean Model Score Improvement\n(with HMO - without HMO)',legend.text=c('Background','Selection','Intersection'))
	#	x=barplot( tmp<-t( data.matrix( sc[idx<-!grepl('HMO',sc$rxn_name),c('mean_MS_background.mean0','mean_MS_background.mean1')] ) ),names.arg=sc$rxn_name[idx],beside=T,las=2,ylab='Background Model Score Improvement\n(with HMO - without HMO)')
	#	error.bar(x,tmp, 1*(t( data.matrix( sc[idx<-!grepl('HMO',sc$rxn_name),c('mean_MS_background.SD0','mean_MS_background.SD1')] ) )))

		barplot( t( data.matrix( sc[idx,c('enrichment','enrichment_inter')] ) ),names.arg=sc$rxn_name[idx],beside=T,las=2,
			ylab='Hypergeometric Enrichment (p-value)',legend.text=c('Selection','Intersection'))
	})
	dev.off()

	# pretty
	Structure_Choice=lapply(DS,function(x){
		mean0=grep('.mean0',colnames(Structure_Choice[[x]]),value=T)
		mean1=grep('.mean1',colnames(Structure_Choice[[x]]),value=T)
		tmp = (Structure_Choice[[x]][,mean1] - Structure_Choice[[x]][,mean0])
		colnames(tmp) = gsub('mean1','mean_diff',colnames(tmp))
		cbind(Structure_Choice[[x]],tmp)
	})
	names(Structure_Choice) = DS

	datas=list(data=c('data1','data2','data'),secretor=c('data1_secretor','data2_secretor','data_secretor'))
	vars = list(prop = c('p_background','p_selection','p_selection_inter'),
		enrich=grep('enrichment',colnames(Structure_Choice[[1]]),value=T),
		means=grep('.mean_diff',colnames(Structure_Choice[[1]]),value=T))
	df = do.call(rbind, lapply(DS,function(x) cbind(ds=x,Structure_Choice[[x]])))[,c('ds','rxn_name',vars$prop,vars$enrich,vars$means)]
	dfm = melt(df,c('ds','rxn_name'))
	dfm$ds = ifelse(grepl('*_inter',dfm$variable)&!grepl('mean_diff',dfm$variable),gsub('1|2','',as.character(dfm$ds)),as.character(dfm$ds))
#	dfm$ds = ifelse(dfm$variable=='p_selection_inter',gsub('1|2','',as.character(dfm$ds)),as.character(dfm$ds))
	dfm$var = paste(dfm$ds,dfm$variable,sep='.')
	idx=!grepl('HMO',df$rxn_name)
	strs = c('^DFLNH','^DFLNT','^DSLNH','^FLNH','^FDSLNH')
	for(d_i in names(datas)){ for(v_i in names(vars)){
		g = grid.arrange(grobs=lapply(strs,function(x){
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[[v_i]],])
			if( v_i=='enrich' ){ data_tmp$value = -log( data_tmp$value+2e-16 , 10 ) }
			ggplot(data=data_tmp , aes( x=rxn_name,y=value,fill=var))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
				ggtitle(gsub('\\^','',x))
		}))
		ggsave(g,filename=paste0('figures/Structure_Choice/Structure_Choice.',d_i,'.',v_i,'.pdf'),height=10,width=10)
	}}


		grbs =lapply(strs,function(x){
			print(x)
			#back = mean(unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][1],])$value)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1],])
			data_tmp$i = gsub('[A-Z]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']],])
			data_tmp_e$sec=ifelse(grepl('secretor',data_tmp_e$ds),'secretor','all')
			data_tmp$enrichment = apply( data_tmp , 1 , function(i){
				if(grepl('inter',i[[3]])){
					mean( data_tmp_e$value[data_tmp_e$ds==i[[1]]&data_tmp_e$rxn_name==i[[2]]&grepl('inter',data_tmp_e$variable)&data_tmp_e$sec==i[[8]]] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==i[[1]]&data_tmp_e$rxn_name==i[[2]]&!grepl('inter',data_tmp_e$variable)&data_tmp_e$sec==i[[8]]] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.01,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+#geom_hline(yintercept=back)+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+facet_grid(~sec)+theme_classic()
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grbs[[1]],grbs[[2]],grbs[[3]],grbs[[4]],grbs[[5]],l)
		ggsave(g,filename=paste0('figures/Structure_Choice.pdf'),height=10,width=15)

		grbs =lapply(strs,function(x){
			print(x)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1] & !grepl('secretor',dfm$ds),])
			data_tmp$i = gsub('[A-Z]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']] & !grepl('secretor',dfm$ds),])
			#data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
			data_tmp$enrichment = apply( data_tmp , 1 , function(i){
				if(grepl('inter',i[3])){
					mean( data_tmp_e$value[data_tmp_e$ds==i[1]&data_tmp_e$rxn_name==i[2]&grepl('inter',data_tmp_e$variable)] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==i[1]&data_tmp_e$rxn_name==i[2]&!grepl('inter',data_tmp_e$variable)] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.01,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+theme_classic()
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grbs[[1]],grbs[[2]],grbs[[3]],grbs[[4]],grbs[[5]],l)
		ggsave(g,filename=paste0('figures/Structure_Choice.all.pdf'),height=10,width=15)
	#}


	## str CBX anyl
	lim=c(1,100,1000,2000,4000,6000)
	out$min.data_inter.secretor_inter = apply( out[,c('data1_intersection','data1_secretor_intersection')] , 1 , min )
	tmp0=do.call( cbind, lapply(lim,function(l) apply(out[,-1],2,function(x) sum( x>=l ,na.rm=T ) ) ))
	tmp1=do.call( cbind, lapply(lim,function(l) apply(out[,-1],2,function(x) sum( x>=l & is.na(out$data1_infeasible),na.rm=T ) ) ))
	tmp2=do.call( cbind, lapply(lim,function(l) apply(out[,-1],2,function(x) sum( x>=l & (is.na(out$data1_infeasible) | out$data1_infeasible<4000) ,na.rm=T) ) ))
	write.csv(cbind(tmp0,tmp1,tmp2),file='data_out/cbx_counts.csv')
	#write.csv(tmp2,file='data_out/cbx_counts.wo_infeasible.csv')
#}


########################
### gene linkage choice

if(run_GC){
lg = list()

cl<-makeCluster(spec = cli)
registerDoParallel(cl = cl)

## calculate relevant proportions
Gene_Choice = lapply(DS,function(dataset){
	# get selected models and scores
	model_score_i = Model_Scores[,DS==dataset]
	infeasible = is.na(model_score_i)

	groups_i = Model_Score_Groups[,DS==dataset]
	selected_models = which( groups_i==1 & !infeasible)
	selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible)

	# gene linkage score and gene inclusion
	ls = attach.big.matrix( Gene_Linkage_Scores[[which(DS==dataset)]] )
	inclusion = attach.big.matrix( Best_Genes[[which(DS==dataset)]] ) 
#	inclusion = t(inclusion)
	#inclusion = inclusion[!grepl('_NA',rownames(inclusion)),]
	save(DS,Best_Genes,Gene_Linkage_Scores,dataset,MAD,invert_ds,file='tmp.rda')

	if(grepl('secretor',dataset)&skip_secretor){
		return(list(0))
	}

	print('calc proportions')
	out<<-data.frame(
		rxn_name = gl_names,
		background_size_complete= dim(inclusion)[1], 
		background_size= dim(inclusion)[2]-sum(infeasible),
		selection_size = length( selected_models ),
		selection_inter_size = length( selected_models_inter ),
		background_success=foreach(idx=1:ncol(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix( paste0('best_gene.',dataset,'.desc') )[!infeasible,idx] ,na.rm=T) },
		selection_success=foreach(idx=1:ncol(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix( paste0('best_gene.',dataset,'.desc') )[selected_models,idx] ,na.rm=T) },
		selection_inter_success=foreach(idx=1:ncol(inclusion),.combine=c,.packages=c('bigmemory')) %dopar% {
			require(bigmemory); sum( attach.big.matrix( paste0('best_gene.',dataset,'.desc') )[selected_models_inter,idx] ,na.rm=T) }
	)	
	out$background_failure=out$background_size-out$background_success
	out$selection_failure=out$selection_size-out$selection_success
	out$selection_inter_failure=out$selection_inter_size-out$selection_inter_success
	print(tail(out))

	print('enrichment')
	out$enrichment = phyper( out$selection_success-1 , m = out$background_success , n = out$background_failure , k = out$selection_size ,lower.tail=F)
	out$enrichment_inter = phyper( out$selection_inter_success-1 , m = out$background_success , n = out$background_failure , k = out$selection_inter_size ,lower.tail=F)
	out$p_background = out$background_success / out$background_size
	out$p_selection = out$selection_success / out$selection_size
	out$p_selection_inter = out$selection_inter_success / out$selection_inter_size
	print(tail(out))

	gc(reset=T)
	print('mean -background')
	#out$t_test_background = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_background = foreach(idx = 1:ncol(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory);load('tmp.rda');  i = attach.big.matrix( paste0('best_gene.',dataset,'.desc') )[,idx] ; 
		model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset] 
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
		#tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
	print(tail(out))
	gc(reset=T)
	print('mean - selected')
#	out$t_test_selected = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_selected = foreach(idx = 1:ncol(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory);load('tmp.rda');  i = attach.big.matrix( paste0('best_gene.',dataset,'.desc')  )[,idx] ; 
		model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1
		i = i[selected_models]; model_score_i = model_score_i[selected_models]
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
#		tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
	print(tail(out))
	gc(reset=T)
	print('mean - selected_inter')
#	out$t_test_selected = do.call(rbind, apply(inclusion,1,function(i){
	out$mean_MS_selected_inter = foreach(idx = 1:ncol(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
		require(bigmemory);load('tmp.rda');  i = attach.big.matrix( paste0('best_gene.',dataset,'.desc')  )[,idx] ; 
		model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1 & (attach.big.matrix('Model_Score_Groups.desc')[,DS==invert_ds(dataset)])==1
		i = i[selected_models]; model_score_i = model_score_i[selected_models]
		idx_down = sample(c(TRUE,FALSE),size=round(length(i)/100),prob=c(.09,.91),replace=T) 
		if( length(table( i[idx_down] ))<2 ) {return(NA)}
		c( mean0 = (mn0<-mean(model_score_i[idx_down & i==0],na.rm=T)) , SD0=sd(model_score_i[idx_down & i==0],na.rm=T),
			mean1 = (mn1<-mean(model_score_i[idx_down & i==1],na.rm=T)) , SD1=sd(model_score_i[idx_down & i==1],na.rm=T) )
#		tmp=NA; try( tmp<-tidy( t.test( model_score_i[idx_down & i==1] , model_score_i[idx_down & i==0]  ,na.action="na.omit") ) ); tmp
	} #) )
	print(tail(out))
#	print('wilcox test')
###	out$t_wilcox_background = do.call(rbind, apply(inclusion,1,function(i){
#	out$t_wilcox_background = foreach(idx = 1:ncol(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
#		require(bigmemory);load('tmp.rda');  i = attach.big.matrix( paste0('best_gene.',dataset,'.desc')  )[,idx] ; 
#		model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
#		if( length(table( sample( i[selected_models] , 1e6)))<2 ){ if( length(table(i[selected_models] ))<2 ) {return(NA)}}
#		tmp=NA; try( tmp<-tidy( wilcox.test( model_score_i[i==1] , model_score_i[i==0]  ,na.action="na.omit") ) ); tmp
#	} #) )
###	out$t_wilcox_selected = do.call(rbind, apply(inclusion,1,function(i){
#	out$t_wilcox_selected = foreach(idx = 1:ncol(inclusion), .combine = rbind,.packages=c('bigmemory','broom'))  %dopar% {
#		require(bigmemory);load('tmp.rda');  i = attach.big.matrix( paste0('best_gene.',dataset,'.desc')  )[,idx] ; 
#		model_score_i=attach.big.matrix('Model_Scores.desc')[,DS==dataset]; 
#		selected_models=(attach.big.matrix('Model_Score_Groups.desc')[,DS==dataset])==1
#		if( length(table( sample( i[selected_models] , 1e6)))<2 ){ if( length(table(i[selected_models] ))<2 ) {return(NA)}}
#		tmp=NA; try( tmp<-tidy( wilcox.test( model_score_i[selected_models & i==1] , model_score_i[selected_models & i==0]  ,na.action="na.omit") ) ); tmp
#	} #) )
#	print(tail(out))
#	out$glm_include_background =    do.call(rbind,apply(inclusion,2,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i ~ i  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))
#	out$glm_include_selected =      do.call(rbind,apply(inclusion,2,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i[selected_models] ~ i[selected_models]  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))
#	out$glm_linkscore_background =  do.call(rbind,biganalytics::apply(ls,2,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i ~ i  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))
#	out$glm_linkscore_selected =    do.call(rbind,biganalytics::apply(ls,2,function(i){
#					if(length(table(i))>1){
#						tmp=tidy( glm( model_score_i[selected_models] ~ i[selected_models]  ,na.action="na.omit") ) 
#						data.frame(intercept=tmp[1,],reaction=tmp[2,])
#					}else{NA}
#				}))

	# select co-occurance

	save(out,file=paste0('data_out/Gene_Choice.',dataset,'.rda'))
	out
})

stopCluster(cl)

if(skip_secretor){
	Gene_Choice[[3]] = Gene_Choice[[1]]
	Gene_Choice[[4]] = Gene_Choice[[2]]
}

names(Gene_Choice) = DS

Gene_Choice=lapply(Gene_Choice, function(sc){
	idx=which( unlist(lapply(sc,is.data.frame) ) | unlist(lapply(sc,is.matrix) ) )
	cbind( sc[,-idx], do.call(cbind, lapply(idx,function(i){ 
		tmp=sc[,i]
		colnames(tmp) = paste(colnames(sc)[i],colnames(tmp),sep='.')
		tmp
	} ) ) )
})

save(Gene_Choice,file='data_out/Gene_Choice.rda')
#lapply(DS,function(ds) write.csv(Structure_Choice[[ds]],file=paste0('data_out/Structure_Choice.',ds,'.csv')))
#library(xlsx)
try(write.xlsx(Gene_Choice,file=paste0('data_out/Gene_Choice.xlsx')))

}






#####################################################
#####################################################
#####################################################
#####################################################
##### end
#####################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
##########################################################################################################
#####################################################
#####################################################
#####################################################
#####################################################


if(F){
#}else{
#}
#if(F){
	load('data_out/Gene_Choice.rda')
#	load('data_out/Gene_Choice_cbx.rda')

	pdf('figures/Gene_Choice.pdf',height=15)
	par(mfrow=c(3,1),mar=c(10,8,2,2))
	lapply(Gene_Choice,function(sc){
		barplot( t( data.matrix( sc[,c('p_background','p_selection','p_selection_inter')] ) ),names.arg=sc$rxn_name,beside=T,las=2,
			ylab='proportion',cex.names=.5,legend.text=c('Background','Selection','Intersection'))

		barplot( t( data.matrix( sc[,grepl('.mean1',colnames(sc))]  - sc[,grepl('.mean0',colnames(sc))] )) ,
			names.arg=sc$rxn_name,beside=T,las=2,ylab='Mean Model Score Improvement\n(with HMO - without HMO)',legend.text=c('Background','Selection','Intersection'))
	#	x=barplot( tmp<-t( data.matrix( sc[idx<-!grepl('HMO',sc$rxn_name),c('mean_MS_background.mean0','mean_MS_background.mean1')] ) ),names.arg=sc$rxn_name[idx],beside=T,las=2,ylab='Background Model Score Improvement\n(with HMO - without HMO)')
	#	error.bar(x,tmp, 1*(t( data.matrix( sc[idx<-!grepl('HMO',sc$rxn_name),c('mean_MS_background.SD0','mean_MS_background.SD1')] ) )))

		barplot( t( data.matrix( sc[,c('enrichment','enrichment_inter')]+1e-10 ) ),names.arg=sc$rxn_name,
			beside=T,las=2,ylab='Hypergeometric Enrichment (p-value)',cex.names=.5,legend.text=c('Selection','Intersection'))
	})
	dev.off()

	# pretty
#	df = do.call(rbind, lapply(DS,function(x) cbind(ds=x,Gene_Choice[[x]])))[,c('ds','rxn_name','p_background','p_selection','p_selection_inter')]
#	dfm = melt(df,c('ds','rxn_name'))
#	dfm$ds = ifelse(dfm$variable=='p_selection_inter',gsub('1|2','',as.character(dfm$ds)),as.character(dfm$ds))
#	dfm$var = paste(dfm$ds,dfm$variable,sep='.')
#	dfm$gene_name = unlist(lapply(strsplit(as.character(dfm$rxn_name),'_'),function(x) x[3]))
#	strs = gsub('_gene','',names(linkage))
#	g = grid.arrange(grobs=lapply(strs[-5],function(x){
#		ggplot(data=unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%c('data1','data2','data'),]) ,aes( x=gene_name,y=value,fill=var))+ # +facet_grid(~ds)
#			geom_bar(stat='identity',position='dodge')+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ggtitle(x)
#	}))
#	ggsave(g,filename='figures/Gene_Choice.data.pdf',height=10,width=15)

	##
	Gene_Choice=lapply(DS,function(x){
		mean0=grep('.mean0',colnames(Gene_Choice[[x]]),value=T)
		mean1=grep('.mean1',colnames(Gene_Choice[[x]]),value=T)
		tmp = (Gene_Choice[[x]][,mean1] - Gene_Choice[[x]][,mean0])
		colnames(tmp) = gsub('mean1','mean_diff',colnames(tmp))
		cbind(Gene_Choice[[x]],tmp)
	})
	names(Gene_Choice) = DS

	datas=list(data=c('data1','data2','data'),secretor=c('data1_secretor','data2_secretor','data_secretor'))
	vars = list(prop = c('p_background','p_selection','p_selection_inter'),
		enrich=grep('enrichment',colnames(Gene_Choice[[1]]),value=T),
		means=grep('.mean_diff',colnames(Gene_Choice[[1]]),value=T))
	df = do.call(rbind, lapply(DS,function(x) cbind(ds=x,Gene_Choice[[x]])))[,c('ds','rxn_name',vars$prop,vars$enrich,vars$means)]
	dfm = melt(df,c('ds','rxn_name'))
##	dfm$ds = ifelse(grepl('*_inter',dfm$variable)&!grepl('mean_diff',dfm$variable),gsub('1|2','',as.character(dfm$ds)),as.character(dfm$ds)) (this is causing the >100% error bc inter_12_1 and inter_12_2 are different because they are based on best genes specific to each dataset
#####	dfm$ds = ifelse(dfm$variable=='p_selection_inter',gsub('1|2','',as.character(dfm$ds)),as.character(dfm$ds))
	dfm$var = paste(dfm$ds,dfm$variable,sep='.')
	dfm$gene_name = unlist(lapply(strsplit(as.character(dfm$rxn_name),'_'),function(x) x[3]))
	strs = gsub('_gene','',names(linkage))
	for(d_i in names(datas)){ for(v_i in names(vars)){
		g = grid.arrange(grobs=lapply(strs[-5],function(x){
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[[v_i]],])
			if( v_i=='enrich' ){ data_tmp$value = -log( data_tmp$value+2e-16 , 10 ) }
			ggplot(data=data_tmp , aes( x=gene_name,y=value,fill=var))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ 
				theme(axis.text.x = element_text(angle = 90, hjust = 1),legend.text=element_text(size=5)) + 
				ggtitle(x)
		}))
		ggsave(g,filename=paste0('figures/Gene_Choice/Gene_Choice.',d_i,'.',v_i,'.pdf'),height=15,width=15)
	}}

		grbs =lapply(strs[-5],function(x){
			print(x)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1],])
#			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['prop']][-1],])
			data_tmp$i = gsub(paste0(x,'_'),'',as.character(data_tmp$rxn_name))
			#data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']],])
			data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
#			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['enrich']],])
			data_tmp$enrichment = apply( data_tmp , 1 , function(x){
				if(grepl('inter',x[3])){
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&grepl('inter',data_tmp_e$variable)] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&!grepl('inter',data_tmp_e$variable)] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.001,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+
#				geom_point()+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+facet_grid(~sec)+theme_classic()+
				theme(axis.text.x = element_text(angle = 90, hjust = 1))
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grobs=grbs)
		ggsave(g,filename=paste0('figures/Gene_Choice.pdf'),height=15,width=15)
		ggsave(l,filename=paste0('figures/Gene_Choice.legend.pdf'),height=5,width=5)

		grbs =lapply(strs[-5],function(x){
			print(x)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1] & !grepl('secretor',dfm$ds),])
#			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['prop']][-1],])
			data_tmp$i = gsub(paste0(x,'_'),'',as.character(data_tmp$rxn_name))
			#data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']]& !grepl('secretor',dfm$ds),])
			#data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
#			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['enrich']],])
			data_tmp$enrichment = apply( data_tmp , 1 , function(x){
				if(grepl('inter',x[3])){
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&grepl('inter',data_tmp_e$variable)] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&!grepl('inter',data_tmp_e$variable)] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.001,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+
#				geom_point()+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+theme_classic()+
				theme(axis.text.x = element_text(angle = 90, hjust = 1))
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grobs=grbs)
		ggsave(g,filename=paste0('figures/Gene_Choice.all.pdf'),height=15,width=15)
		ggsave(l,filename=paste0('figures/Gene_Choice.all.legend.pdf'),height=5,width=5)

		grbs =lapply(strs[-5],function(x){
			print(x)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1] & !grepl('secretor',dfm$ds)& grepl('inter',dfm$ds),])
#			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['prop']][-1],])
			data_tmp$i = gsub(paste0(x,'_'),'',as.character(data_tmp$rxn_name))
			#data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']]& !grepl('secretor',dfm$ds)& grepl('inter',dfm$ds),])
			#data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
#			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['enrich']],])
			data_tmp$enrichment = apply( data_tmp , 1 , function(x){
				if(grepl('inter',x[3])){
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&grepl('inter',data_tmp_e$variable)] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&!grepl('inter',data_tmp_e$variable)] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.001,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+
#				geom_point()+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+theme_classic()+
				theme(axis.text.x = element_text(angle = 90, hjust = 1))
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grobs=grbs)
		ggsave(g,filename=paste0('figures/Gene_Choice.all_inter.pdf'),height=15,width=15)
		ggsave(l,filename=paste0('figures/Gene_Choice.all_inter.legend.pdf'),height=5,width=5)


		grbs =lapply(strs[-5],function(x){
			print(x)
			data_tmp = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['prop']][-1] & grepl('secretor',dfm$ds),])
#			data_tmp = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['prop']][-1],])
			data_tmp$i = gsub(paste0(x,'_'),'',as.character(data_tmp$rxn_name))
			#data_tmp$HMO = gsub('[1-9]|[a-z]|_','',as.character(data_tmp$rxn_name))
			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name)& dfm$variable%in%vars[['enrich']]& grepl('secretor',dfm$ds),])
			#data_tmp$sec=ifelse(grepl('secretor',data_tmp$ds),'secretor','all')
#			data_tmp_e = unique(dfm[grepl(x,dfm$rxn_name) & dfm$ds%in%datas[[d_i]] & dfm$variable%in%vars[['enrich']],])
			data_tmp$enrichment = apply( data_tmp , 1 , function(x){
				if(grepl('inter',x[3])){
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&grepl('inter',data_tmp_e$variable)] )
				}else{
					mean( data_tmp_e$value[data_tmp_e$ds==x[1]&data_tmp_e$rxn_name==x[2]&!grepl('inter',data_tmp_e$variable)] )
				}})
			data_tmp$enriched = ifelse(data_tmp$enrichment>0.05,'',ifelse(data_tmp$enrichment>0.001,'*','**') )
			g=ggplot(data=data_tmp , aes( x=i,y=value,fill=var,label=enriched))+ # +facet_grid(~ds)
				geom_bar(stat='identity',position='dodge')+ggtitle(gsub('\\^','',x))+
#				geom_point()+
				geom_text(aes(y=value+.05),size=4,position=position_dodge(width = 1))+theme_classic()+
				theme(axis.text.x = element_text(angle = 90, hjust = 1))
			l <<- g_legend(g)
			g + theme(legend.position = 'none')
		})
		g = grid.arrange(grobs=grbs)
		ggsave(g,filename=paste0('figures/Gene_Choice.secretor.pdf'),height=15,width=15)
		ggsave(l,filename=paste0('figures/Gene_Choice.secretor.legend.pdf'),height=5,width=5)


#}

##################
lst=list()
lg=list()
if(cooc){
	print('co-occurance')
	cl<-makeCluster(spec = cli)
	registerDoParallel(cl = cl)
	print('background structure co-occurance')
	rxn_full =  annotate_inclusion()
	rxn = !grepl('HMO',rxn_full)

	lst[['background_coocur_count']] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),T,cocount,names=rxn_full[rxn],dim_i=1)
	lst[['background_coocur_jaccard']] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),T,jaccard,names=rxn_full[rxn],dim_i=1)

	for(dataset in DS){
		print(dataset)
		selected_models = Model_Score_Groups[,DS==dataset]==1
		selected_models_inter = Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1

		print('co-structures')
		
		lst[[paste('select_coocur_count',dataset,sep='_')]] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),selected_models,cocount,names=rxn_full[rxn],dim_i=1)
		lst[[paste('select_coocur_jaccard',dataset,sep='_')]] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),selected_models,jaccard,names=rxn_full[rxn],dim_i=1)

		lst[[paste('selectInter_coocur_count',dataset,sep='_')]] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),selected_models_inter,cocount,names=rxn_full[rxn],dim_i=1)
		lst[[paste('selectInter_coocur_jaccard',dataset,sep='_')]] = coocur_par( '../HMO_biosynthesis/inclusion.rxn.desc',which(rxn),selected_models_inter,jaccard,names=rxn_full[rxn],dim_i=1)


		print('co-gene')
		all_links = 1:length(gl_names)
		lg[[paste('background_coocur_count',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,all_links,T,cocount,names=gl_names,dim_i=2)
		lg[[paste('background_coocur_jaccard',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,all_links,T,jaccard,names=gl_names,dim_i=2)

		lg[[paste('select_coocur_count',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,all_links,selected_models,cocount,names=gl_names,dim_i=2)
		lg[[paste('select_coocur_jaccard',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,all_links,selected_models,jaccard,names=gl_names,dim_i=2)

		lg[[paste('selectInter_coocur_count',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,
			all_links,selected_models_inter,cocount,names=gl_names,dim_i=2)
		lg[[paste('selectInter_coocur_jaccard',dataset,sep='_')]] = coocur_par( paste0('best_gene.',dataset,'.desc') ,
			all_links,selected_models_inter,jaccard,names=gl_names,dim_i=2)
	}

	stopCluster(cl)

	print('output - structure')
	pdf('figures/cooccurance_heatmap.structures.pdf',height=10,width=10)
	lapply(names(lst),function(li){
		rn = gsub('_rxn','',rownames(lst[[li]]))
		heatmap.2( lst[[li]],main=li,trace='none',mar=c(10,10),col=topo.colors(15),Colv=F,Rowv=F,
			RowSideColors=rainbow(5)[factor(gsub('[1-9]','',rn))] , ColSideColors=rainbow(5)[factor(gsub('[1-9]','',rn))] , labRow=rn,labCol=rn)
	})
	dev.off()

	pdf('figures/cooccurance_heatmap.structures.cluster.pdf',height=10,width=10)
	lapply(names(lst),function(li){
		rn = gsub('_rxn','',rownames(lst[[li]]))
		heatmap.2( lst[[li]],main=li,trace='none',mar=c(10,10),col=topo.colors(15),
			RowSideColors=rainbow(5)[factor(gsub('[1-9]','',rn))] , ColSideColors=rainbow(5)[factor(gsub('[1-9]','',rn))] , labRow=rn,labCol=rn)
	})
	dev.off()

	print('output - gene')
	pdf('figures/cooccurance_heatmap.genes.pdf',height=10,width=10)
	lapply(names(lg),function(li){
		rn = rownames(lg[[li]])
		heatmap.2( lg[[li]],main=li,trace='none',mar=c(10,10),col=topo.colors(15),Colv=F,Rowv=F,
			RowSideColors=rainbow(10)[factor(unlist(lapply(strsplit(rn,split='_'),function(x) x[1])))] ,
			ColSideColors=rainbow(10)[factor(unlist(lapply(strsplit(rn,split='_'),function(x) x[1])))] , labRow=rn,labCol=rn)
	})
	dev.off()

	pdf('figures/cooccurance_heatmap.genes.cluster.pdf',height=10,width=10)
	lapply(names(lg),function(li){
		rn = rownames(lg[[li]])
		heatmap.2( lg[[li]],main=li,trace='none',mar=c(10,10),col=topo.colors(15),
			RowSideColors=rainbow(10)[factor(unlist(lapply(strsplit(rn,split='_'),function(x) x[1])))] ,
			ColSideColors=rainbow(10)[factor(unlist(lapply(strsplit(rn,split='_'),function(x) x[1])))] , labRow=rn,labCol=rn)
	})
	dev.off()
	save(lg,lst,file='data_out/coocur.rda')
}

if(anyl_CBX){
	########
}

hypergeo_overlap <- function(inter,l1,l2,background){
	if(l1<l2){
		tmp=l1
		l2=l1
		l1=tmp
	}
	phyper(inter-1, l1, background-l2, l2, lower.tail=FALSE)
}
hypergeo_overlap.mat<-function(m){
	m_out=matrix(1,nrow(m),nrow(m))
	for(i in 1:nrow(m)){
		for(j in 1:nrow(m)){
			l1 = which(m[i,]==1)
			l2 = which(m[j,]==1)
			m_out[i,j] = hypergeo_overlap( length(intersect(l1,l2)) , length(l1) , length(l2) , length(union(l1,l2)) )
		}
	}
	m_out
}

if(network){
	library(igraph)

	inclusion = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
	rxns=gsub('_rxn','',annotate_inclusion())
	rxn_sel = !grepl('HMO',rxns)
	
	# all
	model_sel = Model_Score_Groups[,1] & Model_Score_Groups[,2]
	dat = inclusion[rxn_sel,model_sel]
	rownames(dat) = rxns[rxn_sel]
	d = 1-as.matrix(dist( dat,method='binary'))
	d_rand = 1-as.matrix(dist( matrix(sample(as.vector(dat)),nrow(dat),ncol(dat)),method='binary'))
	diag(d) = 0; diag(d_rand)=0
	write.csv(d,file='data_out/cooc.network.data12.csv')
	#d_sig = hypergeo_overlap.mat(dat) ###
	d[d<quantile(d,.9)]= 0
	#d[d_sig>.01] = 0
	d=d[rowSums(d)>0,colSums(d)>0]
	g = graph.adjacency(d,mode='undirected',weighted=T)
	V(g)$color =  rainbow(5)[factor(gsub('[1-9]','', V(g)$name ) ) ]
	#V(g)$weight = sapply(V(g)$name,function(v) sum(inclusion[rxns==v,model_sel],na.rm=T)/sum(model_sel,na.rm=T) )
	#V(g)$weight=exp(-unlist(alpha_centrality(g)))
	pdf('figures/cooc.network.data12.pdf',height=10,width=10)
	plot(g , edge.color="black",edge.width=E(g)$weight*5  , vertex.size=unlist(V(g)$weight)) # ,layout=layout.lgl)
	dev.off()

	# secretor
	model_sel = Model_Score_Groups[,3] & Model_Score_Groups[,4]
	dat = inclusion[rxn_sel,model_sel]
	rownames(dat) = rxns[rxn_sel]
	d = 1-as.matrix(dist( dat,method='binary'))
	diag(d) = 0
	d[d<quantile(d,.9)]= 0
	d=d[rowSums(d)>0,colSums(d)>0]
	g = graph.adjacency(d,mode='undirected',weighted=T)
	V(g)$color =  rainbow(5)[factor(gsub('[1-9]','', V(g)$name ) ) ]
	V(g)$weight = sapply(V(g)$name,function(v) sum(inclusion[rxns==v,model_sel],na.rm=T)/sum(model_sel,na.rm=T) )
	pdf('figures/cooc.network.secretor12.pdf',height=10,width=10)
	plot(g , edge.color="black",edge.width=E(g)$weight*5  , vertex.size=V(g)$weight*10 ) # ,layout=layout.lgl)
	dev.off()
}

if(write_scores){
	print('write_scores')
	for(ds in DS){ 
		print(ds);
		ms = Model_Scores[,DS==ds]
		sel = Model_Score_Groups[,DS==ds]
		idx<-sel==1
		print(sum(idx,na.rm=T))
		tmp = cbind( 	selected_models=which(idx),
				model_scores=ms[idx],
				GeneLinkageScore=attach.big.matrix(paste0('gene_linkage_score.',ds,'.desc'))[idx,] )
		write.csv(tmp,file=paste0('Model_Scores.',ds,'.csv'),quote=F,row.names=F,col.names=T)
	}
	write.big.matrix(Model_Scores,filename='data_out/Model_Scores.csv')
	write.big.matrix(Model_Score_Groups,filename='data_out/Model_Score_Groups.csv')
}

library(gridExtra)
if(anyl_GLSvMS){
	print('GLSvMS')

	div=1e3
	#cl<-makeCluster(spec = 5)
	#registerDoParallel(cl = cl)
	glm_out = list()
	df_b_all = list()
	glm_full_out = list()
	for(ds in DS){
	#jnk=mclapply(DS,function(ds){
	#jnk=foreach(ds=DS) %dopar% {
		#png(paste0( 'figures/gene_linkage_score.vs.model_score.',ds,'.png')) #,height=20,width=20)
		#par(mfrow=c(7,8))
		print(ds)
		ms = Model_Scores[,DS==ds]
		sel = Model_Score_Groups[,DS==ds]
		best = attach.big.matrix(Best_Genes[[which(DS==ds)]])
		#idx<-which(sel==1)
		idx = 1:nrow(Model_Scores)
		idx = sample( idx , round(length(idx)/div))
		GLS = attach.big.matrix(Gene_Linkage_Scores[[which(DS==ds)]]) #[idx,]

		split = strsplit(gl_names,'_')
		genes = unlist(lapply(split,function(x) x[3]))
		links = unlist(lapply(split,function(x) paste(x[1:2],collapse='_')))

		df_b = do.call(rbind,lapply(1:ncol(GLS),function(gls){
			cat(paste(':',gl_names[gls]))
			na.omit( data.frame( gls_i = GLS[idx,gls], ms_i = ms[idx],included = best[idx,gls] ,
				linkage=links[gls], gene=genes[gls] ,selected=sel[idx] ) )
		}))
		df_b$included=as.factor(df_b$included)
		df_b_all[[ds]] = df_b #[df_b$included==1,]

		g=grid.arrange( grobs=lapply(unique(df_b$linkage)[-5],function(l){
			dat = df_b[df_b$linkage==l & df_b$selected==1,]
			ggplot(aes(y=gls_i,x=gene),data=dat)+geom_boxplot()+#facet_wrap(~included,scale='free_x')+
				geom_hline(aes(yintercept=min(-.5,median(gls_i))))+
				#geom_hline(aes(yintercept=median(gls_i)-sd(df_b$gls_i)))+
				#stat_smooth(method='glm')+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
		}) )
		ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score.',ds,'.pdf'),height=10,width=10)

		g=grid.arrange( grobs=lapply(unique(df_b$linkage)[-5],function(l){
			dat = df_b[df_b$linkage==l,]
			dat$keep_gene = dat$gene %in% unique(dat$gene)[sapply( unique(dat$gene) , function(g) median(dat$gls_i[dat$gene==g & dat$selected==1])>(-.5) )]
			ggplot(aes(y=gls_i,x=ms_i-mean(ms_i),color=gene),data=dat[dat$keep_gene,])+geom_jitter(size=.01,alpha=.05)+ #geom_density2d()+#facet_wrap(~included,scale='free_x')+
				stat_smooth(method='glm')+
				ggtitle(l)+xlab('Model Score (Selected Modeles)')+ylab('Gene-Linkage Score')
		}) )
		ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score.vs.model_score.',ds,'.png'),height=12,width=12)

		# univariate regression
		#df_b_all[[ds]] = df_b[df_b$included==1,]
		glm_out[[ds]] = as.data.frame(t(apply( unique(df_b[df_b$included==1,c('linkage','gene')]) , 1 , function(x){
			mod = glm( ms_i ~ gls_i , data= df_b[ df_b$linkage==x[1] & df_b$gene==x[2],] )
			c(x[1],x[2],coef(summary(mod))[1,],coef(summary(mod))[2,],ds)
		})))
		# multivariate regression
		sel_genes =  unique(df_b$gene)[sapply( unique(df_b$gene) , function(g) median(df_b$gls_i[df_b$gene==g & df_b$selected==1])>(-.5) )]
		dat = as.data.frame( GLS[idx,grepl( paste(sel_genes,collapse='|') , gl_names)] )
		colnames(dat) = gl_names[grepl( paste(sel_genes,collapse='|') , gl_names)]
		dat$model_score = ms[idx]
		glm_full_out[[ds]] = stepAIC( glm( model_score ~ . , data=dat ),k=20)
	}# ,mc.cores=5)
	#stopCluster(cl)
	tmp = do.call(rbind,df_b_all)
	tmp$ds = unlist(lapply( strsplit( rownames(tmp) , '\\.') , function(x) x[1] ))
	tmp$sec = ifelse(grepl('secretor',tmp$ds),'secretor','all')
	#tmp = tmp[tmp$selected==1,]
	g=grid.arrange( grobs=lapply(unique(tmp$linkage)[-5],function(l){
			dat = tmp[tmp$linkage==l & tmp$selected==1,]
			ggplot(aes(y=gls_i,x=gene,color=ds),data=dat)+geom_boxplot()+facet_grid(~sec)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score.pdf'),height=10,width=20)
	#tmp_gls = merge(merge( 
	tmp_gls = tmp[tmp$selected==1,] %>% group_by(linkage,gene,ds) %>% summarize(sd_gls=sd(gls_i),median_gls=median(gls_i),cor_glsVms=cor(gls_i,ms_i,method='spearman')) 
			#tmp[tmp$selected==1,] %>% group_by(linkage,gene,ds) %>% summarize(median(gls_i))),
			#tmp[tmp$selected==1,] %>% group_by(linkage,gene,ds) %>% mutate(tidy(cor.test(ms_i,gls_i,method='spearman'))))
	#colnames(tmp_gls) = c('linkage','gene','ds','sd_gls','median_gls','cor_glsVms')
	tmp_gls$sec = ifelse(grepl('secretor',tmp_gls$ds),'secretor','all')
	tmp_gls$set = ifelse(grepl('1',tmp_gls$ds),'1','2')

	g=grid.arrange( grobs=lapply(unique(tmp_gls$linkage)[!grepl('L5_',unique(tmp_gls$linkage))],function(l){
			dat = tmp_gls[tmp_gls$linkage==l,] # & tmp$selected==1,]
			ggplot(aes(x=cor_glsVms,y=median_gls,color=gene),data=dat)+facet_grid(sec~set)+
				geom_pointrange(aes(ymin = median_gls - sd_gls*1.96, ymax = median_gls + sd_gls*1.96)) +
				#geom_errorbarh(aes(xmax = slope + slope_SE*1.96, xmin = slope - slope_SE*1.96), height = 0) +
				geom_hline(yintercept=0)+geom_vline(xintercept=0)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score_vs_cor.pdf'),height=20,width=20)

	g=grid.arrange( grobs=lapply(unique(tmp_gls$linkage)[!grepl('L5_',unique(tmp_gls$linkage))],function(l){
			dat = tmp_gls[tmp_gls$linkage==l & tmp_gls$sec=='all',] # & tmp$selected==1,]
			ggplot(aes(x=cor_glsVms,y=median_gls,color=gene),data=dat)+facet_grid(~set)+
				geom_pointrange(aes(ymin = median_gls - sd_gls*1.96, ymax = median_gls + sd_gls*1.96)) +
				#geom_errorbarh(aes(xmax = slope + slope_SE*1.96, xmin = slope - slope_SE*1.96), height = 0) +
				geom_hline(yintercept=0)+geom_vline(xintercept=0)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score_vs_cor.all.pdf'),height=20,width=20)

	g=grid.arrange( grobs=lapply(unique(tmp_gls$linkage)[!grepl('L5_',unique(tmp_gls$linkage))],function(l){
			dat = tmp_gls[tmp_gls$linkage==l & tmp_gls$sec=='secretor',] # & tmp$selected==1,]
			ggplot(aes(x=cor_glsVms,y=median_gls,color=gene),data=dat)+facet_grid(~set)+
				geom_pointrange(aes(ymin = median_gls - sd_gls*1.96, ymax = median_gls + sd_gls*1.96)) +
				#geom_errorbarh(aes(xmax = slope + slope_SE*1.96, xmin = slope - slope_SE*1.96), height = 0) +
				geom_hline(yintercept=0)+geom_vline(xintercept=0)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score_vs_cor.secretor.pdf'),height=20,width=20)

	tmp = as.data.frame(do.call(rbind,glm_out))
	colnames(tmp) = c('linkage','gene','intercept','intercept_SE','intercept_t','intercept_p','slope','slope_SE','slope_t','slope_p','ds')
	tmp = merge(tmp,tmp_gls,by=c('linkage','gene','ds'))
	tmp$sec = ifelse(grepl('secretor',tmp$ds),'secretor','all')
	tmp$set = ifelse(grepl('1',tmp$ds),'1','2')
	tmp$slope = as.numeric(as.character(tmp$slope))
	tmp$slope_SE = as.numeric(as.character(tmp$slope_SE))
	tmp$intercept = as.numeric(as.character(tmp$intercept))
	tmp$intercept_SE = as.numeric(as.character(tmp$intercept_SE))
	#tmp = tmp[tmp$selected==1,]
	g=grid.arrange( grobs=lapply(unique(tmp$linkage)[!grepl('L5_',unique(tmp$linkage))],function(l){
			dat = tmp[tmp$linkage==l,] # & tmp$selected==1,]
			ggplot(aes(x=slope,y=median_gls,color=gene),data=dat)+facet_grid(sec~set)+
				geom_pointrange(aes(ymin = median_gls - sd_gls*1.96, ymax = median_gls + sd_gls*1.96)) +
				geom_errorbarh(aes(xmax = slope + slope_SE*1.96, xmin = slope - slope_SE*1.96), height = 0) +
				geom_hline(yintercept=0)+geom_vline(xintercept=0)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/gene_linkage_score_vs_slope.pdf'),height=20,width=20)

	g=grid.arrange( grobs=lapply(unique(tmp$linkage)[!grepl('L5_',unique(tmp$linkage))],function(l){
			dat = tmp[tmp$linkage==l,] #& tmp$selected==1,]
			ggplot(aes(x=slope,y=intercept ,color=gene),data=dat)+facet_grid(sec~set)+
				geom_pointrange(aes(ymin = intercept - intercept_SE*1.96, ymax = intercept + intercept_SE*1.96)) +
				geom_errorbarh(aes(xmax = slope + slope_SE*1.96, xmin = slope - slope_SE*1.96), height = 0) +
				#geom_hline(yintercept=0)+geom_vline(xintercept=0)+
				ggtitle(l)+ylab('Gene-Linkage Score (Selected Modeles)')
	}) )
	ggsave(plot=g,filename=paste0( 'figures/intercept_vs_slope.pdf'),height=20,width=20)
}
#dat = do.call(rbind,glm_out)
#colnames(dat) = c("linkage","gene","intercept_Estimate","intercept_Std.Error","intercept_t","intercept_Pr_t","var_Estimate","var_Std.Error","var_t","var_Pr_t","ds")
#ggplot( data=dat , aes(x=as.numeric(intercept_Estimate),y=as.numeric(var_Estimate),label=gene,color=gene))+geom_label()+facet_wrap(~linkage,scale='free')

if(colink){
	print('co-usage: linkages')
	pdf('figures/cousage_linkages.pdf',height=10,width=10)
	lapply(DS,function(li){
		tmp = cor(as.matrix(na.omit(attach.big.matrix(Linkage_Scores[[which(DS==li)]])[sample(1:48e6,10e6),])))
		diag(tmp) = 0
		heatmap.2( tmp,labRow=names(linkage),labCol=names(linkage),
			main=li,trace='none',mar=c(15,14),col=cm.colors(15),Colv=F,Rowv=F)
	})
	dev.off()

	print('co-usage: gene-linkage scores')
	pdf('figures/cousage_gene_linkages.pdf',height=10,width=10)
	lapply(DS,function(li){
		tmp = cor(as.matrix(na.omit(attach.big.matrix(Gene_Linkage_Scores[[which(DS==li)]])[sample(1:48e6,10e6),])))
		diag(tmp) = 0
		heatmap.2( tmp,labRow=gl_names,labCol=gl_names,
			main=li,trace='none',mar=c(10,10),col=cm.colors(15),Colv=F,Rowv=F)
	})
	dev.off()
}

run_lPCA<-function(dat,ks=NULL,ms=NULL,q=FALSE,downsample=10,k=NULL,m=NULL,title='gen'){
#	down_idx = sample(1:nrow(dat),round(nrow(dat)/downsample))
	print(paste('input dim:',dim(dat)))
	dat = dat[,apply(dat,2,sd)>.01]
	print(paste('cleaned dim:',dim(dat)))
	if(is.null(k)){ 
		k_i = foreach(k=ks,.combine='c',.packages='logisticPCA') %dopar% logisticSVD(dat, k = k,quiet=q)$prop_deviance_expl
		k = which.min(abs(k_i-.9))
		print(paste('k selected:',k))
	}
	if(is.null(m)){
		logpca_cv = foreach(m=ms,.combine='c',.packages='logisticPCA') %dopar% as.vector( cv.lpca(dat, ks = k, ms = m,quiet=q) )
		m= which.min(logpca_cv)
		print(paste('m selected:',m))
	}
	try(save(k_i,logpca_cv,file=paste0('pca_tmp.',title,'.rda')))
	print(paste('k,m set:',k,m))
	return(list( logpca_model = logisticPCA(dat, k = k, m = m,quiet=q),
		clogpca_model = convexLogisticPCA(dat, k = k, m = m,quiet=q)))
}

if(runPCA){
	library(logisticPCA)

	rxns = annotate_inclusion()
	out_lPCA = list()
	ks = seq(2,30,1)
	ms = seq(2,30,1)

	inclusion=attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')

	## all samples

	cl<-makeCluster(spec = cli); registerDoParallel(cl = cl)
	dat12=t(inclusion[,idx12<-(Model_Score_Groups[,1]==1 & Model_Score_Groups[,2]==1)])
	out_lPCA$data12 = run_lPCA( dat = dat12 ,ks=ks,ms=ms,title='ds12')
	plot(out_lPCA$data12[[2]],type='scores')+geom_point(color=col,size=.1)
	save(out_lPCA,file='out_lpca.rda')
	stopCluster(cl)

	#col = factor( ifelse( dat12[,rxns[apply(dat12,2,sd)>.01]=='FDSLNH6'] , 'red','black') )
	#plot(out_lPCA$data12[[2]],type='scores')+geom_point(color=col,size=.1)

	## variance distribution
	

	## where sugars
	score12 = 10^((Model_Scores[idx12,1]+Model_Scores[idx12,2])/2)
	hmos = c('^FDSLNH[1-9]','^FLNH[1-9]','^DSLNH[1-9]','^DFLNH[1-9]','^DFLNT[1-9]')
	g=grid.arrange(grobs=c(lapply( hmos , function(hi){
		col_f = factor( apply( dat12[,grepl(hi,rxns)] , 1 , function(x) (rxns[grepl(hi,rxns)])[x==1] ) )
		qplot(x=out_lPCA$data12[[2]]$PCs[,1],y=out_lPCA$data12[[2]]$PCs[,2],col=col_f)+geom_point(size=.1,alpha=.5)+theme_classic()
	}),
		list(qplot(x=out_lPCA$data12[[2]]$PCs[,1],y=out_lPCA$data12[[2]]$PCs[,2],col=score12,alpha=score12)+
			geom_point(size=.01)+theme_classic()+scale_color_gradient(low='white', high='red')
	)))
	ggsave(g,filename='PCA.hmo_location.all.pdf')

	#qplot(out_lPCA$data12[[2]]$PCs[,1] , score12)

	## loadings vis
	loadings = data.frame(out_lPCA$data12[[2]]$U)
	loadings$hmoi = rxns[apply(dat12,2,sd)>.01]
	loadings$hmo = gsub('[1-9]|_rxn','',loadings$hmoi)

	g=grid.arrange(grobs=c(lapply( hmos , function(hi){
		ggplot( data = loadings[grepl(hi,loadings$hmoi),] , aes( x=X1,y=X2,color=hmoi,shape=hmo) ) + 
		geom_segment(xend=0,yend=0,lineend='round') #+ xlim(c(-.2,.2)) + ylim(c(-.2,.2))
	}),
		list(qplot(x=out_lPCA$data12[[2]]$PCs[,1],y=out_lPCA$data12[[2]]$PCs[,2],col=score12,alpha=score12)+
			geom_point(size=.01)+theme_classic()+scale_color_gradient(low='white', high='red')
	)))
	ggsave(g,filename='PCA.hmo_loadings.all.pdf')


	##### secretor samples

	cl<-makeCluster(spec = cli); registerDoParallel(cl = cl)
	dat34=t(inclusion[,idx34<-(Model_Score_Groups[,3]==1 & Model_Score_Groups[,4]==1)])
	out_lPCA$secretor12 = run_lPCA( dat = dat34 ,k=ks,m=ms,title='sec12')
	save(out_lPCA,file='out_lpca.rda')
	stopCluster(cl)

	#col = factor( ifelse( dat34[,rxns[apply(dat34,2,sd)>.01]=='FDSLNH6'] , 'red','black') )
	#plot(out_lPCA$data12[[2]],type='scores')+geom_point(color=col,size=.1)

	# where sugars
	score34 = 10^((Model_Scores[idx34,3]+Model_Scores[idx34,4])/2)
	hmos = c('^FDSLNH[1-9]','^FLNH[1-9]','^DSLNH[1-9]','^DFLNH[1-9]','^DFLNT[1-9]')
	g=grid.arrange(grobs=c(lapply( hmos , function(hi){
		col_f = factor( apply( dat34[,grepl(hi,rxns)] , 1 , function(x) (rxns[grepl(hi,rxns)])[x==1] ) )
		qplot(x=out_lPCA$secretor12[[2]]$PCs[,1],y=out_lPCA$secretor12[[2]]$PCs[,2],col=col_f)+geom_point(size=.1,alpha=.5)+theme_classic()
	}),
		list(qplot(x=out_lPCA$secretor12[[2]]$PCs[,1],y=out_lPCA$secretor12[[2]]$PCs[,2],col=score34,alpha=score34)+
			geom_point(size=.01)+theme_classic()+scale_color_gradient(low='white', high='red')
	)))	
	ggsave(g,filename='PCA.hmo_location.secretor.pdf')

	# loading vis
	loadings = data.frame(out_lPCA$secretor12[[2]]$U)
	loadings$hmoi = rxns[apply(dat34,2,sd)>.01]
	loadings$hmo = gsub('[1-9]|_rxn','',loadings$hmoi)

	g=grid.arrange(grobs=c(lapply( hmos , function(hi){
		ggplot( data = loadings[grepl(hi,loadings$hmoi),] , aes( x=X1,y=X2,color=hmoi,shape=hmo) ) + 
		geom_segment(xend=0,yend=0,lineend='round') #+ xlim(c(-.2,.2)) + ylim(c(-.2,.2))
	}),
		list(qplot(x=out_lPCA$secretor12[[2]]$PCs[,1],y=out_lPCA$secretor12[[2]]$PCs[,2],col=score34,alpha=score34)+
			geom_point(size=.01)+theme_classic()+scale_color_gradient(low='white', high='red')
	)))
	ggsave(g,filename='PCA.hmo_loadings.secretor.pdf')


	library(ProjectionBasedClustering)
	#tsne12 = tSNE(DataOrDists=dat12,k=10,OutputDimension=2,method="binary",Whitening=TRUE,
	#	InitialDimensions=NULL, Iterations=1000,PlotIt=FALSE) #,Cls)

}

if(cluster){
library(factoextra)
	#cl<-makeCluster(spec = 25)
	#registerDoParallel(cl = cl)
	out=list()
	n=1e4
	if(F){
	for(dataset in DS){
		print(dataset)
		# load inclusion
		inclusion_i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
		# select models
		model_score_i = Model_Scores[,DS==dataset]
		infeasible = is.na(model_score_i)
		selected_models = which( Model_Score_Groups[,DS==dataset] & !infeasible)
		selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible)
		#stop()
		#kmm = kmeans(t(Matrix(inclusion_i[,selected_models],sparse=T)),3,nstart = 50,iter.max = 15)
		df=t(Matrix(inclusion_i[,selected_models],sparse=T))
		out[[paste0(dataset,'_selection_wss')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "wss")
		out[[paste0(dataset,'_selection_silhouette')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "silhouette")
		out[[paste0(dataset,'_selection_gap')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "gap_stat")

		df=t(Matrix(inclusion_i[,selected_models_inter],sparse=T))
		out[[paste0(dataset,'_intersection_wss')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "wss")
		out[[paste0(dataset,'_intersection_silhouette')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "silhouette")
		out[[paste0(dataset,'_intersection_gap')]]=fviz_nbclust(as.matrix(df[sample(1:nrow(df),n),]), kmeans, method = "gap_stat")

		gc(reset=T)
		#tmp=NbClust(data = as.matrix(df[sample(1:nrow(df),n),]), distance = "binary", min.nc = 2, max.nc = 10, method = 'kmeans', index = "silhouette", alphaBeale = 0.1)
		#tmp=NbClust(data = as.matrix(df[sample(1:nrow(df),n),]), distance = "binary", min.nc = 2, max.nc = 10, method = 'kmeans', index = "gap", alphaBeale = 0.1)

		#out[[dataset]] = list(selected_wss=wss,selected_sil=silhouette,selected_inter_wss=wss2,selected_inter_sil=silhouette2)
	}
	try( grid.arrange(grobs=lapply(names(out),function(x) out[[x]]+ggtitle(x)),ncol=6) )
	}

	set_k = list( 	data1_selection = 3 , data1_intersection = 5, data2_selection = 2 , data2_intersection = 5,
			data1_secretor_selection = 3 , data1_secretor_intersection = 3, data2_secretor_selection = 2 , data2_secretor_intersection = 3)
	library(mclust)
	out2=list()
	out3=list()
	for(dataset in DS){
		print(dataset)
		# load inclusion
		inclusion_i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
		# select models
		model_score_i = Model_Scores[,DS==dataset]
		infeasible = is.na(model_score_i)
		selected_models = which( Model_Score_Groups[,DS==dataset] & !infeasible)
		selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible)

		df=t(Matrix(inclusion_i[,selected_models],sparse=T))
		km = kmeans( df , centers= set_k[[paste0(dataset,'_selection')]] )
		#km = kmeans( as.matrix(df[sample(1:nrow(df),n),]) , centers= set_k[[paste0(dataset,'_selection')]] )
		out2[[paste0(dataset,'_selection_pred')]] = km #predict(df,km)
		out2[[paste0(dataset,'_selection_mclust')]] = Mclust(df) #predict(df,km)

		if(grepl('1',dataset)){
		 df=t(Matrix(inclusion_i[,selected_models_inter],sparse=T))
		 km = kmeans( df , centers= set_k[[paste0(dataset,'_intersection')]] )
		 #km = kmeans( as.matrix(df[sample(1:nrow(df),n),]) , centers= set_k[[paste0(dataset,'_selection')]] )
		 out2[[paste0(dataset,'_intersection_pred')]] = km #predict(df,km)
		 out2[[paste0(dataset,'_intersection_mclust')]] = Mclust(df) #predict(df,km) 
		 #km = bigkmeans(x, centers, iter.max = 10, nstart = 1, dist = "euclid")
		}
		
		out3[[paste0(dataset,'_selection_idx')]] = selected_models
		out3[[paste0(dataset,'_intersection_idx')]] = selected_models_inter
		
		gc(reset=T)
	}
	save(out,out2,out3,file='data_out/clusters.rda')
	
	library(useful)
	library(MASS)
	out4=list()
	for(dataset in DS){
		print(dataset)
		# load inclusion
		inclusion_i = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
		# select models
		model_score_i = Model_Scores[,DS==dataset]
		infeasible = is.na(model_score_i)
		selected_models = which( Model_Score_Groups[,DS==dataset] & !infeasible)
		selected_models_inter = which( Model_Score_Groups[,DS==dataset]==1 & Model_Score_Groups[,DS==invert_ds(dataset)]==1 & !infeasible)

		df=t(Matrix(inclusion_i[,selected_models],sparse=T))
		km=out2[[paste0(dataset,'_selection_pred')]]
		out3[[paste0(dataset,'_selection_pca')]] = logisticPCA(as.matrix(df))
 
		if(grepl('1',dataset)){
		 df=t(Matrix(inclusion_i[,selected_models_inter],sparse=T))
		 km=out2[[paste0(dataset,'_intersection_pred')]]
		 out3[[paste0(dataset,'_intersection_pca')]] = logisticPCA(as.matrix(df))
		} 
	}
	save(out,out2,out3,out4,file='data_out/clusters.rda')

	out5 = list()
	out6 = list()
	# focused clustering
	dat = inclusion[ ,idx<-which(!is.na(Model_Scores[,1]) & Model_Score_Groups[,1]==1 & Model_Score_Groups[,2]==1 & Model_Score_Groups[,3]==1 & Model_Score_Groups[,4]==1) ]
	scores = Model_Scores[idx,]; colnames(scores) = DS
	out5$data12_secretor12_intersection = hclust(dist(t(dat),method='binary'),method="ward.D2")
	out6$data12_secretor12_intersection = data.frame(model=idx,cut=cutree( out5$data12_secretor12_intersection , k=4 ),model_scores=scores)

	dat = inclusion[ ,idx<-which(!is.na(Model_Scores[,1]) & Model_Score_Groups[,1]==1 & Model_Score_Groups[,2]==1)  ]
	scores = Model_Scores[idx,]; colnames(scores) = DS
	out5$data12_intersection = hclust(dist(t(dat),method='binary'),method="ward.D2")
	out6$data12_intersection = data.frame(model=idx,cut=cutree( out5$data12_intersection , k=5 ),model_scores=scores)

	dat = inclusion[ ,idx<-which(!is.na(Model_Scores[,1]) & Model_Score_Groups[,3]==1 & Model_Score_Groups[,4]==1) ]
	scores = Model_Scores[idx,]; colnames(scores) = DS
	out5$secretor12_intersection = hclust(dist(t(dat),method='binary'),method="ward.D2")
	out6$secretor12_intersection = data.frame(model=idx,cut=cutree( out5$secretor12_intersection , k=c(3,6) ),model_scores=scores)
	colnames( out6$secretor12_intersection )[2] = 'cut'
	lapply(names(out6),function(n){
		write.csv(out6[[n]],file=paste0('data_out/clusters/cluster.',n,'.csv'))
	})
	save(out,out2,out3,out4,out5,out6,file='data_out/clusters.rda')

	# cluster matching
	matching = list()
	cnt=0
	for(i in names(out6)){
		for(j in names(out6)){
		if(i==j){next}
			for(ii in unique(out6[[i]]$cut)){
				for(jj in unique(out6[[j]]$cut)){
					icut = out6[[i]]$model[out6[[i]]$cut==ii]
					jcut = out6[[j]]$model[out6[[j]]$cut==jj]
					matching[[paste(i,ii,j,jj,sep='.')]] = c(ii,jj, length(intersect( icut , jcut )), 
						length(union( icut , jcut)) , min(length(icut) , length(jcut)) )
				}
			}
		}
	}
	tmp = as.data.frame(do.call(rbind,matching))
	colnames(tmp) = c( 'cut1','cut2','intersection','union','lengthMin')
	tmp$similarity = tmp$intersection/tmp$lengthMin
	tmp = tmp[tmp$intersection>0,]
	tmp[order(-tmp$similarity),]

	# cluster proportions
	tmp=do.call(cbind,lapply( names(out6) , function(xn){
		x = out6[[xn]]
		i=sapply(sort(unique(x$cut)) , function(ci){
			sel_i = x$model[x$cut==ci]
			rxn_names = gsub('_rxn','',annotate_inclusion())
			inc_rxns = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx<-!(grepl('HMO',rxn_names)),sel_i]
			rxn_fact  =  factor( gsub('[1-9]','',rxn_names[idx]))
	
			rownames(inc_rxns) = rxn_names[idx]
			apply(inc_rxns,1,function(x) sum(x)/length(x)) 
		})
		colnames(i) = paste(xn,paste0('cut',sort(unique(x$cut))),sep='.')
		i
	}))

	cluster_prop = melt(tmp)
	colnames(cluster_prop) = c('HMOi','X2','proportion')
	cluster_prop$HMO = gsub('[1-9]','',cluster_prop$HMOi)
	cluster_prop$i = gsub('[A-Z]','',cluster_prop$HMOi)
	cluster_prop$compare = unlist(lapply( strsplit( as.character(cluster_prop$X2) , '\\.') , function(x) x[1]))
	cluster_prop$cut = unlist(lapply( strsplit( as.character(cluster_prop$X2) , '\\.') , function(x) gsub('cut','',x[2])))
	g = ggplot(cluster_prop,aes(x=i,y=proportion,fill=cut))+geom_bar(stat='identity',position='dodge')+facet_grid(compare~HMO,scale='free')
	ggsave(g,filename='figures/clusters/proportions.pdf')


}

if(cluster_vis){
	load('data_out/clusters.rda') # out,out2,out3,out4,out5,out6
	sel_clust = out6$data12_intersection
	#el = get_network_edgelist(T); 
	el = get_network_edgelist(F)
	dif = get_rxn()
	el =  data.frame(el,dif,Structure_Choice[[1]]$p_selection_inter,1-Structure_Choice[[1]]$enrichment_inter,
		1-Structure_Choice[[1]]$enrichment,1-Structure_Choice[[2]]$enrichment,
		1-Structure_Choice[[1]]$p_selection,Structure_Choice[[2]]$p_selection)
	colnames(el) = c('reactant_id','product_id','reaction','intersection_proportion','intersection_enrichment','dataset1_enrichment','dataset2_enrichment','dataset1_proportion','dataset2_proportion')
	el_na = na.omit(el)
	weighted_mean_net = aggregate( t( tmp<-(inclusion[,sel_clust$model] * rowMeans(scale(sel_clust[,3:4],center=F))) ) , by=list(sel_clust$cut) , mean)
	# median_net = aggregate( t( tmp<-(inclusion[,sel_clust$model] ) ) , by=list(sel_clust$cut) , median)
	cluster_score = aggregate( sel_clust[,3:4],by=list(sel_clust$cut),mean)
	gscore=ggplot( data=cluster_score , aes(x=model_scores.data1, y=model_scores.data2 ,shape=factor(Group.1),color=factor(Group.1)) )+geom_point()+theme_classic()
	# gscore=ggplot( data=sel_clust , aes(x=model_scores.data1, y=model_scores.data2 ,shape=factor(cut),color=factor(cut)) )+geom_point()+geom_density2d()+theme_classic()
	net_graph = list()
	mets=get_mets()

	# write.csv(mets,file='figures/networks.metabolites')

	# probability and enrichment
	g = graph_from_edgelist(as.matrix(el_na[,1:2]))
	E(g)$propotion_inter = el_na[,4]
	E(g)$enrichment_inter = el_na[,5]
	E(g)$enrichment1 = el_na[,6]
	E(g)$enrichment2 = el_na[,7]
	E(g)$weight = E(g)$propotion_inter +E(g)$propotion_inter*E(g)$enrichment_inter
	E(g)$weight = E(g)$weight/max(E(g)$weight)
	#E(g)$color = rainbow(10)[factor(el_na[,3])]
	E(g)$color = c('gold3','gold3','red','red','red','blue','blue','purple','purple','grey')[factor(el_na[,3])]
	#E(g)$color = apply(data.frame(E(g)$color,E(g)$weight),1,function(x) adjustcolor(x[1],x[2]))
	E(g)$color = apply(data.frame(E(g)$color,2^(E(g)$weight)/(max(2^(E(g)$weight)))),1,function(x) adjustcolor(x[1],x[2]))
	E(g)$lty = c(1,2,1,2,3,1,2,1,2,1)[factor(el_na[,3])]
	#E(g)$weight = ifelse(E(g)$weight==0,NA,E(g)$weight)
	V(g)$weight = strength(g)/degree(g)
	V(g)$name_num = 1:101
	known = data.frame(id=c(97:101,3,4,6,7,5,11,13,18,15,21,31, 28,29,39,57,58,46,48,49,60,56,69,72,79,78,80,81,70,76,77,89,90,91,92,96,95,94,93),
		name=gsub('\\;HMO\\[c\\]','',c(mets[97:101],'2FL','3FL','LNT','LNnT','3SL','LNFPI','LNFPII','LNFPIII','LSTb','LSTc','DSLNT',
			paste0('DFLNT',1:3),paste0('FLNH',1:7),paste0('DFLNH',1:7),paste0('DSLNH',1:3),paste0('FDSLNH',1:7) )))
	V(g)$label = sapply(V(g)$name_num,function(i){ ifelse(i %in% known$id[1:16], as.character(known$name[known$id==i]) , i ) })
	V(g)$label_all = sapply(V(g)$name_num,function(i){ ifelse(i %in% known$id, as.character(known$name[known$id==i]) , i ) })
	tmp = sapply(shortest_paths(g,from=V(g)[1])$vpath,length)
	tmp[1] = 1
	V(g)$depth = tmp
	V(g)$depth_norm = tmp/max(tmp)
	V(g)$weight_norm = V(g)$weight/V(g)$depth_norm


	# path weight distribution
	paths=all_simple_paths(g,from=1,to=known$id[1:16])
	weights = unlist(lapply(paths,function(p) sum(E(g,path=p)$weight)))
	term = 	unlist(lapply(paths,function(p) as.character(rev(p)[1])))
	paths_char = unlist(lapply(paths,function(p) paste(as.character(p),collapse='->') ) )
	df = data.frame(term,weights,percentile=ecdf(weights)(weights),percentile_term=weights,paths=paths_char)	
	for(t in unique(term)){
		if(sum(df$term==t)<=2){
			df$percentile_term[df$term==t] = 1
			# df$percentile_term_q[df$term==t] = 0
		}else{
			df$percentile_term[df$term==t] = (ecdf(df$weights[df$term==t])(df$weights[df$term==t]))
			# df$percentile_term_q[df$term==t] = p.adjust(1-df$percentile_term[df$term==t],'fdr')			
		}
	}
	sel=sort(unique(unlist(strsplit(as.character(df[ df$percentile_term>.95 , 'paths']),'->'))))
	V(g)$frame.color = ifelse(1:length(V(g)) %in% sel,'black',NA)



	#if(!'lay'%in%objects()){
		lay = layout_as_tree(g)
	#	lay[,2] = lay[,2]+runif(nrow(lay))/2
		lay[,1] = lay[,1]*10
	#}

	pdf('figures/networks.prob.enrich.pdf',width=15,height=15)
	par(mfrow=c(2,2),mar=c(0,0,0,0))
	plot(g,layout=lay,
		edge.width=E(g)$propotion_inter*5,edge.arrow.size=0,edge.color='grey',edge.lty=1,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	plot(g,layout=lay,
		edge.width=E(g)$enrichment_inter*5,edge.arrow.size=0,edge.color='grey',edge.lty=1,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	plot(g,layout=lay,
		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color='grey',edge.lty=1,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	plot(g,layout=lay,
		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color=E(g)$color,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	dev.off()

	pdf('figures/networks.prob.enrich.sel.pdf',width=15,height=10)
	par(mfrow=c(1,2),mar=c(0,0,0,0))
	plot(g,layout=lay,
		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color=E(g)$color,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA)
	idx=which(V(g)$frame.color=='black')
	g2=induced_subgraph(g,V(g)[idx])
	plot(g2,layout=layout_as_tree(g2), #lay[idx,],
		edge.width=E(g2)$weight*5,edge.arrow.size=0,edge.color=E(g2)$color,
		vertex.size=10,vertex.label.cex=.8,vertex.label = V(g2)$label_all,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	dev.off()

	### shortest paths between measured sugars
	id=grepl('[A-Z]',V(g2)$label)
	get_d <- function(g,mode,id){ d=distances(g,mode=mode,weights=NA); colnames(d) = rownames(d) = V(g)$label; melt(d[id,id]) }
	d = cbind(get_d(g2,'all',id),get_d(g2,'in',id)[,3],get_d(g2,'out',id)[,3])
	d[,4] = ifelse(is.finite(d[,4]),d[,4],NA)
	d[,5] = ifelse(is.finite(d[,5]),d[,5],NA)
	colnames(d) = c('parent','child','d_all','d_in','d_out')

	# load conc. 
	d12_hmo = read.csv('data/raw/Data12.csv')
	d12_hmo$Dataset = paste0('d',d12_hmo$Dataset)
	colnames(d12_hmo) = gsub('\\.','',colnames(d12_hmo))
	cr = melt(cr_m<-cor(na.omit(d12_hmo[d12_hmo$Dataset=='d1',4:19]),method='spearman'))
	#cr = melt(cor(scale(log(na.omit(d12_hmo[d12_hmo$Dataset=='d1',4:19]))),method='pearson'))
	cr[,1] = gsub('X','',as.character(cr[,1])); 	cr[,2] = gsub('X','',as.character(cr[,2]))
	colnames(cr) = c('parent','child','cor')

	network = list(
		metabolites=data.frame(metabolites=mets,id=1:length(mets),label=V(g)$label_all, in_percentile_95th_path = ifelse(1:length(V(g)) %in% sel,TRUE,FALSE),
			strength.edge_weight_sum=strength(g),degree=degree(g),depth=V(g)$depth),
		reactions=data.frame(reaction=el_na$reaction,
			reactant_id=el_na$reactant_id,reactant = mets[el_na$reactant_id],product_id=el_na$product_id,product=mets[el_na$product_id],
			proportion_dataset1 = el_na$dataset1_proportion,proportion_dataset2 = el_na$dataset2_proportion,proportion_intersection = el_na$intersection_proportion,
			enrichment_dataset1 = el_na$dataset1_enrichment,enrichment_dataset2 = el_na$dataset2_enrichment,enrichment_intersection = el_na$intersection_enrichment,
			weight.enrichment_scaled_proportion = E(g)$weight ),
		paths=df,correlation=cr_m)
	write.xlsx(network,file='figures/network.xlsx')

	m = merge(d,cr)
	m = merge(m,network$metabolites[,c('label','degree','depth','metabolites')],by.x='parent',by.y='label')

	# cooperative regulation (same pathway)
	gp=ggplot(data=droplevels(na.omit(m[m$cor<1,c('cor','d_out','parent','child','depth')])),aes(y=cor,x=d_out,label=child)) + 
		geom_point() + geom_label(aes(fill=parent)) #+ stat_smooth()
	ggsave(gp,filename='figures/network.corVupstream.pdf')
	gp=ggplot(data=droplevels(na.omit(m[m$cor<1,c('cor','d_out','parent','child','depth')])),aes(y=cor,x=factor(depth+1),label=child)) + 
		geom_boxplot() + geom_label(aes(fill=parent))
	ggsave(gp,filename='figures/network.corVdepth.pdf')

	# competative regulation (different pathways)
	library(stringr)
	has<-function(x,sacch){strcount(x,sacch)}
	makeup<-function(x,sacch=c('GN','Ab','Fa','NNa')){sapply(sacch,function(s) str_count(x,s))}
	m2 = merge(m,network$metabolites[,c('label','degree','depth','metabolites')],by.x='child',by.y='label')
	m2$d_compete = ifelse(is.na(m$d_out) | is.na(m$d_in),m$d_all,NA)
	m2$min_depth = apply(m2[,c('depth.x','depth.y')],1,min)
	m2$max_depth = apply(m2[,c('depth.x','depth.y')],1,max)
	m2$mean_depth = apply(m2[,c('depth.x','depth.y')],1,mean)
	m2$from_LNT.y = m2$child %in% m2[!is.na(m2$d_out) & m2$parent=='LNT' , 'child']
	m2$from_LNT.x = m2$parent %in% m2[!is.na(m2$d_out) & m2$parent=='LNT' , 'child']
	#m2 = data.frame(m2,abs(makeup(m2$metabolites.x)*makeup(m2$metabolites.y)))
	m2 = data.frame(m2,abs(makeup(m2$metabolites.x)-makeup(m2$metabolites.y)))


	mod=stepAIC(glm(cor ~ (GN+Ab+Fa+NNa)^2,data=m2[m2$cor<1,]))
	mod=stepAIC(glm(cor ~ GN+Ab+Fa+NNa,data=m2[m2$cor<1,]))

	gp=ggplot(m2,aes(x=factor(Ab),color=factor(Fa),y=cor))+geom_boxplot()+geom_point()
	ggsave(gp,filename='figures/network.competition.pdf')

	#gp=ggplot(m2,aes(x=factor(GN),color=factor(NNa),y=cor))+geom_boxplot()+geom_point()
	#ggsave(gp,filename='figures/network.competition.pdf')


	ggplot(data=droplevels(na.omit(m2[m2$cor<1,c('cor','d_compete','parent','child','mean_depth')])),aes(x=mean_depth,y=cor )) + 
		geom_point() +stat_smooth(method='glm')+ facet_wrap(~d_compete)
	ggplot(data=droplevels(na.omit(m2[m2$cor<1 ,c('cor','d_compete','parent','child','max_depth')])),aes(x=d_compete,y=cor,label=paste(parent,child,sep='-'))) + 
		geom_point() +stat_smooth(method='glm')+ facet_wrap(~max_depth) + geom_label()

	gp=ggplot(data=droplevels(na.omit(m2[m2$cor<1 ,c('cor','d_compete','parent','child','max_depth')])),aes(x=d_compete,y=cor,label=child,fill=child)) + 
		geom_point() +stat_smooth(method='glm')+ facet_wrap(~parent) + geom_label()
	ggsave(gp,filename='figures/network.d_compete.all.pdf',height=20,width=20)

	gp=ggplot(data=droplevels(na.omit(m2[m2$cor<1 & !m2$child%in%c('FLNH','FDSLNH','2FL','3FL','3SL') & !m2$parent%in%c('FLNH','FDSLNH','2FL','3FL','3SL') ,
		c('cor','d_compete','parent','child','max_depth','from_LNT.y')])),aes(x=d_compete,y=cor,label=child,color=from_LNT.y)) + 
		geom_point() +stat_smooth(method='glm')+ facet_wrap(from_LNT.y~parent) + geom_text()
	ggsave(gp,filename='figures/network.d_compete.pdf')
	gp=ggplot(data=droplevels(na.omit(m2[m2$cor<1 & !m2$child%in%c('FLNH','FDSLNH','2FL','3FL','3SL') & !m2$parent%in%c('FLNH','FDSLNH','2FL','3FL','3SL') ,
		c('cor','d_compete','parent','child','max_depth','from_LNT.x','from_LNT.y','depth.y')])),
		aes(y=cor,label=child,x=from_LNT.x,color=from_LNT.y)) + geom_boxplot() + facet_wrap(~depth.y)
	ggsave(gp,filename='figures/network.d_compete.agg.pdf')

	# stats
	corVdepth = droplevels(na.omit(m[m$cor<1,c('cor','d_out','parent','child','depth')]))
	t.test( cor ~ factor(depth) , data=corVdepth[corVdepth$depth<5,])

	# ### try sugar specific comparison

	# segment compare
	names =list(top=c('2','2FL','3FL','3SL'),
		left=c('LNT','LNFPI','LNFPII','LSTb','DFLNT','DSLNH','DSLNT'),
		right=c('LNnT','LNFPIII','LSTc','DFLNH'))
	distr=lapply(names,function(n1){
		lapply(names,function(n2){
				cr[cr[,1]%in%n1 & cr[,2]%in%n2 & cr[,3]<1,]
			})
		})
	distr_m = melt(distr)
	ggplot(distr_m,aes(x=L1,y=value,color=L2))+geom_boxplot()

	# vis correlation graph 
	cri = cr[cr[,3]<1,]
	cri[,1] = sapply(cri[,1],function(h) which(V(g2)$label==h) )
	cri[,2] = sapply(cri[,2],function(h) which(V(g2)$label==h) )

	library(colorRamps)
	#cri$col = cm.colors(8)[cut(cri[,3],breaks=c(-1,-.5,-.25,-.1,0,.1,.25,.5,1))]
	cri$col = matlab.like2(8)[cut(cri[,3],breaks=c(-1,-.5,-.25,-.1,0,.1,.25,.5,1))]

	for(e in 1:nrow(cri)){
		cri[e,1:2] = sort(as.numeric(cri[e,1:2]))
	}
	cri=unique(cri)

	g3 = g2
	for(e in 1:nrow(cri)){
		g3 = g3+edges(c(cri[e,1],cri[e,2]),weight=abs(cri[e,3]),color=cri$col[e],lty=1,arrow.size=0,type=ifelse(cri[e,3]>0,'+','-'))
	}
	V(g3)$label_g3 = ifelse(grepl('[A-Z]',V(g3)$label),V(g3)$label,NA)
	# g3 = simplify(g3)

	pdf('figures/networks.prob.enrich.sel.cor.pdf',width=15,height=15)
	par(mfrow=c(2,2),mar=c(0,0,0,0))
	plot(g,layout=lay,
		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color=E(g)$color,
		vertex.size=5,vertex.label.cex=.5,vertex.label.family='sans',
		vertex.color=NA)
	idx=which(V(g)$frame.color=='black')
	g2=induced_subgraph(g,V(g)[idx])
	lay2 = layout_as_tree(g2)
	plot(g2,layout=lay2, #lay[idx,],
		edge.width=E(g2)$weight*5,edge.arrow.size=0,edge.color=E(g2)$color,
		vertex.size=10,vertex.label.cex=.8,vertex.label = V(g2)$label_all,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	legend(x = "topleft",#inset = 0,
        legend = c('b1,3-Galactose','b1,4-Galactose','a1,2-Fucose','a1,3-Fucose','a1,4-Fucose','b1,3-GlcNAc','b1,6-GlcNAc','a2,3-Sialic Acid','a2,6-Sialic Acid'),
        lty=c(1,2,1,2,3,1,2,1,2),
        col=c('gold3','gold3','red','red','red','blue','blue','purple','purple'), lwd=3) #, cex=.5) #, horiz = TRUE)

	# lay2[,2] = lay2[,2]+runif(nrow(lay2))/2
	lay2 = layout_as_tree(g2)
	idx=c(3,4,5,9,10,11,13,14)
	lay2[idx,2]=lay2[idx,2]+c(0,0,-.25,.25,-.5,-.25,-.25,.25)
	plot(g3,layout=lay2, #lay[idx,],
		edge.width=E(g3)$weight*5,edge.arrow.size=0,edge.color=ifelse(E(g3)$type=='+',E(g3)$color,NA),
		vertex.size=10,vertex.label.cex=.8,vertex.label = V(g3)$label_g3,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	plot(g3,layout=lay2, #lay[idx,],
		edge.width=E(g3)$weight*5,edge.arrow.size=0,edge.color=ifelse(E(g3)$type=='-',E(g3)$color,NA),
		vertex.size=10,vertex.label.cex=.8,vertex.label = V(g3)$label_g3,vertex.label.family='sans',
		vertex.color=NA,vertex.frame.color=NA)
	legend(x = "bottomright",#inset = 0,
        legend = levels(cut(cri[,3],breaks=c(-1,-.5,-.25,-.1,0,.1,.25,.5,1))),
        col=matlab.like2(8), lwd=5) #, cex=.5, horiz = TRUE)
	dev.off()

	# g4 = delete_vertices(g3,which(!is.na(V(g3)$label_g3)))
	# plot(g3,#layout=lay2, #lay[idx,],
	# 	edge.width=E(g3)$weight*5,edge.arrow.size=0,edge.color=ifelse(E(g3)$type %in% c('+','-'),E(g3)$color,NA),
	# 	vertex.size=10,vertex.label.cex=.8,vertex.label = V(g4)$label_g3,vertex.label.family='sans',
	# 	vertex.color=NA,vertex.frame.color=NA)
	gi = graph_from_edgelist(as.matrix(cr[idx<-(cr[,3]<1&abs(cr[,3])>.25),1:2]))
	E(gi)$weight=cr[idx,3]
	E(gi)$color = matlab.like2(8)[cut(cr[idx,3],breaks=c(-1,-.5,-.25,-.1,0,.1,.25,.5,1))]
	V(gi)$weight=strength(gi)
	plot(gi,edge.arrow.size=0,layout=layout_with_fr(gi),#vertex.shape='oval',
		vertex.size=20,vertex.label.cex=.8,vertex.color='grey',vertex.frame.color=NA)



	# tmp =droplevels(na.omit(m2[m2$cor<1 & !m2$child%in%c('FLNH','FDSLNH','2FL','3FL','3SL') & !m2$parent%in%c('FLNH','FDSLNH','2FL','3FL','3SL'),
	# 	c('cor','d_compete','parent','child','max_depth','from_LNT.x','from_LNT.y','depth.y','depth.x')]))
	# anova(from_LNT.x ~ from_LNT.y , data=d_compete.agg)


	# # cluster specific networks
	# pdf('figures/networks.clusters.pdf',height=30,width=15)
	# par(mfrow=c(5,2))
	# for(i in 1:5){
	# 	g = graph_from_edgelist(as.matrix(el_na[,1:2]))
	# 	E(g)$weight = unlist(weighted_mean_net[i,-c(1,1+as.numeric(attr(el_na,'na.action')))])
	# 	E(g)$color = rainbow(10)[factor(el_na[,3])]
	# 	#E(g)$weight = ifelse(E(g)$weight==0,NA,E(g)$weight)
	# 	V(g)$weight = strength(g)/degree(g)
	# 	tmp = sapply(shortest_paths(g,from=V(g)[1])$vpath,length)
	# 	tmp[1] = 1
	# 	V(g)$depth = tmp
	# 	V(g)$depth_norm = tmp/max(tmp)
	# 	V(g)$weight_norm = V(g)$weight/V(g)$depth_norm
	# 	if(!'lay'%in%objects()){
	# 		lay = layout_as_tree(g)
	# 		#lay[,2] = lay[,2]+runif(nrow(lay))/2
	# 		lay[,1] = lay[,1]*5
	# 	}
	# 	plot(g,layout=lay,
	# 		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color='grey',
	# 		vertex.size=5,vertex.label.cex=.5,
	# 		#vertex.size=V(g)$weight,
	# 		vertex.color='white',vertex.frame.color=NA)
	# 		#vertex.label.cex=exp(V(g)$weight)/5)
	# 	plot(g,layout=lay,
	# 		edge.width=E(g)$weight*5,edge.arrow.size=0,edge.color=E(g)$color,
	# 		vertex.size=5,vertex.label.cex=.5,
	# 		#vertex.size=V(g)$weight,
	# 		vertex.color='white',vertex.frame.color=NA)
	# 		#vertex.label.cex=exp(V(g)$weight)/5)

	# 	# E(g)$weight = median_net[i,-c(1,1+as.numeric(attr(el_na,'na.action')))]
	# 	# E(g)$weight = ifelse(E(g)$weight==0,NA,E(g)$weight)
	# 	# plot(g,layout=layout_as_tree,
	# 	# 	edge.width=E(g)$weight,edge.arrow.size=0,
	# 	# 	vertex.label.cex=.2,vertex.color=NA,vertex.frame.color=NA)
	# }
	# dev.off()
}

# padj = apply( Model_Scores , 2, function(x) fdrtool(pnorm(as.vector(scale(na.omit(x))),lower.tail=F),statistic='pvalue' ))


if(new_sels){	
	data12_secretor12 = which(!is.na(Model_Scores[,1]) & Model_Score_Groups[,1]==1 & Model_Score_Groups[,2]==1 &
		Model_Score_Groups[,3]==1 & Model_Score_Groups[,4]==1)
	
	rxn_names = gsub('_rxn','',annotate_inclusion())

	inc_rxns = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')[idx<-!(grepl('HMO',rxn_names)),data12_secretor12]

	rxn_fact  =  factor( gsub('[1-9]','',rxn_names[idx]))
	
	pdf('figures/data12_secretor12/heatmaps.pdf')
	heatmap.2(inc_rxns,trace='none',labRow=rxn_names[idx],mar=c(5,8),RowSideColors=rainbow(5)[rxn_fact],distfun=function(x) dist(x,method='binary'))
	heatmap.2(inc_rxns[idx2<-rowSums(inc_rxns)>0,],trace='none',labRow=(rxn_names[idx])[idx2],mar=c(5,8),RowSideColors=(rainbow(5)[rxn_fact])[idx2],distfun=function(x) dist(x,method='binary'))
	dev.off()
	
	rownames(inc_rxns) = rxn_names[idx]
	m = melt( apply(inc_rxns,1,function(x) sum(x)/length(x)) )
	m$HMOi = rownames(m)
	m$i = gsub('[A-Z]','',m$HMOi)
	m$HMO = factor(gsub('[1-9]','',m$HMOi))

	g=ggplot(data=m,aes(x=i,y=value))+geom_bar(stat='identity') + facet_wrap(~HMO,scale='free_x')
	ggsave(plot=g,filename='figures/data12_secretor12/proportions.pdf')
}

if(gene_v_str){
	sels=list(data1_selection = Model_Score_Groups[,1]==1,
		data2_selection = Model_Score_Groups[,2]==1,
		secretor1_selection = Model_Score_Groups[,3]==1,
		secretor2_selection = Model_Score_Groups[,4]==1)
	sels$data12_inter = sels$data1_selection & sels$data2_selection
	sels$secretor12_inter = sels$secretor1_selection & sels$secretor2_selection
	inclusion = attach.big.matrix('../HMO_biosynthesis/inclusion.rxn.desc')
	idx=!(grepl('HMO',rxn_names <- gsub('_rxn','',annotate_inclusion())))

	out=list()
	out2=list()
	for(n in names(sels)){
		sel = sels[[n]]
		for(gli in 1:2){
			gli=ifelse(grepl('secretor',n),gli+2,gli)
			gene_score = attach.big.matrix( Gene_Linkage_Scores[[gli]])[sel,]
			incl_i = inclusion[idx,sel]
			for(i in 1:ncol(gene_score)){
				for(j in 1:nrow(incl_i)){
					included = incl_i[j,] == 1
					tmp=gene_score[included,i] 
					out[[paste(n,DS[gli],gl_names[i],(rxn_names[idx])[j],sep='.')]] = c(median(tmp),mean(tmp),sd(tmp))
					out2[[paste(n,DS[gli],gl_names[i],(rxn_names[idx])[j],sep='.')]] = tmp
				}
			}
		}
	}
	tmp=as.data.frame( do.call(rbind,out))
	tmp = cbind(tmp, do.call(rbind, strsplit(rownames(tmp),split='\\.') ) )
	colnames(tmp) = c('median','mean','sd','selection','DS','gene','HMOi')
	tmp$links = unlist(lapply(strsplit( as.character(tmp$gene) , '_' ) , function(x) paste(x[1:2],collapse='_')))
	tmp$genei = unlist(lapply(strsplit( as.character(tmp$gene) , '_' ) , function(x) x[3]))
	tmp$HMO = gsub('[1-9]','',as.character(tmp$HMOi))
	tmp$i = gsub('[A-Z]','',as.character(tmp$HMOi))
	tmp$selDS = paste(tmp$DS,tmp$selection,sep='.')
	tmp = tmp[grepl('inter',tmp$selection),]
	g=grid.arrange(grobs=apply( unique(tmp[,c('selection','DS')]) , 1, function(x){
		print(x)
		keep= tmp$selection==x[1] & tmp$DS==x[2] & !grepl('L5',tmp$link)
		ggplot(data=droplevels(tmp[keep,]),aes(x=genei,y=exp(median),color=HMO))+geom_boxplot() +facet_wrap(~links,scale='free') +
			theme(axis.text.x = element_text(angle = 45, hjust = 1,size=)) + ggtitle(paste(x,collapse='.'))
	}))
	ggsave(plot=g,filename='figures/gene_v_structure.pdf',height=15,width=15)

	g=grid.arrange(grobs=apply( unique(tmp[,c('selection','HMO')]) , 1, function(x){
		print(x)
		keep= tmp$selection==x[1] & tmp$HMO==x[2] & !grepl('L5',tmp$link)
		print(sum(keep))
		ggplot(data=droplevels(tmp[keep,]),aes(x=genei,y=exp(median),color=HMOi))+geom_boxplot() + geom_jitter(aes(shape=DS)) +
			facet_wrap(~links,scale='free') +
			theme(axis.text.x = element_text(angle = 45, hjust = 1,size=)) + ggtitle(paste(x,collapse='.'))
	}),ncol=5)
	ggsave(plot=g,filename='figures/gene_v_structure.granular.pdf',height=20,width=40)

#	tmp2=as.data.frame( do.call(rbind,out2))
}

if( inclusion_rank){
	rnk = lapply(DS,function(ds){
		mat = attach.big.matrix(paste0('best_gene_order.',ds,'.desc'))
		tmp = apply( mat, 2, function(x) table(factor(x[Model_Score_Groups[,DS==ds]==1],levels=as.character(1:5))) )		
		colnames(tmp) = gl_names
		tmp = melt(tmp)
		colnames(tmp) = c('rank','gene_link','count')
		tmp$ds = ds
		tmp
	})
	dat = do.call(rbind,rnk)
	dat$link = unlist( lapply(strsplit( as.character(dat$gene_link) , '_') , function(x) paste(x[1:2],collapse='_') ) )
	dat$gene = unlist( lapply(strsplit( as.character(dat$gene_link) , '_') , function(x) x[3] ))
	dat$dsi = ifelse(grepl('1',dat$ds),'ds1','ds2')
	dat$sec = ifelse(grepl('secretor',dat$ds),'secretor','all')
	dat$rank = as.numeric(as.character(dat$rank))
	#dat = dat[dat$count>0,]
	grid.arrange(grobs=lapply(unique(dat$link)[-5],function(x){
		ggplot(data=dat[dat$link==x,],aes(x=rank,y=count,color=gene)) + geom_line() + facet_grid(sec ~ dsi) + scale_y_log10()+ggtitle(x)
	}))
	

	heatmap.2( t(log(tmp+1,10)) , trace='none',mar=c(5,12),Colv=F)

}

}



	min_max_norm = function(x) (x-min(x))/(max(x)-min(x))
	scale_p = function(x) pnorm(scale(x),lower.tail=F)
	load('data_out/Gene_Choice.rda')
	################
	enrich_d1 = Gene_Choice[[1]][,c('rxn_name','selection_inter_size','selection_inter_success','enrichment_inter')]
	rownames(enrich_d1) = enrich_d1$rxn_name
	enrich_d1$PROP = enrich_d1$selection_inter_success / enrich_d1$selection_inter_size
	
	enrich_d2 = Gene_Choice[[2]][,c('rxn_name','selection_inter_size','selection_inter_success','enrichment_inter')]
	rownames(enrich_d2) = enrich_d2$rxn_name
	enrich_d2$PROP = enrich_d2$selection_inter_success / enrich_d2$selection_inter_size
	
	data_all_1 = data.frame( enrich_d1[,c('PROP','enrichment_inter','selection_inter_size','selection_inter_success')] )
	data_all_2 = data.frame( enrich_d2[,c('PROP','enrichment_inter','selection_inter_size','selection_inter_success')])

	data_all_1r = data_all_1
	data_all_2r = data_all_2

	data_all_1b = data.frame( scale(enrich_d1[,c('PROP','enrichment_inter','selection_inter_size','selection_inter_success')] ))
	data_all_2b = data.frame( scale(enrich_d2[,c('PROP','enrichment_inter','selection_inter_size','selection_inter_success')]))

	# Get Commonly Selected Models
	########### load
	Model_Score_Groups = attach.big.matrix('Model_Score_Groups.desc')
	#####################
	common = which( Model_Score_Groups[,1] & Model_Score_Groups[,2] )

	# Gene Linkage Score (min-max normalized median GLS)
	###########load 
	Gene_Linkage_Scores = lapply( DS,function(x) dget( paste0("gene_linkage_score.",x,".desc") ) ) 
	#################
	gls_d1 = attach.big.matrix(Gene_Linkage_Scores[[1]])[common,]
	tmp = apply(gls_d1,2,median,na.rm=T)
	data_all_1r$GLS = tmp
	data_all_1$GLS = min_max_norm(tmp)
	data_all_1b$GLS = scale(tmp)

	gls_d2 = attach.big.matrix(Gene_Linkage_Scores[[2]])[common,]
	tmp = apply(gls_d2,2,median,na.rm=T)
	data_all_2r$GLS = tmp
	data_all_2$GLS = min_max_norm(tmp)
	data_all_2b$GLS = scale(tmp)

	# Model Score Contribution: correlation of model score to GLS in commonly selected models (max-min normalized)
	##### load
	Model_Scores = attach.big.matrix('Model_Scores.desc')[common,1:2]
	##############
	tmp = apply(gls_d1,2,function(x){
		msx = na.omit(cbind(Model_Scores[,1] , x))
		cor( msx[,1],msx[,2], method='pearson')
		})
	data_all_1r$MSC = tmp
	data_all_1$MSC = min_max_norm(tmp)
	data_all_1b$MSC = scale(tmp)

	tmp = apply(gls_d2,2,function(x){
		msx = na.omit(cbind(Model_Scores[,2] , x))
		cor( msx[,1],msx[,2], method='pearson')
		})
	data_all_2r$MSC = tmp
	data_all_2$MSC = min_max_norm(tmp)
	data_all_2b$MSC = scale(tmp)

	# Name
	data_all_1$dataset = 'd1'
	data_all_2$dataset = 'd2'
	data_all_1b$dataset = 'd1'
	data_all_2b$dataset = 'd2'
	data_all_1r$dataset = 'd1'
	data_all_2r$dataset = 'd2'

	labs = data.frame( do.call(rbind,strsplit(gsub('_gene','',rownames(data_all_1)),'_')) )
	colnames(labs) = c('Linkage','Link','Gene')

	data_all = rbind(data_all_1,data_all_2)
	data_all = cbind(data_all,rbind(labs,labs))

	data_allb = rbind(data_all_1b,data_all_2b)
	data_allb = cbind(data_allb,rbind(labs,labs))

	data_allr = rbind(data_all_1r,data_all_2r)
	data_allr = cbind(data_allr,rbind(labs,labs))

	save(data_all,data_allb,data_allr,file='data_out/Gene_Choice_aggstats.rda')
