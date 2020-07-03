# install.packages(c('ggplot2','gridExtra','R.matlab','dplyr','broom','mixtools','diptest','Matrix','cluster','DeducerExtras','bigmemory','bigalgebra',
# 	'biganalytics','dummies','doParallel','foreach','igraph','logisticPCA','ProjectionBasedClustering','factoextra','mclust','useful','MASS',
# 	'stringr','colorRamps','ggrepel','GGally','metap','bigalgebra','gplots','ggExtra'))
#source("https://bioconductor.org/biocLite.R")
#biocLite("multtest")
#biocLite("Biostring")


library(ggplot2)
library(gridExtra)
library(reshape)
library(R.matlab)
library(ggplot2)
library(dplyr)
library(broom)
library(mixtools)
library(diptest)
library(Matrix)
library(cluster)
#library(DeducerExtras)

library(bigmemory)
library(biganalytics)
#library(bigalgebra)
library(dummies)

library(gridExtra)

library(gplots)

library(doParallel)
library(foreach)
library(igraph)

library(ggplot2)
#library(ggExtra)
#library("ggpubr")
#library(xlsx)
library(gridExtra)
library(logisticPCA)
#library(ProjectionBasedClustering)
library(factoextra)
library(mclust)
library(useful)
library(MASS)
library(stringr)
library(colorRamps)
library(ggrepel)
library(GGally)
library(metap)
library(dplyr)
library(ggExtra)
library(ggpubr)
library(Biostrings)
	library(mclust)
	library(useful)
	library(MASS)
	library(stringr)
	library(colorRamps)
	library(ggrepel)
	library(GGally)
	library(gplots)
	library(ggExtra)



#library(parallel)

# 3/26/19 removed fut2 from a13 and a14, rm Fut8 bc it's an a1,6 fucose transferase
load_linkages_bad <- function(){
	list(	L1_b3GnT_gene =	c('B3GNT2','B3GNT3','B3GNT4','B3GNT8','B3GNTL1'),
		L2_a2FucT_gene=	c('POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b'),
		L3_a3FucT_gene=	c('POFUT1b','POFUT2a','POFUT2b','FUT11','FUT3','FUT4','FUT6b'),
		L4_ST3GalT_gene=c('ST3GAL1a','ST3GAL2','ST3GAL3b','ST3GAL4','ST3GAL5b'),
		L5_ST6Gal_gene =c('ST6GAL1b'),
		L6_b3GalT_gene =c('B3GALNT2','B3GALT4','B3GALT5','B3GALT6'),
		L7_b4GalT_gene =c('B4GALNT4','B4GALT1','B4GALT2a','B4GALT3','B4GALT4a','B4GALT6','B4GALT7'),
		L8_b6GnT_gene =c('GCNT1','GCNT2c','GCNT3'),
		L9_a4FucT_gene =c('POFUT1b','POFUT2a','POFUT2b','FUT3','FUT6b'),
		L10_ST6GnT_gene =c('ST6GALNAC1','ST6GALNAC2','ST6GALNAC4a','ST6GALNAC6'))
}

load_linkages<-function(){
	list(	L1_b3GnT_gene =	c('B3GNT2','B3GNT3','B3GNT4','B3GNT8','B3GNTL1'),
		L2_a2FucT_gene=	c('POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'),
		L3_a3FucT_gene=	c('POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'),
		L4_ST3GalT_gene=c('ST3GAL1a','ST3GAL2','ST3GAL3b','ST3GAL4','ST3GAL5b'),
		L5_ST6Gal_gene =c('ST6GAL1b'),
		L6_b3GalT_gene =c('B3GALNT2','B3GALT4','B3GALT5','B3GALT6'),
		L7_b4GalT_gene =c('B4GALNT4','B4GALT1','B4GALT2a','B4GALT3','B4GALT4a','B4GALT6','B4GALT7'),
		L8_b6GnT_gene =c('GCNT1','GCNT2c','GCNT3'),
		L9_a4FucT_gene =c('POFUT1b','POFUT2a','POFUT2b','FUT11','FUT2','FUT3','FUT4','FUT6b','FUT8'),
		L10_ST6GnT_gene =c('ST6GALNAC1','ST6GALNAC2','ST6GALNAC4a','ST6GALNAC6'))
}

# b3gnt6 (zero expression in rnaseq, remove) and b4galt5 failed during flux correlation and were therefore excluded

linkage_update<-function(linkage_orig){
	lapply(names(linkage_orig),function(xn){
	        x = linkage_orig[[xn]]
	        if(grepl('L3_',xn)){
	                rm_genes_i = c(rm_genes,"FUT2","FUT8")
	        }else if(grepl('L2_',xn)){
	                rm_genes_i = c(rm_genes,'FUT11','FUT3','FUT4','FUT6b','FUT8')
	        }else if(grepl('L9_',xn)){
	                rm_genes_i = c(rm_genes,'FUT11',"FUT2",'FUT4','FUT11','FUT6b')
	        }else if(grepl('L1_',xn)){
	                rm_genes_i = c(rm_genes,'B3GNT4','B3GNTL1')
	        }else if(grepl('L4_',xn)){
	                rm_genes_i = c(rm_genes,'ST3GAL2','ST3GAL4')
	        }else if(grepl('L6_',xn)){
	                rm_genes_i = c(rm_genes,'B3GALNT2','B3GALT6') # B3GALT5 is borderline, should it be excluded?
	        }else if(grepl('L7_',xn)){
	                rm_genes_i = c(rm_genes,'B4GALNT4','B4GALT6','B4GALT7')
	        }else if(grepl('L10_',xn)){
	                rm_genes_i = c(rm_genes,'ST6GALNAC1')
	        }else{
	                rm_genes_i = rm_genes
	        }
	        x[!x%in%rm_genes_i]
	} )
}


### write multiple sheets


## cooccurance heatmap

jaccard = function(x,y=NULL){
	if(is.null(y)){ sum(x,na.omit=T)/length(na.omit(x)) }
	length(intersect(which(x==1),which(y==1)))/length(union(which(x==1),which(y==1)))
}
cocount = function(x,y=NULL){
	if(is.null(y)){ sum(x,na.rm=T) }
	sum(x+y==2,na.rm=T)
}
coocur = function(X,f){
	m = matrix(0,nrow(X),nrow(X),dimnames=list(rownames(X),rownames(X)))
	for(i in rownames(X)){
		for(j in rownames(X)){
			if(i==j){
				m[i,j] = f(X[i,])
			}else{
				m[i,j] = f(X[i,],X[j,])
			}
		}
	}
	m
}

get_idx = function(X,i,dim_i=1){
	if(dim_i==1){
		X[i,]
	}else if(dim_i==2){
		X[,i]
	}
}

coocur_par = function(mat,idx,sub_idx=T,f=jaccard,names=NULL,dim_i=1){
	save(mat,get_idx,idx,f,names,dim_i,file='tmp.rda')
	m = foreach(i=idx,.combine=rbind,.packages=c('bigmemory')) %dopar% {
		load('tmp.rda')
		X = attach.big.matrix(mat)
		sapply(idx,function(j) {
			if(i==j){
				f(get_idx(X,i,dim_i)[sub_idx])
			}else{
				f(get_idx(X,i,dim_i)[sub_idx],get_idx(X,j,dim_i)[sub_idx])
			}
		})
	}
	dimnames(m) = list(names,names)
	m
}

###

select_save<-function(data,plot,filename,localdir,outputdir='output',...){
	ggsave(plot=plot,filename=file.path(localdir,paste(filename,'pdf',sep='.')),...)
	ggsave(plot=plot,filename=file.path(outputdir,paste(filename,'pdf',sep='.')),...)
	save(data,	 file=file.path(outputdir,paste(filename,'rda',sep='.')))
}

##
t_enriched_group = function( group , value , ... ){
	data.frame( gene=unique(group) , do.call(rbind,lapply( unique(group) , function(i){ 
		#if( sd(value[group==i]) ==0 | sum(group==i)==1){
#			data.frame(estimate=((e1<-mean(value[group==i]))-(e2<-mean(value[group!=i])
#			tidy(t.test(value[group!=i],mu=mean(value[group==i]) ,alternative ='less',...))
#			tidy(t.test(value[group==i],mu=mean(value[group!=i]) ,alternative ='greater',...)) 
		#}else{
#			tidy(t.test(value~(group!=i) ,alternative ='greater',...))# currently using this on ##### this one is currently used
#			tidy(t.test(value[group==i],value[group!=i] ,alternative ='greater',var.equal=T,...))
#			tidy(t.test(value[group==i],mu=mean(value[group!=i]) ,alternative ='greater',...))
			c(estimate=(mu1<-mean(value[group==i]))-(mu2<-mean(value[group!=i])),
				statistic=(z<-(mu1-mean(value))/sd(value)), p.value=pnorm(z,lower.tail=F)) 
		#}
		} ) ) )
}
##

error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
stop("vectors must be same length")
arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
}

### get residuals
#k = 1000
get_sigmoid_bounds<-function(cluster_score,k,n=1,col=heat.colors(3),roll=10){
	#tmpi = unique( tmpi[,c('k','cluster_score','cut')] )
	#df = data.frame( score = sort(tmpi$cluster_score) , i = 1:nrow(tmpi))
	df = data.frame( score = sort(cluster_score) , i = 1:length(cluster_score))  
	if(nrow(df)>k){stop('mismatched rows')}
	#df = data.frame( score = c(rnorm(round(k/2)),rnorm(round(k/4),mean=-1),rnorm(round(k/4),mean=1)) , i = 1:k)
	# extract central scores 
	md = median(df$score); sd = sd(df$score)
	keep= df$score>(md-n*sd) & df$score<(md+n*sd)
	# construct model based on central scores
	mod = glm( score ~ i , data=df[keep,])
	# predict all scores using model
	df$pred = predict( mod , df )
	# calculate residuals and standard residuals based on model residuals
	df$res = df$score - df$pred
	df$std.res = (df$res - 0)/sd(mod$residuals)
	# determine percentiles for scores
	d = density(df$score)
	d$y = d$y/sum(d$y) # normalize density
	df$perc = sapply(df$score , function(x) sum(d$y[1:which.min(abs(x-d$x))]))
	group_tmp = ifelse( df$std.res < -1.96 & df$score < md , "Low" ,ifelse( df$std.res > 1.96 & df$score > md, "High", "Average" ))
	rng = range(df$score[group_tmp=='Average'])
	df$group = ifelse( df$score < min(rng) , "Low" ,ifelse( df$score >max(rng), "High", "Average" ))
#	for(i in roll:(nrow(df)-roll)){
#		tab = table(df$group[(i-roll):(i+roll)])
#		if( length(tab)>1 ) {
#		  max_lab = names(tab)[which.max(tab)]
#		  if((max_lab == 'High'    & df$score[i] > (md+n*sd)) |
#		     (max_lab == 'Low'     & df$score[i] < (md-n*sd)) | 
#		     (max_lab == 'Average' & df$score[i] > (md-(2*n)*sd)) & df$score[i] < (md+(2*n)*sd) )
#		    df$group[i] = max_lab
#		}
#	}
		plot(df$pred,df$score,col=col[factor(df$group)])
		rg=range( c(df$pred,df$score) )
		lines(rg,rg)
	# get percentile range for middle group
	#ret= range(df$perc[df$group=='Average'])
	return(df$group)
}
###

load_inclusion<-function(file='data/models/Network_Red.expanded.mat'){
	n=211; m=15e4
	# load models
	if(is.na(file)){
		Matrix( sample(0:1,n*m,replace=T) , n,m,sparse=T)
	}else{
		Matrix(readMat(file)[[1]][[2]],sparse=T)
		#Matrix(as.matrix(readMat(file)[[1]][[2]]),sparse=T)
	}
}

annotate_inclusion<-function(inclusion=NULL){
	# load network
	net = readMat('data/models/Network_Red.expanded.mat')
#	net = readMat('data/models/Network_Red.mat')
	red_net = as.matrix(net[[1]][[12]])
	el = t(apply(red_net,2,function(x) c(ifelse(-1%in%x,which(x==-1),NA),ifelse(1%in%x,which(x==1),NA))))
	rownames(el) = unlist(net[[1]][[4]])
	mets = unlist(net[[1]][[1]])
	el[,1] = mets[as.numeric(el[,1])]
	el[,2] = mets[as.numeric(el[,2])]
	# combine
	if(is.null(inclusion)){return(rownames(el))}
	rownames(inclusion) = rownames(el)
	colnames(inclusion) = paste0('model',1:ncol(inclusion))
	return(inclusion)
}

get_mets<-function(){
	net = readMat('data/models/Network_Red.expanded.mat')
	mets = unlist(net[[1]][[1]])
	return(mets)
}

get_network_edgelist<-function(names=TRUE){
		# load network
	net = readMat('data/models/Network_Red.expanded.mat')
#	net = readMat('data/models/Network_Red.mat')
	red_net = as.matrix(net[[1]][[12]])
	el = t(apply(red_net,2,function(x) c(ifelse(-1%in%x,which(x==-1),NA),ifelse(1%in%x,which(x==1),NA))))
	rownames(el) = unlist(net[[1]][[4]])
	if(names){
		mets = unlist(net[[1]][[1]])
		el[,1] = mets[as.numeric(el[,1])]
		el[,2] = mets[as.numeric(el[,2])]		
	}
	return(el)
}

get_rxn <- function(){
	library(Biostrings)
	net = readMat('data/models/Network_Red.expanded.mat')
#	net = readMat('data/models/Network_Red.mat')
	red_net = as.matrix(net[[1]][[12]])
	el = t(apply(red_net,2,function(x) c(ifelse(-1%in%x,which(x==-1),NA),ifelse(1%in%x,which(x==1),NA))))
	# name vertexes
	rownames(el) = unlist(net[[1]][[4]])
	mets = unlist(net[[1]][[1]])
	el[,1] = mets[as.numeric(el[,1])]
	el[,2] = mets[as.numeric(el[,2])]
	rownames(el) = unlist(net[[1]][[4]])
	# align metabolite names
	el1 = gsub(')|\\(','',unlist(lapply(strsplit(el[,1],';'),function(x) x[1])))
	el2 = gsub(')|\\(','',unlist(lapply(strsplit(el[,2],';'),function(x) x[1])))
	removeRS <- function(str){
		strc = rle(strsplit(str, "")[[1]])
		strc$lengths[which(strc$values=='-')] = 1
		paste(inverse.rle(strc),collapse="")
	}
	dif = apply(cbind(el1,el2,rownames(el)),1,function(x){
		if(any(is.na(x))){return(NA)}
		if(!grepl('Syn_',x[3])){return('out')}
		tmp=pairwiseAlignment(x[1],x[2],type='local',gapOpening=.5, gapExtension=0)
		gsub(gsub('-','|',removeRS(as.character(tmp@pattern))),'',x[2])
		})
	m=cbind(el1,el2,dif)
	return(dif)
}

get_sel<-function(prefix,Model_Score_Groups){
	l=list(
		d1 = which(Model_Score_Groups[,grepl(paste0(prefix,'_data1'),colnames(Model_Score_Groups))]=='High'),
		d2 = which(Model_Score_Groups[,grepl(paste0(prefix,'_data2'),colnames(Model_Score_Groups))]=='High'),
		inter=intersect(d1,d2)
	)
	maxlen = max(unlist(lapply(l,length)))
	l = lapply(l,function(x) c(x,rep(NA,maxlen-length(x))))
	sel = data.frame(common=l$inter,data1=l$d1,data2=l$d2)
}

invert_ds<-function(dataset) ifelse( grepl('1',dataset), gsub('1','2',dataset),gsub('2','1',dataset))

MAD = function(x,md=NULL) median(abs(x-ifelse(is.null(md),median(x,na.rm=T),md)),na.rm=T)

Iglewicz.Hoaglin.outlier = function(x) {
	data.frame( modifiedz=(tmp<-(x-median(x))/MAD(x)) , scaled_mz = (tmp<-0.6745*tmp) , outlier = abs(tmp)>3.5 )
}

g_legend<-function(a.gplot){
    tmp <- ggplot_gtable(ggplot_build(a.gplot))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    legend
}
