check_inclusion<-function(i1,i2){
	if(!all(i1==i2)){stop('The inclusion matrix has been changed')}
}

check_cluster<-function(hc){
	K_rec=unique(apply(cut<-cutree(hc,h=seq(0,1,0.05)),2,function(x) length(unique(x)) ))
	K = max(K_rec)
	ci = cutree(hc,k=K)

	tmp=sapply(unique(ci),function(j){
	if(sum(ci==j)>1){
		apply( inclusion[select_reactions,ci==j] ,1,function(x) all(x)|(!any(x)))
	}else{
		rep(T,sum(select_reactions))
	}
	})
	if(!all(tmp)){stop('Some clusters contain more than one structure combination')}
}

check_correlation_score<-function(cs){
	if( any(is.na(cs)) | any(!is.finite(cs)) ) {stop('there are na or infinate values in ri')}
}
