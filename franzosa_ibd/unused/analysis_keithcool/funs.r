# other way around:
gg2taxtable <- function(gg, ids=NULL){
	sp <- lapply(X=gg, FUN=function(x){unlist(strsplit(x, split=";"))})
	# check all have same number of levels
	lens <- sapply(X=sp, FUN=length )
	if(! all(lens == lens[1])){
		stop("Not all tax strings have same number of levels.")
	}
	# format sp into matrix
	sp <- t(simplify2array(sp))
	# check column coherence, assign column labels
	cn <- rep("", ncol(sp))
	for(j in 1:ncol(sp)){
		prfxs <- substr(sp[,j], 1, 1)
		if(all(prfxs == prfxs[1])){
			cn[j] <- prfxs[1]
		}else{
			stop(paste("Prefix mismatch in column", j))
		}
		# get rid of prefixes
		sp[,j] <- sub(pattern=paste0("^", prfxs[1], "__"), x=sp[,j], replacement="")
	}
	colnames(sp) <- cn
	# add in is column if provided
	if(!is.null(ids)){
		if( length(ids) != nrow(sp)){
			stop("Provided ids are wrong length")
		}else{
			sp <- cbind(id=ids, sp)
		}
		
	}
	return(sp)
}

