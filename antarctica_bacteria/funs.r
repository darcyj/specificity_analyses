# other way around:
gg2taxtable <- function(gg, ids=NULL, nlvl=0, delim="; "){

	sp <- lapply(X=gg, FUN=function(x){unlist(strsplit(x, split=delim))})
	# check all have same number of levels
	lens <- sapply(X=sp, FUN=length )
	# estimate nlvl if nlvl=0 (auto)
	if(nlvl == 0){
		# an overly verbose way to just get the most common value!
		nlvl <- as.numeric(names(rev(sort(table(lens))))[1])
	}

	if(sum(lens == nlvl) == 0){
		stop("No tax strings had specified nlvl.")
	}else if(! all(lens == nlvl)){
		n_mism <- sum(lens != nlvl)
		warning(paste0(n_mism, "/", length(lens), " tax strings had wrong nlvl. Changing them to unknown."))
		# make unknown
		unk <- sp[[which(lens == nlvl)[1]]]
		unk <- paste0(substr(unk, 1, 1), "__unknown")
		# replace nlvl mismatches with unk - now all will have same nlvl
		sp <- lapply(X=1:length(lens), FUN=function(i){
			if(lens[i] != nlvl){
				return(unk)
			}else{
				return(sp[[i]])
			}
		})
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
	# replace blank entries with "unknown"
	sp[sp == ""] <- "unknown"
	return(sp)
}


agg_spec_by_tax <- function(spec, tax, lvl, min_species=0){
	# make unknown in case we need it
	unk <- tax[1, ,drop=FALSE]
	unk[1,] <- rep("unknown", nrow(unk))

	# fetch tax line for each species in spec, make 'em into a matrix
	tax2 <- t(sapply(
		X=rownames(spec),
		FUN=function(rn){
			if(sum(tax[,1] == rn) == 1){
				return(tax[tax[,1] == rn])
			}else if(! rn %in% tax[,1]){
				warning(paste0("species \"", rn, "\" not found in tax, using \"unknown\"."))
				return(unk)
			}else{
				stop(paste0("species \"", rn, "\" in tax multiple times!"))
			}
		}
	))
	colnames(tax2) <- colnames(tax)

	# get the column we're interested in from tax2
	if(sum(colnames(tax2) == lvl) == 1){ # if lvl is in tax2's colnames exactly once:
		colind <- which(colnames(tax2) == lvl)
	}else if(sum(startsWith(colnames(tax2), lvl)) == 1){ # if exactly one tax2 colname starts with lvl:
		colind <- which(startsWith(colnames(tax2), "c"))
	}else{
		stop(paste0("Taxonomic lvl \"", lvl, "\" not in tax colnames."))
	}
	lvl_vec <- tax2[,colind]
	lvl_uniq <- unique(lvl_vec)

	# do aggregation
	output <- lapply(
		X=lvl_uniq,
		FUN=function(x){ spec[lvl_vec == x, , drop=FALSE] }
	)
	names(output) <- lvl_uniq

	# get number of species per taxon
	nspec <- lapply(X=output, FUN=nrow)

	# remove taxa with too few species
	output <- output[nspec >= min_species]

	return(output)
}

	replace_tip_with_polytomy <- function(tree, tipname, newtips){
		require(ape)
		# which node number in the main tree? node number is just the index of tree$tip.label.
		node_ind <- which(tree$tip.label == tipname)
		# polytomous tree for newtips
		subtree_i <- read.tree(text=paste(
			"(", 
			paste( paste(newtips, ":0", sep=""), collapse=","),
			"):0;"
		))
		# add them to tree
		return( bind.tree(tree, subtree_i, where=node_ind) )
	}