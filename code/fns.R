theme_set(theme_bw())
knitr::opts_chunk$set(fig.height = 8)

### My random keys

G2M_acro = 23.92278
G2M_acro_sd = 1.488121
G2M_acro_lowlim = 23.92278 - 1.488121
G2M_acro_upplim = (23.92278 + 1.488121) - 24

S_acro =16.14223
S_acro_sd = 0.390664
S_acro_lowlim = 16.14223 - 0.390664
S_acro_upplim = (16.14223 + 0.390664)

## Set 2
# G2M_acro = 19.76729
# G2M_acro_sd = 1.488121
# G2M_acro1 = 19.76729 - 1.488121
# G2M_acro2 = 19.76729 + 1.488121
#              
# S_acro = 12.11769
# S_acro_sd = 0.390664
# S_acro1 = 12.11769 - 0.390664
# S_acro2 = 12.11769 + 0.390664
###

colkey <- viridis::viridis(20)
# tod_to_zt <- data.frame(time_of_day=seq(0,23,1), zt=c(seq(17,23,1),seq(0,16,1)))
# ct_to_zt <- data.frame(ct=seq(0,95,1), zt=rep(c(seq(16,23,1),seq(0,15,1)), rep=4))

ct2zt <- function(x){
  if(max(x) <=24){
    
    ifelse(x < 8, x+16, x-8)
    
  } else if(max(x) > 24){
    
    tmp <- data.frame(ct=0:239, zt=rep(c(16:23, 0:15), 10))
    sapply(1:length(x), function(y) tmp[tmp$ct == x[y],]$zt)
    
  }
}

ct2zt_e3e4 <- function(x){
  if(max(x) <=24){
    
    ifelse(x < 12, x+12, x-12)
    
  } else if(max(x) > 24){
    
    tmp <- data.frame(ct=0:239, zt=rep(c(12:23, 0:11), 10))
    sapply(1:length(x), function(y) tmp[tmp$ct == x[y],]$zt)
    
  }
}

sincos2acro <- function(sin,cos, per=24) (atan2(sin,cos)/2/pi*per+per)%%per
sincos2amp <- function(sin,cos) 2*sqrt(sin^2+cos^2)

### Fetch pctm
fetch_pctm_pipelinerunid <- function(res, raw.table){
  res.table = res
  idx = res$modid
  raw.table = raw.table
  
  coltmp <- raw.table[raw.table@modid %in% idx,]
  
  if(is.annmatrix(coltmp)==TRUE){
    ss <- colanns(coltmp) %>%
      mutate(pipelinerunid = as.character(pipelinerunid))
    
    df1 <- cbind(rowanns(coltmp), data.table(coltmp))
    
    df1 <- melt(df1, id.vars = c(colnames(rowanns(coltmp)))) %>%
      rename(pipelinerunid = variable, pctm = value) %>%
      select(modid, pipelinerunid, pctm) %>%
      mutate(pipelinerunid = as.character(pipelinerunid)) %>%
      left_join(., ss, by="pipelinerunid") %>%
      left_join(., res.table, by="modid")
    
  }else{
    tmp.res <- rowanns(raw.table)
    tmp.res <- tmp.res[tmp.res$modid %in% idx,]
    df1 <- data.frame(pipelinerunid = names(coltmp), pctm = coltmp, tmp.res) %>%
      mutate(pipelinerunid = as.character(pipelinerunid)) %>%
      left_join(., colanns(raw.table), by="pipelinerunid")
  }
  
  return(df1)
}

fetch_pctm <- function(res, raw.table){
  res.table = res
  idx = res$modid
  raw.table = raw.table
  
  coltmp <- raw.table[raw.table@modid %in% idx,]
  
  if(is.annmatrix(coltmp)==TRUE){
    ss <- colanns(coltmp) %>%
      mutate(sampleid = as.character(sampleid))
    
    df1 <- cbind(rowanns(coltmp), data.table(coltmp))
    
    df1 <- melt(df1, id.vars = c(colnames(rowanns(coltmp)))) %>%
      rename(sampleid = variable, pctm = value) %>%
      select(modid, sampleid, pctm) %>%
      mutate(sampleid = as.character(sampleid)) %>%
      left_join(., ss, by="sampleid") %>%
      left_join(., res.table, by="modid")
    
  }else{
    tmp.res <- rowanns(raw.table)
    tmp.res <- tmp.res[tmp.res$modid %in% idx,]
    df1 <- data.frame(sampleid = names(coltmp), pctm = coltmp, tmp.res) %>%
      mutate(sampleid = as.character(sampleid)) %>%
      left_join(., colanns(raw.table), by="sampleid")
  }
  
  return(df1)
}

### fetch expression
fetch_expression <- function(res, raw.table){
  res.table = res
  idx = res$binid
  raw.table = raw.table
  
  coltmp <- raw.table[raw.table@binid %in% idx,]
  if(is.annmatrix(coltmp)==TRUE){
    ss <- colanns(coltmp) %>%
      mutate(sampleid = as.character(sampleid))
    
    df1 <- cbind(rowanns(coltmp), data.table(coltmp))
    
    df1 <- melt(df1, id.vars = c(colnames(rowanns(coltmp)))) %>%
      rename(sampleid = variable, rpm = value) %>%
      select(binid, sampleid, rpm) %>%
      mutate(sampleid = as.character(sampleid)) %>%
      left_join(., ss, by="sampleid") %>%
      left_join(., res.table, by="binid")
    
  }else{
    tmp.res <- rowanns(raw.table)
    tmp.res <- tmp.res[tmp.res$binid %in% idx,]
    df1 <- data.frame(sampleid = names(coltmp), rpm = coltmp, tmp.res) %>%
      mutate(sampleid = as.character(sampleid)) %>%
      left_join(., colanns(raw.table), by="sampleid")
  }
  
  return(df1)
}

### Harmonic stuff
getSlope <- function(zt,acrophase,amplitude){
  amplitude*cos(zt/24*2*pi-(acrophase-6)/24*2*pi)*pi*2/24
}

### Chi sqr test res grapher
chi_grapher <- function(toplot, two_by_two = FALSE){
  
  if(two_by_two == FALSE){
    xlabel = unlist(strsplit(toplot$data.name, "[(]|[)]|, "))[2]
    ylabel = unlist(strsplit(toplot$data.name, "[(]|[)]|, "))[3]
    residual <- data.frame(toplot$residuals) %>%
      rename(grp1 = 1, grp2 = 2) 
    
    contribution <- data.frame(100*toplot$residuals^2/toplot$statistic) %>%
      rename(grp1 = 1, grp2 = 2) 
    
    compounded <- data.frame(100*toplot$residuals^2/toplot$statistic) %>%
      mutate(Freq = Freq*sign(data.frame(toplot$residuals)$Freq)) %>%
      rename(grp1 = 1, grp2 = 2)
    
    count_table <- data.frame(toplot$observed) %>%
      rename(grp1 = 1, grp2 = 2)
    
  } else {
    xlabel = "timelapse\nacrophase"
    ylabel = "osc-modC\nacrophase"
    residual <- data.frame(toplot$residuals) %>%
      separate(Var1, c("grp1","grp2"), by=":") 
    
    contribution <- data.frame(100*toplot$residuals^2/toplot$statistic) %>%
      separate(Var1, c("grp1","grp2"), by=":") 
    
    compounded <- data.frame(100*toplot$residuals^2/toplot$statistic) %>%
      mutate(Freq = Freq*sign(data.frame(toplot$residuals)$Freq)) %>%
      separate(Var1, c("grp1","grp2"), by=":")
    
    count_table <- data.frame(toplot$observed) %>%
      separate(Var1, c("grp1","grp2"), by=":")
  }
  
  
  residual <- residual %>%
    ggplot(aes(x=grp1, y= grp2, fill=Freq)) +
    viridis::scale_fill_viridis()+
    geom_tile()  +
    xlab(xlabel) + ylab(ylabel) +
    geom_text(aes(grp1, grp2, label=format(Freq, digits=1), col=Freq<0)) +
    scale_color_manual(values = c("black","white")) +
    theme(text = element_text(size=12),
          title = element_text(size=8)) + 
    ggtitle(paste(xlabel,"vs." ,ylabel, "\nresiduals")) +
    theme(legend.position = "none")
  
  contribution <- contribution %>%
    mutate(Freq2 = paste0(format(round(Freq, 1), digits=2),"%")) %>%
    ggplot(aes(x=grp1, y= grp2, fill=Freq2)) +
    viridis::scale_fill_viridis()+
    geom_tile()  +
    xlab(xlabel) + ylab(ylabel) +
    geom_text(aes(grp1, grp2, label=Freq2, col=Freq<mean(Freq))) +
    scale_color_manual(values = c("black","white")) +
    theme(text = element_text(size=12),
          title = element_text(size=8)) + 
    ggtitle(paste(xlabel,"vs." ,ylabel, "\ncontributions")) +
    theme(legend.position = "none")
  
  compounded <- compounded %>%
    ggplot(aes(x=grp1, y= grp2, fill=Freq)) +
    viridis::scale_fill_viridis()+
    geom_tile()  +
    xlab(xlabel) + ylab(ylabel) +
    geom_text(aes(grp1, grp2, label=format(round(Freq, 1), digits=2), col=Freq<mean(Freq))) +
    scale_color_manual(values = c("black","white")) +
    theme(text = element_text(size=12),
          title = element_text(size=8)) + 
    ggtitle(paste(xlabel,"vs." ,ylabel, "\ncompounded")) +
    theme(legend.position = "none")
  
  count_table <- count_table %>%
    ggplot(aes(x=grp1, y= grp2)) +
    # viridis::scale_fill_viridis()+
    geom_tile(fill="grey50", col="white")  +
    xlab(xlabel) + ylab(ylabel) +
    geom_text(aes(grp1, grp2, label=format(round(Freq, 1), digits=2), col="white")) +
    scale_color_manual(values = c("black","white")) +
    theme(text = element_text(size=12),
          title = element_text(size=8)) + 
    ggtitle(paste(xlabel,"vs." ,ylabel, "\n count")) +
    theme(legend.position = "none")
  
  return(list(residual, contribution, compounded, count_table))
}

summaryPvalTable <- function(p,as_kable=FALSE) {
  ps <- list()
  ps[[1]] <- p < 0.05
  ps[[2]] <- p < 0.01
  ps[[3]] <- p.adjust(p, method="fdr") < 0.05
  ps[[4]] <- p.adjust(p, method="bonferroni") < 0.05
  
  levs <- c("p < 0.05", "p < 0.01", "FDR", "Bonferroni")
  perc <- sapply(ps, function(p) round(mean(p, na.rm=TRUE), 4)) * 100
  sums <- sapply(ps, function(p) sum(p, na.rm=TRUE))
  
  ret <- data.frame(level=levs, count=sums, percent=paste0(perc, "%"),
                    stringsAsFactors=FALSE
  )
  
  if(as_kable) ret <- knitr::kable(ret)
  return(ret)
}

### Extending rsrc functions to groups
# getOutliersPCA
groupWiseOutliersPCA <- function(dat,addseqrun=F,...){
	
	if(addseqrun){
		grp <- paste0(dat$age,"_",dat$media,"_",dat$seqrunid)
	}else{
		grp <- paste0(dat$age,"_",dat$media)
	}
	pcares <- lapply(csplit(dat,grp),function(x){
		tmpdat <- x[apply(x,1,function(y) min(y)!=max(y)),]
		tmp <- getOutliersPCA(tmpdat,npcs=4,...)
		PCscores <- tmp$PCscores
		isoutlier <- tmp$isoutlier
		colnames(isoutlier) <- rownames(PCscores)
		merge(melt(isoutlier),melt(PCscores),by.x=c("Var1","Var2"),by.y=c("Var2","Var1"))
	}) %>%  melt(id.var=c("Var1","Var2","value.x","value.y")) %>% 
		mutate(isoutlier=value.x)
}
groupWisePCAvar <- function(dat,addseqrun=F,...){
	if(addseqrun){
		grp <- paste0(dat$age,"_",dat$media,"_",dat$seqrunid)
	}else{
		grp <- paste0(dat$age,"_",dat$media)
	}
	lapply(csplit(dat,grp),function(x){
		tmpdat <- x[apply(x,1,function(y) min(y)!=max(y)),]
		pca <- prcomp(t(tmpdat))
		data.frame(pctVar=pca$sdev/sum(pca$sdev),"PC"=factor(1:length(pca$sdev),levels = 1:length(pca$sdev)))
	}) %>% data.table::rbindlist(idcol="grp")
}
plot.groupWiseOutliersPCA <- function(pcares,dat){
	pcares %>% 	
		merge(colanns(dat),by.x="Var2",by.y="pipelinerunid") %>%
		ggplot(aes(Var1,value.y,color=isoutlier))+
		ylab("PC Score")+
		geom_point()+
		facet_wrap(~L1)
}
# getPCAEffects
groupWisePCAEffects <- function(dat,addseqrun=F,npc=10){
	
	if(addseqrun){
		grp <- paste0(dat$age,"_",dat$media,"_",dat$seqrunid)
	}else{
		grp <- paste0(dat$age,"_",dat$media)
	}
	
	pcaefflist <- lapply(csplit(dat,grp),function(x){
		tmpdat <- x[apply(x,1,function(y) min(y)!=max(y)),]
		meta <- colanns(tmpdat)
		if(addseqrun) meta <- colanns(tmpdat) %>% select(-seqrunid)
		getPCAEffects(tmpdat,meta,npc = npc)
	})
	
	melt(pcaefflist) %>%
		mutate(fdr=p.adjust(value,method="fdr")) %>% 
		ggplot(aes(Var2,Var1,size=-log10(fdr),color=value<0.05,shape=fdr<0.05))+
		geom_point()+
		facet_wrap(~L1)
}

# Stats
fastOLS <- function(X,Y,robust=F){
	if(robust){
		fits <- lm.fit(X,t(Y))
		betas <- t(fits$coefficients)
		residuals <- t(fits$residuals)
	}else{
		betas <- solve(t(X) %*% X) %*% t(X) %*% t(Y)
		fits <- t(X %*% betas)
		residuals <- fits-as.matrix(Y)
	}
	return(list(residuals=residuals,betas=betas,fit=fits))
}


fastFtest <- function(fitres, 
											nullres,
											p_fit, # number of paramters in fit
											p_null, # number of parameters in null
											retF=F
											){
	SStot <- rowSums(nullres^2)
	SSres <- rowSums(fitres^2)
	Rsqs  <- 1 - (SSres/SStot)
	SSmod <- SStot - SSres
	DFres <- ncol(fitres) - p_fit
	DFmod <- p_fit - p_null 
	MSres <- SSres / DFres
	MSmod <- SSmod / DFmod
	Fstatistic <- MSmod / MSres
	pval <- pf(Fstatistic, DFmod, DFres, lower.tail=FALSE)
	if(retF){
		return(Fstatistic)
	}else{
		return(pval)
	}
}

# Use custom OLS on fit and null model, then test significance using custom F test

OLSwFtest <- function(y,mod,mod0,ver=F,robust=F,retF=F){
	if(ncol(y)<=ncol(mod)){
		warning("Cannot fit with ncol(y)<=ncol(mod)")
		return(rep(NA,nrow(y)))
	}

	if(nrow(y)==0){
		warning("No rows of data input")
		return(NULL)
	}
	
	# use base lm
	if(robust){
		pval <- apply(y,1,function(x){
			fit <- lm(x ~ mod)
			null <- lm(x ~ mod0)
			anova(fit,null)[2,6]
		})
		return(pval)
	} else {
		
		# Fit the fast model
		fit <- fastOLS(mod,y,T)
		null <- fastOLS(mod0,y,T)
		pval <- fastFtest(fit$residuals,null$residuals,ncol(mod),ncol(mod0),retF=retF)
		
		if(identical(mod,mod0)){
			warning("Fit and null models are the same")	
		}else{
			# warning, fast version is not always robust for complex models 
			# 	The model matrix is not optimized for degrees of freedom
			# Checking this on the first row each time
			x=y[1,]
			tmpfit <- lm(x ~ mod)
			tmpnull <- lm(x ~ mod0)
			rpv <- anova(tmpfit,tmpnull)[2,6]
			fastpv <- pval[1]
			if(rpv!=fastpv){warning("fastOLS/fastFtest does not equal lm, use robust=T")}
		}
		
		if(ver){
			return(list("pvalue"=pval,"fit"=fit,"null"=null))
		} else{
			return(pval)
		}
	}
}

# OLS and Ftest using a single annmatrix and formula
annOLSwFtest <- function(x,f1,f2,ver=F,robust=F,retF=F){
	mod <- model.matrix(f1,data = colanns(x))
	mod0 <- model.matrix(f2,data = colanns(x))
	OLSwFtest(x,mod,mod0,ver,robust,retF)
}

#' annmatrix combine
#'
#' Combine multiple annmatrix objects into one.
#'
#' Returns objects combined by row \code{rbind} or column \code{cbind}.
#'
#' @name combine
#'
#' @param ... annmatrix objects.
#'
#' @return combined annmatrix object.
#'
#' @examples
#'   cbinded <- cbind(annMatExample, annMatExample)
#'   rbinded <- rbind(annMatExample, annMatExample)
#'
#' @author Karolis Koncevicius \email{karolis@koncevicius.lt}
#' @export
#' TODO keep cann
rbind.annmatrix <- function(...) {
	allAnnMat <- list(...)
	if(!all(sapply(allAnnMat, is.annmatrix))) {
		tmp <- lapply(allAnnMat,function(x){
			if(is.annmatrix(x)){
				as.matrix(x)	
			} else {
				x	
			}
		})
		do.call(rbind,tmp)
		warning("Some elements are not of the class annMatrix, coercing all annmatrix to matrix")
	}
	mat <- as.matrix(allAnnMat[[1]])
	rann <- rowanns(allAnnMat[[1]])
	if(length(allAnnMat) > 1) {
		for(i in 2:length(allAnnMat)) {
			mat <- cbind(mat, as.matrix(allAnnMat[[i]]))
			rann <- rbind(rann, rowanns(allAnnMat[[i]]))
		}
	}
	annmatrix(x = mat,rann = rann)
}

#' @rdname combine
#' @export
#' TODO keep rann
cbind.annmatrix <- function(...) {
	allAnnMat <- list(...)
	if(!all(sapply(allAnnMat, is.annmatrix))) {
		tmp <- lapply(allAnnMat,function(x){
			if(is.annmatrix(x)){
			  as.matrix(x)	
			} else {
			  x	
			}
		})
		do.call(cbind,tmp)
		warning("Some elements are not of the class annMatrix, coercing all annmatrix to matrix")
	}
	mat <- as.matrix(allAnnMat[[1]])
	cann <- colanns(allAnnMat[[1]])
	if(length(allAnnMat) > 1) {
		for(i in 2:length(allAnnMat)) {
			mat <- cbind(mat, as.matrix(allAnnMat[[i]]))
			cann <- rbind(cann, colanns(allAnnMat[[i]]))
		}
	}
	annmatrix(x = mat,cann = cann)
}

# Split annmatrix by cols and lapply a function
annapply <- function(x,cols,fun,comb=T){
	grp <- apply(colanns(x)[,cols,drop=F],1,paste,collapse="_")
	res <- lapply(csplit(x,grp),fun) 
	
	# Tries to generate a data.frame
	if(comb){
		if("data.frame" %in% class(res[[1]])){
			return(res %>% data.table::rbindlist(idcol="grp"))
		} else if(is.vector(res[[1]])) {
			return(lapply(res,function(x) data.frame("value"=x)) %>% data.table::rbindlist(idcol="grp"))
		}else {
			warning("Could not aggregate results, returning a list")
			return(res)
		}
	} else {
		return(res)
	}
}

annapply_fill <- function(x,cols,fun,comb=T){
  grp <- apply(colanns(x)[,cols,drop=F],1,paste,collapse="_")
  res <- lapply(csplit(x,grp),fun) 
  
  # Tries to generate a data.frame
  if(comb){
    if("data.frame" %in% class(res[[1]])){
      return(res %>% data.table::rbindlist(idcol="grp"))
    } else if(is.vector(res[[1]])) {
      return(lapply(res,function(x) data.frame("value"=x)) %>% data.table::rbindlist(idcol="grp", fill=TRUE))
    }else {
      warning("Could not aggregate results, returning a list")
      return(res)
    }
  } else {
    return(res)
  }
}

# wrapper for PCA with annmatrix
annPCA <- function(x){
	npcs <- min(10,ncol(x))
	pca <- prcomp(t(x),scale. = T)
	mat <- t(pca$x[,1:npcs])
	cann <- colanns(x)
  varE <- data.frame(t(summary(pca)$importance))[1:npcs,]
	return(annmatrix(mat,cann=cann,rann=data.frame(PC=rownames(mat),varE)))
}

# Take median of biotechnical replicates
# By default discards samples with only 1 replicate
# Input: annmatrix
avgBioTech <- function(dat,minreps=2){
	cat(sprintf("Merging technical replicates, removing samples with Nreps < %i\n",minreps))
	tmp <- annapply(dat,"sampleid",function(x){
		if(ncol(x)<minreps){
			NULL
		} else{
			matrixStats::rowMedians(x)
		}
	},comb=F) %>% do.call(cbind,.)
	
	return(
		annmatrix(tmp,
							rann = rowanns(dat),
							cann = colanns(dat)[match(colnames(tmp),colanns(dat)$sampleid),])
	)
}

# Regress out effects and then average technical replicates
regEffavgReps <- function(mydat,resgrp,fres,minreps=2){
  
	# Remove samples without enough reps prior to regressing	
	multisamps <- names(which(table(mydat$sampleid)>=minreps))
	mydat2 <- mydat[,mydat$sampleid %in% multisamps]
	
	# Regress out fres
	cordat <- annapply(mydat2,resgrp,function(x)	annOLSwFtest(x,fres,~1,ver=T)$fit$residuals,comb = F) %>% do.call(cbind,.)
	
	# Average biotechnical replicates
	corAdat <- avgBioTech(cordat,minreps)

	return(corAdat)
}	

##### Plotting ##### 

annOLSwFtest_permSinCos <- function(x,f1,f2,N){
  lapply(1:N,function(i){
    set.seed(i)
    tmpdat <- x
    colanns(tmpdat) <- colanns(tmpdat) %>% 
      mutate(shuf=sample(1:length(sin))) %>% 
      mutate(sin=sin[shuf],cos=cos[shuf])
    data.frame(pvalue=annOLSwFtest(tmpdat,f1,f2))
  }) %>% data.table::rbindlist(idcol="iter")
}

permANOVAwithin_SinCos <- function(dat,f1,f2,N,shufvar="sampleid"){
  permdf <- annapply(dat,c("age","media"),function(x){
    annOLSwFtest_permSinCos(x,f1,f2,N)
  })
  
  realdf <- annapply(dat,c("age","media"),function(x){
    res <- annOLSwFtest(x,f1,f2)
    data.frame(pvalue=res,iter="real")
  })
  
  # p <- plotPerms(realdf,permdf)
  
  # print(p)
  
  return(list("perms"=permdf,"real"=realdf
              # p
              ))
}

permANOVAwithin_SinCos_media <- function(dat,f1,f2,N,shufvar="sampleid"){
  permdf <- annapply(dat,c("media"),function(x){
    annOLSwFtest_permSinCos(x,f1,f2,N)
  })
  
  realdf <- annapply(dat,c("media"),function(x){
    res <- annOLSwFtest(x,f1,f2)
    data.frame(pvalue=res,iter="real")
  })
  
  # p <- plotPerms(realdf,permdf)
  
  # print(p)
  
  return(list("perms"=permdf,"real"=realdf
              # p
  ))
}
