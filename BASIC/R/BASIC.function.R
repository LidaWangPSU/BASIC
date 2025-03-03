#' Title
#' @param x A sets of p values from 0 to 1
#' @return
#' @export
#'
#' @examples cauchy_byrow(x)
cauchy_byrow <- function(x) {
  x <- x[!is.na(x)]
  if(length(x) == 0) return(NA)  # all NA
  p <- ACAT(x)
  return(p)
}

#' Title
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param se This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes. Column 1: Gene-SNP pair identifiers (should match the effect size matrix). Column 2: Bulk effect size standard error. Columns 3+: Standard errors for the cell type-specific eQTL effect sizes.
#' @param nPC Number fo cell types
#' @return
#' @export
#'
#' @examples meta_regression_fast(beta,se,nPC)
meta_regression_fast<-function(beta,se,nPC){



  cat(paste0("Step 1-1 - calculate MDS \n"))

  beta[is.na(beta)]<-0
  for(c in 3:ncol(se)){
    se[is.na(se[,c]),c]<-mean(se[,c],na.rm=T)
  }
  eqtl<-beta[,3:ncol(beta)]
  eqtl<-t(eqtl)
  d <- dist(eqtl, method = "euclidean")

  ####pca
  k<-(ncol(beta)-3)
  MDS <- cmdscale(d, k)
  MDS<-as.data.frame(MDS)
  MDS$study <- gsub("BETA_","",rownames(MDS))
  colnames(MDS) <- c(paste0("PC",1:k),"Cell_Type")

  cat(paste0("Step 1-2 - meta-regression \n"))
  a<-Sys.time()
  gamma_all<-data.frame()
  X<-as.matrix(MDS[,1:k])

  #X_1<-cbind(1,X)

  X_1<-as.matrix(cbind(1,X))
  X_1<-as.matrix(X_1[,1:(nPC+1)])

  X_o<-gramSchmidt(X_1)
  X_1<-X_o$Q


  Y_se<-as.numeric(se[1,3:ncol(beta)])
  R<-diag(1/Y_se^2)
  se_new<-diag(ginv(t(X_1)%*%R%*%X_1))^0.5

  #se_new<-diag(ginv(R))^0.5


  Y<-beta[,3:ncol(beta)]
  gamma_tmp<-ginv(t(X_1)%*%R%*%X_1)%*%t(X_1)%*%R
  gamma<-as.matrix(Y)%*%t(gamma_tmp)
  dim(gamma)

  se<-matrix(se_new,ncol=nrow(gamma),nrow=ncol(gamma))
  se<-t(se)

  pair<-beta[,1]
  gamma_all<-cbind(pair,gamma,se)


  b<-Sys.time()
  b-a
  print(b-a)
  colnames(gamma_all)<-c("pair",paste0("PC",0:nPC,"_effect_size"),paste0("PC",0:nPC,"_se"))
  gamma_all<-as.data.frame(gamma_all)
  gamma_all$gene_id<-sapply(strsplit(gamma_all$pair, "-"), function(x){as.character(x[1])})
  gamma_all$snp<-sapply(strsplit(gamma_all$pair, "-"), function(x){as.character(x[2])})
  gamma_all<-gamma_all[,c("pair","gene_id","snp",paste0("PC",0:nPC,"_effect_size"),paste0("PC",0:nPC,"_se"))]

  dim(gamma)
  dim(X_1)
  gamma_cells<-gamma%*%t(X_1)

  return(list(X_1,gamma_all))
}

#' Title
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param se This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes. Column 1: Gene-SNP pair identifiers (should match the effect size matrix). Column 2: Bulk effect size standard error. Columns 3+: Standard errors for the cell type-specific eQTL effect sizes.
#' @return
#' @export
#' @param nPC Number of PCs to use, (usually from 1 to number fo cell types - 1)
#' @examples nnls.weights(beta,se)
nnls.weights<-function(beta,se){
  df<-as.data.frame(beta)
  df_s<-as.data.frame(se)
  ###check NA
  #cat(paste0("Check NA in sc-eqtl \n"))
  na.list<-list()
  for(r in 3:ncol(df_s)){
    snps<-df[which(is.na(df[,r])),1]
    names<-colnames(df)[r]
    na.list[[names]]<-snps
    #cat(paste0(length(snps)," NA in ",names," \n"))
  }


  df[is.na(df)]<-0
  #df_s[is.na(df_s)]<-1
  for(c in 3:ncol(df_s)){
    df_s[is.na(df_s[,c]),c]<-mean(df_s[,c],na.rm=T)
  }

  ###check if all sc-eqtl effects are 0
  x_old<-df[,3:ncol(df)]
  #cat(paste0("Check if all sc-eqtl effects are 0 \n"))
  #cat(paste0(nrow(x_old), " gene-snp pairs in total \n"))
  sum1<-rowSums(abs(x_old))
  df<-df[which(sum1!=0),]
  #cat(paste0(nrow(df), " gene-snp pairs passed the filtering \n"))
  df_s<-df_s[which(sum1!=0),]

  ########estimate frac by NNLS
  pair<-df[,1]
  Y<-df[,2]
  X<-as.matrix(df[,3:ncol(df)])

  cat(paste0("Step1: Estimating cell type weights using NNLS \n"))
  library(pracma)
  Aeq <- matrix(rep(1, ncol(X)), nrow= 1)
  beq <- c(1)

  # Lower and upper bounds of the parameters, i.e [0, 1]
  lb <- rep(0, ncol(X))
  ub <- rep(1, ncol(X))

  # And solve:
  frac<-lsqlincon(X, Y, Aeq= Aeq, beq= beq, lb= lb, ub= ub)


  frac<-t(as.data.frame(frac))

  colnames(frac)<-colnames(X)

  weight<-frac

  weight<-as.numeric(weight)

  return(weight)
}


#' Title
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param se This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes. Column 1: Gene-SNP pair identifiers (should match the effect size matrix). Column 2: Bulk effect size standard error. Columns 3+: Standard errors for the cell type-specific eQTL effect sizes.
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param gamma_all Output from meta_regression_fast
#' @param MDS_ALL Output from meta_regression_fast
#' @param weight Weights for each cell types
#'
#' @return
#' @export
#'
#' @examples basic_function(beta,se,gamma_all,MDS_ALL,weight,nPC)
basic_function<-function(beta,se,gamma_all,MDS_ALL,weight,nPC){

  logl <- function( x0, beta_sc,se_sc,beta_bk,se_bk,prop,MDS){
    0.5*(sum((beta_sc-x0%*%t(MDS))^2/se_sc^2)+(rowSums(t(t(x0%*%t(MDS))*prop))-beta_bk)^2/se_bk^2)
  }


  pairs<-gamma_all[,1]
  beta<-beta[match(pairs,beta[,1]),]
  se<-se[match(pairs,se[,1]),]

  MDS<-as.matrix(MDS_ALL[,1:(1+nPC)])
  gamma_m<-gamma_all[,c(1,4:(4+nPC))]
  gamma_m[,2:ncol(gamma_m)]<-apply(as.matrix(gamma_m[,2:ncol(gamma_m)]),2,as.numeric)

  gamma_input<-gamma_m

  prop<-weight

  df_pc<-data.frame()
  df_factor<-data.frame()

  genes<-sapply(strsplit(gamma_all[,1], "-"), function(x){as.character(x[1])})
  gene_list<-unique(genes)
  cat(paste0(length(gene_list)," genes in total \n"))
  for(g in 1:length(gene_list)){
    cat(paste0("Analyzing gene: ",gene_list[g],"\n"))
    gamma_sub<-gamma_all[which(genes==gene_list[g]),]

    pair_sub<-intersect(beta[,1],gamma_sub[,1])

    beta_sub<-beta[match(pair_sub,beta[,1]),]
    se_sub<-se[match(pair_sub,se[,1]),]
    gamma_sub<-gamma_sub[match(pair_sub,gamma_sub[,1]),]
    gamma_sub<-as.data.frame(gamma_sub)

    beta_sub[is.na(beta_sub)]<-0
    for(c in 3:ncol(se_sub)){
      se_sub[is.na(se_sub[,c]),c]<-mean(se_sub[,c],na.rm=T)
    }

    sums<-colSums(beta_sub[,3:9])
    index<-which(sums!=0)
    index

    if(length(index)==0){
      next
    }

    if(length(index)>1){

      num<-min(length(index),(nPC+1))
      X_1<-as.matrix(MDS[index,1:num])
      X_o<-gramSchmidt(X_1)
      X_1<-X_o$Q

      df<-data.frame()
      p_val<-c()
      for(i in 1:nrow(beta_sub)){
        #print(i)
        beta_sc<-as.numeric(beta_sub[i,3:ncol(beta)])[index]
        se_sc<-as.numeric(se_sub[i,3:ncol(beta)])[index]
        beta_bk<-as.numeric(beta_sub[i,2])
        se_bk<-as.numeric(se_sub[i,2])

        x0<-as.matrix(gamma_sub[i,c(paste0("PC",0:(num-1),"_effect_size"))])

        s<-nlminb(start = x0, logl,beta_sc=beta_sc,se_sc=se_sc,beta_bk=beta_bk,se_bk=se_bk,prop=prop[index]/sum(prop[index]),MDS=as.matrix(X_1))
        logli<-s$objective
        lognull<-logl( x0=rep(0, length(x0)), beta_sc=beta_sc,se_sc=se_sc,beta_bk=beta_bk,se_bk=se_bk,prop=prop[index]/sum(prop[index]),MDS=as.matrix(X_1))
        delta<-2*(lognull-logli)
        p_val[i]<-pchisq(delta,df=length(x0),lower.tail=FALSE)

        df<-rbind(df,s$par)
      }

      head(df)

      se_sc<-as.data.frame(se_sub[,3:ncol(beta)])
      se_sc<-se_sc[,index]

      se_x<-data.frame()
      for(l in 1:num){
        theta_l<-X_1[,l]*prop[index]
        V1 = rowSums(as.matrix(1/se_sc^2)%*%diag(X_1[,l]^2))
        V2 = sum(theta_l)^2/(se_sub[,2]^2)
        V=1/(V1+V2)
        gamma_se<-V^0.5
        se_x<-rbind(se_x,gamma_se)
      }
      se_x<-t(se_x)
      head(se_x)

      ####pc eqtls
      beta_PC<-matrix(NA,nrow=nrow(gamma_sub),ncol=(nPC+1))
      se_PC<-matrix(NA,nrow=nrow(gamma_sub),ncol=(nPC+1))

      beta_PC[,1:num]<-as.matrix(df)
      se_PC[,1:num]<-se_x

      output_PC<-cbind(gamma_sub[,1:3],beta_PC,se_PC)
      colnames(output_PC)<-c("pair","gene_id","snp",paste0("PC",0:nPC,"_effect_size"),paste0("PC",0:nPC,"_se"))

      #output_PC$LR_p_val<-p_val

      df_pc<-rbind(df_pc,output_PC)
    }else
    {

      beta_sc<-as.data.frame(beta_sub[,3:ncol(beta)])
      se_sc<-as.data.frame(se_sub[,3:ncol(beta)])
      beta_sc<-beta_sc[,index]
      se_sc<-se_sc[,index]

      beta_PC<-matrix(0,nrow=nrow(gamma_sub),ncol=(nPC+1))
      se_PC<-matrix(1,nrow=nrow(gamma_sub),ncol=(nPC+1))

      beta_PC[,1]<-beta_sc
      se_PC[,1]<-se_sc

      output_PC<-cbind(gamma_sub[,1:3],beta_PC,se_PC)
      colnames(output_PC)<-c("pair","gene_id","snp",paste0("PC",0:nPC,"_effect_size"),paste0("PC",0:nPC,"_se"))

      p_val<-pchisq((beta_sc/se_sc)^2,df=1,lower.tail=FALSE);
      #output_PC$LR_p_val<-p_val

      df_pc<-rbind(df_pc,output_PC)

    }


  }

  return(df_pc)
}

#' Title
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param se This matrix has the same dimensions as the effect size matrix and represents the standard errors associated with the eQTL effect sizes. Column 1: Gene-SNP pair identifiers (should match the effect size matrix). Column 2: Bulk effect size standard error. Columns 3+: Standard errors for the cell type-specific eQTL effect sizes.
#' @param beta This matrix contains eQTL effect sizes for both bulk and single-cell data. Column 1: Gene-SNP pair identifiers. Column 2: Bulk effect size. Columns 3+: Cell type-specific eQTL effect sizes.
#' @param gamma_all Output from meta_regression_fast
#' @param MDS_ALL Output from meta_regression_fast
#' @param weight Weights for each cell types
#'
#' @return
#' @export
#'
#' @examples basic_function(beta,se,gamma_all,MDS_ALL,weight,nPC)
basic_axisQTL<-function(beta,se,gamma_all,MDS_ALL,weight,nPC){

  pairs<-gamma_all[,1]
  beta<-beta[match(pairs,beta[,1]),]
  se<-se[match(pairs,se[,1]),]

  cauchy_pvalue<-gamma_all[,1:3]

  cat(paste0("Fit ",nPC," PCs BASIC model \n"))
  basic_list<-basic_function(beta,se,gamma_all,MDS_ALL,weight,nPC=nPC)
  cauchy_z<-basic_list[,paste0("PC",0:nPC,"_effect_size")]/basic_list[,paste0("PC",0:nPC,"_se")]
  cauchy_z<-apply(cauchy_z,2,as.numeric)
  cauchy_p<-pchisq(cauchy_z^2,df=1,lower.tail = F)
  cauchy_pvalue<-cbind(cauchy_pvalue,cauchy_p)

  colnames(cauchy_pvalue)<-c(colnames(cauchy_pvalue)[1:3],paste0("PC",0:nPC,"_pval"))

  cauchy_pvalue$p_cauchy <- unlist(apply(cauchy_pvalue[,4:ncol(cauchy_pvalue)], 1, cauchy_byrow))

  return(list("pvalue"=cauchy_pvalue,"axis-QTL"=basic_list))
}


