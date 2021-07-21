require(tibble)
require(dplyr)


lasso_feature_selector <- function(df,feat,flag,feat_only=T,max.feat.num=30) {
  require(glmnet)
  while(length(feat) >= max.feat.num) {
    feat.old <- feat
    lasso.out <- list(x=df %>% select(all_of(feat)) %>% as.matrix(), y=df[[flag]]) %>% 
      { glmnet::cv.glmnet(.$x, .$y, family = 'binomial', alpha = 1) }
    result <- lasso.out %>% coef(s = 'lambda.min') %>% as.matrix() %>% as.data.frame() %>% 
      tibble::rownames_to_column("VALUE") %>% setNames(c('var','weight')) %>% 
      mutate(weight = abs(weight)) %>% dplyr::slice(2:nrow(.))
    feat <- result %>% filter(weight > 0) %>% pull('var')
    if(length(feat.old) - length(feat) == 0) break
  }
  if(feat_only) feat else result
}


stepwise_feature_selector <- function(df,feat,flag){
  m <- glm(paste0(flag,'~.'),family='binomial',
           data=df %>% select(c(all_of(feat),all_of(flag))))
  result.stepaic <- step(m,direction="both",trace = 0)  
  result.stepaic$coefficients %>% names() %>% `[`(2:length(.)) %>% return()
}


vif_feature_selector <- function(df,feat,flag,feat_only=T,threshold=2){
  require(scorecard)
  while(length(feat) > 0){
    m <- glm(paste0(flag,'~.'),family='binomial',data=df.woe %>% 
      select(c(all_of(feat),all_of(flag))))
    result <- scorecard::vif(m)
    feat.rm <- result %>% filter(gvif > threshold) %>% arrange(desc(gvif))
    if(nrow(feat.rm) == 0) break
    feat.rm <- feat.rm %>% head(1) %>% pull(variable)
    feat <- setdiff(feat, feat.rm)
  }
  if(feat_only) feat else m
}


corr_feature_selector <- function(df,feat,threshold=0.6){
  mat.corr <- cor(df %>% select(all_of(feat)),use = "na.or.complete") %>% as.matrix()
  cor.pair <- list()
  cor.sum <- list()
  
  for(i in 1:(nrow(mat.corr))){
    cor.sum[[feat[i]]] <- 0
    for(j in 1:nrow(mat.corr)){
      if(abs(mat.corr[i,j]) > threshold && i!=j)
        cor.pair[[feat[i]]] <- ifelse(feat[i] %in% names(cor.pair),
                                      c(cor.pair[[feat[i]]],feat[[j]]),feat[[j]])
      if(i!=j && abs(mat.corr[i,j])<=threshold)
        cor.sum[[feat[i]]] <- abs(mat.corr[i,j]) + cor.sum[[feat[i]]]
    }
  }
  features.ignore <- c()
  while(length(cor.pair)>0){
    cnt.max <- 0
    cnt.max.key <- NULL
    # rule 1: rm key of max cnt
    for(i in 1:length(cor.pair)){
      if(length(cor.pair[[i]]) > cnt.max){
        cnt.max <- length(cor.pair[[i]])
        cnt.max.key <- names(cor.pair)[i]
      } else if(length(cor.pair[[i]]) == cnt.max){
        cnt.max.key <- c(cnt.max.key, names(cor.pair)[i])
      }
    }
    # rule 2: rm key with higher corr
    if(length(cnt.max.key) > 0){
      cnt.max.2 <- 0
      cnt.max.key.2 <- NULL
      for(i in 1:length(cnt.max.key)){
        if(cor.sum[[cnt.max.key[i]]] > cnt.max.2){
          cnt.max.2 <- cor.sum[cnt.max.key[i]]
          cnt.max.key.2 <- cnt.max.key[i]
        }
      }
      cnt.max.key <- cnt.max.key.2
    }
    cor.pair[cnt.max.key] <- NULL
    features.ignore <- c(features.ignore, cnt.max.key)
    if(length(cor.pair) == 0) break
    i <- 1
    while(i<=length(cor.pair)){
      cor.pair[[i]] <- cor.pair[[i]][!cnt.max.key %in% cor.pair[[i]]]
      if(length(cor.pair[[i]]) == 0)
        cor.pair[[i]] <- NULL
      else
        i <- i + 1
    }
  }
  return(setdiff(feat,features.ignore))
}


# TESTME 
# search optimal feature set to keep with no corr higher than threshold
corr_feature_selector.opt <- function(df,feat,threshold=0.6){
  mat.corr <- cor(df %>% select(all_of(feat)),use = "na.or.complete") %>% as.matrix()
  cor.pair <- list()
  for(i in 1:(nrow(mat.corr))){
  	for(j in 1:nrow(mat.corr)){
  		if(abs(mat.corr[i,j]) > threshold && i!=j)
  			cor.pair[[i]] <- ifelse(i %in% names(cor.pair),c(cor.pair[[i]],j),j)
  	}
  }
  state <- c(rep(-1,length(cor.pair)))  # 0=no,1=yes
  search <- function(vert,choice){
    if(choice == 0){
      state[vert] <- 0
      for(v in cor.pair[vert]){
        if(state[v] == -1){
          # if any 1 is adjacent to v
          adj.1 <- F
          for( w in cor.pair[v]){
            if(w != vert && state[w] == 1)
              adj.1 <- T
          }
          if(!adj.1) search(v,0)
          else {
            search(v,0)
            search(v,1)
          }
        }
      }
    } else {
      state[vert] <- 1
      for(v in cor.pair[vert]){
        if(state[v] == -1){
          search(v,0)
        }
      } 
    }
    cnt.vis <- 0
    for(i in state){
      if(state[i] == 0) cnt.vis <- cnt.vis + 1
    }
    if(cnt.vis == length(cor.pair) & cnt < cnt.min){
      cnt.min <- cnt
      state.min <- state
    }
    state[vert] <- -1
  }
  cnt.min <- Inf
  search(1,1)
  cnt.min.yes <- cnt.min
  state.min.yes <- state.min
  state.min <- NULL
  cnt.min <- Inf
  search(1,0)
  feat.rm <- data.frame(key = feat, 
                        hit = ifelse(state.min<state.min.yes,state.min,state.min.yes)) %>% 
    filter(hit == 1) %>% pull(key)
  return(setdiff(feat, feat.rm))
}


pval_feature_selector <- function(df,feat,flag,threshold=0.001,feat_only=T){
  feat.orgin <- feat
  while(length(feat) > 0){
    m <- glm(paste0(flag,'~.'),
             family='binomial',
             data=df.woe %>% select(c(all_of(feat),all_of(flag))))
    df.coef <- summary(m)$coefficients %>% as.data.frame()
    df.coef <- df.coef %>% tibble::rownames_to_column("var") %>% 
      filter(`Pr(>|z|)`>0.001) %>% arrange(desc(`Pr(>|z|)`))
    if(nrow(df.coef) == 0) break
    feat.ignore <- df.coef %>% head(1) %>% pull('var')
    feat <- setdiff(feat,feat.ignore)
  }
  if(feat_only) feat
  else if(length(setdiff(feat.orgin,feat))>0)
    list(coef=summary(m)$coefficients %>% as.data.frame() %>% 
      tibble::rownames_to_column("var"),feat=feat)
  else NULL
}


## TODO check features against following properties
# 1. gvif
# 2. p-value
# 3. iv
# 4. psi
# 5. lasso importance
# 6. feature AIC
