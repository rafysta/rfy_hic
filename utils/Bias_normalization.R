# Hi-Cのbiasを処理する

suppressPackageStartupMessages(library("optparse"))
option_list <- list(  
  make_option(c("-i", "--in"),help="Distance normalized matrix file"),
  make_option(c("-o", "--out"),help="Normalized matrix file"),
  make_option(c("--inter"), default="NA", help="read for inter-chromosome"),
  make_option(c("--times"), default="30", help="how many times apply normalization")
)
opt <- parse_args(OptionParser(option_list=option_list))

FILE_matrix <- as.character(opt["in"])
FILE_out <- as.character(opt["out"])
FILE_inter <- as.character(opt["inter"])

FILE_object <- sub(".matrix", ".rds", FILE_matrix)
FILE_log <- sub(".matrix", ".log", FILE_out)
if(!file.exists(FILE_object)){
  if(!file.exists(FILE_matrix)){
    map <- as.matrix(read.table(FILE_matrix, header=TRUE, check.names = FALSE))
  }else{
    cat("There is not", FILE_matrix, "file\n")
    q()
  }
}else{
  map <- readRDS(FILE_object)
}
r <- rownames(map)

#------------------------------------------------
# 全体を100に分割し、0readの区画が60%以上を占めるbinは計算にいれない
#------------------------------------------------
if(nrow(map) > 100){
  # 全体を約100に分割
  n <- cut(1:nrow(map), breaks=c(seq(0,nrow(map),by=round(nrow(map)/100)), nrow(map)+1))
  # 各区間について合計を計算
  SUM_sec <- aggregate(map, by=list(n), sum)
  # 0である区間の合計を調べる
  NUM0 <- apply(SUM_sec[,2:ncol(SUM_sec)]==0, 2, sum)
  
  
  #------------------------------------------------
  # 残りのbinのうちBottom 2%は除く
  #------------------------------------------------
  v1 <- apply(map,1,sum)
  v1_sort <- sort(v1[NUM0 < 60])
  Threshold <- v1_sort[length(v1_sort)*0.02]
  
  ### 除くbin
  index_remove <- names(which(v1 < Threshold))
  
  map[index_remove, ] <- NA
  map[,index_remove] <- NA
}

#------------------------------------------------
# inter-chromosomeの情報がある場合に読む
#------------------------------------------------
if(FILE_inter != "NA"){
  D <- read.table(FILE_inter, sep="\t", header=F, row.names = 1)
  INTER <- D[,1]
  names(INTER) <- rownames(D)
}

Single <- function(m, b){
  Db <- apply(m, 2, sum, na.rm=TRUE)
  
  # inter-chromosomeの値があるときにはbias修正
  if(FILE_inter != "NA"){
    Db <- Db + INTER / b / mean(b[b!=0], na.rm=T)
  }
  
  Db <- Db / mean(Db[Db != 0])
  Db <- ifelse(Db == 0, 1, Db)
  
  d <- nrow(m)
  m <- m / matrix(rep(Db, times=d), nrow=d)
  m <- m / matrix(rep(Db, each=d), nrow=d)

  list(map=m, bias=Db)
}


multi <- function(m, times){
  Variances <- c()
  B <- rep(1, length(r))
  
  # 初期のvarianceを求める
  initial_var <- apply(m, 2, sum, na.rm=TRUE)
  if(FILE_inter != "NA"){
    initial_var <- initial_var + INTER
  }
  initial_var <- initial_var / mean(initial_var[initial_var != 0], na.rm=T)
  initial_var <- ifelse(initial_var == 0, 1, initial_var)
  initial_var <- var(initial_var)
  cat("Variance at 0 times:\t", initial_var, "\n", sep="", file = FILE_log, append = TRUE)
  
  for(i in 1:times){
    S <- Single(m, B)
    m <- S$map
    B <- B * S$bias
    cat("Variance at ", i, " times:\t", var(S$bias), "\n", sep="", file = FILE_log, append = TRUE)
  }
  m
}

options(warn=-1)
TIMES_apply <- as.numeric(as.character(opt["times"]))
map <- multi(map, TIMES_apply)


colnames(map) <- r
rownames(map) <- r
write.table(map, file=FILE_out, append=TRUE, quote=FALSE, sep="\t", eol="\n", row.names=TRUE, col.names=NA)



