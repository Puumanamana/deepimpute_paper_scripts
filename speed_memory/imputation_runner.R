loadScData <- function(path) {
    data <- read.csv(path,row.names=1,header=TRUE,check.names=F)
    return(data)
}


runSAVER <- function(path,ncores=1) {
    library(SAVER)
    data <- loadScData(path)
    imputed <- saver(data, ncores=ncores)$estimate
    return(imputed)
}

runScImpute <- function(path,ncores=1) {
    library(scImpute)
    scimpute(path,labeled=FALSE,Kcluster=1,ncores=ncores)
}

runDrImpute <- function(path) {
    library(DrImpute)
    
    data <- as.matrix(loadScData(path))
    
    data <- preprocessSC(data)
    size.factor <- apply(data, 2, mean)
    
    data.preproc <- log1p(
        t( t(data) / size.factor )
    )
    imputed <- DrImpute(data.preproc)
    return(imputed)
}

runVIPER <- function(path,ncores=1) {
    library(VIPER)

    data <- loadScData(path)
    norm.factor <- 10^6 / colSums(data)
    data.norm <- sweep(data.pos,MARGIN=2,FUN="*",STATS=norm.factor)

    imputed <- VIPER(data.norm, num=5000, outdir=NULL)$imputed
    return(imputed)
}

run <- function(path,method,ncores=1,writeback=TRUE) {
    method.l = method.tolower()
    t0 = proc.time()
    
    if(method.l=='saver'){
        imputed <- runSAVER(path,ncores=ncores)
    } else if(method.l=='scimpute'){
        runScImpute(path,ncores=ncores)
    } else if(method.l=='drimpute'){
        imputed <- runDrImpute(path)
    } else if(method.l=='viper'){
        imputed <- runVIPER(path,ncores=ncores)
    } else {
        stop("Unknown method. Aborting")
    }

    if(method.l!='scimpute' && writeback){
        write.csv(imputed,'imputed.csv')
    }
    seconds = (proc.time() - ptm)[[3]]
    return(seconds)
}

options <- commandArgs(trailingOnly = TRUE)
path <- options[1]
method <- options[2]
trial <- option[3]
threads <- options[4]

dt <- run(path,method,ncores)

write(
    sprintf("{} {} {}",method,trial,dt),
    sprintf("time_{}_{}.txt",method,trial)
)
