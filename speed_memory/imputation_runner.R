loadScData <- function(path) {
    data <- read.csv(path,row.names=1,header=TRUE,check.names=F)
    return(data)
}

runSAVER <- function(path,ncores) {
    library(SAVER)
    data <- loadScData(path)
    imputed <- saver(data, ncores=ncores)$estimate
    return(imputed)
}

runScImpute <- function(path,ncores) {
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

runVIPER <- function(path,ncores) {
    library(VIPER)

    data <- loadScData(path)
    norm.factor <- 10^6 / colSums(data)
    data.norm <- sweep(data.pos,MARGIN=2,FUN="*",STATS=norm.factor)

    imputed <- VIPER(data.norm, num=5000, outdir=NULL)$imputed
    return(imputed)
}

run <- function(path,method,ncores,writeback=TRUE) {
    method.l = tolower(method)
    t0 = proc.time()
    
    if(method.l=='saver'){
        imputed <- runSAVER(path,ncores)
    } else if(method.l=='scimpute'){
        runScImpute(path,ncores)
    } else if(method.l=='drimpute'){
        imputed <- runDrImpute(path)
    } else if(method.l=='viper'){
        imputed <- runVIPER(path,ncores)
    } else {
        stop("Unknown method. Aborting")
    }

    if(method.l!='scimpute' && writeback){
        write.csv(imputed,'imputed.csv')
    }
    seconds = (proc.time() - t0)[[3]]
    return(seconds)
}

options <- commandArgs(trailingOnly = TRUE)
path <- options[1]
method <- options[2]
ncells <- options[3]
trial <- options[4]
threads <- as.numeric(options[5])

dt <- run(path,method,threads)

write(
    sprintf("%s,%s,%s,%s,runTime",method,ncells,trial,dt),
    sprintf("time_%s_%s_%s.txt",method,ncells,trial)
)
