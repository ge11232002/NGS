### -----------------------------------------------------------------
### format the writeLines
### Not Exported!
my.write = function(..., sep="", collapse=" ", con=stdout()){
  args = list(...)  ## see function message
  #args = sapply(args, print)# as.character)
  text = paste(sapply(args, paste, collapse=collapse), collapse=sep)
  writeLines(text, con=con)
}

### -----------------------------------------------------------------
### Some functions to monitor the running of jobs
### Not exported!
my.time <- function(){
  format(Sys.time(), "%Y-%m-%d--%H-%M-%S")
}

my.jobStart <- function(name){
  my.write(my.time(), " ", name, " start")
  list(name=name, start=proc.time())
}

my.writeElapsed <- function(job, status="done"){
  my.write(my.time(), " ", job$name, " ", status, ": ", getElapsed(job$start))
}

getElapsed <- function(x){
  paste(signif((proc.time() - x)[ "elapsed"]/60, digits=4), "min")
}

### -----------------------------------------------------------------
### grep for a patternList with the combination of "or" or "and"
### Exported!
grepPatternList <- function(patternList, x, combine=c("or", "and"), ...){
  combine <- match.arg(combine)
  result <- logical(length(x))
  if(combine == "or"){
    for(pattern in patternList){
      idx <- grep(pattern, x, ...)
      result[idx] <- TRUE
    }
  }else if(combine == "and"){
    result[] <- TRUE
    for(pattern in patternList){
      idx <- grepl(pattern, x, ...)
      result[!idx] <- FALSE
    }
  }
  result
}

### -----------------------------------------------------------------
### return the number of threads to use
### Not Exported!
getThreads <- function(){
  ## The environmental variables of allocated cores from SGE
  nslots <- as.integer(Sys.getenv("NSLOTS"))
  
  ## Then we use the option: cores
  if (is.na(nslots)){
    if (!is.null(getOption("cores"))){
      return(getOption("cores"))
    }
    ## Use the detected cores by parallel::detectCores
    nslots <- detectCores(logical=TRUE)
  }
  return(nslots)
}

### -----------------------------------------------------------------
### More friendly mclapply
### Not Exported!
my.mclapply <- function(x, FUN, ..., mc.preschedule=TRUE, 
                        mc.set.seed=TRUE, mc.silent=FALSE, 
                        mc.cores=min(length(x), getThreads()), 
                        errorLogFile=""){
  
  if (mc.cores == 1){
    return(lapply(x, FUN, ...))
  }
  result <- mclapply(x, FUN, ..., mc.preschedule=mc.preschedule, 
                     mc.set.seed=mc.set.seed,
                     mc.silent=mc.silent, mc.cores=min(mc.cores, length(x)))
  isError <- sapply(result, function(x){any(grepPatternList("error", class(x), 
                                                    ignore.case=TRUE))})
  if (any(isError)){
    sapply(result[isError], print)
    stop("mclapply failed")
  }
  return(result)
}

### -----------------------------------------------------------------
### shrinkToRange
### Not Exported!
shrinkToRange <- function(x, theRange){
  x[ x > theRange[2]] = theRange[2]
  x[ x < theRange[1]] = theRange[1]
  return(x)
}

