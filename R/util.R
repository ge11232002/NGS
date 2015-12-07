### -----------------------------------------------------------------
### grep for a patternList with the combination of "or" or "and"
### Exported!
grepPatternList <- function(patternList, x, combine=c("or", "and"), ...){
  combine <- match.arg(combine)
  result <- logical(length(x))
  if (combine == "or"){
    for (pattern in patternList){
      idx <- grep(pattern, x, ...)
      result[idx] <- TRUE
    }
  } else if (combine == "and"){
    result[] <- TRUE
    for (pattern in patternList){
      idx <- grep(pattern, x, ...)
      result[-idx] <- FALSE
    }
  }
  result
}



