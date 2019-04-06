###produce the first column 
FirstColumn <- function(data){
  data <- data.frame(data)
  Taxa <- rownames(data)
  rownames(data) <- NULL
  data <- cbind(Taxa,data)
  return(data)
}
