####Scatterplot function 
Scatterplot<- function(mydata){
  panel.cor <- function(x, y, digits=4, prefix="", cex.cor) 
  {
    usr <- par("usr"); on.exit(par(usr)) 
    par(usr = c(0, 1, 0, 1)) 
    r <- abs(cor(x, y, use="pairwise")) 
    txt <- format(c(r, 0.123456789), digits=digits)[1] 
    txt <- paste(prefix, txt, sep="") 
    if(missing(cex.cor)) cex <- 0.9/strwidth(txt) 
    
    test <- cor.test(x, y, use="pairwise")
    # borrowed from printCoefmat
    Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                     cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                     symbols = c("***", "**", "*", ".", " ")) 
    text(0.5, 0.5, txt, cex = cex * r) 
    text(.5, .8, Signif, cex = cex, col=2)
  }
  pairs(mydata, lower.panel=panel.smooth, upper.panel=panel.cor)
}
