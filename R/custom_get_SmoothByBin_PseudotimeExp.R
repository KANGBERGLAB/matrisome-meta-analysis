function (seurat_object, FC = TRUE, Bin = 50, method = "Pseudotime", 
    FcType = "Q95") 
{
    FcType1 <- FcType
    ByBin1 <- c("Equal.Pseudotime")
    pbmc1 <- seurat_object
    if (method == "Pseudotime") {
        TotalBin1 <- Bin
        Exp1 <- custom_smoothByBin(as.matrix(pbmc1@assays$integrated@data), 
            pbmc1@meta.data, SmoothLength1 = TotalBin1, ByBin1 = ByBin1[1])
    }
    if (method == "State") {
        TotalBin1 <- round(Bin/length(levels(as.factor(pbmc1@meta.data$State))))
        Exp1 <- smoothByState(pbmc1, pbmc1@meta.data, TotalBin1)
    }
    if (FC == TRUE) {
        if (FcType1 == "Q90") {
            Prob1 <- c(0.1, 0.9)
        }
        else {
            Prob1 <- c(0.05, 0.95)
        }
        Fc1 <- apply(Exp1, 1, function(x1) {
            x12 <- quantile(x1, probs = Prob1)
            x2 <- x12[2] - x12[1]
            return(x2)
        })
        Exp2 <- cbind(Fc1, Exp1)
        colnames(Exp2)[1:ncol(Exp2)] <- c(paste0("FoldChange", 
            FcType1), paste0("SmExp", 1:ncol(Exp1)))
        var1 <- as.data.frame(Exp2)
        return(var1)
    }
    else {
        colnames(Exp1)[1:ncol(Exp1)] <- c(paste0("SmExp", 1:ncol(Exp1)))
        var1 <- as.data.frame(Exp1)
        return(var1)
    }
}
