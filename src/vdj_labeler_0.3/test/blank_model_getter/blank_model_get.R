# Getting model for V-genes identities.
df <- read.csv("V_id/V_id.csv", sep='\t', header=F)
probs <- read.csv("V_id/V_probs.csv", sep='\t', header=F)
probs <- as.vector(unlist(probs))
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="V_id/V_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# Getting model for D-genes identities.
df <- read.csv("D_id/D_id.csv", sep='\t', header=F)
probs <- read.csv("D_id/D_probs.csv", sep='\t', header=F)
probs <- rowSums(probs)
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="D_id/D_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# Getting model for J-genes identities.
df <- read.csv("J_id/J_id.csv", sep='\t', header=F)
probs <- read.csv("J_id/J_probs.csv", sep='\t', header=F)
probs <- colSums(probs)
df$V3 <- format(probs[df$V2], scientific=FALSE)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, c(1,3)], file="J_id/J_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# Getting model for all nucleotides is straight forward from MATLAB source

# V_dels
df <- read.csv("V_pal_del/V_id.csv", sep='\t', header=F)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
probs <- as.matrix(t(read.csv("V_pal_del/V_probs.csv", sep='\t', header=F)))
row.names(probs) <- NULL
probs <- apply(probs, c(1, 2), function(x) format(x, scientific=FALSE, digits = 15))
df <- cbind(df, probs[df$V2,])
V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, -2], file="V_pal_del/V_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# J_dels
df <- read.csv("J_pal_del/J_id.csv", sep='\t', header=F)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
probs <- as.matrix(t(read.csv("J_pal_del/J_probs.csv", sep='\t', header=F)))
row.names(probs) <- NULL
probs <- apply(probs, c(1, 2), function(x) format(x, scientific=FALSE, digits = 15))
df <- cbind(df, probs[df$V2,])
V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
write.table(df[, -2], file="J_pal_del/J_model.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE, )

# D_dels
df <- read.csv("D_pal_del/D_id.csv", sep='\t', header=F)
df$V1 <- sapply(df$V1, function(x) gsub(pattern = "\'", replacement = '', x=x))
temp <- list.files(path = "D_pal_del/matlab_matrices_for_each_gen/", pattern="*.csv")
D.pal.del.files <- lapply(temp, function(x) { read.csv(paste("D_pal_del/matlab_matrices_for_each_gen/", x, sep = ''), header=F, sep=',') } )
df$V3 <- D.pal.del.files[df$V2]
df$Left <- lapply(df$V3, rowSums)
df$Right <- lapply(df$V3, colSums)
df$Left <- lapply(df$Left, function(x) format(x, scientific = FALSE))
df$Right <- lapply(df$Right, function(x) format(x, scientific = FALSE))

write.table(cbind(df$V1, t(simplify2array(df$Left))), file="D_pal_del/D_model_Left.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(cbind(df$V1, t(simplify2array(df$Right))), file="D_pal_del/D_model_Right.csv", sep=",", quote = FALSE, row.names = FALSE, col.names = FALSE)
