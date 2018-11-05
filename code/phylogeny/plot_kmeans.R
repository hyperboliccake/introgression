a=read.table('strain_gene_introgressed_matrix.tsv', sep='\t',header=T,row.names=1)
wss <- (nrow(a)-1)*sum(apply(a,2,var))
for (i in 1:15) {

    wss[i] = sum(kmeans(a, centers=i)$withinss)
    
}
print(wss)
pdf('kmeans.pdf', width = 7, height = 7)
plot(1:15, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
dev.off()
