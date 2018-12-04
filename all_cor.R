setwd("/Volumes/Mac_Workplace/Mac_workspace/leptin/all")
library(DESeq2)
library(ggplot2)
library(reshape2)
library(Rtsne)
#library(scatterplot3d)
#library(rgl)
library(stringr)
library(ggsci)
library(data.table)
library(pheatmap)
library(dplyr)
#library(WGCNA)
library(qpcR)
########################################
################# pre-define functions
########################################
gran_theme <- theme_classic() +
  theme(legend.title =element_blank(),legend.text = element_text(size = 12, face = "bold"),axis.title =element_text(size=12,face = "bold") ,axis.text=element_text(size=12))

depict_vol_plot <- function(deseq2_res,pass_title){
  deseq2_res <- na.omit(as.data.frame(deseq2_res)) ###Exclude lines with 'NA' qvalue.
  significant <- c()
  for (i in row.names(deseq2_res)) 
    if (deseq2_res[i,"log2FoldChange"] >log2(1.5) && deseq2_res[i,"padj"] < 0.05) {
      significant = c(significant,"up")} else if (deseq2_res[i,"log2FoldChange"] < -log2(1.5) && deseq2_res[i,"padj"] < 0.05){
        significant = c(significant,"down")} else{
          significant = c(significant,"non")}
  depict_matrix <- cbind(deseq2_res,significant)
  depict_matrix <- depict_matrix[which(abs(depict_matrix$log2FoldChange) <=  10,),]
  depict_matrix <- depict_matrix[which(abs(depict_matrix$padj) >= 1e-50,),]
  volcano <- ggplot(depict_matrix, aes(x= log2FoldChange, y= -1*log10(padj)))
  volcano+geom_point(aes(color=significant))+ 
    scale_color_manual(values =c("#4DBBD5B2","grey","#E64B35B2"))+
    labs(title=pass_title,x="log2(fold change)", y="-log10(adj.pvalue)")+
    geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-0.585,0.585),linetype=3)+
    gran_theme + theme(legend.position='none',plot.title = element_text(hjust = 0.5,size=12,face = "bold")) 
}

depict_point_plot <- function(expr_cols,high_list,low_list,pass_title){ ###Notice: expr_cols - WT_KO
  significant <- rep(0,nrow(expr_cols))
  depict_matrix <- cbind(expr_cols,significant)
  depict_matrix[high_list,]$significant <- "high"
  depict_matrix[low_list,]$significant <- "low"
  depict_matrix[which(depict_matrix$significant == "0"),]$significant <- "not_sign"
  depict_matrix$significant <- factor(depict_matrix$significant,levels = c("low","not_sign","high"))
  point <- ggplot(depict_matrix, aes(x=depict_matrix[,1], y= depict_matrix[,2]))
  point + geom_point(aes(color=significant))+ 
    scale_color_manual(values =c("green3","grey","firebrick2"))+
    labs(title=pass_title,x="Normalized counts of WT", y="Normalized counts of KO")+
    gran_theme + theme(legend.position='none',plot.title = element_text(hjust = 0.5,size=14,face = "bold")) 
}

david2gp <- function(csv,title_def) {
  csv <- read.csv(paste0("./all_go/",csv),header = T,col.names = c("Cate","Desc.","Pvalue","LogPvalue","dir"))
  p <- ggplot(csv) + geom_bar(aes(x=reorder(Desc.,LogPvalue),y=LogPvalue),stat = "identity",fill="#4DBBD5B2")+coord_flip()
  p <- p + theme(panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.title = element_text(),
                 axis.text = element_text(size=14),
                 axis.line.x = element_line(colour = "black"),
                 axis.title.x = element_text(size = 14,  face = "bold"),
                 axis.title.y = element_text(size = 14, face = "bold"),
                 axis.ticks.y = element_blank()) +scale_y_continuous(expand = c(0, 0))+ xlab("Terms")+
    labs(title=title_def)+ theme(plot.title = element_text(hjust = 0.5,size=14,face = "bold")) 
  return(p)
}

########################################
#################### Set-up 
########################################
inflam_geneset <- as.character(read.csv("positive_regulation_of_inflammatory_response.txt",col.names = "V1")$V1)
#fat_geneset <- as.character(read.csv("fat_geneset_LP.csv",col.names = "V1")$V1) Fig2F
fat_geneset <- as.character(read.csv("fatty_acid_biosynthetic_process.txt",col.names = "V1")$V1)
rawmat1 <- read.csv("4_8_16_merge_geneCounts_featureCounts_hisat2.txt",sep="\t",header=T)
rawmat1 <- rawmat1[!duplicated(rawmat1[,1]),]
rownames(rawmat1) <- rawmat1$X.id
rawmat1 <- rawmat1[,-1]
rawmat2 <- read.csv("new_rat_counts.txt",sep="\t",header=T)
rawmat2 <- rawmat2[!duplicated(rawmat2[,1]),]
rownames(rawmat2) <- rawmat2$X.id
rawmat2 <- rawmat2[,-1]
rawmat <- cbind(rawmat1,rawmat2)


coldata <- read.csv("coldata.csv",row.names = 1,header = T,strip.white = T)
all(rownames(coldata) %in% colnames(rawmat))
all(rownames(coldata) == colnames(rawmat))
rawmat <- rawmat[,rownames(coldata)]
all(rownames(coldata) == colnames(rawmat))

####################Deseq
dds <- DESeqDataSetFromMatrix(countData = rawmat,colData = coldata,design = ~ condition)
keep <- rowSums(counts(dds)) >= 10 
dds <- dds[keep,] 
#dds$condition <- relevel(dds$condition, ref = "untreated") 
dds <- DESeq(dds)
res <- results(dds) #FDR=0.1 by default
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > 1, na.rm = TRUE)
rlog <- rlog(dds,blind=F)

##################Normalized matrix
normat <- counts(dds,normalized=T)
normat <- log2(counts(dds,normalized=T)+1)
write.csv(normat,"rat_all_normat.csv",quote=F)
#normat <- read.csv("rat_all_normat.txt",header=T,row.names = 1)

########################################
#################### quality control
########################################

#########PCA
pca<- prcomp(t(normat))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
pca_result<-as.data.frame(pca$x)
pca_result<-cbind(pca_result,coldata)
pca_result$treat <- factor(pca_result$treat,levels = c("treated","untreated"))
pca_result$type <- factor(pca_result$type,levels = c("W4","W8","W16","W32","W48"))
p<-ggplot(pca_result)+geom_point(aes(x=pca_result[,1],y=pca_result[,2],color=pca_result$type,shape = pca_result$treat),size=5)
p<-p+gran_theme+theme(legend.title =element_blank())+labs(x="PC1: 21 % variance",y="PC2: 18 % variance")
p + scale_color_manual(values=c("LightGray","Cyan","Blue","magenta","red")) 
ggsave("newdata_all_PCA_newcol.pdf",dpi=300)

# library(scatterplot3d)
# scatterplot3d(x=pca_result$PC1,y=pca_result$PC2,z=pca_result$PC3,
#               color=c("red","yellow","green","blue","black")[coldata$type])

# library(rgl)
# plot3d(x=pca_result$PC1,y=pca_result$PC2,z=pca_result$PC3,
#        col=c("red","yellow","green","blue","black")[coldata$type],
#        type="s",radius=0.5)


#########tSNE
iris_unique <- unique(normat) # Remove duplicates
iris_matrix <- t(as.matrix(iris_unique))
for(i in (5:9)){
  tsne_out <- Rtsne(iris_matrix,perplexity=i,max_iter = 8000) 
  p<-ggplot(data.frame(tsne_out$Y))+
    geom_point(aes(x=tsne_out$Y[,1],y=tsne_out$Y[,2],color=pca_result$type,shape = pca_result$treat),size=5)
  p<-p+labs(x="TSNE 1",y="TSNE 2")  + gran_theme +  scale_color_manual(values=c("LightGray","LightCyan","Cyan","DeepSkyBlue","Green1"))  
  ggsave(paste('TSNE_',i,'_all_newcolor.pdf'),dpi = 300) 
}
 
######################################
############ Equalization for Reps
######################################

mean_normat <- data.frame(gene=rownames(normat))
rep_class <- levels(factor(str_sub(colnames(normat),1,-2)))
for(i in rep_class){
  mean_normat <- cbind(mean_normat,apply(normat[,grep(paste0(i,"."),colnames(normat))],1,mean))
}
colnames(mean_normat) <- c('gene',rep_class)
mean_normat <- mean_normat[,-1][,c("W4_WT","W8_WT","W16_WT","W32_WT","W48_WT","W4_KO","W8_KO","W16_KO","W32_KO","W48_KO")]





################################################
################ Modelling on inflam & fat dev
################################################

#######################################
##### Scenario1: Intersect then merge
####################################### 

############ Seperated heatmap. 

inflam_normat <- mean_normat[which(rownames(mean_normat) %in% de_genes),]
inflam_normat <- inflam_normat[which(rownames(inflam_normat) %in% inflam_geneset),]
inflam_normat <- t(apply(inflam_normat,1, function(x) (x-mean(x))/sd(x)))
a = pheatmap(inflam_normat,cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,fontsize_row =6,fontsize_col =8) 
temp = inflam_normat[,c(6:10)]-inflam_normat[,c(1:5)]
#a = pheatmap(temp,cluster_cols = F,border_color = NA,legend = F,fontsize_row =6,fontsize_col =8)
#a=pheatmap(temp[which(cutree(a$tree_row,2) == "1",),],cluster_cols = F,border_color = NA,fontsize_row =6,fontsize_col =8)
ggsave("inflam_geneset_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 4)
#ggsave("zhu_want_fig3a_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 4)

fat_normat <- mean_normat[which(rownames(mean_normat) %in% de_genes),]
fat_normat <- fat_normat[which(rownames(fat_normat) %in% fat_geneset),]
#fat_normat <- mean_normat  ##### For specific genes picked by LU. fig2f
#fat_normat <- mean_normat[which(rownames(mean_normat) %in% fat_geneset),] ## fig2f
fat_normat <- t(apply(fat_normat,1, function(x) (x-mean(x))/sd(x)))
a=pheatmap(fat_normat,cluster_cols = F,cluster_rows = F,border_color = NA,legend = F,fontsize_row =6,fontsize_col =8)
a=pheatmap(fat_normat[,c(6:10)]-fat_normat[,c(1:5)],cluster_cols = F,border_color = NA,fontsize_row =6,fontsize_col =8) #fig2f
ggsave("fat_geneset_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 5)
ggsave("zhu_want_fig2f_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 3)
############# Data cleaning.

two_geneset <- as.data.frame(rbind(inflam_normat,fat_normat))
two_geneset$cate <- rep(c("inflam","fat"),c(dim(inflam_normat)[1],dim(fat_normat)[1]))
two_geneset$gene <- rownames(two_geneset) 


two_geneset<-melt(
  two_geneset,                                      #待转换的数据集名称
  id.vars=c("cate","gene"),  #要保留的主字段
  variable.name="name",         #转换后的分类字段名称（维度）
  value.name="expr"            #转换后的度量值名称
)
two_geneset$time <- factor(str_split_fixed(two_geneset$name,"_",2)[,1],levels = c("W4","W8","W16","W32","W48"))
two_geneset$treat <- str_sub(str_split_fixed(two_geneset$name,"_",2)[,2],1,2)


###############################################
################################ Mean_dataframe
###############################################
inflam_mean <- c(apply(inflam_normat,2, mean))
fat_mean <- apply(fat_normat,2, mean)

two_geneset_mean <- data.frame(rbind(inflam_mean,fat_mean))
two_geneset_mean$cate <- c("inflam","fat")
two_geneset_mean<-melt(
  two_geneset_mean,                                      #待转换的数据集名称
  id.vars=c("cate"),  #要保留的主字段
  variable.name="name",         #转换后的分类字段名称（维度）
  value.name="expr"            #转换后的度量值名称
)
two_geneset_mean$time <- factor(str_split_fixed(two_geneset_mean$name,"_",2)[,1],levels = c("W4","W8","W16","W32","W48"))
two_geneset_mean$treat <- str_sub(str_split_fixed(two_geneset_mean$name,"_",2)[,2],1,2)

############# Timescale line.
ggplot(data=two_geneset)+geom_jitter(aes(x=time,y=expr,group=treat,color=treat),alpha = 0.3 )+
  geom_line(data = two_geneset_mean , aes(x=time,y=expr,group=treat,color=treat),alpha = .7, size = 2)+
  facet_grid(. ~ cate) + gran_theme + theme(strip.text.x = element_text(size = 12)) +  scale_color_npg() 
ggsave("2geneset_timescale_final_LP.pdf",dpi=300,width = 5,height = 4)
ggsave("test_fat_inflam_timescale.pdf",dpi=300)

##############Wilcox.test of each stages.

wilcox.test(fat_normat[,1],fat_normat[,6],paired = T,alternative = "less")  #W4:0.02 W8:0.05  W48:0.005
wilcox.test(inflam_normat[,3],inflam_normat[,8],paired = T,alternative = "less") #W16:0.04 W32:0.02 W48: 0.02

 
# #########################################
# ##### Scenario2: Intersect at each stage
# ######################################### 
# 
# ############ Seperated heatmap. 
# 
# inflam_normat <- inflam_normat[which(rownames(inflam_normat) %in% inflam_geneset),]
# inflam_normat <- t(apply(inflam_normat,1, function(x) (x-mean(x))/sd(x)))
# a = pheatmap(inflam_normat,cluster_cols = F,cluster_rows = F,border_color = NA,legend = F) 
# ggsave("inflam_geneset_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 6)
# 
# fat_normat <- fat_normat[which(rownames(fat_normat) %in% fat_geneset),]
# fat_normat <- t(apply(fat_normat,1, function(x) (x-mean(x))/sd(x)))
# a=pheatmap(fat_normat,cluster_cols = F,cluster_rows = F,border_color = NA,legend = F)
# #ggsave("fat_geneset_heatmap.png"",a$gtable,dpi=300,width = 3,height = 6)
# ggsave("fat_geneset_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 6)
# 
# ############# Data cleaning.
# two_geneset <- as.data.frame(rbind(inflam_normat,fat_normat))
# two_geneset$cate <- rep(c("inflam","fat"),c(dim(inflam_normat)[1],dim(fat_normat)[1]))
# two_geneset$gene <- rownames(two_geneset) 
# 
# 
# two_geneset<-melt(
#   two_geneset,                                      #待转换的数据集名称
#   id.vars=c("cate","gene"),  #要保留的主字段
#   variable.name="name",         #转换后的分类字段名称（维度）
#   value.name="expr"            #转换后的度量值名称
# )
# two_geneset$time <- factor(str_split_fixed(two_geneset$name,"_",2)[,1],levels = c("W4","W8","W16","W32","W48"))
# two_geneset$treat <- str_sub(str_split_fixed(two_geneset$name,"_",2)[,2],1,2)
# 
# ################ Note: We intersect it with DE genes SEPERATELY.
# two_geneset <- rbind(
# two_geneset[which(two_geneset$time == "W4"),][which(two_geneset[which(two_geneset$time == "W4"),]$gene %in% diff_w4$gene),],
# two_geneset[which(two_geneset$time == "W8"),][which(two_geneset[which(two_geneset$time == "W8"),]$gene %in% diff_w8$gene),],
# two_geneset[which(two_geneset$time == "W16"),][which(two_geneset[which(two_geneset$time == "W16"),]$gene %in% diff_w16$gene),],
# two_geneset[which(two_geneset$time == "W32"),][which(two_geneset[which(two_geneset$time == "W32"),]$gene %in% diff_w32$gene),],
# two_geneset[which(two_geneset$time == "W48"),][which(two_geneset[which(two_geneset$time == "W48"),]$gene %in% diff_w48$gene),]
# )
# 
# 
# ###############################################
# ################################ Mean_dataframe
# ###############################################
# 
# two_geneset_mean <- two_geneset %>% group_by(name,cate,time,treat) %>% summarise(mean=mean(expr))
# 
#   
#   
#   
# 
# ############# Timescale line.
# ggplot(data=two_geneset)+geom_jitter(aes(x=time,y=expr,group=treat,color=treat),alpha = 0.4 )+
#   geom_line(data = two_geneset_mean , aes(x=time,y=mean,group=treat,color=treat),alpha = .8, size = 3)+
#   facet_grid(. ~ cate) + gran_theme + theme(strip.text.x = element_text(size = 16)) +  scale_color_npg() 
# ggsave("2geneset_timescale_final_LP.pdf",dpi=300)
# ggsave("test_fat_inflam_timescale.pdf",dpi=300)
# 
# ##############Wilcox.test of each stages.
# 
# wilcox.test(fat_normat[,5],fat_normat[,10],paired = TRUE) #W4~W48 wilcox=0.04,0.09,0.98,0.93,0.01
# wilcox.test(inflam_normat[,5],inflam_normat[,10],paired = TRUE) #W4~W48  wilcox=0.44,0.87,0.07,0.04,0.03
# 
# 






################################################
#################### Differential Expression
################################################

res <- results(dds,contrast=c("condition","W4_KO","W4_WT")) 
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, na.rm = TRUE)
diff_w4_up <- rownames(res[which((res$log2FoldChange) > log2(1.5) & res$padj < 0.05),])
diff_w4_down <- rownames(res[which((res$log2FoldChange) < log2(1/1.5)  & res$padj < 0.05),])
depict_vol_plot(res,"Genes differentially expressed in Week4")
ggsave("Volcano_w4.pdf",dpi=300,width = 4,height = 3)

res <- results(dds,contrast=c("condition","W8_KO","W8_WT")) 
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, na.rm = TRUE)
diff_w8_up <- rownames(res[which((res$log2FoldChange) > log2(1.5) & res$padj < 0.05),])
diff_w8_down <- rownames(res[which((res$log2FoldChange) < log2(1/1.5)  & res$padj < 0.05),])
depict_vol_plot(res,"Genes differentially expressed in Week8")
ggsave("Volcano_w8.pdf",dpi=300,width = 4,height = 3)

res <- results(dds,contrast=c("condition","W16_KO","W16_WT")) 
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, na.rm = TRUE)
diff_w16_up <- rownames(res[which((res$log2FoldChange) > log2(1.5) & res$padj < 0.05),])
diff_w16_down <- rownames(res[which((res$log2FoldChange) < log2(1/1.5)  & res$padj < 0.05),])
depict_vol_plot(res,"Genes differentially expressed in Week16")
ggsave("Volcano_w16.pdf",dpi=300,width = 4,height = 3)

res <- results(dds,contrast=c("condition","W32_KO","W32_WT")) 
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, na.rm = TRUE)
diff_w32_up <- rownames(res[which((res$log2FoldChange) > log2(1.5) & res$padj < 0.05),])
diff_w32_down <- rownames(res[which((res$log2FoldChange) < log2(1/1.5)  & res$padj < 0.05),])
depict_vol_plot(res,"Genes differentially expressed in Week32")
ggsave("Volcano_w32.pdf",dpi=300,width = 4,height = 3)
# diff_w32 <- data.frame(gene=c(diff_w32_up,diff_w32_down),property=rep(c("up","down"),c(length(diff_w32_up),length(diff_w32_down))))
# write.csv(diff_w32,"diff_w32.csv",quote = F)

res <- results(dds,contrast=c("condition","W48_KO","W48_WT")) 
resOrdered <- res[order(res$padj),] 
sum(abs(res$log2FoldChange) > log2(1.5) & res$padj < 0.05, na.rm = TRUE)
diff_w48_up <- rownames(res[which((res$log2FoldChange) > log2(1.5) & res$padj < 0.05),])
diff_w48_down <- rownames(res[which((res$log2FoldChange) < log2(1/1.5)  & res$padj < 0.05),])
depict_vol_plot(res,"Genes differentially expressed in Week48")
ggsave("Volcano_w48.pdf",dpi=300,width = 4,height = 3)
#diff_w48 <- data.frame(gene=c(diff_w48_up,diff_w48_down),property=rep(c("up","down"),c(length(diff_w48_up),length(diff_w48_down))))
#write.csv(diff_w48,"diff_w48.csv",quote = F)

depict_point_plot(mean_normat[,c(1,6)],diff_w4_up,diff_w4_down,"Week 4 DE genes")
ggsave("Point_w4.pdf",dpi=300,width = 5,height = 4)
depict_point_plot(mean_normat[,c(2,7)],diff_w8_up,diff_w8_down,"Week 8 DE genes")
ggsave("Point_w8.pdf",dpi=300,width = 5,height = 4)
depict_point_plot(mean_normat[,c(3,8)],diff_w16_up,diff_w16_down,"Week 16 DE genes")
ggsave("Point_w16.pdf",dpi=300,width = 5,height = 4)
depict_point_plot(mean_normat[,c(4,9)],diff_w32_up,diff_w32_down,"Week 32 DE genes")
ggsave("Point_w32.pdf",dpi=300,width = 5,height = 4)
depict_point_plot(mean_normat[,c(5,10)],diff_w48_up,diff_w48_down,"Week 48 DE genes")
ggsave("Point_w48.pdf",dpi=300,width = 5,height = 4)


diff_up <- unique(c(diff_w4_up,diff_w8_up,diff_w16_up,diff_w32_up,diff_w48_up)) 
diff_down <- unique(c(diff_w4_down,diff_w8_down,diff_w16_down,diff_w32_down,diff_w48_down)) 
de_genes <- unique(c(diff_up,diff_down))

write.csv(diff_up,"rat_diff_up.csv",quote = F,row.names = F)
write.csv(diff_down,"rat_diff_down.csv",quote = F,row.names = F)

##############For evolution analysis
library(qpcR)
all_diff <- qpcR:::cbind.na(diff_w4_up,diff_w4_down,diff_w8_up,diff_w8_down,diff_w16_up,diff_w16_down,
                            diff_w32_up,diff_w32_down,diff_w48_up,diff_w48_down)
write.csv(all_diff,"rat_all_diff.csv",quote = F,row.names = F)

write.csv(sort(c(diff_w4_up,diff_w4_down)),"rat_diff_w4.csv",quote = F,row.names = F)
write.csv(sort(c(diff_w8_up,diff_w8_down)),"rat_diff_w8.csv",quote = F,row.names = F)
write.csv(sort(c(diff_w16_up,diff_w16_down)),"rat_diff_w16.csv",quote = F,row.names = F)
write.csv(sort(c(diff_w32_up,diff_w32_down)),"rat_diff_w32.csv",quote = F,row.names = F)
write.csv(sort(c(diff_w48_up,diff_w48_down)),"rat_diff_w48.csv",quote = F,row.names = F)

#########################################################
############################### DAVID GO enrichment
#########################################################

#############################
######## Simple DAVID GO
#############################
write.csv(diff_w8_up,"diff_w8_up.csv",quote = F,row.names = F)
write.csv(diff_w16_up,"diff_w16_up.csv",quote = F,row.names = F)
write.csv(diff_w32_up,"diff_w32_up.csv",quote = F,row.names = F)
write.csv(diff_w48_up,"diff_w48_up.csv",quote = F,row.names = F)

david2gp("w8up_go_有脂肪积累.csv","Week8_KO_up")
ggsave("w8up_go_有脂肪积累.pdf",dpi=300,width = 8,height = 3)
david2gp("w16up_炎症.csv","Week16_KO_up")
ggsave("w16up_炎症.pdf",dpi=300,width = 8,height = 3)
david2gp("w32_免疫衰老.csv","Week32_KO_up")
ggsave("w32_免疫衰老.pdf",dpi=300,width = 8,height = 3)
david2gp("w48_脂代谢免疫衰老.csv","Week48_KO_up")
ggsave("w48_免疫衰老.pdf",dpi=300,width = 8,height = 3)
















###################################################################
################################# Archived Codes
##################################################################

# ################################################
# #################### WGCNA
# ################################################
# 
# ##########Input: expr_data & trait_data
# datExpr0 <- t(normat[c(diff_up,diff_down),])
# datTrait=data.frame(
#   sample = rownames(datExpr0),
#   sample_class =str_sub(rownames(datExpr0),1,-2),
#   stringsAsFactors = F)
# 
# ##########Test and h-clust
# gsg = goodSamplesGenes(datExpr0, verbose = 3);
# gsg$allOK
# sampleTree = hclust(dist(datExpr0), method = "average")
# par(cex = 0.6)
# par(mar = c(0,4,2,0))
# plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
#      cex.axis = 2, cex.main = 2,cex=2)
# abline(h = 105, col = "red");
# clust = cutreeStatic(sampleTree, cutHeight = 80000, minSize = 10)
# table(clust)
# keepSamples = (clust==1)
# datExpr = datExpr0[keepSamples, ]
# nGenes = ncol(datExpr)
# nSamples = nrow(datExpr)
# 
# ##########Choose threhold
# powers = c(c(1:10), seq(from = 12, to=20, by=2))
# sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
# 
# par(mfrow = c(1,2));
# cex1 = 0.9;
# plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
#      main = paste("Scale independence"));
# text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
#      labels=powers,cex=cex1,col="red");
# abline(h=0.85,col="red")
# plot(sft$fitIndices[,1], sft$fitIndices[,5],
#      xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
#      main = paste("Mean connectivity"))
# text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
# 
# ################One-step network construction and module detection
# net = blockwiseModules(datExpr, power = 10, maxBlockSize = 6000,
#                        TOMType = "signed", minModuleSize = 30,
#                        reassignThreshold = 0, mergeCutHeight = 0.25,
#                        numericLabels = TRUE, pamRespectsDendro = FALSE,
#                        saveTOMs = TRUE,
#                        saveTOMFileBase = "AS-green-FPKM-TOM",
#                        verbose = 3)
# table(net$colors)
# 
# ############### Plot the dendrogram and the module colors underneath
# sizeGrWindow(12, 9)
# mergedColors = labels2colors(net$colors)
# plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
#                     "Module colors",
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05)
# 
# 
# ############### Extract genes from colorful modules
# design=model.matrix(~0+datTrait$sample_class)
# rownames(design)=datTrait$sample
# moduleColors <- labels2colors(net$colors)
# 
# MEs0=moduleEigengenes(datExpr,moduleColors)$eigengenes
# MEs=orderMEs(MEs0)
# moduleTraitCor = cor(MEs, design , use = "p")
# moduleTraitCor <- moduleTraitCor[,c(6,10,2,4,8,5,9,1,3,7)] ### <-- This is the order of trait cor heatmap.
# moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# 
# sizeGrWindow(10,6)
# textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
#                    signif(moduleTraitPvalue, 1), ")", sep = "")
# dim(textMatrix) = dim(moduleTraitCor)
# par(mar = c(6, 8.5, 3, 3))
# 
# 
# 
# labeledHeatmap(Matrix = moduleTraitCor,
#                xLabels = str_sub(colnames(moduleTraitCor),22,-1),
#                yLabels = names(MEs),
#                ySymbols = names(MEs),
#                colorLabels = FALSE,
#                colors = greenWhiteRed(50),
#                textMatrix = textMatrix,
#                setStdMargins = FALSE,
#                cex.text =0.9,
#                zlim = c(-1,1),
#                main = paste("Module-trait relationships"))
# 
# ############### Extract genes from colorful modules
# relation_matrix=moduleTraitCor
# #colnames(relation_matrix)=rownames(datExpr)
# 
# #modules=apply(relation_matrix,2,function(x) rownames(relation_matrix)[which.max(x)])
# #modules_matrix=data.frame(names(modules),modules,stringsAsFactors = F)
# #sample_num=length(rownames(modules_matrix))   <---选最大可以相详细参考刘文举的code
# 
# gene_names=colnames(datExpr)
# gene_color <- data.frame(name = gene_names,color = moduleColors)
# color2gene <- function(col){
#   as.character(gene_color[which(gene_color$color == col),]$name)
# }
# high_cor <- c(color2gene("turquoise"),color2gene("purple"),color2gene("pink"),color2gene("yellow"))
# 
# 
# 
# mean_scale_normat <- t(apply(mean_normat,1,function(x) (x-mean(x))/sd(x) ))
# pheatmap(mean_scale_normat[high_cor,grep("KO",colnames(mean_scale_normat))],cluster_cols = F,show_rownames = F)
# ######切一切
# a = data.frame(mean_scale_normat[high_cor,grep("KO",colnames(mean_scale_normat))])
# write.csv(rownames(a[a$W32_KO-a$W16_KO >= 1,]),"W32ko>W16ko.csv",quote = F,row.names = F)
# 
# ####
# clust=pheatmap(mean_scale_normat[high_cor,grep("KO",colnames(mean_scale_normat))],cluster_cols = F,show_rownames = F)
# ann_row = data.frame(
#   ClassGene = factor(paste0('Cluster',cutree(clust$tree_row,3)))
# )
# ann_row = cbind(rownames(mean_scale_normat[high_cor,grep("KO",colnames(mean_scale_normat))]),ann_row)
# colnames(ann_row) <- c("gene","cluster")
# #write.csv(as.character(ann_row[which(ann_row$cluster == "Cluster1"),]$gene),"ttt.csv")
# #pheatmap(mean_scale_normat[as.character(ann_row[which(ann_row$cluster == "Cluster2"),]$gene),grep("KO",colnames(mean_scale_normat))],cluster_cols = F,show_rownames = F)
# pheatmap(mean_scale_normat[high_cor,grep("WT",colnames(mean_scale_normat))],cluster_cols = F,show_rownames = F)
# 

# ############# Accumulated 
# two_geneset_mean_accu <- data.table(two_geneset_mean,class =paste0(two_geneset_mean$treat,"_",two_geneset_mean$cate))
# fat_accu <- data.frame(two_geneset_mean_accu[,cumsum(expr),by=class])
# fat_accu <- fat_accu[grep("fat$",fat_accu$class),]
# colnames(fat_accu) <- c('class','value')
# rownames(fat_accu) <- c(1:10)
# fat_accu$time <- factor(rep(c("W4","W8","W16","W32","W48"),2),levels = c("W4","W8","W16","W32","W48"))
# fat_ko_accu <- data.frame(value = fat_accu[6:10,]$value - fat_accu[1:5,]$value,
#                           time = factor(c("W4","W8","W16","W32","W48"),levels = c("W4","W8","W16","W32","W48")),
#                           group = rep(1,5))
# 
# 
# inflam_index <- data.frame(value = inflam_mean[6:10] - inflam_mean[1:5],row.names = c(1:5))
# inflam_index$time <- factor(c("W4","W8","W16","W32","W48"),levels = c("W4","W8","W16","W32","W48"))
# inflam_index$group <- rep(1,5)
# ggplot(fat_ko_accu) + geom_line(aes(x=time,y=value,group=group)) +
#   geom_line(data=inflam_index,aes(x=time,y=value,group=group),alpha=.3) + gran_theme + scale_color_npg() + scale_fill_npg() 
# 


################## Vennplot for intersection.
#install.packages("VennDiagram")
library(grid)
library(VennDiagram)
A = c(diff_w4_up)
B = c(diff_w8_up)
C = c(diff_w16_up)
D = c(diff_w32_up)
E = c(diff_w48_up)
D1<-venn.diagram(list(w4_up=A,w8_up=B,w16_up=C,w32_up=D,w48_up=E),filename=NULL,lwd=1,lty=1,col=c("LightGray","Cyan","Blue","magenta","red"),
                 fill=c("LightGray","Cyan","Blue","magenta","red"),rotation.degree=0,
                 main = "Up-regulated genes within timescale")
grid.draw(D1)

################## Matrix for intersection.
result_mat  <- matrix(rep(0,25),nrow=5,ncol=5)
a <- qpcR:::cbind.na(diff_w4_up=A,diff_w8_up=B,diff_w16_up=C,diff_w32_up=D,diff_w48_up=E)

for(i in 1:ncol(a)){
  for(j in 1:ncol(a)){
    result_mat[i,j] = length(intersect(as.character(a[,i]),as.character(a[,j])))
  }
}
colnames(result_mat) <- rownames(result_mat) <- colnames(a)
write.csv(result_mat,"up_venn_to_matrix.csv",quote=F)



A = c(diff_w4_down)
B = c(diff_w8_down)
C = c(diff_w16_down)
D = c(diff_w32_down)
E = c(diff_w48_down)
D1<-venn.diagram(list(w4=A,w8=B,w16=C,w32=D,w48=E),filename=NULL,lwd=1,lty=1,col=c("LightGray","Cyan","Blue","magenta","red"),
                 fill=c("LightGray","Cyan","Blue","magenta","red"),rotation.degree=0,
                 main = "down-regulated genes within timescale")
grid.draw(D1)


################## Matrix for intersection.
result_mat  <- matrix(rep(0,25),nrow=5,ncol=5)
a <- qpcR:::cbind.na(diff_w4_down=A,diff_w8_down=B,diff_w16_down=C,diff_w32_down=D,diff_w48_down=E)

for(i in 1:ncol(a)){
  for(j in 1:ncol(a)){
    result_mat[i,j] = length(intersect(as.character(a[,i]),as.character(a[,j])))
  }
}
colnames(result_mat) <- rownames(result_mat) <- colnames(a)
write.csv(result_mat,"down_venn_to_matrix.csv",quote=F)



############################################################
############################### Related to clinical target
############################################################
clinical_target <- read.csv("target_clinical.csv",header=F)
clinical_target_mean <- na.omit(mean_normat[as.character(clinical_target$V1),])
temp = clinical_target_mean[,c(6:10)]-clinical_target_mean[,c(1:5)]
a = pheatmap(temp,breaks = c(seq(-2, 0, length.out = 50) ,seq(0.1, 5.5, length.out = 50) ),color= colorRampPalette(c("blue", "white", "red"))(100),cluster_cols = F,border_color = NA,fontsize_row =6,fontsize_col =8)
ggsave("clinical_heatmap.pdf",a$gtable,dpi=300,width = 3,height = 4)


######################################
############### Individial codes
######################################
#"LightGray","Cyan","Blue","magenta","red"
david2gp <- function(csv,title_def) {
  csv <- read.csv(paste0("./all_go/",csv),header = T,col.names = c("Cate","Desc.","Pvalue","LogPvalue","dir"))
  csv$Desc. <- str_split_fixed(csv$Desc.,"~",2)[,2]
  p <- ggplot(csv,aes(x=reorder(Desc.,LogPvalue),y=LogPvalue)) + 
    geom_bar(stat = "identity",fill="LightGray", width = 0.6)+coord_flip()
  p <- p + theme(panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 text = element_text(),
                 legend.text = element_text(size=1),
                 axis.text = element_text(size=10),
                 axis.line.x = element_line(colour = "black"),
                 axis.title.x = element_text(size = 10,  face = "bold"),
                 axis.title.y = element_text(size = 10, face = "bold"),
                 axis.ticks.y = element_blank()) + scale_y_continuous(expand = c(0, 0))+ xlab("Terms")+
    labs(title=title_def)+ theme(plot.title = element_text(hjust = 0.5,size=10,face = "bold")) 
  return(p)
}

david2gp("rat_w4_up_go-zhu.csv","W4_up_genes_go")
ggsave("rat_w4_up_go.pdf",dpi=300,width = 8,height = 2)
david2gp("rat_w8_up_go-zhu.csv","W8_up_genes_go")
ggsave("rat_w8_up_go.pdf",dpi=300,width = 8,height = 2.2)
david2gp("rat_w16_up_go-zhu.csv","W16_up_genes_go")
ggsave("rat_w16_up_go.pdf",dpi=300,width = 8,height = 2.4)
david2gp("rat_w32_up_go-zhu.csv","W32_up_genes_go")
ggsave("rat_w32_up_go.pdf",dpi=300,width = 8,height = 2.2)
david2gp("rat_w48_up_go-zhu.csv","W48_up_genes_go")
ggsave("rat_w48_up_go.pdf",dpi=300,width = 8,height = 2.4)


david2gp("rat_w4_down_go-zhu.csv","W4_down_genes_go")
ggsave("rat_w4_down_go.pdf",dpi=300,width = 8,height = 1.6)
david2gp("rat_w8_down_go-zhu.csv","W8_down_genes_go")
ggsave("rat_w8_down_go.pdf",dpi=300,width = 8,height = 2)
david2gp("rat_w16_down_go-zhu.csv","W16_down_genes_go")
ggsave("rat_w16_down_go.pdf",dpi=300,width = 8,height = 2.2)
david2gp("rat_w32_down_go-zhu.csv","W32_down_genes_go")
ggsave("rat_w32_down_go.pdf",dpi=300,width = 8,height = 2)
david2gp("rat_w48_down_go-zhu.csv","W48_down_genes_go")
ggsave("rat_w48_down_go.pdf",dpi=300,width = 8,height = 2)

diff_w4 <- data.frame(gene=c(diff_w4_up,diff_w4_down),property=rep(c("up","down"),c(length(diff_w4_up),length(diff_w4_down))))
write.csv(diff_w4,"diff_w4.csv",quote = F)
diff_w8 <- data.frame(gene=c(diff_w8_up,diff_w8_down),property=rep(c("up","down"),c(length(diff_w8_up),length(diff_w8_down))))
write.csv(diff_w8,"diff_w8.csv",quote = F)
diff_w16 <- data.frame(gene=c(diff_w16_up,diff_w16_down),property=rep(c("up","down"),c(length(diff_w16_up),length(diff_w16_down))))
write.csv(diff_w16,"diff_w16.csv",quote = F)
diff_w32 <- data.frame(gene=c(diff_w32_up,diff_w32_down),property=rep(c("up","down"),c(length(diff_w32_up),length(diff_w32_down))))
write.csv(diff_w32,"diff_w32.csv",quote = F)
diff_w48 <- data.frame(gene=c(diff_w48_up,diff_w48_down),property=rep(c("up","down"),c(length(diff_w48_up),length(diff_w48_down))))
write.csv(diff_w48,"diff_w48.csv",quote = F)

