### Volcano plot ###

library(ggplot2)
library(ggrepel)

adjust.pvalue<-as.numeric(0.05)
fold.change<-as.numeric(2)


#======================================================================================================
# Simplified version
VolcanoPlot1<-function(){
  volcano=as.data.frame(res)
  volcano$padj=ifelse(is.na(volcano$padj), 1, volcano$padj)
  volcano$Group=rep(0)

  # seperate groups (significantly up-/down-regulated or non-significant)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & volcano$log2FoldChange <= -log2(fold.change)), "Group"]=rep(1)  # down-regulated
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & volcano$log2FoldChange >= log2(fold.change)), "Group"]=rep(2)  # up-regulated
  volcano[which(-log10(volcano$padj) <= -log10(adjust.pvalue) | (-log10(volcano$padj) > -log10(adjust.pvalue) & abs(volcano$log2FoldChange) < log2(fold.change))), "Group"]=rep(3)  # non-significant
  volcano$Group<-as.factor(volcano$Group)
  volcano.label<-volcano[which(volcano$Group!=3),]  # groups to show gene labels

  # plot
  volcano.p=ggplot(volcano, aes(log2FoldChange, -log10(padj)))+geom_point(aes(color=Group), alpha=0.7)+scale_color_manual(name="Significance", breaks=c("2", "1", "3"), labels=c("Up", "Down", "Non"), values=c("1"="royalblue3", "2"="firebrick3", "3"="grey55"))+geom_text_repel(data=volcano.label, aes(label=rownames(volcano.label)), size=2, color="black")

  volcano.p=volcano.p+labs(title={volcano.title}, x="log2foldchange", y="-log10(padj)")+geom_hline(yintercept=-log10(adjust.pvalue), linetype=2, color="black", size=0.3)+geom_vline(xintercept=c(-log2(fold.change), log2(fold.change)), linetype=2, color="black", size=0.3)+theme_bw()+theme(text=element_text(size=12, color="black"), panel.grid=element_blank(), panel.border=element_rect(size=0.5), axis.line=element_line(color="black", size=0.2), axis.ticks=element_line(color="black"), axis.text=element_text(color="black"))

  volcano.p
}
VolcanoPlot1(res, "treatment vs control")


#======================================================================================================
# Detailed version
VolcanoPlot2<-function(res, x.min, x.max, y.max, volcano.title){
  volcano=as.data.frame(res)
  volcano$padj=ifelse(is.na(volcano$padj), 1, volcano$padj)
  volcano$Group=rep(0)

  # seperate groups
  # group1: down-regulated within normal fc and padj range (blue dots)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange >= x.min & volcano$log2FoldChange <= -log2(fold.change)), "Group"]=rep(1)

  # group2: up-regulated within normal fc and padj range (red dots)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange >= log2(fold.change) & volcano$log2FoldChange <= x.max), "Group"]=rep(2)

  # group3: non-significant within normal fc or padj range (grey dots)
  volcano[which((-log10(volcano$padj) <= -log10(adjust.pvalue) & abs(volcano$log2FoldChange) <= x.max) | (-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) < y.max & abs(volcano$log2FoldChange) <= log2(fold.change))), "Group"]=rep(3)

  # group4: non-significant outside of normal fc or padj range (grey trangle)
  volcano[which(-log10(volcano$padj) <= -log10(adjust.pvalue) & volcano$log2FoldChange < x.min), c("log2FoldChange", "Group")]=list(log2FoldChange = x.min, Group=4)
  volcano[which(-log10(volcano$padj) <= -log10(adjust.pvalue) & volcano$log2FoldChange > x.max), c("log2FoldChange", "Group")]=list(log2FoldChange = x.max, Group=4)
  volcano[which(-log10(volcano$padj) > y.max & abs(volcano$log2FoldChange) < log2(fold.change)), c("padj", "Group")]=list(padj = 10^(-y.max), Group=4)

  # group5: down-regulated outside of normal fc or padj range (blue trangle)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange < x.min), c("log2FoldChange", "Group")]=list(log2FoldChange = x.min, Group=5)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange < x.min), c("log2FoldChange", "padj", "Group")]=list(log2FoldChange = x.min, padj = 10^(-y.max), Group=5)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange <= -log2(fold.change) & volcano$log2FoldChange >= x.min), c("padj", "Group")]=list(padj = 10^(-y.max), Group=5)

  # group6: up-regulated outside of normal fc or padj range (red trangle)
  volcano[which(-log10(volcano$padj) > -log10(adjust.pvalue) & -log10(volcano$padj) <= y.max & volcano$log2FoldChange > x.max), c("log2FoldChange", "Group")]=list(log2FoldChange = x.max, Group=6)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange > x.max), c("log2FoldChange", "padj", "Group")]=list(log2FoldChange = x.max, padj = 10^(-y.max), Group=6)
  volcano[which(-log10(volcano$padj) > y.max & volcano$log2FoldChange >= log2(fold.change) & volcano$log2FoldChange <= x.max), c("padj", "Group")]=list(padj = 10^(-y.max), Group=6)

  volcano$Group<-as.factor(volcano$Group)

  # which groups to show labels
  volcano.label=volcano[which(volcano$Group!=3 & volcano$Group!=4),]

  # plot
  volcano.p=ggplot(volcano, aes(log2FoldChange, -log10(padj), color=Group))+geom_point(data=volcano[which(volcano$Group==1 | volcano$Group==2 | volcano$Group==3),], aes(color=Group), alpha=0.7)+scale_color_manual(name="Significance", breaks=c("2", "1", "3"), labels=c("Up", "Down", "Non"), values=c("1"="royalblue3", "2"="firebrick3", "3"="grey55"))+geom_point(data=volcano[which(volcano$Group==4),], shape=2, color="gray55", alpha=0.7)+geom_point(data=volcano[which(volcano$Group==5),], shape=2, color="royalblue3", alpha=0.7)+geom_point(data=volcano[which(volcano$Group==6),], shape=2, color="firebrick3", alpha=0.7)+geom_text_repel(data=volcano.label, aes(label=rownames(volcano.label)), size=2, color="black")

  volcano.p=volcano.p+labs(title={volcano.title}, x="log2foldchange", y="-log10(padj)")+geom_hline(yintercept=-log10(adjust.pvalue), linetype=2, color="black", size=0.3)+geom_vline(xintercept=c(-log2(fold.change), log2(fold.change)), linetype=2, color="black", size=0.3)+theme_bw()+theme(text=element_text(size=12, color="black"), panel.grid=element_blank(), panel.border=element_rect(size=0.5), axis.line=element_line(color="black", size=0.2), axis.ticks=element_line(color="black"), axis.text=element_text(color="black"))+xlim({x.min}, {x.max})+ylim(0, {y.max})

  volcano.p
}

VolcanoPlot2(res, -5, 5, 200, "treatment vs control")


#======================================================================================================
# Explanation:
# res=results(dds, alpha=adjust.pvalue), result from DESeq2.
# adjust.pvalue, adjusted pvalue for defining significant DEGs used by DESeq2.
# fold.change, fold change for defining significant DEGs used by DESeq2.
# x.min, x-axis left limit.
# x.max, x-axis right limit.
# y.max, y-axis up/maximum limit.
# volcano.title

