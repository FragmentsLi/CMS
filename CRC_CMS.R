###############################
#
#构建以无重复BAM的样本Barcode为主的信息表格
#
###############################
#读取微生物数据/RNA-Seq数据/CMS数据
#构建表格
library(CMSclassifier)

rna<-read.csv('/data3/data/JN/CRC/TCGAbiolinks/TCGAbiolinks_combat_data_geneid_tumor.csv',row.names=1,check.names=F)
cms_label<-read.delim('/data3/data/JN/CRC/CMS/cms_labels_public_all.txt')
cms_label<-subset(cms_label,cms_label$dataset=='tcga')
rownames(cms_label)<-cms_label$sample

#构建表格,第一列为肿瘤样本Barcode，共594个
#第二列为对应患者ID，
#第三列为CMS金标准预测结果，
#第四列为使用CMS的自己预测结果
library(stringr)
#patinet_id的顺序和barcode一致，也为594个
patient_id<-unlist(lapply(colnames(rna),function(x){substr(x,1,12)}))
cms_class<-unlist(lapply(patient_id,function(x)
  if(x %in% cms_label$sample)return(cms_label[x,5])
  else return('')))

#按RF(暂选择)
Rfcms <- CMSclassifier::classifyCMS(rna,method="RF")[[3]]
#按barcode顺序取预测值和最近值
cms_predict<-unlist(lapply(colnames(rna),function(x){Rfcms[x,'RF.predictedCMS']}))
cms_nearest<-unlist(lapply(colnames(rna),function(x){Rfcms[x,'RF.nearestCMS']}))

#按SSP
SScms <- CMSclassifier::classifyCMS(data,method="SSP")[[3]]
#按barcode顺序取预测值和最近值
cms_predict<-unlist(lapply(colnames(rna),function(x){SScms[x,'SSP.predictedCMS']}))
cms_nearest<-unlist(lapply(colnames(rna),function(x){SScms[x,'SSP.nearestCMS']}))

#金标准优先，如果没有加入预测值
cms_merge<-unlist(lapply(patient_id,function(x)
  if(x %in% cms_label$sample)return(cms_label[x,5])
  else return(Rfcms[x,'RF.predictedCMS'])))
cms_merge[which(is.na(cms_merge))]='NOLBL'

info<-data.frame(barcode=barcode,patient_id=patient_id,cms_class=cms_class,
                 cms_predict=cms_predict,cms_nearest=cms_nearest,cms_merge=cms_merge)
write.csv(info,'/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=F)



###############################
#
#按照CMS的分型做生存分析
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
#选RF
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
#选SSP
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_ssp.csv',row.names=1)

#考虑平台
meta_RNA<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1)
info<-info[rownames(subset(meta_RNA,meta_RNA$platform=='Illumina HiSeq')),]


info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
survival_info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_survival_nmf.csv',row.names=1)

meta<-as.data.frame(matrix(ncol=4,nrow=length(unique(info$patient_id))))
colnames(meta)<-c('sample_id','time','event','group')
rownames(meta)<-unique(info$patient_id)
#3个患者有2个replicate，所得CMS一致
for(i in 1:nrow(meta)){
  meta[i,1]=unique(info$patient_id)[i]
  meta[i,2]=survival_info[meta[i,1],"OS.time"]
  meta[i,3]=survival_info[meta[i,1],"OS"]
  meta[i,4]=subset(info,info$patient_id==meta[i,1])[1,'cms_merge']
}

#meta<-subset(meta,meta$group!='Mixed-CMS')
meta$time<-as.numeric(meta$time)
meta$event<-as.numeric(meta$event)


library(survival)
library(survminer)
sfit <- survfit(Surv(time, event) ~ group,
                data = meta)

res<-pairwise_survdiff(Surv(time, event) ~ group,
                          data = meta, p.adjust.method = "BH",rho = 0)
res

pdf('survival_cms_rf_all_hiseq.pdf',width = 6,height =6,onefile = F)
ggsurvplot(sfit,pval = T,palette = "jco")
dev.off() 

#考虑单个stage
survival_info$ajcc_pathologic_tumor_stage[survival_info$ajcc_pathologic_tumor_stage=='Stage IA']='Stage I'
survival_info$ajcc_pathologic_tumor_stage[survival_info$ajcc_pathologic_tumor_stage=='Stage IIA'|survival_info$ajcc_pathologic_tumor_stage=='Stage IIB'|survival_info$ajcc_pathologic_tumor_stage=='Stage IIC']='Stage II'
survival_info$ajcc_pathologic_tumor_stage[survival_info$ajcc_pathologic_tumor_stage=='Stage IIIA'|survival_info$ajcc_pathologic_tumor_stage=='Stage IIIB'|survival_info$ajcc_pathologic_tumor_stage=='Stage IIIC']='Stage III'
survival_info$ajcc_pathologic_tumor_stage[survival_info$ajcc_pathologic_tumor_stage=='Stage IVA'|survival_info$ajcc_pathologic_tumor_stage=='Stage IVB']='Stage IV'

stage_patient<-rownames(subset(survival_info,survival_info$ajcc_pathologic_tumor_stage=='Stage III'))
meta<-as.data.frame(matrix(ncol=4,nrow=length(intersect(unique(info$patient_id),stage_patient))))
colnames(meta)<-c('sample_id','time','event','group')
rownames(meta)<-intersect(unique(info$patient_id),stage_patient)
#3个患者有2个replicate，所得CMS一致
for(i in 1:nrow(meta)){
  meta[i,1]=unique(info$patient_id)[i]
  meta[i,2]=survival_info[meta[i,1],27]
  meta[i,3]=survival_info[meta[i,1],26]
  meta[i,4]=subset(info,info$patient_id==meta[i,1])[1,2]
}

meta<-subset(meta,meta$group!='Mixed-CMS')
meta<-subset(meta,meta$group!='')

meta$time<-as.numeric(meta$time)
meta$event<-as.numeric(meta$event)


library(survival)
library(survminer)
sfit <- survfit(Surv(time, event) ~ group,
                data = meta)

res<-pairwise_survdiff(Surv(time, event) ~ group,
                       data = meta, p.adjust.method = "BH",rho = 0)
res

pdf('survival_cms_gs_nomix_stageIII.pdf',width = 6,height =6,onefile = F)
ggsurvplot(sfit,pval = T,palette = "jco")
dev.off() 


###############################
#
#按照CMS的分型做菌种α多样性分析
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
#microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)
survival_info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_survival_nmf.csv',row.names=1)

#考虑平台
meta_RNA<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1)
#microbes<-microbes[rownames(subset(meta_RNA,meta_RNA$platform=='Illumina GA')),]

#使用NMF的矫正负值方法矫正VOOM-SNM负值问题,加上最小值使得所有数为正值
library(NMF)
microbes<-nneg(as.matrix(microbes),method='min')
write.csv(microbes,'/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode_nng.csv',row.names=T)


#物种丰度抽平分析（暂无法做）
library(vegan)
head(rowSums(microbes))
set.seed(315)
otu_Flattening = as.data.frame(t(vegan::rrarefy(microbes, min(rowSums(microbes)))))
head(colSums(otu_Flattening))

#Chao1/ACE指数大小和生存曲线(暂无法做)
#Chao1  <- estimateR(t(otu_Flattening))[2,]
#ACE  <- estimateR(t(otu_Flattening))[4,]

#Shannon/Simpson指数大小和生存曲线
#Shannon <- diversity(t(otu_Flattening), index = 'shannon', base = exp(1))
#Gini_simpson  <- diversity(t(otu_Flattening), index = 'simpson')

Shannon <- diversity(microbes, index = 'shannon', base = exp(1))
Gini_simpson  <- diversity(microbes, index = 'simpson')

#构建meta
meta<-as.data.frame(matrix(ncol=6,nrow=length(Shannon)))
#colnames(meta)<-c('sample_id','time','event','CMS','Chao1','ACE','Shannon','Gini_simpson')
colnames(meta)<-c('sample_id','time','event','CMS','Shannon','Gini_simpson')
rownames(meta)<-names(Shannon)
for(i in 1:nrow(meta)){
  meta[i,1]=substr(names(Shannon)[i],1,12)
  meta[i,2]=survival_info[meta[i,1],"OS.time"]
  meta[i,3]=survival_info[meta[i,1],"OS"]
  meta[i,4]=info[rownames(meta)[i],"cms_merge"]
  #meta[i,5]=Chao1[i]
  #meta[i,6]=ACE[i]
  meta[i,5]=Shannon[i]
  meta[i,6]=Gini_simpson[i]
}

#绘制不同指数的boxplot
library(ggpubr)
library(RColorBrewer)
p <- ggplot(meta, aes(x = CMS , y = Gini_simpson)) +
  geom_boxplot(aes(fill = CMS), show.legend = FALSE, width = 0.4) +  #绘制箱线图
  scale_fill_manual(values = brewer.pal(n=5,name='Set2')) +  #箱线图的填充色
  geom_point(size = 1,alpha = 0.5) +  #绘制样本点
  #geom_line(aes(group = Individual), color = 'gray50', lwd = 0.5) +  #绘制配对样本间连线
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major=element_line(colour=NA), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12, color = 'black',face="bold")) +
  labs(x = '', y = 'Gini_simpson')

my_comparisons <- list(c("CMS1","CMS2"),c("CMS1","CMS3"),c("CMS1","CMS4"),c("CMS1","Mixed-CMS"),
                       c("CMS2","CMS3"),c("CMS2","CMS4"),c("CMS2","Mixed-CMS"),c("CMS3","CMS4"),
                       c("CMS3","Mixed-CMS"),c("CMS4","Mixed-CMS"))
p<-p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",vjust=0.1,size=2)
#p<-p +facet_wrap(~variable,nrow = 1)



ggsave(p, file='type_comparison_Gini_simpson_vs_display.pdf', width=6, height=6)

#只留下CMS4和其他差异的展示用图
p <- ggplot(meta, aes(x = CMS , y = Gini_simpson)) +
  geom_boxplot(aes(fill = CMS), show.legend = FALSE, width = 0.4) +  #绘制箱线图
  scale_fill_manual(values = brewer.pal(n=5,name='Set2')) +  #箱线图的填充色
  geom_point(size = 1,alpha = 0.2) +  #绘制样本点
  #geom_line(aes(group = Individual), color = 'gray50', lwd = 0.5) +  #绘制配对样本间连线
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major=element_line(colour=NA), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12, color = 'black',face="bold")) +
  labs(x = '', y = 'Gini Simpson')
  #ylim(c(0.75,1.5))

my_comparisons <- list(c("CMS1","CMS4"),c("CMS2","CMS4"),c("CMS3","CMS4"),c("CMS4","NOLBL"))
p<-p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",vjust=0.1,size=3)
#p<-p +facet_wrap(~variable,nrow = 1)

ggsave(p, file='type_comparison_Gini_simpson_CMS4.pdf', width=5, height=5)

#按照CMS的分型做菌种α多样性分析（全部基于金标准数据/全部基于自己预测数据）
#构建meta
meta<-as.data.frame(matrix(ncol=8,nrow=length(Chao1)))
colnames(meta)<-c('sample_id','time','event','CMS','Chao1','ACE','Shannon','Gini_simpson')
rownames(meta)<-names(Chao1)
for(i in 1:nrow(meta)){
  meta[i,1]=substr(names(Chao1)[i],1,12)
  meta[i,2]=survival_info[meta[i,1],"OS.time"]
  meta[i,3]=survival_info[meta[i,1],"OS"]
  #meta[i,4]=info[rownames(meta)[i],"cms_class"]#金标准
  meta[i,4]=info[rownames(meta)[i],"cms_predict"]#自己预测
  meta[i,5]=Chao1[i]
  meta[i,6]=ACE[i]
  meta[i,7]=Shannon[i]
  meta[i,8]=Gini_simpson[i]
}
#meta<-subset(meta,meta$CMS!='')#金标准
meta$CMS[which(is.na(meta$CMS))]<-'NOLBL'#自己预测

#绘制不同指数的boxplot
library(ggpubr)
library(RColorBrewer)
p <- ggplot(meta, aes(x = CMS , y = Shannon)) +
  geom_boxplot(aes(fill = CMS), show.legend = FALSE, width = 0.4) +  #绘制箱线图
  scale_fill_manual(values = brewer.pal(n=5,name='Set3')) +  #箱线图的填充色
  geom_point(size = 1,alpha = 0.5) +  #绘制样本点
  #geom_line(aes(group = Individual), color = 'gray50', lwd = 0.5) +  #绘制配对样本间连线
  theme_bw() + 
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major=element_line(colour=NA), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12, color = 'black',face="bold")) +
  labs(x = '', y = 'Shannon')

my_comparisons <- list(c("CMS1","CMS2"),c("CMS1","CMS3"),c("CMS1","CMS4"),c("CMS1","NOLBL"),
                       c("CMS2","CMS3"),c("CMS2","CMS4"),c("CMS2","NOLBL"),c("CMS3","CMS4"),
                       c("CMS3","NOLBL"),c("CMS4","NOLBL"))
p<-p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",vjust=0.1,size=2)
#p<-p +facet_wrap(~variable,nrow = 1)

ggsave(p, file='type_comparison_Shannon_pr.pdf', width=10, height=16)


###############################
#
#按照CMS的分型做菌种β多样性分析（PCoA）
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
#microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)
#读取非负微生物丰度矩阵
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode_nng.csv',row.names=1)
meta_RNA<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1)
#去掉Mixed-CMS
info<-subset(info,info$cms_merge!='Mixed-CMS')
microbes<-microbes[rownames(info),]


#基于VOOM-SNM的数据
library(vegan)
library(ggplot2)
library(ggrepel)
library(ade4)
distance<-vegdist(microbes,method='bray')

pcoa<-cmdscale(distance,k=(nrow(microbes)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]
names(plot_data)<-c('PCoA1','PCoA2')
eig<-pcoa$eig
data<-plot_data[match(rownames(info),rownames(plot_data)),]
data<-data.frame(info,plot_data,meta_RNA[rownames(info),])
head(data)
colnames(data)[which(colnames(data)=='cms_merge')]<-'Type'
colnames(data)[which(colnames(data)=='platform')]<-'Platform'


library(RColorBrewer)
#按CMS分型
p<-ggplot(data,aes(x=PCoA1,y=PCoA2,shape=Type,color=Type))+
  geom_point(alpha=1,size=1)+stat_ellipse(level=0.95,size=1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digit=4),"%)",sep=""),
       y=paste("PCoA2(",format(100*eig[2]/sum(eig),digit=4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  scale_color_manual(values = brewer.pal(n=5,name='Set2'))+
  theme(panel.background=element_rect(fill='white',color='black'),
        axis.title=element_text(size = 6, color = 'black',face="bold"),
        axis.text=element_text(size=6),
        legend.text=element_text(size=6),legend.title=element_text(size=6))+
  annotate('text',label='ANOSIM analysis: R=0.021, P-value=0.102',
           x=-0.05,y=0.05,size=2,hjust=0)+

#ANOSIM:9999
png('beta_diversity_vs_rf_nng_CMS.png',width = 1200,height = 900,res=300,pointsize=6)
p
dev.off() 

ggsave(p, file='beta_diversity_vs_rf_nng_CMS.pdf', width=7, height=6)


#按平台
p<-ggplot(data,aes(x=PCoA1,y=PCoA2,shape=Platform,color=Platform))+
  geom_point(alpha=1,size=3)+stat_ellipse(level=0.95,size=1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digit=4),"%)",sep=""),
       y=paste("PCoA1(",format(100*eig[2]/sum(eig),digit=4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  scale_color_manual(values = brewer.pal(n=5,name='Set2'))+
  theme(panel.background=element_rect(fill='white',color='black'),
        axis.title.x=element_text(size = 12, color = 'black',face="bold"),axis.title.y=element_text(size = 12, color = 'black',face="bold"),
        legend.text=element_text(size=10),legend.title=element_text(size=12))
ggsave(p, file='beta_diversity_vs_rf_nng_platform.pdf', width=7.5, height=6)


#基于原始数据和抽平后的数据
microbes2<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
meta_RNA<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1)
#物种丰度抽平分析
library(vegan)
library(ggplot2)
library(RColorBrewer)
head(rowSums(microbes2))
set.seed(315)
otu_Flattening = as.data.frame(vegan::rrarefy(microbes2, min(rowSums(microbes2))))
head(rowSums(otu_Flattening))

distance<-vegdist(microbes2,method='bray')
#distance<-vegdist(otu_Flattening,method='bray')
pcoa<-cmdscale(distance,k=(nrow(microbes2)-1),eig=TRUE)
#pcoa<-cmdscale(distance,k=(nrow(otu_Flattening)-1),eig=TRUE)
plot_data<-data.frame({pcoa$point})[1:2]
names(plot_data)<-c('PCoA1','PCoA2')
eig=pcoa$eig
data<-plot_data[match(rownames(info),rownames(plot_data)),]
data<-data.frame(info,plot_data,meta_RNA[rownames(info),])
colnames(data)[which(colnames(data)=='cms_merge')]<-'Type'
colnames(data)[which(colnames(data)=='platform')]<-'Platform'

library(RColorBrewer)
#按CMS分型
p<-ggplot(data,aes(x=PCoA1,y=PCoA2,shape=Type,color=Type))+
  geom_point(alpha=1,size=3)+stat_ellipse(level=0.95,size=1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digit=4),"%)",sep=""),
       y=paste("PCoA1(",format(100*eig[2]/sum(eig),digit=4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  scale_color_manual(values = brewer.pal(n=5,name='Set2'))+
  theme(panel.background=element_rect(fill='white',color='black'),
        axis.title.x=element_text(size = 12, color = 'black',face="bold"),axis.title.y=element_text(size = 12, color = 'black',face="bold"),
        legend.text=element_text(size=10),legend.title=element_text(size=12))


ggsave(p, file='beta_diversity_vs_rf_CMS.pdf', width=7, height=6)


#按平台
p<-ggplot(data,aes(x=PCoA1,y=PCoA2,shape=Platform,color=Platform))+
  geom_point(alpha=1,size=3)+stat_ellipse(level=0.95,size=1)+
  labs(x=paste("PCoA1(",format(100*eig[1]/sum(eig),digit=4),"%)",sep=""),
       y=paste("PCoA1(",format(100*eig[2]/sum(eig),digit=4),"%)",sep=""))+
  geom_vline(aes(xintercept=0),linetype='dotted')+
  geom_hline(aes(yintercept=0),linetype='dotted')+
  scale_color_manual(values = brewer.pal(n=5,name='Set2'))+
  theme(panel.background=element_rect(fill='white',color='black'),
        axis.title.x=element_text(size = 12, color = 'black',face="bold"),axis.title.y=element_text(size = 12, color = 'black',face="bold"),
        legend.text=element_text(size=10),legend.title=element_text(size=12))
ggsave(p, file='beta_diversity_vs_rf_platform.pdf', width=7.5, height=6)

###############################
#
#按照CMS的分型做菌种β多样性分析（NMDS）
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)

#基于raw data
library(vegan)
distance<-vegdist(microbes,method='bray')
nmds_dis<-metaMDS(distance,k=2)
#应力函数值一般<0.2
nmds_dis$stress

nmds_dis_site<-data.frame(nmds_dis$points)
nmds_dis_species<-wascores(nmds_dis$point,microbes)

library(ggplot2)
library(RColorBrewer)
abundance<-apply(microbes,2,sum)
abundance_top10<-names(abundance[order(abundance,decreasing=TRUE)][1:10])
species_top10<-data.frame(nmds_dis_species[abundance_top10,1:2])
species_top10$name<-unlist(lapply(rownames(species_top10),function(x){return(strsplit(x,'\\.')[[1]][length(strsplit(x,'\\.')[[1]])])}))

#分组信息
nmds_dis_site$name<-rownames(nmds_dis_site)
merged<-merge(nmds_dis_site,info,by="row.names",all.x=TRUE)

p <- ggplot(data = merged, aes(MDS1, MDS2)) +
  geom_point(size=2,aes(color =  cms_merge,shape = cms_merge)) +
  stat_ellipse(aes(fill = cms_merge), geom = 'polygon', level = 0.95, alpha = 0.1, show.legend = FALSE) +   #添加置信椭圆，注意不是聚类
  scale_color_manual(values =brewer.pal(n=5,name='Set3')[1:length(unique(info$cms_merge))]) +
  scale_fill_manual(values = brewer.pal(n=5,name='Set3')[1:length(unique(info$cms_merge))]) +
  theme(panel.grid.major = element_line(color = 'gray', size = 0.2), panel.background = element_rect(color = 'black', fill = 'transparent'), 
        plot.title = element_text(hjust = 0.5),legend.title = element_blank()) +
  #, legend.position = 'none'
  geom_vline(xintercept = 0, color = 'gray', size = 0.5) +
  geom_hline(yintercept = 0, color = 'gray', size = 0.5)+
  geom_text(data = species_top10, aes(label = name), color ="royalblue4", size = 4)
  #geom_text(data =merged, aes(label = Row.names,x =MDS1, y = MDS2), size=4, check_overlap = TRUE)
ggsave(p, file='beta_diversity_nmds.pdf', width=10, height=10)


###############################
#
#ANOSIM分析
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
#读取非负微生物丰度矩阵
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode_nng.csv',row.names=1)
meta_RNA<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1)
#去掉Mixed-CMS
info<-subset(info,info$cms_merge!='Mixed-CMS')
microbes<-microbes[rownames(info),]

library(vegan)
anosim=anosim(microbes, info[rownames(microbes),]$cms_merge, permutations=9999)
summary(anosim)
anosim_plot<-data.frame(dis=anosim$dis.rank,CMS=anosim$class.vec)

library(ggpubr)
library(RColorBrewer)
p <- ggplot(anosim_plot, aes(x = CMS , y = dis)) +
  geom_boxplot(aes(fill = CMS), show.legend = FALSE, width = 0.4) +  #绘制箱线图
  scale_fill_manual(values = brewer.pal(n=6,name='Set2')[c(6,1:5)]) +  #箱线图的填充色
  geom_point(size = 1,alpha = 0.5) +  #绘制样本点
  theme_bw() + 
  annotate('text',label='Anosim:\n
           R=0.01855:\n
           p-value=0.106',
           x=0.1,y=200000,size=4,hjust=0)+
  theme(panel.grid = element_blank(), panel.background = element_rect(colour = "black",size = 1),
        panel.grid.major=element_line(colour=NA), panel.grid.minor=element_blank(),
        plot.title = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 12, color = 'black'), axis.title = element_text(size = 12, color = 'black',face="bold")) +
  labs(x = '', y = 'Bray-Curtis Rank')+
  

#my_comparisons <- list(c("CMS1","CMS2"),c("CMS1","CMS3"),c("CMS1","CMS4"),c("CMS1","Mixed-CMS"),
#                       c("CMS2","CMS3"),c("CMS2","CMS4"),c("CMS2","Mixed-CMS"),c("CMS3","CMS4"),
#                       c("CMS3","Mixed-CMS"),c("CMS4","Mixed-CMS"))
#p<-p+stat_compare_means(comparisons = my_comparisons,method = "wilcox.test",vjust=0.1,size=2)
#p<-p +facet_wrap(~variable,nrow = 1)

ggsave(p, file='anosim_rf_vs.pdf', width=9, height=6)


###############################
#
#门(Phylum)和属（genus）层面上的微生物组分分析
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
#microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)
#读取原始数据
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
#读取矫正后归非负的数据
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode_nng.csv',row.names=1)
#取不同平台
annot<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1,check.names = 1)
microbes<-microbes[rownames(subset(annot,annot$platform=='Illumina HiSeq')),]
microbes<-microbes[rownames(subset(annot,annot$platform=='Illumina GA')),]


#取genus层面的名称
library(stringr)
names<-lapply(colnames(microbes),function(x){strsplit(x,split='\\.')})
names<-lapply(names,function(x){if(length(str_subset(x[[1]],"g__"))==0){return('g__nogenus')}else{return(str_subset(x[[1]],"g__"))}})
names<-unlist(lapply(names,function(x){strsplit(x[[1]],split='g__')[[1]][2]}))
colnames(microbes)<-names

#计算5种CMS中，所属样本在各genus中的平均丰度
genus_cms<-data.frame(matrix(ncol=ncol(microbes),nrow=length(unique(info$cms_merge))))
rownames(genus_cms)<-unique(info$cms_merge)
colnames(genus_cms)<-colnames(microbes)
for(i in 1:length(unique(info$cms_merge))){
  genus_cms[i,]<-colMeans(subset(microbes,rownames(microbes) %in% rownames(subset(info,info$cms_merge==unique(info$cms_merge)[i]))))
}
#绘制堆叠图，取5种CMS分型的平均genus表达累计前十的菌种
genus_top10<-names(sort(colSums(genus_cms),decreasing=T))[1:10]
genus_other<-rowSums(genus_cms[,setdiff(colnames(microbes),genus_top10)])
genus_new<-data.frame(genus_cms[,genus_top10],Other=genus_other)
genus_pro<-t(apply(genus_new,1,function(x){return(x/sum(x))}))
rownames(genus_pro)<-rownames(genus_new)
colnames(genus_pro)<-colnames(genus_new)

library(reshape2)
genus_pro<-data.frame('CMS'=rownames(genus_pro),genus_pro)
genus_pro<-genus_pro %>%
  melt(id.vars='CMS',variable.name='genus',value.name='proportion')

#画堆叠图
library(ggplot2)
library(RColorBrewer)
p<-ggplot(genus_pro,aes(x=CMS,y=100*proportion,fill=genus))+
  geom_col(position='stack',width=0.6)+
  #scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= brewer.pal(n=11,name='Set3'))+
  labs(x='',y='Relative Abundance(%)')+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) + #设置主题背景，根据自己的需求定制
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.title = element_blank(), legend.text=element_text(size=10))

ggsave(p, file='proportion_top10_genus_vs_rf.pdf', width=8, height=6)

#microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode_nng.csv',row.names=1)
#病毒合并，剩1153种
library(stringr)
virus<-unlist(lapply(colnames(microbes),function(x){grepl('Viruses',x)}))
virus_sum<-rowSums(microbes[,virus])
microbes<-data.frame(microbes[,!virus],p__Virus=virus_sum)
#取phylum层面的名称,没有phylum层面注释的独立为一种
names<-lapply(colnames(microbes),function(x){strsplit(x,split='\\.')})
names<-lapply(names,function(x){if(length(str_subset(x[[1]],"p__"))==0){return('p__nophylm')}else{return(str_subset(x[[1]],"p__"))}})
names<-unlist(lapply(names,function(x){strsplit(x[[1]],split='p__')[[1]][2]}))
#phylum层面加和，共34个phylum,1个virus，结果为每个样本中各phylum的总和abundance
phylum<-as.data.frame(matrix(ncol=length(unique(names)),nrow=nrow(microbes)))
colnames(phylum)<-unique(names)
rownames(phylum)<-rownames(microbes)
for(i in 1:ncol(phylum)){
  if(length(which(names==colnames(phylum)[i]))!=1)
    phylum[,i]=rowSums(microbes[,which(names==colnames(phylum)[i])])
  else
    phylum[,i]=microbes[,which(names==colnames(phylum)[i])]
}

phylum_raw_relative<-t(apply(phylum,1,function(x){return(x/sum(x))}))
colnames(phylum_raw_relative)<-colnames(phylum)
rownames(phylum_raw_relative)<-rownames(phylum)
write.csv(phylum_raw_relative,'/data3/data/JN/CRC/TCGA_microbes/phylum_vs_rf_relative.csv',row.names=T)
#计算5种CMS中，所属样本在各phylum中的平均丰度
phylum_cms<-data.frame(matrix(ncol=ncol(phylum),nrow=length(unique(info$cms_merge))))
rownames(phylum_cms)<-unique(info$cms_merge)
colnames(phylum_cms)<-colnames(phylum)
for(i in 1:length(unique(info$cms_merge))){
  phylum_cms[i,]<-colMeans(subset(phylum,rownames(phylum) %in% rownames(subset(info,info$cms_merge==unique(info$cms_merge)[i]))))
}

#绘制堆叠图，取5种CMS分型的平均phylum表达累计前十的菌种

phylum_sum<-colSums(phylum_cms)
phylum_top10<-names(sort(colSums(phylum_cms),decreasing=T))[1:10]
phylum_other<-rowSums(phylum_cms[,setdiff(colnames(phylum),phylum_top10)])
phylum_new<-data.frame(phylum_cms[,phylum_top10],Other=phylum_other)
phylum_pro<-t(apply(phylum_new,1,function(x){return(x/sum(x))}))
rownames(phylum_pro)<-rownames(phylum_new)
colnames(phylum_pro)<-colnames(phylum_new)
write.csv(phylum_pro,'/data3/data/JN/CRC/TCGA_microbes/phylum_vs_rf_relative_top10.csv',row.names=T)

library(reshape2)
phylum_pro<-data.frame('CMS'=rownames(phylum_pro),phylum_pro)
phylum_pro<-phylum_pro %>%
  melt(id.vars='CMS',variable.name='phylum',value.name='proportion')

#画堆叠图
library(ggplot2)
library(RColorBrewer)
p<-ggplot(phylum_pro,aes(x=CMS,y=100*proportion,fill=phylum))+
  geom_col(position='stack',width=0.6)+
  #scale_y_continuous(expand=c(0,0))+
  scale_fill_manual(values= brewer.pal(n=11,name='Set3'))+
  labs(x='',y='Relative Abundance(%)')+
  theme(panel.grid = element_blank(), panel.background = element_rect(color = 'black', fill = 'transparent'), strip.text = element_text(size = 12)) + #设置主题背景，根据自己的需求定制
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.title = element_blank(), legend.text=element_text(size=10))

ggsave(p, file='proportion_top10_phylum_vs_rf_GA.pdf', width=8, height=6)
phylum_cms<-t(apply(phylum_cms,1,function(x){return(x/sum(x))}))
write.csv(phylum_cms[sort(rownames(phylum_cms)),],'CMS_relative_abundances_vs_rf_phylum.csv',row.names=T)




###############################
#
#差异菌分析(wilcox)
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
microbes_raw<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)
#microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)

#计算相对丰度(暂不用)
#microbes<-t(apply(microbes_raw,1,function(x){return(x/sum(x))}))
#rownames(microbes)<-rownames(microbes_raw)
#colnames(microbes)<-colnames(microbes_raw)

#去除prevalance<20%的菌种,剩余797种
microbes<-microbes[,which(unlist(apply(microbes_raw,2,function(x){sum(x>0)/length(x)}))>0.2)]

#CMS1：193;CMS2:113;CMS3:53;CMS4:209;NOLBL:33
#CMS1：293;CMS2:153;CMS3:66;CMS4:109;NOLBL:38
#单个分型vs其他
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'

info$cms_merge[which(info$cms_merge!='Mixed-CMS')]='OTHER'
wilcox_pvalue<-apply(t(microbes),1,function(x){wilcox.test(x~info[rownames(microbes),]$cms_merge,paired=FALSE)$p.value})
cmsmix_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge=='Mixed-CMS'))])
other_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge=='OTHER'))])
means_diff<-cmsmix_means-other_means
diff<-data.frame(cmsmix_means=cmsmix_means,other_means=other_means,means_diff=means_diff,pvalue=wilcox_pvalue)
rownames(diff)<-colnames(microbes)
diff<-subset(diff,diff$pvalue<0.05)
rcmsmix<-diff

write.csv(diff,'/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_Mixed-CMSvsOTHER.csv',quote=F)

#生成sparseCC需要的矩阵
rcms1<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS1vsOTHER.csv',row.names=1)
rcms2<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS2vsOTHER.csv',row.names=1)
rcms3<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS3vsOTHER.csv',row.names=1)
rcms4<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS4vsOTHER.csv',row.names=1)

library(readr)
microbes_raw<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes_raw<-microbes_raw[,which(unlist(apply(microbes_raw,2,function(x){sum(x>0)/length(x)}))>0.2)]
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'


data<-data.frame(otu_id=rownames(rcms4),
                 t(microbes_raw[rownames(subset(info,info$cms_merge=='CMS4')),rownames(rcms4)]),check.names = F)
colnames(data)[1]<-'#OTU ID'
dim(data)

write_tsv(data,'/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc_rf_vs_n/pre-sparcc_rf_vs_n_cms4.tsv')


#CMS1：14;CMS2:0;CMS3:0;CMS4:7;NOLBL:33
#ANOVA+后验
library(reshape2)
library(DescTools)
aov_pvalue<-c()
cms2_dun_pvalue<-c()
cms3_dun_pvalue<-c()
cms4_dun_pvalue<-c()
nolbl_dun_pvalue<-c()
cms2_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge=='CMS2'))])
other_means<-rowMeans(t(microbes)[,rownames(subset(info,info$cms_merge!='CMS2'))])
for(i in 1:ncol(microbes)){
  data=data.frame(values=microbes[,i],group=info[rownames(microbes),]$cms_merge)
  data$group<-as.factor(data$group)
  aov_pvalue[i] <- summary(aov(values~group, data=data))[[1]][1,5]
  dun<-DunnettTest(x=data$values, g=data$group, control ="CMS2")$CMS2
  cms2_dun_pvalue[i]<-dun[1,4]
  cms3_dun_pvalue[i]<-dun[2,4]
  cms4_dun_pvalue[i]<-dun[3,4]
  nolbl_dun_pvalue[i]<-dun[4,4]
  
}
diff<-data.frame(aov_pvalue=aov_pvalue,cms1_dun_pvalue=cms1_dun_pvalue,cms2_dun_pvalue=cms2_dun_pvalue,
                 cms3_dun_pvalue=cms3_dun_pvalue,nolbl_dun_pvalue=nolbl_dun_pvalue)
rownames(diff)<-colnames(microbes)
diff<-subset(diff,((diff$aov_pvalue<0.05)&((diff$cms1_dun_pvalue<0.05)&(diff$cms2_dun_pvalue<0.05))
                   |(diff$cms3_dun_pvalue<0.05)&(diff$nolbl_dun_pvalue<0.05)))
dim(diff)

diff<-subset(diff,((diff$aov_pvalue<0.05)&(diff$cms2_dun_pvalue<0.05)&(diff$cms3_dun_pvalue<0.05)&(diff$cms4_dun_pvalue<0.05)&(diff$nolbl_dun_pvalue<0.05)))
diff<-data.frame(cms1_means=cms1_means[rownames(diff)],other_means=other_means[rownames(diff)],diff)

write.csv(diff,'/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_aov_CMS4.csv',quote=F)


#venn图
library(ggvenn)
library(ggtext)
library(RColorBrewer)
library(ggVennDiagram)
rcms1<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS1vsOTHER.csv',row.names=1)
rcms2<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS2vsOTHER.csv',row.names=1)
rcms3<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS3vsOTHER.csv',row.names=1)
rcms4<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS4vsOTHER.csv',row.names=1)
rcmsno<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_Mixed-CMSvsOTHER.csv',row.names=1)


dat <-list(
  CMS1=rownames(rcms1),
  CMS2=rownames(rcms2),
  CMS3=rownames(rcms3),
  CMS4=rownames(rcms4),
  Mixed.CMS=rownames(rcmsno)
)
brewer.pal(n=5,name='Set2')
color1<- alpha( "#f8766d",0.9)
p<-ggVennDiagram(dat, label_alpha=0,label_size =4,edge_size= 0.5,label = "count",category.names = c("CMS1","CMS2","CMS3", "CMS4","Mixed-CMS")) +
  scale_color_brewer(palette = "Set2")+
  scale_fill_gradient(low= "white",high = alpha( "#f8766d",0.9),guide= "none")
  

ggsave(p, file='Venn_CMS_vs_rf.pdf', width=8, height=5)

dat <-list(
  CMS1=rownames(rcms1),
  CMS2=rownames(rcms2),
  CMS3=rownames(rcms3),
  CMS4=rownames(rcms4)
)
brewer.pal(n=4,name='Set2')

library(VennDiagram)
p<-ggvenn(dat,fill_color=brewer.pal(n=4,name='Set2'),show_percentage=F,stroke_color='white',fill_alpha = 0.5,set_name_size = 6,text_size = 6)

png(file = "Venn_CMS_vs_rf.png",width = 3000,height = 2000,bg='white',res=300)
p
dev.off()


library(readr)
#write.csv(data.frame(otu_id=cms4,t(microbes_raw[,cms4]),check.names = F),'/data3/data/JN/CRC/TCGA_microbes/CMS/pre-sparcc_cms4.csv',row.names=F,quote=F)
#write.table(data.frame(otu_id=cms4,t(microbes_raw[,cms4]),check.names = F),'/data3/data/JN/CRC/TCGA_microbes/CMS/pre-sparcc_cms4.txt',row.names=F,quote=F,sep='\t')
write_tsv(data,'/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/pre-sparcc_cms4.tsv')

#重新构架correlation pvalue关系
fastsapr<-read.table('/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/cms4_pvalues.tsv',row.names=1)
colnames(fastsapr)<-rownames(fastsapr)
microbe_pair<-c()
pvalue<-c()
n=0
for(i in 1:(nrow(fastsapr)-1)){
  for(j in 1:(nrow(fastsapr)-i)){
    if(fastsapr[i,j]<0.05){
      microbe_pair[n]<-paste(rownames(fastsapr)[i],colnames(fastsapr)[j],sep='_and_')
      pvalue[n]<-fastsapr[i,j]
      n=n+1
    }
  }
}

head(data.frame(microbe_pair,pvalue))
write.csv(data.frame(microbe_pair,pvalue),'/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/cms4_p0.05_100.csv',row.names=F)

pairs<-c()
for(i in 1:nrow(res2)){
  pairs[i]<-paste(res2[i,1],res2[i,2],sep='_and_')
}

#重新构建correlation matrix关系
#重新构架correlation pvalue关系
fastsapr_p<-read.table('/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/cms4_pvalues_test1.tsv',row.names=1)
fastsapr_cor<-read.table('/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/median_correlation.tsv',row.names=1)
colnames(fastsapr_p)<-rownames(fastsapr_p)
colnames(fastsapr_cor)<-rownames(fastsapr_cor)
library(reshape2)

fastsapr_p<-as.matrix(fastsapr_p)
fastsapr_p[upper.tri(fastsapr_p)]<-NA
diag(fastsapr_p)<-NA
fastsapr_p2<-reshape::melt(fastsapr_p,na.rm=T)
fastsapr_p2<-subset(fastsapr_p2,!is.na(fastsapr_p2$value))

fastsapr_cor<-as.matrix(fastsapr_cor)
fastsapr_cor[upper.tri(fastsapr_cor)]<-NA
diag(fastsapr_cor)<-NA
fastsapr_cor2<-reshape::melt(fastsapr_cor,na.rm=T)
fastsapr_cor2<-subset(fastsapr_cor2,!is.na(fastsapr_cor2$value))

cor_table<-data.frame(fastsapr_cor2,fastsapr_p2$value)
colnames(cor_table)<-c('microbe1','microbe2','cor','pvalue')
cor_table$microbe1<-as.character(sapply(as.character(cor_table$microbe1),function(x){strsplit(x,'\\.')[[1]][length(strsplit(x,'\\.')[[1]])]}))
cor_table$microbe2<-as.character(sapply(as.character(cor_table$microbe2),function(x){strsplit(x,'\\.')[[1]][length(strsplit(x,'\\.')[[1]])]}))

test1<-apply(as.matrix(subset(cor_table,(abs(cor_table$cor)>0.2)&(cor_table$pvalue<0.05))),1,function(x){paste(x,collapse=';')})
data<-subset(cor_table,(abs(cor_table$cor)>0.5)&(cor_table$pvalue<0.05))
write.table(data,'/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc/cms4_cor_graph.txt',row.names=F,quote=F,sep='\t')


#门水平上差异
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='NOLBL'
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/phylum_raw_relative.csv',row.names=1)

#wilcox检验
#CMS4 vs CMS1:14
#CMS4 vs CMS2:11
#CMS4 vs CMS3:7
#CMS4 vs NOLBL:1,Fibrobacteres

cms<-subset(info,info$cms_merge=='CMS4'|info$cms_merge=='NOLBL')
wilcox_pvalue<-apply(t(microbes[rownames(cms),]),1,function(x){wilcox.test(x~cms$cms_merge,paired=FALSE)$p.value})
cms4_means<-rowMeans(t(microbes)[,rownames(subset(cms,cms$cms_merge=='CMS4'))])
nolbl_means<-rowMeans(t(microbes)[,rownames(subset(cms,cms$cms_merge=='NOLBL'))])
FC<-cms4_means/nolbl_means
diff<-data.frame(cms4_means=cms4_means,nolbl_means=nolbl_means,FC=FC,pvalue=wilcox_pvalue)
rownames(diff)<-colnames(microbes)
diff<-subset(diff,diff$pvalue<0.05)

write.csv(diff,'/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_r_CMS4vsNOLBL.csv',quote=F)

dat <-list(
  CMS4vsCMS1=rownames(rcms41),
  CMS4vsCMS2=rownames(rcms42),
  CMS4vsCMS3=rownames(rcms43),
  CMS4vsNOLBL=rownames(rcms4no)
)

p<-ggvenn(dat,show_percentage = T,show_elements = F,label_sep = ",",
          digits = 1,stroke_color = "white",
          fill_color = c("#E41A1C", "#1E90FF", "#FF8C00",
                         "#4DAF4A"))

ggsave(p, file='venn_diff_r_cms4_phylm_raw.pdf', width=10, height=6)

###############################
#
#co-occurance的菌种对
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS/sparcc_rf_vs_n')
cms1<-read.csv('cms1_cor_graph_0.7.csv')



###############################
#
#差异菌分析(limma)
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
microbes_raw<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_raw_barcode.csv',row.names=1)
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)

#去除prevalance<20%的菌种,剩余797种
microbes<-microbes[,which(unlist(apply(microbes_raw,2,function(x){sum(x>0)/length(x)}))>0.2)]

#CMS1：260;CMS2:129;CMS3:68;CMS4:101;NOLBL:31
#单个分型vs其他
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'

info$cms_merge[which(info$cms_merge!='CMS4')]='Other'
library(limma)
condition<-factor(info$cms_merge,levels=c("Other","CMS4"))
design <- model.matrix(~condition)

fit <- lmFit(t(microbes), design)
fit <- eBayes(fit)
diff <- topTable(fit, coef="conditionCMS4",n = Inf, sort = "p")
diff<-subset(diff,diff$P.Value<0.05)

write.csv(diff,'/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_limma_CMS4vsOTHER.csv',quote=F)


###############################
#
#差异菌分析热图展示(py38)
#
###############################
setwd('/data3/data/JN/CRC/TCGA_microbes/CMS')
#取微生物数据
rcms1<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS1vsOTHER.csv',row.names=1)
rcms2<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS2vsOTHER.csv',row.names=1)
rcms3<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS3vsOTHER.csv',row.names=1)
rcms4<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS4vsOTHER.csv',row.names=1)
rcmsno<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_Mixed-CMSvsOTHER.csv',row.names=1)
microbes<-read.csv('/data3/data/JN/CRC/TCGA_microbes/TCGA_RNA_cri3_barcode.csv',row.names=1)


#归一化maxmin转化
#maxmin <- function(x) (x - min(x))/(max(x)-min(x))
#microbes2<-apply(microbes, 2, maxmin)
#rna<-apply(rna, 2, maxmin)

#秩转换+maxmin转化
microbes<-apply(microbes,2,rank)
maxmin <- function(x) (x - min(x))/(max(x)-min(x))
microbes<-apply(microbes, 2, maxmin)

#Z-score
#microbes2<-apply(microbes,2,function(x){scale(x)})
#rownames(microbes2)<-rownames(microbes)
#microbes<-microbes2

#ocms1<-setdiff(rownames(rcms1),c(rownames(rcms2),rownames(rcms3),rownames(rcms4),rownames(rcmsno)))
#ocms2<-setdiff(rownames(rcms2),c(rownames(rcms1),rownames(rcms3),rownames(rcms4),rownames(rcmsno)))
#ocms3<-setdiff(rownames(rcms3),c(rownames(rcms1),rownames(rcms2),rownames(rcms4),rownames(rcmsno)))
#ocms4<-setdiff(rownames(rcms4),c(rownames(rcms1),rownames(rcms2),rownames(rcms3),rownames(rcmsno)))
#ocmsno<-setdiff(rownames(rcmsno),c(rownames(rcms1),rownames(rcms2),rownames(rcms3),rownames(rcms4)))

#取genus
ocms1<-rownames(subset(rcms1,abs(rcms1$means_diff)>0.5))
ocms2<-rownames(subset(rcms2,abs(rcms2$means_diff)>0.5))
ocms3<-rownames(subset(rcms3,abs(rcms3$means_diff)>0.5))
ocms4<-rownames(subset(rcms4,abs(rcms4$means_diff)>0.5))
ocmsno<-rownames(subset(rcmsno,abs(rcmsno$means_diff)>0.5))

ocms1_up<-rownames(subset(rcms1,rcms1$means_diff>0.5))
ocms1_down<-rownames(subset(rcms1,rcms1$means_diff<(-0.5)))
ocms2_up<-rownames(subset(rcms2,rcms2$means_diff>0.5))
ocms2_down<-rownames(subset(rcms2,rcms2$means_diff<(-0.5)))
ocms3_up<-rownames(subset(rcms3,rcms3$means_diff>0.5))
ocms3_down<-rownames(subset(rcms3,rcms3$means_diff<(-0.5)))
ocms4_up<-rownames(subset(rcms4,rcms4$means_diff>0.5))
ocms4_down<-rownames(subset(rcms4,rcms4$means_diff<(-0.5)))
#ocmsno_up<-rownames(subset(rcmsno,rcmsno$means_diff>0.5))
ocmsno_down<-rownames(subset(rcmsno,rcmsno$means_diff<(-0.5)))


#取rna
#ocms1<-rownames(rcms1[order(rcms1$pvalue),])[1:100]
#ocms2<-rownames(rcms2[order(rcms2$pvalue),])[1:100]
#ocms3<-rownames(rcms3[order(rcms3$pvalue),])[1:100]
#ocms4<-rownames(rcms4[order(rcms4$pvalue),])[1:100]
#ocmsno<-rownames(rcmsno[order(rcmsno$pvalue),])[1:100]

ocms1<-rownames(subset(rcms1,rcms1$means_FC<0.25|rcms1$means_FC>4))
ocms2<-rownames(subset(rcms2,rcms2$means_FC<0.25|rcms2$means_FC>4))
ocms3<-rownames(subset(rcms3,rcms3$means_FC<0.25|rcms3$means_FC>4))
ocms4<-rownames(subset(rcms4,rcms4$means_FC<0.25|rcms4$means_FC>4))
ocmsno<-rownames(subset(rcmsno,rcmsno$means_FC<0.25|rcmsno$means_FC>4))


#读取CMS分型注释
info<-read.csv('/data3/data/JN/CRC/TCGA_microbes/sample_info_no_multibam_tumor_rf.csv',row.names=1)
info$cms_merge[which(is.na(info$cms_merge))]='Mixed-CMS'
info$cms_merge[which(info$cms_merge=='NOLBL')]='Mixed-CMS'
#读取meta注释
annot<-read.csv('/data3/data/JN/CRC/TCGA_microbes/meta_RNA_with_submitter_id_barcode.csv',row.names=1,check.names = 1)


#Genus
#按CMS分布
microbes<-rbind(microbes[rownames(subset(info,info$cms_merge=='CMS1')),],
                microbes[rownames(subset(info,info$cms_merge=='CMS2')),],
                microbes[rownames(subset(info,info$cms_merge=='CMS3')),],
                microbes[rownames(subset(info,info$cms_merge=='CMS4')),],
                microbes[rownames(subset(info,info$cms_merge=='Mixed-CMS')),])
microbes<-t(microbes)


#取各类CMS的特征菌种合并(594样本，601重复genus)
microbes_cms<-rbind(microbes[ocms1,],
                    microbes[ocms2,],
                    microbes[ocms3,],
                    microbes[ocms4,],
                    microbes[ocmsno,])

microbes_cms<-rbind(microbes[c(ocms1_up,ocms1_down),],
                    microbes[c(ocms2_up,ocms2_down),],
                    microbes[c(ocms3_up,ocms3_down),],
                    microbes[c(ocms4_up,ocms4_down),],
                    microbes[c(ocmsno_down),])

#差异基因
#按CMS分布
rna<-rbind(rna[rownames(subset(info,info$cms_merge=='CMS1')),],
           rna[rownames(subset(info,info$cms_merge=='CMS2')),],
           rna[rownames(subset(info,info$cms_merge=='CMS3')),],
           rna[rownames(subset(info,info$cms_merge=='CMS4')),],
           rna[rownames(subset(info,info$cms_merge=='Mixed-CMS')),])
rna<-t(rna)


#取各类CMS的特征菌种合并(594样本，601重复genus)
rna_cms<-rbind(rna[ocms1,],
               rna[ocms2,],
               rna[ocms3,],
               rna[ocms4,],
               rna[ocmsno,])



library(RColorBrewer)
library(circlize)
library(ComplexHeatmap)
#样本的CMS注释

cms_sample_split<-c(rep('CMS1',table(info$cms_merge)['CMS1']),rep('CMS2',table(info$cms_merge)['CMS2']),
                    rep('CMS3',table(info$cms_merge)['CMS3']),rep('CMS4',table(info$cms_merge)['CMS4']),
                    rep('Mixed.CMS',table(info$cms_merge)['Mixed-CMS']))
#cms_taxa_split<-c(rep('CMS1',length(ocms1)),rep('CMS2',length(ocms2)),rep('CMS3',length(ocms3)),
#                  rep('CMS4',length(ocms4)),rep('Mixed.CMS',length(ocmsno)))

cms_taxa_split<-c(rep('CMS1_up',length(ocms1_up)),rep('CMS1_down',length(ocms1_down)),rep('CMS2_up',length(ocms2_up)),rep('CMS2_down',length(ocms2_down)),
                  rep('CMS3_up',length(ocms3_up)),rep('CMS3_down',length(ocms3_down)),rep('CMS4_up',length(ocms4_up)),rep('CMS4_down',length(ocms4_down)),
                  rep('Mixed.CMS_dowm',length(ocmsno_down)))


#所有注释按microbe_cms列名排列!
info<-info[colnames(microbes_cms),]
annot<-annot[colnames(microbes_cms),]

#所有注释按rna_cms列名排列!
info<-info[colnames(rna_cms),]
annot<-annot[colnames(rna_cms),]


#增加CMS注释
cms_col<-brewer.pal(n=5,name='Set2')[1:5]
names(cms_col)<-names(table(info$cms_merge))

#增加性别注释
annot$gender<-unlist(lapply(annot$gender,tolower))
gender_col<-brewer.pal(n=5,name='Set1')[1:2]
names(gender_col)<-names(table(annot$gender))

#增加年龄注释
annot$age[which(annot$age<=50)]<-"<50"
annot$age[which((annot$age>50)&(annot$age<=70))]<-"50-70"
annot$age[which(annot$age>70)]<-">70"

age_col<-brewer.pal(n=5,name='Set1')[1:3]
names(age_col)<-c("<50","50-70",">70")

#增加样本点位置注释
clin_coad<-read.csv('/data3/data/JN/CRC/TCGAbiolinks/TCGAbiolinks_coad_clin.csv',row.names=1,check.names=F)
clin_read<-read.csv('/data3/data/JN/CRC/TCGAbiolinks/TCGAbiolinks_read_clin.csv',row.names=1,check.names=F)

clin_sub<-rbind(clin_coad[,c("ajcc_pathologic_m","tissue_or_organ_of_origin")],clin_read[,c("ajcc_pathologic_m","tissue_or_organ_of_origin")])
clin_sub['submitter']<-rownames(clin_sub)
clin_sub<-merge(clin_sub,data.frame(sample=rownames(annot),submitter=unlist(lapply(rownames(annot),function(x){substr(x,1,12)}))),all.Y=TRUE)
rownames(clin_sub)<-clin_sub$sample
#所有注释按microbe_cms列名排列!
clin_sub<-clin_sub[colnames(microbes_cms),]

clin_sub$tissue_or_organ_of_origin[which(clin_sub$tissue_or_organ_of_origin=='Colon, NOS')]<-"Colon"
clin_sub$tissue_or_organ_of_origin[which(clin_sub$tissue_or_organ_of_origin=='Connective, subcutaneous and other soft tissues of abdomen')]<-"Abdomen"
clin_sub$tissue_or_organ_of_origin[which(clin_sub$tissue_or_organ_of_origin=='Hepatic flexure of colon')]<-"Hepatic flexure"
clin_sub$tissue_or_organ_of_origin[which(clin_sub$tissue_or_organ_of_origin=='Rectum, NOS')]<-"Rectum"
clin_sub$tissue_or_organ_of_origin[which(clin_sub$tissue_or_organ_of_origin=='Splenic flexure of colon')]<-"Splenic flexure"

tissue_or_organ_of_origin_col<-c(brewer.pal(n=8,name='Set1')[1:8],brewer.pal(n=5,name='Set3')[c(1,3,5)],'grey')
names(tissue_or_organ_of_origin_col)<-names(table(clin_sub$tissue_or_organ_of_origin))

#增加T期注释
annot$pathologic_t_label[which(annot$pathologic_t_label=='T4a'|annot$pathologic_t_label=='T4b')]<-"T4"
pathologic_t_label_col<-c('grey',brewer.pal(n=5,name='Set1')[1:5])
names(pathologic_t_label_col)<-names(table(annot$pathologic_t_label))

#增加N期注释
annot$pathologic_n_label[which(annot$pathologic_n_label=='N1a'|annot$pathologic_n_label=='N1b'|annot$pathologic_n_label=='N1c')]<-"N1"
annot$pathologic_n_label[which(annot$pathologic_n_label=='N2a'|annot$pathologic_n_label=='N2b')]<-"N2"

pathologic_n_label_col<-c(brewer.pal(n=5,name='Set1')[1:3],'grey',brewer.pal(n=5,name='Set1')[4])
names(pathologic_n_label_col)<-names(table(annot$pathologic_n_label))

#增加M期注释
clin_sub$ajcc_pathologic_m[which(clin_sub$ajcc_pathologic_m=='M1a'|clin_sub$ajcc_pathologic_m=='M1b')]<-"M1"

ajcc_pathologic_m_col<-c(brewer.pal(n=5,name='Set1')[1:3])
names(ajcc_pathologic_m_col)<-names(table(clin_sub$ajcc_pathologic_m))

#增加stage注释
annot$pathologic_stage_label[(annot$pathologic_stage_label=='Stage IA')|(annot$pathologic_stage_label=='Stage I')]='I'
annot$pathologic_stage_label[(annot$pathologic_stage_label=='Stage IIA')|(annot$pathologic_stage_label=='Stage IIB')|(annot$pathologic_stage_label=='Stage IIC')|(annot$pathologic_stage_label=='Stage II')]='II'
annot$pathologic_stage_label[(annot$pathologic_stage_label=='Stage IIIA')|(annot$pathologic_stage_label=='Stage IIIB')|(annot$pathologic_stage_label=='Stage IIIC')|(annot$pathologic_stage_label=='Stage III')]='III'
annot$pathologic_stage_label[(annot$pathologic_stage_label=='Stage IVA')|(annot$pathologic_stage_label=='Stage IVB')|(annot$pathologic_stage_label=='Stage IV')]='IV'

pathologic_stage_label_col<-c(brewer.pal(n=5,name='Set1')[1:4],'grey')
names(pathologic_stage_label_col)<-names(table(annot$pathologic_stage_label))


#创建行注释

ha_row<-rowAnnotation(foo=c(rep('CMS1',length(c(ocms1_up,ocms1_down))),rep('CMS2',length(c(ocms2_up,ocms2_down))),rep('CMS3',length(c(ocms3_up,ocms3_down))),
                            rep('CMS4',length(c(ocms4_up,ocms4_down))),rep('Mixed-CMS',length(ocmsno_down))),
                      col=list(foo=cms_col),
                      show_legend=FALSE,
                      show_annotation_name=FALSE)

ha_row<-rowAnnotation(foo=c(rep('CMS1',length(ocms1)),rep('CMS2',length(ocms2)),rep('CMS3',length(ocms3)),rep('CMS4',length(ocms4)),rep('Mixed-CMS',length(ocmsno))),
                      col=list(foo=cms_col),
                      show_legend=FALSE,
                      show_annotation_name=FALSE)

#创建列注释
df<-data.frame(CMS=info$cms_merge,
               Gender=annot$gender,
               Age=annot$age,
               tissue_or_organ_of_origin=clin_sub$tissue_or_organ_of_origin,
               pathologic_t_label=annot$pathologic_t_label,
               pathologic_n_label=annot$pathologic_n_label,
               ajcc_pathologic_m=clin_sub$ajcc_pathologic_m,
               pathologic_stage_label=annot$pathologic_stage_label
               )
ha_column<-HeatmapAnnotation(df=df,
                             col=list(CMS=cms_col,
                                      Gender=gender_col,
                                      Age=age_col,
                                      tissue_or_organ_of_origin=tissue_or_organ_of_origin_col,
                                      pathologic_t_label=pathologic_t_label_col,
                                      pathologic_n_label=pathologic_n_label_col,
                                      ajcc_pathologic_m=ajcc_pathologic_m_col,
                                      pathologic_stage_label=pathologic_stage_label_col),
                             show_legend=FALSE,
                             annotation_label=c('CMS','Gender','Age','Tissue Site','AJCC Tumor Category','AJCC LN Category',
                                                'AJCC Metastasis Category','Tumor Stage'))
lgd_list = list(
  Legend(labels = names(cms_col), title = "CMS", type = "grid",
         legend_gp = gpar(fill=cms_col)),
  Legend(labels = names(gender_col), title = "Gender", type = "grid",
         legend_gp = gpar(fill=gender_col)),
  Legend(labels = names(age_col), title = "Age", type = "grid",
         legend_gp = gpar(fill=age_col)),
  Legend(labels = names(pathologic_t_label_col[c(6,2,3,4,5)]), title = "AJCC Tumor Category", type = "grid",
         legend_gp = gpar(fill=pathologic_t_label_col[c(6,2,3,4,5)])),
  Legend(labels = names(pathologic_n_label_col[c(1,2,3,5)]), title = "AJCC LN Category", type = "grid",
         legend_gp = gpar(fill=pathologic_n_label_col[c(1,2,3,5)])),
  Legend(labels = names(ajcc_pathologic_m_col), title = "AJCC Metastasis Category", type = "grid",
         legend_gp = gpar(fill=ajcc_pathologic_m_col)),
  Legend(labels = names(pathologic_stage_label_col), title = "Tumor Stage", type = "grid",
         legend_gp = gpar(fill=pathologic_stage_label_col)),
  Legend(labels = names(tissue_or_organ_of_origin_col), title = "Tissue Site", type = "grid",
         legend_gp = gpar(fill=tissue_or_organ_of_origin_col))
         )
#genus画图
p<-Heatmap(
  name="mat",
  microbes_cms, 
  # 添加上方注释
  top_annotation = ha_column,
  left_annotation =ha_row,
  column_split=cms_sample_split,
  row_split=ordered(cms_taxa_split,levels = c("CMS1_up", "CMS1_down", "CMS2_up","CMS2_down",
                                      "CMS3_up", "CMS3_down", "CMS4_up", "CMS4_down",
                                      "Mixed.CMS_dowm")),
  #clustering_distance_rows="kendall",
  #clustering_distance_columns="kendall",
  show_row_names=F,
  show_column_names=F,
  column_title = "Sample",
  row_title='Genus',
  heatmap_legend_param = list(
    #at = c(0:2),
    #labels = c("low", "high"),
    title = "Microbial Abundance(Normalized)",
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot"
  )
)

#rna画图
p<-Heatmap(
  name="mat",
  rna_cms, 
  # 添加上方注释
  top_annotation = ha_column,
  left_annotation =ha_row,
  column_split=cms_sample_split,
  row_split=cms_taxa_split,
  show_row_names=F,
  show_column_names=F,
  column_title = "Sample",
  row_title='Genus',
  heatmap_legend_param = list(
    #at = c(0:2),
    #labels = c("low", "high"),
    title = "Microbial Abundance(Normalized)",
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot"
  )
)


#rp<-c(row_order(p)$CMS1,row_order(p)$CMS2,row_order(p)$CMS3,row_order(p)$CMS4,row_order(p)$Mixed.CMS)
#46+24+3+12
rp<-c(row_order(p)$CMS1_up,row_order(p)$CMS1_down,row_order(p)$CMS2_up,row_order(p)$CMS2_down,
      row_order(p)$CMS3_up,row_order(p)$CMS3_down,row_order(p)$CMS4_up,row_order(p)$CMS4_down,
      row_order(p)$Mixed.CMS_dowm)

cp<-c(column_order(p)$CMS1,column_order(p)$CMS2,column_order(p)$CMS3,column_order(p)$CMS4,column_order(p)$Mixed.CMS)

#microbes_cms[microbes_cms>=3]=3
#microbes_cms[microbes_cms<=-3]=-3
#genus画图
p2<-Heatmap(
  name="mat",
  microbes_cms, 
  col<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  # 添加上方注释
  top_annotation = ha_column,
  left_annotation =ha_row,
  row_order=rp,
  column_order=cp,
  column_split=cms_sample_split,
  row_split=ordered(cms_taxa_split,levels = c("CMS1_up", "CMS1_down", "CMS2_up","CMS2_down",
                                              "CMS3_up", "CMS3_down", "CMS4_up", "CMS4_down",
                                              "Mixed.CMS_dowm")),
  show_row_names=F,
  show_column_names=F,
  cluster_columns=F,
  cluster_rows=F,
  column_title = "Sample",
  row_title='Genus',
  heatmap_legend_param = list(
    #at = c(0:2),
    #labels = c("low", "high"),
    title = "Microbial Abundance(Normalized)",
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot"
  )
)

#RNA画图
rna_cms[rna_cms>=5]=5
rna_cms[rna_cms<=-5]=-5
p2<-Heatmap(
  name="mat",
  rna_cms, 
  col<-colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
  # 添加上方注释
  top_annotation = ha_column,
  left_annotation =ha_row,
  row_order=rp,
  column_order=cp,
  #column_split=cms_sample_split,
  #row_split=cms_taxa_split,
  show_row_names=F,
  show_column_names=F,
  cluster_columns=F,
  cluster_rows=F,
  column_title = "Sample",
  row_title='Genus',
  heatmap_legend_param = list(
    #at = c(0:2),
    #labels = c("low", "high"),
    title = "Microbial Abundance(Normalized)",
    legend_height = unit(4, "cm"),
    title_position = "leftcenter-rot"
  )
)


#colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)

pdf('CMS_heatmap_rf_vs3.pdf',width = 12,height = 13,onefile = F)
draw(p2,ht_gap = unit(7, "mm"), row_km = 1, annotation_legend_list = lgd_list,align_heatmap_legend='heatmap_center')
dev.off() 

###############################
#
#差异菌分析所属phylum分析
#
###############################
cms1<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS1vsOTHER.csv',row.names=1)
cms2<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS2vsOTHER.csv',row.names=1)
cms3<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS3vsOTHER.csv',row.names=1)
cms4<-read.csv('/data3/data/JN/CRC/TCGA_microbes/CMS/microbes_diff_vs_rf_CMS4vsOTHER.csv',row.names=1)
#病毒合并，共1153种
library(stringr)
microbes_name<-rownames(cms4)

names<-lapply(microbes_name,function(x){strsplit(x,split='\\.')})
names<-lapply(names,function(x){if(length(str_subset(x[[1]],"p__"))==0){
    return(str_subset(x[[1]],"k__"))}
  else{
    return(str_subset(x[[1]],"p__"))}})
names<-unlist(lapply(names,function(x){if(x[[1]]=='k__Viruses'){
  return('Virus')}
  else{
  return(strsplit(x[[1]],split='p__')[[1]][2])}}))
cms4_phylum<-names
sort(table(cms4_phylum),decreasing=T)

#CMS1:293 Proteobacteria:101(34.47%) Firmicutes:58(19.80%)
#CMS2:153 Proteobacteria:44(28.76%) Firmicutes:35(22.88%)
#CMS3:66 Proteobacteria:17(25.76%) Firmicutes:10(15.15%)
#CMS4:109 Proteobacteria:54(49.54%) Firmicutes:11(10.09%)

names<-lapply(rownames(cms4),function(x){strsplit(x,split='\\.')})
names<-lapply(names,function(x){str_subset(x[[1]],"g__")})
names<-unlist(lapply(names,function(x){strsplit(x[[1]],split='g__')[[1]][2]}))
cms4_genus<-names
