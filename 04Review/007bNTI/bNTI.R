library(picante)
library(ape)
comm<-read.table("out_table.txt",sep='\t',row.names=1,header=T)#读入OUT文件
phy<-read.tree("rep_set.tre")#读入树文件
prune_tree<-prune.sample(comm,phy)#使树文件和OTU表文件对齐
phydist<-cophenetic(prune_tree)#计算每个OTU之间距离矩阵
mntd<-ses.mntd(comm,phydist,null.model="taxa.labels",abundance.weighted=T, runs=999)#计算MNTD

# 结果列名中、英文解释
# ntaxa Number of taxa in community，群落中OTU的数目；
# mntd.obs Observed MNTD in community，观测到的群落MNTD值；
# mntd.rand.mean Mean MNTD in null communities，随机群落的MNTD均值；
# mntd.rand.sd Standard deviation of MNTD in null communities，随机群落MNTD值的标准差；
# mntd.obs.rank Rank of observed MNTD vs. null communities，观察到的群落的MNTD值与随机群落MNTD值的排序；
# mntd.obs.z Standardized effect size of MNTD vs. null communities (=(mntd.obs-mntd.rand.mean)/mntd.rand.sd, equivalent to -NTI) 观察到的群落的MNTD与随机群落的MNTD值的标准化数值大小；
# mntd.obs.p P-value (quantile) of observed MNTD vs. null communities (= mntd.obs.rank/runs + 1) 观察到的群落的MNTD与随机群落的MNTD值比较的p值；
# runs Number of randomizations，随机化的次数。

comdist.resultS12<-comdistnt(comm,phydist)#计算βMNTD值

#以下为βMNTDnull的计算过程
f<-function(){
    g<-randomizeMatrix(comm,null.model=c("frequency","richness","independentswap","trialswap"),iterations=1000)#针对两组样方实现随机群落构建
    fc<-comdist.result<-comdistnt(g,phydist)#针对随机构建的两个样方进行群落间βMNTD的计算
}
mt<-mean(replicate(999, f()))#上述f计算过程重复进行999次，得到βMNTD的平均值
mt
mt.sd<-sd(replicate(999, f()))#计算上述999次计算得到的βMNTD值的标准差
mt.sd
βNTI=(comdist.resultS12-mt)/mt.sd#如果样品对较多，该步骤需在计算完每组样品对的βMNTD、mt和mt.sd值后，手动计算得出βNTI值
βNTI


