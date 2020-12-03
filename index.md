## 6.1 RNAeditor
1. 理解chr1.editingSites.vcf和chr1.editingSites.gvf的含义
* chr1.editingSites.vcf :对所有的编辑位点的一个粗略的展示。
> chrom 代表位于那条染色体上
> POS 代表variant位于染色体的哪一个位点，如果是INDEL则代表其第一个碱基所在的位置。
> REF 代表在参考基因组上在POS上的碱基是什么
> ALT 代表alter的碱基
> QUAL 是以phred格式的质量值，其值越高代表该位点存在variant的概率越大。
> FILTER 代表过滤结果是否可靠，PASS代表PASS，若不可靠，则该项不为"PASS"或“.”
> INFO 为variant 的详细信息
* chr1.editingSites.gvf ：Gene Variation File 对==每个==编辑位点的更信息的信息注释，包括geneID、name、segment等
![[chr1.editingSites.gvf.png]]
	**例如第一行，代表着，这个Variation 位于染色体1上的GeneID=ENSG00000225159，Name=NPM1P39，SEGMEMT=noncoding-exon,这段基因位于染色体的27206930-27207796，该Variation位于27206932，参考基因组上其为A，而Alter为G，其质量分数为66.28，共有7条reads，其中有6条发生了edited**
2. 根据chr1.editingSites.gvf文件，统计RNA编辑位点在基因组上的分布
```R
chr1.editingSites<-read_tsv("chr1.editingSites.gvf")
# 读取文件
chr1.editingSites%>%rename(Gene_ID=`#Gene_ID`)
# 对列名重命名
editingSites_summary<-chr1.editingSites%>%select(Gene_ID,Name,SEGMENT,Edited_Reads)
# 选取需要的4列
editingSites_summary<-editingSites_summary%>%group_by(Gene_ID,Name,SEGMENT)%>%summarise(number=sum(Edited_Reads))%>%ungroup()
# 将位于同一基因不同位点上的variant进行合并
editingSites_summary<-editingSites_summary%>%spread(SEGMENT,number)
# 将长数据转化为宽数据
editingSites_summary[is.na(editingSites_summary)]<-0
# 将NA值转化为0
```
![[chr1.editorSites.png]]
```R
plot_editingSites<-chr1.editingSites%>%select(SEGMENT,Edited_Reads)%>%group_by(SEGMENT)%>%summarise(number=sum(Edited_Reads))%>%ungroup()

plot_editingSites%>%ggplot(aes(x=SEGMENT,y=number))+geom_bar(stat="identity",aes(fill=SEGMENT))+
theme_bw()+
ggtitle("Summary of editorSites region")+
theme(axis.text.x=element_text(colour="black",family="Times",size=14),axis.text.y=element_text(family="Times",size=14,face="plain"), axis.title.y=element_text(family="Times",size = 14,face="plain"), axis.title.x=element_text(family="Times",size = 14,face="plain"),plot.title = element_text(family="Times",size=25,face="bold",hjust = 0.5))+
theme(legend.title = element_text(size = 20),legend.text = element_text(size = 15),axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text.x = element_text(size = 14,color="black"),axis.text.y = element_text(size = 14,color="black"))
```
![[条形图.png]]

## 6.2 APA(Alternative Polyadenylation)Detection
1. 产生Generate region annotation 
2. Identify APA events

### 对Alternative Polyadenylation (APA) background 的补充学习和理解
[Alternative Polyadenylation: Methods, Findings, and Impacts](https://www.sciencedirect.com/science/article/pii/S1672022917301341?via%3Dihub)
* 聚腺苷酸化发生在mRNA成熟的过程，对mRNA翻译效率、稳定性和定位有着重要的作用
*  真核生物基因拥有不止一个的PolyA位点，被称作alternative PolyA(APA),导致相同的基因产生不同的mRNA亚型
*  通过pA位点的定位，APA能够被分为两大类
	> CR-APA(coding region-APA)：pA sites 位于exons或者introns内部，影响翻译区并且导致拥有不同C 末端的蛋白质亚型。
	> UTR-APA：APA位于3‘ untranslated region，最终产生具有相同的coding region但是具有不同的3’ UTR。更长的3‘UTR能够有更多的microRNA的结合位点，更多RBP识别位点，从而影响RNA二级结构
### DaPars 工作原理
* DaPars  利用两个样本的摆动比对文件通过提取的远端Pa位点来推断近端的Pa位点
### homework
![[DaPars_results.png]]
> Gene: 代表基因名称及其信息
> predicted Proximal polyadenylation: 预测的近端的PolyA所在的位置
> loci: 基因组信息
> exp：A,B两个样本、长短两组的表达量
> Group_A_Mean、Group_B_Mean: PDUI作为Mean值
> PDUI_Group_diff:  等于Group_A_Mean-Group_B_Mean
1. 解释PDUI的含义
**代表一个远端的PolyA位点使用占比，数值在0-1，可以利用这个数值来评价APA事件发生的比例，如果PDUI接近于1代表这个基因更多的存在长的3’UTR，如果PDUI接近于0代表这个基因更多存在短的3‘UTR**
2. 写脚本过滤adjusted.P_val<=0.05,PDUI_Group_diff>=0.5, PDUI_fold_change>=0.59的作为diff-APA events，和Pass_filter为“Y“筛选出来的diff-APA events做比较。

```linux
# 按照adjusted.P_val<=0.05,PDUI_Group_diff>=0.5,PDUI_fold_change>=0.59 过滤数据
awk -F'\t' 'NR!=1{if($13>=0.5 && $15<=0.05 && $11/$12>0.59) print $0}' DaPars_Test_data_All_Prediction_Results.txt
```
![[filter_diff_APA.png]]

```linux
# 按照Pass_filter=="Y"筛选的diff-APA events
awk -F'\t' 'NR!=1{if($16=="Y") print $0}' DaPars_Test_data_All_Prediction_Results.txt
```
![[PASS_filter.png]]

两者筛选出来APA_diff一致
3.  思考
**我自己理解的是A,B两个样本之间通过PDUI的差值，来说明这个基因上的不同PolyA位点之间的距离差距。通过FDR、PDUI、Fold_change来判断这样距离的差距是否是明显的。假如A,B两个样本的PDUI值非常相近，则可能是同一个PolyA位点，而非diff-APA。**
**因此我觉得，假如软件中对A,B没有control和test之间的区别的话，应该以PDUI_Group_diff的绝对值>=某个值作为筛选的标准，因为无论是A>B,B>A应该是以A,B之间的距离来判断两个PolyA之间的距离**


## 6.3.Ribo-seq
### Backgroud
Ribo-seq是==细胞内蛋白翻译图谱的新型二代测序技术==，用来描述全基因组水平蛋白质的翻译情况。主要是==选择性捕捉80S核糖体及其结合的RNA片段而定位核糖体所位于的RNA位置。==
>  在细胞裂解物中富集多聚核糖体(polysome)
>  将多聚体核糖体用核酸酶(RNA nuclease)消化为单核糖体(monosome)
>  选择性的收集和富集80S核糖体并经纯化得到80S核糖体所保护的RNA片段。
>  在此过程中，将80S核糖体保护的RNA片段进行下一步构建文库和测序
