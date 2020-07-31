## 软件安装和数据库部署
	
	# 软件管理器miniconda2
	# 下载和安装
    wget -c https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    bash Miniconda2-latest-Linux-x86_64.sh
	# 查看并记录核心软件版本
    conda -V # 4.8.3
    python --version # 2.7.15

	# 配置conda环境
    # 添加生物频道，才能找到生物学软件
    conda config --add channels bioconda
    conda config --add channels conda-forge

    # 创建虚拟环境，防止污染环境变量
    conda create -n qc2 python=2.7
	# 启动虚拟环境
	conda activate qc2

	# 安装软件和查看版本
    conda install fastqc -y
    fastqc -v # FastQC v0.11.9
    conda install multiqc -y
    multiqc --version # multiqc, version 1.7
    conda install kneaddata -y
    kneaddata --version # 0.6.1
	conda install parallel -y
    parallel --version # 20200522

	# 查看可用的索引，有人类、小鼠基因组、人类转录组和核糖体数据库
	kneaddata_database
	# 下载索引的人类基因组Bowtie 2索引，多文件推荐建立并指定文件夹
    mkdir -p db
    kneaddata_database --download human_genome bowtie2 db/

	# 准备输入数据
	# (可选) 下载公共数据
	# 多文件中新建文件夹保存，-p允许建立多级目录、多个目录且目录已存在也不报错
	mkdir -p seq
	# wget下载单个样本，-c为支持断点续传，-O指定保存位置并可重命名，每个双端样本下载两个文件
	wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR117732/CRR117732_f1.fq.gz \
	  -O seq/C2_1.fq.gz
	wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR117732/CRR117732_r2.fq.gz \
	  -O seq/C2_2.fq.gz
	# 结合for循环再下载3个样本，seq命令产生连续序列，$i替换命令中可变部分，结尾加\保证变量名正确识别
	for i in `seq 3 5`;do
	  wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR11773$i/CRR11773$i\_f1.fq.gz \
		-O seq/C$i\_1.fq.gz
	  wget -c ftp://download.big.ac.cn/gsa/CRA002355/CRR11773$i/CRR11773$i\_r2.fq.gz \
		-O seq/C$i\_2.fq.gz
	done

# 实验步骤

	# 加载软件环境
	conda activate qc2

## 1.	FastQC测序数据质量评估

	# 用*.fq.gz代表所以有.fq.gz结尾的文件，即当前所有测序数据；-t 3指定最多同时使用3个线程
	fastqc seq/*.fq.gz -t 3

## 2.	MultiQC汇总质量评估报告

	# -d指定输入目录，-o指定输出目录
    multiqc -d seq/ -o ./

	# (可选)确定接头文件
	# 确定trimmomatic软件位置
	type trimmomatic
	# 根据上面环境路径+share/trimmomatic/adapters目录匹配接头序列的文件
	grep 'ATCGGAAGAGCACACGTCTGAAC' ~/.conda/envs/qc2/share/trimmomatic/adapters/*

## 3. (可选)检查测序双端序列标签是否唯一

	# 查看样本是否标签有重复
	zcat seq/C2_1.fq.gz|head
	zcat seq/C2_2.fq.gz|head

	# 解压、对序列的左、右端标题行分别添\1、\2
	gunzip seq/*.gz
    sed -i '1~4 s/$/\\1/g' seq/*_1.fq
    sed -i '1~4 s/$/\\2/g' seq/*_2.fq
	# 核对样本是否标签有重复
	head seq/C2_1.fq
	head seq/C2_2.fq
	# 压缩节省空间
    gzip seq/*.fq

## 4. KneadData质量控制和去宿主

    kneaddata -h # 显示帮助
	# 单位个样本质控和去宿主
    kneaddata -i seq/C2_1.fq.gz -i seq/C2_2.fq.gz \
      -o qc/ -v -t 8 --remove-intermediate-output \
      --trimmomatic ~/.conda/envs/qc2/share/trimmomatic \
      --trimmomatic-options 'ILLUMINACLIP:~/.conda/envs/qc2/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' -db db/Homo_sapiens

	# parallel管理队列，允许2个任务同时进行
    parallel -j 2 --xapply \
      "kneaddata -i seq/C{1}_1.fq.gz \
      -i seq/C{1}_2.fq.gz \
      -o qc/ -v -t 8 --remove-intermediate-output \
      --trimmomatic ~/.conda/envs/qc2/share/trimmomatic \
      --trimmomatic-options 'ILLUMINACLIP:~/.conda/envs/qc2/share/trimmomatic/adapters/TruSeq3-PE.fa:2:40:15 SLIDINGWINDOW:4:20 MINLEN:50' \
      --bowtie2-options '--very-sensitive --dovetail' -db db/Homo_sapiens" \
      ::: `seq 3 5`

	# 样本质控汇总表
	kneaddata_read_count_table \
	  --input qc \
	  --output kneaddata_sum.txt
	# 提取raw、trim和final关键结果
	cut -f 1-5,12-13 kneaddata_sum.txt | sed 's/_1_kneaddata//;s/pair//g' \
	  > kneaddata_report.txt
	cat kneaddata_report.txt

## 5.	(可选)质控后质量再评估 

    fastqc qc/*_1_kneaddata_paired_*.fastq -t 2
    multiqc -d qc/ -o ./


# 常见问题

## 1. 软件下载慢或无法下载

	# 添加清华conda镜像
	site=https://mirrors.tuna.tsinghua.edu.cn/anaconda
	conda config --add channels $site/pkgs/free/ 
	conda config --add channels $site/pkgs/main/
	conda config --add channels $site/cloud/conda-forge/
	conda config --add channels $site/pkgs/r/
	conda config --add channels $site/cloud/bioconda/

## 2. 数据库下载慢或无法下载

	# 很多国外发布数据库下载缓慢，甚至托管于Google或Dropbox等国内无法访问的站点。宏基因组公众号团队建立本领域常用数据库下载的国内备份链接和百度云链接，方便国内同行下载，详见：
	https://github.com/YongxinLiu/MicrobiomeStatPlot/blob/master/Data/BigDataDownlaodList.md

## 3. 非默认物种参考基因组下载和建索引，以拟南芥为例

	wget -c ftp://ftp.ensemblgenomes.org/pub/plants/release-47/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz \
		-O ath.fa.gz
	bowtie2-build ath.fa.gz ath.bt2


## 4. 根据质量评估报告确定接头序列
	# 确定trimmomatic软件位置
	type trimmomatic
	# 根据上面环境路径+share/trimmomatic/adapters目录匹配接头序列的文件
	grep 'ATCGGAAGAGCACACGTCTGAAC' ~/.conda/envs/qc2/share/trimmomatic/adapters/*


## 5. 检查测序双端序列标签是否唯一

	# 查看样本是否标签有重复
	zcat seq/C2_1.fq.gz|head
	zcat seq/C2_2.fq.gz|head

	# 解压、对序列的左、右端标题行分别添\1、\2
	gunzip seq/*.gz
    sed -i '1~4 s/$/\\1/g' seq/*_1.fq
    sed -i '1~4 s/$/\\2/g' seq/*_2.fq
	# 核对样本是否标签有重复
	head seq/C2_1.fq
	head seq/C2_2.fq
	# 压缩节省空间
    gzip seq/*.fq