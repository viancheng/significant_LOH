使用方法（work.sh是一个示例）：
python sig_loh.py -i input/input.list -o output -p 0.01 -fdr 0.2 -all gistic/all_lesions.conf_95.TumorName.txt -amp gistic/amp_genes.conf_95.TumorName.txt -del gistic/del_genes.conf_95.TumorName.txt -score gistic/scores.TumorName.gistic -b hg19

输入参数说明：
-i input.list——输入文件列表，该文件含有两列，对一列为CNV文件路径，第二列为maf文件路径，每行为一个样本；
-o outdir——输出目录，不用提前建立；
-p float——p-value显著性阈值；
-fdr float——fdr显著性阈值；
-all file——gistic的输出文件all_lesions.conf_95.*.txt文件；
-amp file——gistic的输出文件amp_genes.conf_95.*.txt文件；
-del file——gistic的输出文件del_genes.conf_95.*.txt文件；
-score file——gistic的输出文件scores.*.gistic文件
-b str——参考基因组版本，支持输入hg18、hg19、hg38；
-freq float——基因在人群中的最小频率阈值；

输出文件说明（在outdir下共有5个文件）：
individual_sig.txt——个体水平显著LOH；
population_sig.txt——群体水平显著LOH；
sigLOH_gene.txt——群体水平显著LOH包含的显著cosmic gene以及该基因在基因组上的坐标位置；
LOH_stats.txt——LOH统计结果，文件四列分别是样本名称、原始LOH数目、个体显著LOH数目、群体显著LOH数目；
sig_cnv_loh.pdf——显著LOH与显著CNV图形化结果；

注意事项：
1.使用python2运行；
2.画图的R脚本为/zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/sig_LOH/gisticplot.r，如果需要修改图片，请拷贝后修改副本再运行；

