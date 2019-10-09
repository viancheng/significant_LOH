from __future__ import division
import operator
from scipy import stats
from os import path,makedirs,system
import argparse
import sys
import pandas as pd
import statsmodels.stats.weightstats as st
from collections import Counter

parser = argparse.ArgumentParser()
parser.add_argument('-i','--inputlist',type=str,required=True)
parser.add_argument('-o','--outdir',type=str,required=True)
parser.add_argument('-p','--pvalue',type=float,required=False,default=0.05)
parser.add_argument('-fdr','--fdr',type=float,required=False,default=0.2)
parser.add_argument('-freq','--freq',type=float,required=False,default=0.05)
parser.add_argument('-b','--build',type=str,required=False,default='hg19')
parser.add_argument('-all','--all_cnv',type=str,required=True)
parser.add_argument('-amp','--amps',type=str,required=True)
parser.add_argument('-del','--dels',type=str,required=True)
parser.add_argument('-score','--score',type=str,required=True)
args = parser.parse_args(sys.argv[1:])

if not path.exists(args.outdir):
	makedirs(args.outdir)

def bh_qvalue(pv):
	if pv == []:
		return []
	m=len(pv)
	args,pv=zip(*sorted(enumerate(pv),None,operator.itemgetter(1)))
	if pv[0] < 0 or pv[-1] > 1:
		raise ValueError("p-value must between 0 and 1")
	qvalue=m*[0]
	mincoeff=pv[-1]
	qvalue[args[-1]]=mincoeff
	for j in xrange(m-2,-1,-1):
		coeff=m*pv[j]/float(j+1)
		if coeff < mincoeff:
			mincoeff=coeff
		qvalue[args[j]]=mincoeff
	return qvalue

def p_fisher(lst):
	p_value=stats.fisher_exact(lst)[0]
	return p_value

def t_test(lst1,lst2):
	samp1=pd.Series(lst1)
	samp2=pd.Series(lst2)
	samp1_var=samp1.var()
	samp2_var=samp2.var()
	df1=len(lst1)-1
	df2=len(lst2)-1
	if samp1_var > samp2_var:
		F=samp1_var/samp2_var
		p_value=stats.f.sf(F,df1,df2)
	else:
		F=samp2_var/samp1_var
		p_value=stats.f.sf(F,df2,df1)
	if p_value<args.pvalue:
		t,p_2tailed,df=st.ttest_ind(samp1,samp2,usevar='unequal')
	else:
		t,p_2tailed,df=st.ttest_ind(samp1,samp2,usevar='pooled')
	return p_2tailed

def get_idx(name,lst):
	idx=[]
	for i in name:
		idx.append(lst.index(i))
	return idx

def rm_idx(name,glist,*lst):
	idx=[]
	for i,g in enumerate(glist):
		if g == name:
			idx.append(i)
	for i in sorted(idx,reverse=True):
		del glist[i]
		for l in lst:
			del l[i]

def get_matrix(idx,lst):
	a=lst[idx[0]]
	b=lst[idx[1]]
	c=lst[idx[2]]
	d=lst[idx[3]]
	return [[a,b],[c,d]]

def overlap(lst1,lst2):
	for g in lst1:
		if g in lst2:
			return True
	return False

def getpos(build,chromosome,start,end):
	pos=[0,0]
	chrlen=[]
	if build == 'hg19':
		chrlen=[249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566]
	elif build == 'hg18':
		chrlen=[247249719, 242951149, 199501827, 191273063,180857866, 170899992, 158821424, 146274826, 140273252,135374737, 134452384, 132349534, 114142980, 106368585,100338915, 88827254, 78774742, 76117153, 63811651,62435964, 46944323, 49691432, 154913754, 57772954]
	elif build == 'hg38':
		chrlen=[248956422, 242193529, 198295559, 190214555,181538259, 170805979, 159345973, 145138636, 138394717,133797422, 135086622, 133275309, 114364328, 107043718,101991189, 90338345, 83257441, 80373285, 58617616,64444167, 46709983, 50818468, 156040895, 57227415]
	chrname=map(str,range(1,23))
	chrname.append('X')
	chrname.append('Y')
	if chromosome in chrname:
		idx=chrname.index(str(chromosome))
		for i,length in enumerate(chrlen):
			if i<idx:
				pos[0] += length
				pos[1] +=length
		pos[0] += int(start)
		pos[1] += int(end)
	return pos
	
def updata_dic(dic1,dic2):
	for k2,v2 in dic2.iteritems():
		if k2 not in dic1.keys():
			dic1[k2]=v2
		else:
			dic1[k2].extend(v2)
	return dic1

def detect_share(gene,dic):
	count=0
	cutoff=len(dic)*args.freq
	for samp,infos in dic.iteritems():
		for info in infos:
			if gene in info.split('\t')[10].split(';'):
				count+=1
				break
	if count > cutoff:
		return True
	else:
		return False
			
def process_main(lst,build):
	genelist={}
	genepos={}
	out={}
	mark=0
	for i in open(lst,'r'):
		i=i.strip()
		cnv,maf=i.split('\t')
		loh={}
		count=[0,0,0]
		sample=path.basename(cnv).split('.')[0]
		for line in open(cnv,'r').readlines():
			line=line.strip()
			lst=line.split('\t')
			if line.startswith('chrom') and mark == 0:
				individual_out.writelines("sample\t"+line+'\n')
				population_out.writelines("sample\t"+line+'\tsigLOH-gene\n')
				mark=1
				continue
			if 'LOH' in lst[8]:
				chromosome=lst[0]
				count[0]+=1
				loh.setdefault(chromosome,[]).append(line)
		flag=0
		pv=[]
		pv_info=[]
		pv_gene=[]
		head=['n_alt_count','n_ref_count','t_alt_count','t_ref_count']
		headidx=[]
		genelist_tmp={}
		for line in open(maf,'r').readlines():
			line=line.strip()
			lst=line.split('\t')
			if flag==0:
				headidx=get_idx(head,lst)
				flag=1
				continue
			if line.startswith('Hugo_Symbol'):
				continue
			if lst[8] in ["3'UTR","5'Flank","5'UTR","IGR","Intron","lincRNA","RNA","Silent"]:
				continue
			chromosome,start,end=lst[4:7]
			gene=lst[0]
			genepos.setdefault(gene,[]).append(getpos(build,chromosome,start,end)[0])
			vaf=float(lst[headidx[2]])/float(lst[headidx[2]]+lst[headidx[3]])
			genelist_tmp.setdefault(gene,[]).append(vaf)
			for ranges in loh.setdefault(chromosome,[]):
				myranges=ranges.split('\t')
				if gene in myranges[10].split(';'):
					pv.append(p_fisher(get_matrix(headidx,lst)))
					pv_info.append(ranges)
					pv_gene.append(gene)
					
		for g in genelist_tmp.keys():
			if g not in pv_gene:
				del genelist_tmp[g]
		q=bh_qvalue(pv)
		tmplst=[]
		for i,n in sorted(enumerate(q),None,operator.itemgetter(0),reverse=True):
			if not (n <= 0.01 and pv[i] <= 0.01):
				if genelist_tmp.get(pv_gene[i]):
					del genelist_tmp[pv_gene[i]]
				del q[i],pv[i],pv_info[i],pv_gene[i]
		count_result=Counter(pv_gene)
		for g,i in count_result.iteritems():
			if i < 1:
				if genelist_tmp.get(g):
					del genelist_tmp[g]
				rm_idx(g,pv_gene,q,pv,pv_info)
		for i,n in enumerate(q):
			if (sample+'\t'+str(pv_info[i])) not in tmplst:
				count[1]+=1
				tmplst.append(sample+'\t'+str(pv_info[i]))
		out.setdefault(sample,tmplst)
		stats_out[sample]=count
		individual_out.writelines('\n'.join(tmplst))
		genelist=updata_dic(genelist,genelist_tmp)
	individual_out.close()
	return out,genelist,genepos

stats_out={}
individual_out=open(args.outdir+'/individual_sig.txt','w')
population_out=open(args.outdir+'/population_sig.txt','w')
gene_out=open(args.outdir+'/sigLOH_gene.txt','w')
sigloh_out=open(args.outdir+'/LOH_stats.txt','w')
gene_out.writelines("gene\tpos\n")
sigloh_out.writelines("sample\torg-loh\tind-loh\tpop-loh\n")
sample_info,genelist,genepos=process_main(args.inputlist,args.build)
p_Ttest=[]
p_Ttest_info=[]
for gene,vaflst in genelist.iteritems():
	lst1=vaflst
	lst2=[]
	for gene2,vaflst2 in genelist.iteritems():
		if gene2 != gene:
			lst2.extend(vaflst2)
	if lst1 != [] and lst2 != []:
		p_Ttest.append(t_test(lst1,lst2))
		p_Ttest_info.append(gene)

q=bh_qvalue(p_Ttest)
tmplst=[]
tmp1={}
for i,n in enumerate(q):
	if n <= args.fdr and p_Ttest[i] <= args.pvalue:
		if not detect_share(p_Ttest_info[i],sample_info):
			continue
		for samp,infos in sample_info.iteritems():
			tmp={}
			for info in infos:
				if p_Ttest_info[i] in info.split('\t')[10].split(';'):
					tmplst.append(p_Ttest_info[i]+'\t'+str(genepos[p_Ttest_info[i]][0]))
					tmp.setdefault(info,[]).append(samp)
					tmp1.setdefault(info,[]).append(p_Ttest_info[i])
			stats_out[samp][2]=len(tmp)
for k,v in tmp1.iteritems():
	population_out.writelines(k+'\t'+';'.join({}.fromkeys(v).keys())+'\n')
population_out.close()

for i in {}.fromkeys(tmplst).keys():
	gene_out.writelines(i+'\n')
gene_out.close()

for s,c in stats_out.iteritems():
	sigloh_out.writelines(s+'\t'+'\t'.join(map(str,c))+'\n')
sigloh_out.close()

cmd="/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/bin/Rscript /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/sig_LOH/gisticplot.r %s %s %s %s %s/sigLOH_gene.txt %s %s/sig_cnv_loh.pdf" % (args.all_cnv,args.amps,args.dels,args.score,args.outdir,args.build,args.outdir)
system("export R_LIBS=/zfssz5/BC_PS/chengyuanfang/lib/R:/hwfssz4/BC_COM_FP/bc_tumor/Pipeline/Cancer_Genomics_Interpretation_Analysis/CGIA_v2.0/lib/R_LIBS:/zfssz5/BC_PS/chengyuanfang/lib/R:$R_LIBS")
system(cmd)
