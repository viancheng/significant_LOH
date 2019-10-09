#!/bin/bash
set -e
set -o pipefail
echo hostname: `hostname`
echo ==========start at : `date` ==========

python sig_loh.py -i input/input.list -o output -p 0.01 -fdr 0.01 -all gistic/all_lesions.conf_95.TumorName.txt -amp gistic/amp_genes.conf_95.TumorName.txt -del gistic/del_genes.conf_95.TumorName.txt -score gistic/scores.TumorName.gistic -b hg19 -freq 0.2

if [ -e /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/sig_LOH/work.sh.e* ]; then ls /zfssz5/BC_PS/chengyuanfang/work/CGIA_3.0_module/sig_LOH/work.sh.e* | while read dd ; do jobid=`echo $dd | sed -r 's/^.*\.e([0-9]*)$/\1/g'`; echo $jobid; qstat -j $jobid | grep 'usage'; done ; fi && \
echo ==========end at : `date` ==========

