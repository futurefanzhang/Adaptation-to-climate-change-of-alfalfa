/public/agis/zhouyongfeng_group/caoshuo/tools/plink/plink --vcf merge_712_allele2_meanDP3-62_miss0.2_addid.plink.vcf.vcf --keep 51_sample_fastsimcoal.txt --allow-extra-chr --recode vcf --out merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink --threads 1 ##select used individual
java -Xmx10g -jar /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/beagle.22Jul22.46e.jar nthreads=10 gt=merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink.vcf out=merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink_beagle.vcf #use beagle for impute
gunzip merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink_beagle.vcf.vcf.gz
###refer https://github.com/xhchauvet/superSFS for superSFS.py
python /public/home/zhangtianhao/soft/superSFS/superSFS.py 1 10fastsimcoal_outgroup.txt 3 merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink_beagle.vcf.vcf merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink_beagle_reverse.vcf  
###function 1，speculate ancestory allele of vcf file by outgroup；10fastsimcoal_outgroup.txt is outgroup，3 means if there are 3 alleles of outgroup (total 20 alleles among 10 outgroup), the allele will be transfer by outgroup (from ref to alter or from alter to ref)，
python /public/home/zhangtianhao/soft/superSFS/superSFS.py 2 51_sample_group_fastsimcoal.txt merge_51_sample_caerulea_falcata_diploid_outgroup_allele2_meanDP3-62_miss0.2_addid_plink_beagle_reverse.vcf ./count ##function 2, generate count value for fastsimcoal

python 2psfs.py ##change count to SFS file
/public/home/zhangtianhao/soft/fsc27_linux64/fsc2709 -t QT1.tpl -e QT1.est -d  -n 100000 -L 50 -s 0 -M -c 8 -C 10 -B 40
