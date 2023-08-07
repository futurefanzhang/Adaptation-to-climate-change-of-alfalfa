##step 1, add header
cut -f1-3  all_399_indel_plink_addid.vcf > all_399_indel_position.txt ##extract pos info
perl /public/agis/zhouyongfeng_group/zhangfan02/vcf_file/match_indel_position.pl all_399_indel_position.txt ./result_lfmm/lfmm_pvalue_maf0.1_new.LFMM_BIO2.env.txt ./result_lfmm/add_pos_indel_maf0.1_bio2.txt
