##use LDBlockShow, ref: https://github.com/BGI-shenzhen/LDBlockShow
/mnt/d/Population_643/all_761_info/population_analysis/LDBlockShow-1.40/bin/LDBlockShow -InVCF MYB5_457
_2KB_use_20_maf0.05_variance.vcf -OutPut out -InGWAS gwas_pvalue.txt -InGFF In.gff.txt -Region 7:101700479-101705734 #生成基本图形
/mnt/d/Population_643/all_761_info/population_analysis/LDBlockShow-1.40/bin/ShowLDSVG -InPreFix out -OutPut out.svg -InGWAS gwas_pvalue.txt -Cutline 19.4 -InGFF In.gff.txt -crGene yellow:lightblue:pink:orange -shownum -OutPdf -SpeSNPName Spe.snp.txt  -ShowGWASSpeSNP -XYLabFontSizeRatio 2 ##调整图形颜色，字体等
