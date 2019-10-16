# 群体过滤脚本使用说明

## 示例脚本
perl mis_het_filter_vcf_V7.pl \
  -vcf vcf_file \
  -out outpath \
  -gq_all 30 \
  -dp_min1 4 \
  -dp_max1 100 \
  -DP_min 300 \
  -DP_max 8000000 \
  -gq_single 10 \
  -miss 0.99 \
  -maf 0.01 \
  -num 2034 \
  -dp_pos 2 \
  -gq_pos 3 \
  -name merge.Chr15.g.vcf.geno.vcf \
  -het 0.2
