rm(list=ls())

setwd("/Users/MacProMatt/Desktop/pileup_results_LP005/")

Sample1_1 = read.csv("Sample1_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample2_1 = read.csv("Sample2_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample3_1 = read.csv("Sample3_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample5_1 = read.csv("Sample5_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample5_2 = read.csv("Sample5_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample5_3 = read.csv("Sample5_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample6_1 = read.csv("Sample6_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample6_2 = read.csv("Sample6_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample6_3 = read.csv("Sample6_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample7_1 = read.csv("Sample7_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample7_2 = read.csv("Sample7_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample7_3 = read.csv("Sample7_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample8_1 = read.csv("Sample8_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample8_2 = read.csv("Sample8_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample8_3 = read.csv("Sample8_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample8_4 = read.csv("Sample8_4.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_1 = read.csv("Sample9_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_2 = read.csv("Sample9_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_3 = read.csv("Sample9_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_4 = read.csv("Sample9_4.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample10_1 = read.csv("Sample10_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample10_2 = read.csv("Sample10_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample10_3 = read.csv("Sample10_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample10_4 = read.csv("Sample10_4.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)

grouping_list_names = c('Sample1_1', 'Sample2_1', 'Sample3_1', 'Sample5_1', 'Sample5_2', 'Sample5_3', 'Sample6_1', 'Sample6_2', 'Sample6_3', 'Sample7_1', 'Sample7_2', 'Sample7_3', 'Sample8_1', 'Sample8_2', 'Sample8_3', 'Sample8_4', 'Sample9_1', 'Sample9_2', 'Sample9_3', 'Sample9_4', 'Sample10_1', 'Sample10_2', 'Sample10_3', 'Sample10_4')
sample_counter = 1
grouping_list = list(Sample1_1, Sample2_1, Sample3_1, Sample5_1, Sample5_2, Sample5_3, Sample6_1, Sample6_2, Sample6_3, Sample7_1, Sample7_2, Sample7_3, Sample8_1, Sample8_2, Sample8_3, Sample8_4, Sample9_1, Sample9_2, Sample9_3, Sample9_4, Sample10_1, Sample10_2, Sample10_3, Sample10_4)
chr_list = c('M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2', 'NEP', 'M2')
chr_counter = 1

group_df = data.frame(row.names = grouping_list_names)
chr_snp_df = data.frame(matrix(ncol = 24, nrow = 9) ,row.names = chr_list)
colnames(chr_snp_df) = grouping_list_names

cov = c()
snps = c()
snp_freq = c()
indels = c()
indel_freq = c()
counter = 1

for (i in grouping_list){
  cov = c(cov, sum(i$cov))
  snps = c(snps, sum(i$snps))
  snp_freq = c(snp_freq, sum(i$snps)/sum(i$cov))
  indels = c(indels, sum(i$indels))
  indel_freq = c(indel_freq, sum(i$indels)/sum(i$cov))
  for ( j in chr_list){
    chr_file = subset(i, i$row.names == j)
    chr_snp_df[chr_list[chr_counter], grouping_list_names[sample_counter]] = (sum(chr_file$snps)/(sum(chr_file$cov)))
    chr_counter = chr_counter + 1
  }
  chr_counter = 1
  sample_counter = sample_counter + 1
}

group_df$cov = cov
group_df$SNPs = snps
group_df$SNP_freq = snp_freq
group_df$indels = indels
group_df$indel_freq = indel_freq

subset_samples = c('Sample5_3', 'Sample6_3', 'Sample7_3', 'Sample8_3', 'Sample9_3', 'Sample10_3')
new_chr_frame_24 = chr_snp_df[subset_samples]
new_chr_frame_24 = t(new_chr_frame_24)
row.names(new_chr_frame_24) = c('Sample5_3_24', 'Sample6_3_24', 'Sample7_3_24', 'Sample8_3_24', 'Sample9_3_24', 'Sample10_3_24')

group_df_24 = group_df

boxplot(new_chr_frame_24)
