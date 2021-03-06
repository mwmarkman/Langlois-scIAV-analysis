rm(list=ls())
#The plots at the bottom rely on LP005 and LP009. They can be correctly generated by runninng LP005_individual.R followed by LP009_individual.R (this file)

setwd("/Users/MacProMatt/Desktop/pileup_results_LP009/")

Sample1_1 = read.csv("Sample1_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample2_1 = read.csv("Sample2_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample3_1 = read.csv("Sample3_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample4_1 = read.csv("Sample4_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample4_2 = read.csv("Sample4_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample4_3 = read.csv("Sample4_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
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
Sample9_1 = read.csv("Sample9_1.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_2 = read.csv("Sample9_2.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)
Sample9_3 = read.csv("Sample9_3.csv", header = T, na.strings = "na", col.names = c("pos", "con", "cov", "snps", "indels", "freq" ), row.names = NULL)

grouping_list_names = c('Sample1_1', 'Sample2_1', 'Sample3_1', 'Sample4_1', 'Sample4_2', 'Sample4_3', 'Sample5_1', 'Sample5_2', 'Sample5_3', 'Sample6_1', 'Sample6_2', 'Sample6_3', 'Sample7_1', 'Sample7_2', 'Sample7_3', 'Sample8_1', 'Sample8_2', 'Sample8_3', 'Sample9_1', 'Sample9_2', 'Sample9_3')
sample_counter = 1
grouping_list = list(Sample1_1, Sample2_1, Sample3_1, Sample4_1, Sample4_2, Sample4_3, Sample5_1, Sample5_2, Sample5_3, Sample6_1, Sample6_2, Sample6_3, Sample7_1, Sample7_2, Sample7_3, Sample8_1, Sample8_2, Sample8_3, Sample9_1, Sample9_2, Sample9_3)
chr_list = c('M1', 'NA', 'NP', 'NS1', 'PA', 'PB1', 'PB2', 'NEP', 'M2')
chr_counter = 1

group_df = data.frame(row.names = grouping_list_names)
chr_snp_df = data.frame(matrix(ncol = 21, nrow = 9) ,row.names = chr_list)
colnames(chr_snp_df) = grouping_list_names

cov = c()
snps = c()
snp_freq = c()
indels = c()
indel_freq = c()

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

subset_samples = c('Sample4_3', 'Sample5_3', 'Sample6_3', 'Sample7_3', 'Sample8_3', 'Sample9_3')
subset_samples2 = c('Sample4_2', 'Sample5_2', 'Sample6_2', 'Sample7_2', 'Sample8_2', 'Sample9_2')
new_chr_frame_12_2 = chr_snp_df[subset_samples2]
new_chr_frame_12 = chr_snp_df[subset_samples]
new_chr_frame_12 = t(new_chr_frame_12)
new_chr_frame_12_2 = t(new_chr_frame_12_2)
row.names(new_chr_frame_12) = c('Sample4_3_12', 'Sample5_3_12', 'Sample6_3_12', 'Sample7_3_12', 'Sample8_3_12', 'Sample9_3_12')

total_df = rbind(new_chr_frame_12, new_chr_frame_24)

group_df_12 = group_df
group_df_both = rbind(group_df_12, group_df_24)

snp_neg_single_12 = c(0.002461162, 0.002543369, 0.002868077)
snp_neg_seq_12 = c(0.002687196, 0.003122100, 0.003049064)
snp_neg_single_24 = c(0.002454447, 0.002629865, 0.002718284)
snp_neg_seq_24 = c(0.002742065, 0.002533088, 0.002641329)

snp_low_single_12 = c(0.003828631, 0.003903105, 0.003643967)
snp_low_seq_12 = c(0.004248529, 0.004079616, 0.004066927)
snp_low_single_24 = c(0.003677945, 0.004005888, 0.003679047)
snp_low_seq_24 = c(0.005159600, 0.004241783, 0.004945732)

snp_high_single_12 = c(0.011470548, 0.011517281, 0.010207231)
snp_high_seq_12 = c(0.011747194, 0.011584108, 0.011400116)
snp_high_single_24 = c(0.010559125, 0.010817415, 0.010482063)
snp_high_seq_24 = c(0.009734451, 0.010949643, 0.010201426)


neg_array_indel = c(0.0004548349, 0.0004174526, 0.0003851417, 0.0004852798, 0.0005645812, 0.0005511848, 0.0005532906, 0.0006757857, 0.0005895550, 0.0006918636, 0.0005908179, 0.0006000639)
low_array_indel = c(0.0005242524, 0.0005424169, 0.0005363677, 0.0006380003, 0.0006328010, 0.0006180884, 0.0008318118, 0.0008395115, 0.0007636351, 0.0011424306, 0.0009524809, 0.0010354586)
high_array_indel = c(0.0016085174, 0.0016340791, 0.0012312646, 0.0015615634, 0.0016024338, 0.0015679428, 0.0019167680, 0.0021086525, 0.0021638305, 0.0020478714, 0.0019718784, 0.0018754347)

#par(mfrow = c(1,1))
#layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))



new_chr_frame_12 = new_chr_frame_12[-c(4, 5, 6), ]
write.csv(new_chr_frame_12, file = "/Users/MacProMatt/Desktop/Chr.csv")
file = read.csv("/Users/MacProMatt/Desktop/Chr.csv")

jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
chr_plot = boxplot(new_chr_frame_12, main = "mCherry High, n = 3", xlab = "IAV Gene Segment", col = 'lightcoral', boxwex = 0.75, las = 2)
dev.off()

jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
chr_plot_low = boxplot(new_chr_frame_12_2, main = "mCherry Low, n = 3", xlab = "IAV Gene Segment", col = 'lightcoral', boxwex = 0.75, las = 2)
dev.off()


jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
low_high_plot_snp = boxplot(snp_neg_single_12, snp_low_single_12, snp_high_single_12, boxwex = 0.75, col = 'lightcoral', boxwex = 0.75, main = "Single Base Variations", xlab = "n = 3", ylab = "Per Base Mutation Frequency", names = c("Neg", "Low","High"), log = 'y')
dev.off()

jpeg(filename = '/Users/MacProMatt/Desktop/all_grey.jpeg', height = 4, width = 5, units = "in", res = 800, pointsize = 12)
low_high_plot_indel = boxplot(neg_array_indel, low_array_indel, high_array_indel, log = 'y', boxwex = 0.75, names = c("Negative", "Low", "High"), col = 'lightcoral', main = 'Insertions/Deletions', xlab = 'n = 12', ylab = "Per Base Mutation Frequency")
dev.off()



















hourhpi_plot_high = boxplot(c(0.011470548, 0.011517281, 0.010207231), c(0.011747194, 0.011584108, 0.011400116), c(0.010559125, 0.010817415, 0.010482063), c(0.009734451, 0.010949643, 0.010201426), names = c("12 hpi", "12 hpi", "24 hpi", "24 hpi"), xlab = "n = 3", ylab = 'Mean Per Base Mutation Frequency', main = "mCherry High", col = c('lightcoral', 'light green', 'light coral', 'light green'), boxwex = 0.70, ylim = c(0.0098, 0.0133))
legend('topright' , legend = c('Single Infection', 'Sequential Infection'), pch = 15, col =  c('lightcoral', 'light green', 'light coral', 'light green'))


#first two (NS)
t.test(c(0.011470548, 0.011517281, 0.010207231), c(0.011747194, 0.011584108, 0.011400116), var.equal = T)

#Last two (NS)
t.test(c(0.010559125, 0.010817415, 0.010482063), c(0.009734451, 0.010949643, 0.010201426), var.equal = T)

#Reds 
t.test(c(0.011470548, 0.011517281, 0.010207231), c(0.010559125, 0.010817415, 0.010482063), var.equal = T)

#Greens (0.02523)
t.test(c(0.011747194, 0.011584108, 0.011400116), c(0.009734451, 0.010949643, 0.010201426), var.equal = T)


hourhpi_plot_low = boxplot(c(0.003828631, 0.003903105, 0.003643967), c(0.004248529, 0.004079616, 0.004066927), c(0.003677945, 0.004005888, 0.003679047), c(0.005159600, 0.004241783, 0.004945732), names = c("12 hpi", "12 hpi", "24 hpi", "24 hpi"), xlab = "n = 3", ylab = 'Mean Per Base Mutation Frequency', main = "mCherry Low", col = c('lightcoral', 'light green', 'light coral', 'light green'), boxwex = 0.70, ylim = c(0.0036, 0.00555))
legend('topleft' , legend = c('Single Infection', 'Sequential Infection'), pch = 15, col =  c('lightcoral', 'light green', 'light coral', 'light green'))

#first two 
t.test(c(0.003828631, 0.003903105, 0.003643967), c(0.004248529, 0.004079616, 0.004066927), var.equal = T)

#last two
t.test(c(0.003677945, 0.004005888, 0.003679047), c(0.005159600, 0.004241783, 0.004945732), var.equal = T)

#Reds (NOT SIGNIFICANT)
t.test(c(0.003828631, 0.003903105, 0.003643967),c(0.003677945, 0.004005888, 0.003679047), var.equal = T)

#Greens ()
t.test(c(0.004248529, 0.004079616, 0.004066927), c(0.005159600, 0.004241783, 0.004945732), var.equal = T)

df = as.data.frame(new_chr_frame_12)
df = t(df)
