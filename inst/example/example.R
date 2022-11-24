demo_cn <- as.data.frame(matrix(nrow = 1, ncol = 5, data = c(1,50000,1000000,2,1)))
colnames(demo_cn) <- c("chr","startpos","endpos","nMaj1_A","nMin1_A")
demo_cn$startpos <- as.numeric(demo_cn$startpos)
demo_cn$endpos <- as.numeric(demo_cn$endpos)
demo_cn$nMaj1_A <- as.numeric(demo_cn$nMaj1_A)
demo_cn$nMin1_A <- as.numeric(demo_cn$nMin1_A)

demo_mult <- as.data.frame(matrix(nrow = 20, ncol = 3, data = c(rep("1",20),57682:57701,rep(2,5),rep(1,15))))
colnames(demo_mult) <- c("chr","end","no.chrs.bearing.mut")
demo_mult$end <- as.numeric(demo_mult$end)
demo_mult$no.chrs.bearing.mut <- as.numeric(demo_mult$no.chrs.bearing.mut)

demo_muts <- as.data.frame(matrix(nrow = 20, ncol = 5, data = c(rep("1",20),57682:57701,57682:57701,rep("C",15),rep("A",5),rep("T",15),rep("G",5))))
colnames(demo_muts) <- c("chr", "start", "end", "ref", "alt")
demo_muts$start <- as.numeric(demo_muts$start)
demo_muts$end <- as.numeric(demo_muts$end)

time_amplification(
  cn_data = demo_cn,
  multiplicity_data = demo_mult,
  mutation_data = demo_muts,
  muts_type = "All",
  sample_id = "test_sample",
  amplification_chrom = 1,
  amplification_start = 50100,
  amplification_stop = 1500000,
  is_WGD = TRUE,
  genome = "hg19"
)