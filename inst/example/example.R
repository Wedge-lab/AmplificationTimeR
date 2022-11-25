library(BSgenome.Hsapiens.UCSC.hg19)
data("demo_cn", "demo_mult", "demo_muts")

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
