
import Snp_windows

main = undefined

testDat_Snps = [
    SNP "chrom_1" 1 0.424
  , SNP "chrom_1" 6 0.123
  , SNP "chrom_1" 24 0.13
  , SNP "chrom_1" 49 0.1023
  , SNP "chrom_1" 104 0.943
  , SNP "chrom_1" 599 0.904123
  , SNP "chrom_1" 100043 0.332
  , SNP "chrom_1" 1000434 0.123
  , SNP "chrom_2" 43 0.23
  , SNP "chrom_2" 43 1.0
  , SNP "chrom_2" 43 0
  , SNP "chrom_2" 43 0.2
  ]

testConfig = WindowConfig 100 10 
testState = initWindows (windowSize testConfig) (windowStep testConfig) "chrom_1"


