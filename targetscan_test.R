library(testthat)

output_mmu = read_tsv('results/targets/no_reannotation/mmu_chrY.contextpp.tsv')

context('targetscan')

test_that('test targetscan values match those on the website', {
        expect_equal(output_mmu$`context++ score`[output_mmu$`Gene ID` == 'ENSMUST00000055032.13'],-0.18,0.018) # 10% tolerance

