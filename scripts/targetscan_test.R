library(testthat)

output_mmu = readr::read_tsv('../results/targets/mmu/oocyte_chrY.contextpp.tsv')
output_mmu = dplyr::filter(output_mmu, `Mirbase ID` == 'mmu-miR-188-5p') 

print(output_mmu)

context('targetscan')

test_that('test filtar targetscan values match those on the official targetscan website', {
        expect_equal(output_mmu$`context++ score`[output_mmu$`Gene ID` == 'ENSMUST00000189888.7'],-0.25,0.025)}
)
