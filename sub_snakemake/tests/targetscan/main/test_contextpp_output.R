library(tidyverse)

# Athough I did not write the targetscan code myself, this test is nonetheless a test for the correctness of the environment (i.e. targetscan dependencies) and the wrapper around targetscan for it to be executed correctly

output = read_tsv('hsa_mock_contextpp.tsv')
print(output$`context++ score`)

test_that('ENST00000393577.3 has the correct context++ score', {
	expect_equal(output$`context++ score`[output$`Gene ID` == 'ENST00000393577.3'],-0.298,0.0298) # 10% tolerance
})

test_that('ENST00000354636.3 has the correct context++ score', {
        expect_equal(output$`context++ score`[output$`Gene ID` == 'ENST00000354636.3'],-0.1760,0.0176) # 10% tolerance
})
