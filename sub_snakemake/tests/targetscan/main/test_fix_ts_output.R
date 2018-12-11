library(plyr)
library(tidyverse)
library(filtar)

#output = fix_ts_output('hsa_chrY_msa.sites.new_method.tsv') #%>% as.data.frame()
output_canonical=read_tsv('hsa_chrY_msa.sites.tsv') # this is output from the unpatched script
output = read_tsv('tmp.tsv') # this is output from patched script

#write.table(output, 'tmp.tsv',sep='\t',quote=FALSE,row.names=FALSE)



#test1 = map(output$Group_type, function (x) {strsplit(x,'+',fixed=TRUE) %>% unlist %>% duplicated})

#test1 = unlist(test1)

#test_that('Group type column contains no duplicate site types',{
#	expect_match(test1 %>% as.character, 'FALSE')
#})

test_that('Canonical and patched script output is identical',{
        expect_equal(output, output_canonical)
})
