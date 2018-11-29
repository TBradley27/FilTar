shuf -n 1000 Homo_sapiens.GRCh38.92.gtf > hsa_mock.gtf
shuf -n 1000 Mus_musculus.GRCm38.92.gtf > mmu_mock.gtf

gtf_to_bed.sh hsa_mock.gtf 3UTR 1 hsa_mock_3UTR.bed
gtf_to_bed.sh hsa_mock.gtf CDS 1 hsa_mock_CDS.bed
gtf_to_bed.sh mmu_mock.gtf 3UTR 1 mmu_mock_3UTR.bed
gtf_to_bed.sh mmu_mock.gtf CDS 1 mmu_mock_CDS.bed
