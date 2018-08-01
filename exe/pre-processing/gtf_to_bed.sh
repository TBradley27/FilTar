if [[ "$1" == '-h' ]] || [[ "$1" == '' ]]; then
  echo 'Arguments:  tx_region   chromosome_number  species_triplet_code  
  --  e.g.  'three_prime_utr'    2   hsa biopython'; exit 1
fi

tx_region="$1"; chromosome="$2"; species="$3"; format="$4" # aliases

if [[ "$species" == 'hsa' ]]; then
  gtf_file=data/Homo_sapiens.GRCh38.92.gtf
elif [[ "$species" == 'mmu' ]]; then
  gtf_file=data/Mus_musculus.GRCm38.92.gtf
else
  echo "Not a valid species code - "$species" "; exit 1
fi

if [[ "$tx_region" == '3UTR' ]]; then
  path_var_one=three_UTR; path_var_two=3UTR; tx_region_long=three_prime_utr
elif [[ "$tx_region" == 'CDS' ]]; then
  path_var_one=CDS; path_var_two=CDS; tx_region_long=CDS
else
  echo "Please enter a valid genomic feature"; exit 1
fi

grep "$tx_region_long" "$gtf_file" |
grep -P "${tx_region_long}\t"|         # Braces needed for correct search
grep '^'$chromosome'\s' > tmp"$chromosome"  # CCDS match confounds CDS search

if [[ "$?" -ne 0 ]] || [[ ! -s tmp"$chromosome"  ]]; then
  echo "Unable to complete initial processing gg" >&2; exit 1
fi

if [[ "$format" == 'biopython' ]]; then
  awk '{print $1,$4,$5,$7,$14$16}' tmp"$chromosome" |
    sed 's/+/1/g'   |
    sed 's/-/-1/g'  |
    sed 's/\"//g'   |
    sed 's/;/./g'   |
    sed 's/\.//2'   |
    tr ' '  \\t     |
    #awk '$2!=$3'    |  # the next line converts from 1 to 0-based indexing
   awk '{ OFS="\t" }{print $1,$2-1,$3,$4,$5,$6}'   > results/"$species"_"$path_var_two"_biopython.chr"$chromosome".bed
    
    if [[ "$?" -ne 0 ]]; then
      err "Unable to complete biopython formatting" >&2
    fi
  
elif [[ "$format" == 'bedtools' ]]; then
  awk '{print $1,$4,$5,$14$16,"1",$7}' tmp"$chromosome" |
    sed 's/\"//g'  |
    sed 's/;/./g'  |
    sed 's/\.//2'   |
    tr ' '  \\t    |
    #awk '$2!=$3' | # the next line converts from 1 to 0-based indexing
    awk '{ OFS="\t" }{print $1,$2-1,$3,$4,$5,$6}'  > results/"$species"_"$path_var_two"_chr"$chromosome".bed

    if [[ "$?" -ne 0 ]]; then
      err "Unable to complete bedtools formatting" >&2
    fi
    
else
  echo "Please enter a valid desired bed file output format"; exit 1
fi

rm tmp"$chromosome"
