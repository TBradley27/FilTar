#       FilTar: Integrating RNA-Seq data to improve microRNA target prediction accuracy in animals
#       Copyright (C) 2019 Thomas Bradley
#       License summary: https://github.com/TBradley27/FilTar/LICENSE_NOTICE

convert_fa_to_tsv () {

awk '{print $1}' $1 | sed 's/--//g' | sed '/^$/d' | sed 'N;s/\n/ /' | sed 's/>//g' | awk '{ print $1, substr($2,2,7) }' | awk '{OFS="\t"}{print $1,$2}' | grep -E "^[a-z][a-z][a-z]-(miR|let|lin)-?[0-9]{1,5}[a-z]?(-[0-9]-)?((-3|5)p)?"

}
