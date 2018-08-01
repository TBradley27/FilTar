#!/bin/bash

err() {
  echo "[$(date +'%Y-%m-%dT%H:%M:%S%z')]: $@" >&2; exit 1
}

bedtools getfasta -name -s -fi $1 -fo $3 -bed $2

cat $3 | sed 's/:.*$//g' > tmp"$4" && mv tmp"$4" $3
