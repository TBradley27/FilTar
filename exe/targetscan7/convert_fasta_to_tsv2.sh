#!/bin/bash

cat $1 | sed 's/(.*//g' | tr '\n' '\t' | tr '>' '\n' | sed "s/\t/\t$2\t/1" | sed '/^$/d' 
