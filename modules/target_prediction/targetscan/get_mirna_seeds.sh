#!/bin/bash

awk '{print $1}' $1 | sed 's/--//g' | sed '/^$/d' | sed 'N;s/\n/ /' | sed 's/>//g' | awk '{ print $1, substr($2,2,7) }' | awk '{ print $1,$2}' 
