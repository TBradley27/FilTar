#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call
import pandas

data = pandas.read_table(snakemake.input['data'])

data = data[['Name','avg']]

data['tissue'] = snakemake.wildcards['tissue']

data.columns = ['mrnas_id','TPM','experiments_id']
data.index.names = ['id']

data['mrnas_id'].replace(regex=True,inplace=True,to_replace=r'\..*',value=r'')


conn = sqlite3.connect('filtar.db')
c = conn.cursor()

#c.execute("DELETE FROM {};".format(sys.argv[2]))
#conn.commit()
#c.execute("PRAGMA foreign_keys = ON")
#conn.commit()

#call("sqlite3 filtar.db -separator '\t' '.import {} {}'".format(sys.argv[1], sys.argv[2]), shell=True)

data.to_sql('expression_profiles', con=conn, if_exists='append')

conn.commit()
conn.close()
call(['touch', snakemake.output[0]])
