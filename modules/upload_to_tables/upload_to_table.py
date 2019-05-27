#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call
import pandas

print (snakemake.params['columns'])

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

#c.execute("DELETE FROM {};".format(snakemake.params['table']))
conn.commit()
c.execute("PRAGMA foreign_keys = ON")
conn.commit()

data = pandas.read_csv(snakemake.input['data'], sep="\t")

print(data)

data.columns = snakemake.params['columns']
data.index.names = ['id']

#df.to_sql('mRNA', conn, if_exists='append', index=False)

#call("sqlite3 filtar.db -separator '\t' '.import {} {}'".format(sys.argv[1], sys.argv[2]), shell=True)

data.to_sql('{}'.format(snakemake.params['table']), con=conn, if_exists='append')

conn.commit()
conn.close()
call(['touch', "{}".format(snakemake.output)])
