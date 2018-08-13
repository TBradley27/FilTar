#!/usr/bin/env python

import sqlite3
import sys
import os
from subprocess import call
import pandas
from snakemake.utils import report


data = pandas.read_table("{}".format(snakemake.input))

data.columns=['mrna_id','tissue_id','utr_length']

data.index.names = ['id']

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("DELETE FROM {};".format(snakemake.params['table']))
conn.commit()
#c.execute("PRAGMA foreign_keys = ON")
#conn.commit()

#call("sqlite3 filtar.db -separator '\t' '.import {} {}'".format(sys.argv[1], sys.argv[2]), shell=True)

data.to_sql(snakemake.params['table'], con=conn, if_exists='append')

conn.commit()
conn.close()
call(['touch', "{}".format(snakemake.output)])
