#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call
import pandas

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("DELETE FROM {};".format(sys.argv[2]))
conn.commit()
c.execute("PRAGMA foreign_keys = ON")
conn.commit()

#df = pandas.read_csv(sys.argv[1])
#df.to_sql('mRNA', conn, if_exists='append', index=False)

call("sqlite3 filtar.db -separator '\t' '.import {} {}'".format(sys.argv[1], sys.argv[2]), shell=True)

conn.commit()
conn.close()
call(['touch', "{}".format(sys.argv[3])])
