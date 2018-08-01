#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='mRNA';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE mRNA (
                mRNA_ID VARCHAR(20) DEFAULT NULL,
		annotator VARCHAR(30) DEFAULT NULL,
                annotation_version VARCHAR(20) DEFAULT NULL,
		Gene_ID VARCHAR(20) DEFAULT NULL,
		id INTEGER PRIMARY KEY,
                CONSTRAINT UNIQUE_mRNA_ID UNIQUE (`mRNA_ID`)
                );''')

	c.execute("CREATE INDEX 'key_Gene_Name' ON mRNA (Gene_Name);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(sys.argv[1])])        

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE mRNA")
	create_table()
	call(['touch', "{}".format(sys.argv[1])])

