#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='expression_profiles';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `expression_profiles` (
	`TPM` decimal(10,2) DEFAULT NULL,
	`experiments_id` varchar(20) DEFAULT NULL,
	`mrnas_id` varchar(20) NOT NULL,
	id INTEGER PRIMARY KEY,
	FOREIGN KEY(mrnas_id) REFERENCES mRNA (mRNA_ID)
	)''')

	c.execute("CREATE INDEX 'key_mRNA_expression_profiles' ON expression_profiles (mrnas_id);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(sys.argv[1])])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE expression_profiles")
	create_table()
