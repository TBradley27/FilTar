#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='contextpp';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `contextpp` (
	id INTEGER PRIMARY KEY,
	mirna_id VARCHAR(15) DEFAULT NULL,
	mrna_id VARCHAR(20) DEFAULT NULL,
        Species VARCHAR(20) DEFAULT NULL,
	UTR_START SMALLINT(6) DEFAULT NULL,
	UTR_END SMALLINT(6) DEFAULT NULL,
	Site_Type VARCHAR(7) DEFAULT NULL,
	`score` decimal(4,3) DEFAULT NULL,
	`weighted_score` decimal(4,3) DEFAULT NULL,
	CONSTRAINT UNIQUE_site UNIQUE (`mirna_id`,`mrna_id`,`UTR_START`,`UTR_END`),
	FOREIGN KEY(Species) REFERENCES species (taxonomic_ID)
	)''')

	c.execute("CREATE INDEX 'key_contextpp_mrna' ON contextpp (mrna_id);")
	c.execute("CREATE INDEX 'key_contextpp_species' ON contextpp (Species);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE contextpp")
	create_table()
