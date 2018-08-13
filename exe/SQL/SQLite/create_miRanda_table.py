#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='miRanda';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE miRanda (
	mirna_id VARCHAR(15) DEFAULT NULL,
	mrna_id VARCHAR(25) DEFAULT NULL,
	Species VARCHAR(5) DEFAULT NULL,
	UTR_START SMALLINT(6) DEFAULT NULL,
	UTR_END SMALLINT(6) DEFAULT NULL,
	score DECIMAL(5,2) DEFAULT NULL,
	id INTEGER PRIMARY KEY,
	CONSTRAINT UNIQUE_KEY UNIQUE (`mirna_id`,`mrna_id`,`UTR_START`,`UTR_END`),
	FOREIGN KEY(Species) REFERENCES species (taxonomic_id)
	);''')

	c.execute("CREATE INDEX 'key_Species' ON miRanda (Species);")
	c.execute("CREATE INDEX 'key_mRNA' ON miRanda (mrna_id);")
	c.execute("CREATE INDEX 'key_mRNA_species' ON miRanda (Species, mRNA_id);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(sys.argv[1])])        

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE miRanda")
	create_table()
	call(['touch', "{}".format(sys.argv[1])])

