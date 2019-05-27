#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='miRNA';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE miRNA (
	name VARCHAR(20) DEFAULT NULL,
	for_inline_id INT(11) DEFAULT NULL,
	test_id INT(11) DEFAULT NULL,
	species VARCHAR(20) DEFAULT NULL,
	id INTEGER PRIMARY KEY,
	FOREIGN KEY(for_inline_id) REFERENCES miRNA (id),
	FOREIGN KEY(test_id) REFERENCES miRNA (id),
	FOREIGN KEY(species) REFERENCES species (taxonomic_ID)
	);''')

	c.execute("CREATE INDEX 'key_for_inline_id' ON miRNA (for_inline_id);")
	c.execute("CREATE INDEX 'key_test_id' ON miRNA (test_id);")
	c.execute("CREATE INDEX 'key_species_mirna' ON miRNA (species);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(sys.argv[1])])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE miRNA")
	create_table()
