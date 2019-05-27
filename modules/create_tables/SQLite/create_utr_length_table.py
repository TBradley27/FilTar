#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='utr_length';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `utr_length` (
	`mrna_id` VARCHAR(50) NOT NULL,
	id INTEGER PRIMARY KEY,
	`tissue_id` VARCHAR(50) DEFAULT NULL,
	`utr_length` INT(11) DEFAULT NULL,
	FOREIGN KEY(mrna_id) REFERENCES mRNA (mRNA_ID),
	FOREIGN KEY(tissue_id) REFERENCES Species (taxonomic_ID)
	)''')

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE utr_length")
	create_table()
