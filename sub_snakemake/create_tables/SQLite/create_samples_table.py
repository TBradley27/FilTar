#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Samples';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `Samples` (
	`name` VARCHAR(50) NOT NULL,
	id INTEGER PRIMARY KEY,
	`tissue_id` VARCHAR(20) DEFAULT NULL,
	`species_id` VARCHAR(20) DEFAULT NULL,
	CONSTRAINT UNIQUE_name UNIQUE (`name`),
	FOREIGN KEY(tissue_id) REFERENCES Tissues (name),
	FOREIGN KEY(species_id) REFERENCES species (taxonomic_ID)
	)''')

	c.execute("CREATE INDEX 'key_Samples_tissues' ON Samples (tissue_id);")
	c.execute("CREATE INDEX 'key_Samples_species' ON Samples (species_id);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE Samples")
	create_table()
