#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Gene_species';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `Gene_species` (
	`gene_id` VARCHAR(50) NOT NULL,
	`species_id` VARCHAR(20) DEFAULT NULL,
	PRIMARY KEY (`gene_id`,`species_id`)
	FOREIGN KEY(gene_id) REFERENCES Gene (name),
	FOREIGN KEY(species_id) REFERENCES species (taxonomic_ID)
	)''')

	c.execute("CREATE INDEX 'key_species_Gene_species' ON Gene_species (species_id);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE Gene_species")
	create_table()
