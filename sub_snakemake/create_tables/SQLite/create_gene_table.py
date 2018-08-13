#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Gene';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `Gene` (
	`name` VARCHAR(50) NOT NULL,
	id INTEGER PRIMARY KEY,
	`test_id` INT(11) DEFAULT NULL,
	`for_inline_id` INT(11) DEFAULT NULL,
	CONSTRAINT UNIQUE_name UNIQUE (`name`),
	FOREIGN KEY(for_inline_id) REFERENCES Gene (id),
	FOREIGN KEY(test_id) REFERENCES Gene (id)
	)''')

	c.execute("CREATE INDEX 'key_for_inline_id_gene' ON Gene (for_inline_id);")
	c.execute("CREATE INDEX 'key_test_id_gene' ON Gene (test_id);")

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE Gene")
	create_table()
