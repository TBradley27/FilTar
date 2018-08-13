#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

# I need to find a way of separating out the SQL auery from the boilerplate python code

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Runs';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `Runs` (
	`name` VARCHAR(50) NOT NULL,
	id INTEGER PRIMARY KEY,
	`sample` VARCHAR(20) DEFAULT NULL,
	CONSTRAINT UNIQUE_name UNIQUE (`name`),
	FOREIGN KEY(sample) REFERENCES Samples (name)
	)''')

	conn.commit()
	conn.close()
	call(['touch', "{}".format(snakemake.output)])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE Runs")
	create_table()
