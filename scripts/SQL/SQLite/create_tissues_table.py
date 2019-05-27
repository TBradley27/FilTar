#!/bin/env python3

import sqlite3
import sys
import os
from subprocess import call

conn = sqlite3.connect('filtar.db')
c = conn.cursor()

c.execute("SELECT name FROM sqlite_master WHERE type='table' AND name='Tissues';")

test_object = c.fetchone()

def create_table():
	c.execute('''CREATE TABLE `Tissues` (
	`name` VARCHAR(50) NOT NULL,
	id INTEGER PRIMARY KEY,
	`taxonomic_ID` VARCHAR(20) DEFAULT NULL,
	`test_id` INT(11) DEFAULT NULL,
	`for_inline_id` INT(11) DEFAULT NULL,
	CONSTRAINT UNIQUE_name_taxonomic_ID UNIQUE (`name`,`taxonomic_ID`)
	)''')

	conn.commit()
	conn.close()
	call(['touch', "{}".format(sys.argv[1])])

	return()

if test_object == None:
	create_table()
else:
	c.execute("DROP TABLE Tissues")
	create_table()
