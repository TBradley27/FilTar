CREATE TABLE miRanda (
  mirna_id VARCHAR(15) DEFAULT NULL,
  mrna_id VARCHAR(25) DEFAULT NULL,
  Species VARCHAR(5) DEFAULT NULL,
  UTR_START SMALLINT(6) DEFAULT NULL,
  UTR_END SMALLINT(6) DEFAULT NULL,
  score DECIMAL(5,2) DEFAULT NULL,
  id INTEGER PRIMARY KEY,
  CONSTRAINT UNIQUE_KEY UNIQUE (`mirna_id`,`mrna_id`,`UTR_START`,`UTR_END`),
  FOREIGN KEY(Species) REFERENCES species (taxonomic_id)
);

CREATE INDEX 'key_Species' ON miRanda (Species);
CREATE INDEX 'key_mRNA' ON miRanda (mrna_id);
CREATE INDEX 'key_mRNA_species' ON miRanda (Species, mRNA_id);
