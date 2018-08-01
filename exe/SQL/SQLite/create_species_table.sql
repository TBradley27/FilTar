CREATE TABLE species (
  taxonomic_ID VARCHAR(15) DEFAULT NULL,
  species_name CHAR(20) DEFAULT NULL,
  genome_build CHAR(20) DEFAULT NULL,
  common_name VARCHAR(50) DEFAULT NULL,
  id INTEGER PRIMARY KEY,
  CONSTRAINT UNIQUE_taxonomic_ID UNIQUE (`taxonomic_ID`),
  CONSTRAINT UNIQUE_common_name UNIQUE (`common_name`),
  CONSTRAINT UNIQUE_species_name UNIQUE (`species_name`)
);

INSERT INTO species (taxonomic_ID,species_name,genome_build,common_name) VALUES (9606,'Homo sapiens','hg38','Human');
INSERT INTO species (taxonomic_ID,species_name,genome_build,common_name) VALUES (10090,'Mus musculus','mm10','Mouse');
