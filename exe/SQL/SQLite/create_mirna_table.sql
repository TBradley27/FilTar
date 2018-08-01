CREATE TABLE miRNA (
  name VARCHAR(20) DEFAULT NULL,
  for_inline_id INT(11) DEFAULT NULL,
  test_id INT(11) DEFAULT NULL,
  species VARCHAR(20) DEFAULT NULL,
  id INTEGER PRIMARY KEY,
  FOREIGN KEY(for_inline_id) REFERENCES miRNA (id),
  FOREIGN KEY(test_id) REFERENCES miRNA (id),
  FOREIGN KEY(species) REFERENCES species (taxonomic_ID)
);

CREATE INDEX 'key_for_inline_id' ON miRNA (for_inline_id);
CREATE INDEX 'key_test_id' ON miRNA (test_id);
CREATE INDEX 'key_species_mirna' ON miRNA (species);

