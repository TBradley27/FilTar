CREATE TABLE `miRanda` (
  `mirna_id` varchar(15) DEFAULT NULL,
  `mrna_id` varchar(25) DEFAULT NULL,
  `Species` varchar(5) DEFAULT NULL,
  `UTR_START` smallint(6) DEFAULT NULL,
  `UTR_END` smallint(6) DEFAULT NULL,
  `score` decimal(5,2) DEFAULT NULL,
  `id` int(11) NOT NULL AUTO_INCREMENT,
  PRIMARY KEY (`id`),
  UNIQUE KEY `unique_index` (`mirna_id`,`mrna_id`,`UTR_START`,`UTR_END`),
  KEY `fk_Species` (`Species`),
  KEY `fk_miRanda_mrna` (`mrna_id`),
  KEY `mrna_species_idx` (`mrna_id`,`Species`),
  CONSTRAINT `fk_Species` FOREIGN KEY (`Species`) REFERENCES `species` (`taxonomic_ID`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1
