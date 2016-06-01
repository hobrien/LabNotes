DROP TABLE IF EXISTS `CufflinksGTF`;
CREATE TABLE `CufflinksGTF` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `frac` float DEFAULT NULL,
  `full_read_support` varchar(255) DEFAULT NULL,
  `cov` float DEFAULT NULL,
  `frame` int(11) DEFAULT NULL,
  `conf_lo` float DEFAULT NULL,
  `FPKM` float DEFAULT NULL,
  `feature` varchar(255) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `transcript_id` varchar(255) DEFAULT NULL,
  `gene_id` varchar(255) DEFAULT NULL,
  `conf_hi` float DEFAULT NULL,
  `seqid` varchar(255) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `score` float DEFAULT NULL,
  `strand` int(11) DEFAULT NULL,
  `source` varchar(255) DEFAULT NULL,
  `sample_id` varchar(255) DEFAULT NULL,
  `exon_number` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;

LOAD DATA INFILE '/Users/heo3/BTSync/FetalRNAseq/Cufflinks/AllGTF.csv'
INTO TABLE CufflinksGTF  
    FIELDS TERMINATED BY ',' 
           OPTIONALLY ENCLOSED BY '"'
    LINES  TERMINATED BY '\n' -- or \r\n
    (frac, 
     full_read_support, 
     exon_number, 
     cov, 
     conf_lo, 
     FPKM, 
     feature, 
     start, 
     transcript_id, 
     gene_id, 
     conf_hi, 
     seqid, 
     end, 
     score, 
     strand, 
     source, 
     sample_id
     );

ALTER TABLE `CufflinksGTF` ADD INDEX (`sample_id`);
ALTER TABLE `CufflinksGTF` ADD INDEX (`transcript_id`);
ALTER TABLE `CufflinksGTF` ADD INDEX (`gene_id`);
ALTER TABLE `CufflinksGTF` ADD INDEX (`seqid`);
ALTER TABLE `CufflinksGTF` ADD INDEX (`start`);
ALTER TABLE `CufflinksGTF` ADD INDEX (`end`)   