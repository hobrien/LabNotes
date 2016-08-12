DROP TABLE IF EXISTS `GencodeGTF`;
CREATE TABLE `GencodeGTF` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `seqid` varchar(255) DEFAULT NULL,
  `source` varchar(255) DEFAULT NULL,
  `feature` varchar(255) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `score` float DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `frame` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`)
) ENGINE=InnoDB AUTO_INCREMENT=1 DEFAULT CHARSET=utf8;

LOAD DATA INFILE '/Users/heo3/BTSync/FetalRNAseq/Reference/genes.csv'
INTO TABLE GencodeGTF  
    FIELDS TERMINATED BY ',' 
           OPTIONALLY ENCLOSED BY '"'
    LINES  TERMINATED BY '\n' -- or \r\n
    (id,
     seqid, 
     source, 
     feature, 
     start, 
     end, 
     score, 
     strand, 
     frame
     );
ALTER TABLE `GencodeGTF` ADD INDEX (`seqid`, `start`, `end`);
ALTER TABLE `GencodeGTF` ADD INDEX (`feature`);


DROP TABLE IF EXISTS `GencodeFeatures`;
CREATE TABLE `GencodeFeatures` (
  `id` int(11) unsigned NOT NULL,
  `feature` varchar(255) NOT NULL,
  `value` varchar(255) DEFAULT NULL,
  PRIMARY KEY (`id`, `feature`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;

LOAD DATA INFILE '/Users/heo3/BTSync/FetalRNAseq/Reference/features.csv'
INTO TABLE `GencodeFeatures`  
    FIELDS TERMINATED BY ',' 
           OPTIONALLY ENCLOSED BY '"'
    LINES  TERMINATED BY '\n' -- or \r\n
    (`id`,
     `feature`, 
     `value`
     );
ALTER TABLE `GencodeFeatures` ADD INDEX (`value`);
ALTER TABLE `GencodeFeatures` ADD INDEX (`feature`);


