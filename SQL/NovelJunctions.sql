DROP TABLE IF EXISTS NovelJunctions;
CREATE TABLE `NovelJunctions` (
  `id` int(11) unsigned NOT NULL AUTO_INCREMENT,
  `junc_id` int(11) DEFAULT NULL,
  `sample_id` varchar(255) DEFAULT NULL,
  `start` int(11) DEFAULT NULL,
  `end` int(11) DEFAULT NULL,
  `strand` char(1) DEFAULT NULL,
  `chr` varchar(255) DEFAULT NULL,
  `status` varchar(255) DEFAULT NULL,
  `score` int(11) DEFAULT NULL,
  `int_start` int(11) DEFAULT NULL,
  `int_end` int(11) DEFAULT NULL,
  PRIMARY KEY (`id`),
  INDEX position (`chr`,`int_start`,`int_end`),
  UNIQUE KEY `junction` (`junc_id`,`sample_id`)
) ENGINE=InnoDB DEFAULT CHARSET=utf8;