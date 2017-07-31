
use mysql;
drop table annotation;

CREATE TABLE annotation ( ann_id INT(10) PRIMARY KEY,
   gene VARCHAR(20),
   genomic VARCHAR(20),
   id INT(10),
   protein VARCHAR(20),
   corrected VARCHAR(20),
   eng VARCHAR(20),
   cosm_base VARCHAR(20),
   COSMIC VARCHAR(20),
   count INT(10)
);


LOAD DATA LOCAL INFILE  
'/home/zjones/COSMIC/mutations_2.csv'
INTO TABLE annotation
FIELDS TERMINATED BY '\t' 
LINES TERMINATED BY '\n'
IGNORE 1 LINES
(ann_id,gene,genomic,id,protein,corrected,eng,cosm_base,COSMIC,count);
