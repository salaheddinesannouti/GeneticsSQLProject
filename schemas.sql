CREATE DATABASE IF NOT EXISTS CancerBioInfo;

USE CancerBioInfo;

CREATE TABLE Gene(
gene_symbol VARCHAR(20) PRIMARY KEY,
gene_name VARCHAR(100) NOT NULL, 
chromosome VARCHAR(5) NOT NULL,
chr_start_pos INT NOT NULL,
chr_end_pos INT NOT NULL,
strand CHAR(1) NOT NULL,
CONSTRAINT check_strand CHECK (strand IN ('+', '-'))
);

CREATE TABLE Protein(
accession_number VARCHAR(50) PRIMARY KEY,
gene_symbol VARCHAR(20),
chain_length INT NOT NULL,
molecular_mass INT NOT NULL,
protein_function TEXT,
FOREIGN KEY (gene_symbol) REFERENCES Gene(gene_symbol)
);


CREATE TABLE Domain(
    domain_name VARCHAR(50) PRIMARY KEY,
    family VARCHAR(50),
    domain_description VARCHAR(200)
);

CREATE TABLE ProteinContainsDomain(
accession_number VARCHAR(50),
domain_name VARCHAR(50),
PRIMARY KEY (accession_number, domain_name),
FOREIGN KEY (accession_number) REFERENCES Protein(accession_number),
FOREIGN KEY (domain_name) REFERENCES Domain(domain_name)
);

