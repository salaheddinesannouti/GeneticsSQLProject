USE cancerbioinfo;

# 1. Donner les 3 protéines avec les chaines d'acide aminé les plus longues
SELECT * FROM protein
ORDER BY chain_length DESC
LIMIT 3;

# 2. Donner le nombre de gènes par chromosome, en classant par nombre de gène  décroissant
SELECT chromosome, count(*) AS gene_count
FROM Gene 
GROUP BY chromosome
ORDER BY gene_count DESC;


# 3. Donner le nombre de protéines qui ont une chaine de longueur entre 300 et 500 aminés
# et dont les gènes codants sont situés sur l'un de ces chromosomes: 18, 19, 20, X ou Y
SELECT count(*) AS protein_number
FROM Protein 
JOIN Gene ON Gene.gene_symbol=Protein.gene_symbol 
WHERE chain_length BETWEEN 300 AND 500 AND chromosome IN ('18', '19','20','X','Y');

# 4. Donner les 10 gènes avec le plus de domaines associés 
SELECT Gene.gene_symbol, count(*) AS domain_count
FROM ProteinContainsDomain pcd
JOIN Protein p
ON pcd.accession_number = p.accession_number
JOIN Gene
ON Gene.gene_symbol = p.gene_symbol
GROUP BY Gene.gene_symbol
ORDER BY domain_count DESC
LIMIT 10;




# 5. Identifier les protéines anormalement longues 
# mais qui ont moins de deux domaines
SELECT p.accession_number, p.chain_length, COUNT(pcd.domain_name) AS domain_count
FROM Protein p
JOIN ProteinContainsDomain pcd ON p.accession_number = pcd.accession_number
GROUP BY p.accession_number, p.chain_length
HAVING p.chain_length > 1200 AND domain_count <= 2;


# 6. Donner la liste des domaines présents dans au moins dix protéines différentes.
SELECT Domain.domain_name
FROM Domain 
JOIN ProteinContainsDomain pcd
  ON Domain.domain_name = pcd.domain_name
GROUP BY Domain.domain_name
HAVING COUNT(DISTINCT pcd.accession_number) >= 10;

# 7. Lister les domaines présents dans des protéines de plus de quatre chromosomes et les chromosomes sur lesquels ils sont présent.
SELECT pcd.domain_name, GROUP_CONCAT(DISTINCT g.chromosome) AS chromosomes
FROM ProteinContainsDomain pcd
JOIN Protein p ON pcd.accession_number = p.accession_number
JOIN Gene g ON p.gene_symbol = g.gene_symbol
GROUP BY pcd.domain_name
HAVING COUNT(DISTINCT g.chromosome) > 4;

# 8. Identifier les gènes qui ont des domaines spécifiques(produits uniquement par ces gènes)
WITH  domains_count AS(
SELECT domain_name, count(DISTINCT accession_number) AS genes_number
FROM ProteinContainsDomain pcd
GROUP BY domain_name
HAVING genes_number = 1
)
SELECT Protein.gene_symbol, domains_count.domain_name FROM domains_count
JOIN ProteinContainsDomain pcd 
ON domains_count.domain_name = pcd.domain_name
JOIN Protein
ON Protein.accession_number = pcd.accession_number
;

# 9.Trouver le chromsome "le plus lourd". C'est-à-dire le chromosome 
# dont la masse moléculaire totale des protéines qu'il produit est la plus haute
SELECT chromosome, total_mass
FROM (
    SELECT 
        g.chromosome, 
        SUM(p.molecular_mass) AS total_mass
    FROM Gene g
    JOIN Protein p ON g.gene_symbol = p.gene_symbol
    GROUP BY g.chromosome
) AS ChromosomeAggregates
ORDER BY total_mass DESC
LIMIT 1; 

# 10. Vérifier si le plus long gène code pour la plus grande protéine
WITH longest_gene AS (
    SELECT gene_symbol, (chr_end_pos - chr_start_pos) AS gene_length
    FROM Gene
    ORDER BY gene_length DESC
    LIMIT 1
),
longest_protein AS (
    SELECT gene_symbol, accession_number, chain_length
    FROM Protein
    ORDER BY chain_length DESC
    LIMIT 1
)
SELECT 
    CASE 
        WHEN lg.gene_symbol = lp.gene_symbol THEN
            CONCAT('Le gène le plus long', lg.gene_symbol, 'de longueur', lg.gene_length,
            'pour la protéine la plus longue.')
        ELSE
            CONCAT('Le gène le plus long ', lg.gene_symbol, 'de longueur', 
            lg.gene_length, 'nt ne code pas pour la protéine', lp.accession_number, 
            'la plus longue(', lp.chain_length, 'aa).')
    END AS result
FROM longest_gene lg
CROSS JOIN longest_protein lp;

SELECT count(*) FROM proteincontainsdomain;