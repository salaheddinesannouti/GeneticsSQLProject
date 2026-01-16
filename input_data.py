import re

import pandas as pd
import numpy as np



PATTERN = re.compile(r"(#=GF)\s(ID|TP|DE|AC)\s+(.+)")

all_domains = []
temp_data = {}

with open("data/domains_data.dat") as domain_file:
    for i, line in enumerate(domain_file):

        # This header indicates the start of the data about a domain 
        if line.startswith("# STOCKHOLM 1.0"):

            # If temp_data is full, we need to add its data to the df 
            # and reset it
            if temp_data != {}:
                all_domains.append(temp_data)
                temp_data = {}
            continue

        # In this case, we are dealing with a line that has data the domaine
        else:
            data = PATTERN.match(line)
            if data is not None:
                # Get type of info about d
                data_type = data.group(2)
                data_entry = data.group(3)
                
                # Decipher what data the line is about
                match data_type:
                    case "AC":
                        stable_version = data_entry.split(".")[0]
                        temp_data["domaine_id"] = stable_version
                    case "ID":
                        temp_data["domain_name"] = data_entry
                    case "DE":
                        temp_data["domain_description"] = data_entry
                    case "TP":
                        temp_data["family"] = data_entry
                        
    
df_domains = pd.DataFrame(all_domains)

genes = pd.read_csv("data/gene_data.txt")
proteins = pd.read_csv("data/proteins_data.tsv", sep="\t")

# Reformat the gene description so that it only has the long name of the gene
genes["Gene description"] = genes["Gene description"].str.split("[").str[0].str.strip()

# Rename the gene dataset so that it matches the column names of our project
genes.rename(columns={"Gene description": "gene_name",
                      "Gene name": "gene_symbol",
                      "Chromosome/scaffold name": "chromosome",
                      "Gene start (bp)": "chr_start_pos",
                      "Gene end (bp)": "chr_end_pos",
                      "Strand": "strand"}, inplace=True)

# Format the strand on the required format
genes["strand"]  = np.where(genes["strand"] == 1, "+", "-")


# Rename the protein dataset so that it matches the column names of our project
proteins.rename(columns={"Entry": "accession_number",
                 "Gene Names (primary)": "gene_symbol",
                 "Function [CC]": "protein_function",
                 "Length": "chain_length",
                 "Mass": "molecular_mass"}, inplace=True)


# Keep only proteins which are coded by genes in our cancer dataset 
# and which have known domains
cancer_prot_filter = ~(proteins["Pfam"].isna()) & (proteins["gene_symbol"].isin(genes["gene_symbol"]))
cancer_prot_with_domains = proteins[cancer_prot_filter]


# Map all proteins to their domains
protein_to_domain_map = (
    cancer_prot_with_domains[["accession_number", "Pfam"]]
    .assign(Pfam=lambda df: df["Pfam"].str.strip(";").str.split(";"))
    .explode("Pfam")
)

print(protein_to_domain_map)

# Keep all used domains
used_domains = protein_to_domain_map["Pfam"].unique()
df_domains = df_domains[df_domains["domaine_id"].isin(used_domains)]

# Map protein to domain name
protein_to_domain = protein_to_domain_map.merge(
    df_domains[["domaine_id", "domain_name"]],
    left_on="Pfam",
    right_on="domaine_id",
    how="inner"
)[["accession_number", "domain_name"]]


# Set family only for Pfams that are of type Family AND are used
df_domains["family"] = np.where(
    (df_domains["family"] == "Family"),
    df_domains["domain_name"],
    np.nan
)



# Drop irrelevant columns for our assignement
df_domains.drop(columns="domaine_id", inplace=True)
cancer_prot_with_domains.drop(columns=["Organism", "Entry Name", "InterPro", "Protein names", "Pfam"], inplace=True)

# Save to csv files
df_domains.to_csv("new_data/domains.csv", index=False)
cancer_prot_with_domains.to_csv("new_data/proteins.csv", index= False)
genes.to_csv("new_data/genes.csv", index= False)
protein_to_domain.to_csv("new_data/protein_domain_mappings.csv", index=False)
