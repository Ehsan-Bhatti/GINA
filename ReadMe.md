# GINA 
GEO Integrated Network Analysis

## Requirements 
GINA currently makes use of the following Python Packages:
- pandas 
- numpy 
- scipy
- scikit-learn 
- drugstone
- goatools
- mygene
- rpy2
- setuptools-git 

For differential expression analysis and GEO data acquisition GINA uses the following R packages: 
- GEOquery
- limma 
- openxlsx
- edgeR



## Installation

To install use the following command 
```bash 
pip install git+https://github.com/Ehsan-Bhatti/GINA
```

If your system meets the following requirements download the DOMINO pacakge from https://github.com/Shamir-Lab/DOMINO 

Otherwise use their web server for active module identification available from https://domino.cs.tau.ac.il/ 

## Documentation
Read the docs https://gina.readthedocs.io/en/latest/

## Usage 
Examples of use cases using this package.

The APID network can be downloaded from http://apid.dep.usal.es:8080/APID/init.action


Differential Expression Analysis: 
``` Python 
from gina import differnetial_expression_analysis as dea

# Download GSE14520 to a directory called "data/GSE14520"
dea.get_data_from_geo('GSE14520', directory_to_save_file='data/GSE14520')


# Perform differential expression analysis on GPL571 dataset of GSE14520
gpl571 = dea.read_geo_file('data/GSE14520/GSE14520-GPL571_series_matrix.txt')
pheonotype_data = dea.phenotype_data_from_geo_data(gpl571)
dea.print_all_columns(phenotype_data)

# We can see all columns present in the phenotype data. 
# The phenotype column is Tissue:ch1

# See phenotypes present in data 
dea.print_unique_values_in_geo_column(phenotype_data, 'Tissue:ch1')

# Create label groups based on phenotypes present 
label_groups = [['Liver Tumor Tissue'], ['Liver Non-Tumor Tissue', 'Liver tissue of six healthy donors', 'Liver tissue of six healthy liver donors']] 

# Run differential expression analysis
dea.differential_expression_analysis(geo_data=gpl571, phenotype_col_name='Tissue:ch1', list_of_phenotype_groups=label_groups, save_file_name='data/GSE14520/gpl571_analysis.xlsx')

```

Active Module Identification if the DOMINO package has been downloaded: 
``` Python
from gina import active_module_identification as ami
from gina import differnetial_expression_analysis as dea 
import pandas as pd 

# Get the differentially expressed genes from our prior analysis
feature_df = pd.read_excel('data/GSE14520/gpl571_analysis.xlsx', sheet_name='feature_df', index_col=0)
diff_dict = dea.list_of_differentially_expressed_genes(feature_df, gene_name_col='Gene_Symbol', differential_analysis_col='group_1-group_2', sep_rows_string='///')
list_of_genes = diff_dict.get('list')

# Using the APID dataframe of human protein-protein interactions as our base network
# Create files for DOMINO
APID_df = pd.read_csv('9606.mitab', sep='\t')
ami.create_files_for_domino_input(network_df=APID_df, network_sif_file_name='data/network/APID.sif', network_node_1_col='interactor_a', network_node_2_col='interactor_b', nodes_of_interest_list=list_of_genes, nodes_of_interest_txt_file_name='data/GSE14520/gpl571_noi.txt')
ami.create_slices_file(network_sif_file_path='data/network/APID.sif', output_file_path='data/network/APID.slices.txt')

# Run DOMINO
ami.run_domino(network_sif_file_path='data/network/APID.sif', slices_file_path='data/network/APID.slices.txt', active_genes_file_path='data/GSE14520/gpl571_noi.txt', output_folder_path='data/GSE14520/gpl571_domino') 

```

Gene Ontology Network Analysis 

``` Python 
from gina import go_analysis as ga 
import pandas as pd

# Download OBO and GAF files
ga.download_basic_go_obo()
ga.download_gaf(species='goa_human')

# Get all genes present in used network 
network_df = pd.read_csv('data/network/APID.sif', sep='\t')
all_genes = ga.all_nodes_from_df(network_df, list_of_cols=['interactor_a', 'interactor_b']

# GO Analysis
ga.go_enrichment_analysis_modules(path_to_obo_file='data/GO/go-basic.obo', path_to_gaf_file='data/GO/go_human.gaf', output_directory='data/GSE14520/gpl571_domino/GO_analysis', path_to_modules_file='data/GSE14520/gpl571_domino/modules.out', all_genes_in_network=all_genes)   

```

PoLoBag Analysis and drug search of active genes:

``` Python 
from gina import polobag_and_drug_search as pbds
import pandas as pd 

# Get feature and expression data 
path_to_xl_file = 'data/GSE14520/gpl571_analysis.xlsx'

feature_df = pd.read_excel(path_to_xl_file, sheet_name='feature_df', index_col=0)
expression_df = pd.read_excel(path_to_xl_file, sheet_name='expression_df', index_col=0)

# Run PoLoBag on modules 
pbds.run_polobag_on_modules(modules='data/GSE14520/gpl571_domino/modules.out', expression_df=expression_df, output_directory='data/GSE14520/gpl571_polobag', feature_df=feature_df, feature_gene_id_col='Gene Symbol', feature_gene_sep_string='///', n22=2)

# Analysis of PoLoBag files and drug search the active genes in the files
apid_network = pd.read_csv('data/network/APID.sif', sep='\t')

pbds.polobag_analysis_multiple_files(path_to_polobag_output_directory='data/GSE14520/gpl571_polobag', network_df=apid_network, active_genes='data/GSE14520/gpl571_noi.txt', analysis_output_directory='data/GSE14520/gpl571_polobag_analysis', active_gene_drug_search_output_directory='data/GSE14520/gpl571_polobag_drug_search') 

```

## Creation
GINA was created on MacOS (Sonoma) and Ubuntu 22.04 (VM for Windows).

  
## Future Work
Other examples will be added. 
