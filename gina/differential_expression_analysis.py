import pandas as pd
import os
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage as STAP
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


r_file_path = os.path.dirname(__file__)
r_file_path = os.path.join(r_file_path, 'util_files/diff_analysis.R')

with open(r_file_path, 'r') as f:
    rf = f.read()

r_func = STAP(rf, 'r_func')

decide_adj_test_algorithms = ['none', 'BH', 'fdr', 'BY', 'holm']
normalising_method_algorithms = ['TMM', 'RLE', 'upperquartile', 'none']


def helper_convert_r_to_python(li):
    """
    Helper function to convert between R and Python

    :param li: R list (dictionary) to be converted.
    :return: A python dictionary of the converted items.
    """

    output = dict(zip(li.names, list(li)))
    with (ro.default_converter + pandas2ri.converter).context():
        return {k: ro.conversion.get_conversion().rpy2py(v) for k, v in output.items()}


def get_data_from_geo(data_id, directory_to_save_file=None):
    """
    Retrieves data from the Gene Expression Omnibus using the accession ID

    :param data_id: A string, the accession ID of the desired dataset e.g. 'GSE1234'
    :param directory_to_save_file: Path to the location to save the downloaded data. Defaults to the current working directory if no path is supplied.
    :return: Downloads data nad saves it to the directory.
    """

    if directory_to_save_file is None:
        directory_to_save_file = os.getcwd()

    r_func.get_data_from_geo(data_id, directory_to_save_file)


def read_geo_file(path_to_file):
    """
    Load a GEO dataset and reads it into a dataframe
    :param path_to_file: String path to the dataset
    :return A dataframe of the GEO data
    """
    return r_func.read_geo_file(path_to_file)


def phenotype_data_from_geo_data(geo_data, save_file_name=None):
    """
    Retrieve the phenotype data from a GEO dataframe. Useful for finding the col labelling samples and the unique labels
    :param geo_data: The GEO dataframe most likely accessed using the read_geo_file method
    :param save_file_name: File name to save the phenotype data as a csv file. Must end in '.csv'
    :return: A dataframe of the phenotypic data
    """
    if save_file_name is None:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.pheno_data(geo_data))

    else:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.pheno_data(geo_data, save_file_name))


def expression_data_from_geo_data(geo_data, save_file_name=None):
    """
    Retrieve the expression data from a GEO dataframe
    :param geo_data: The GEO dataframe most likely accessed using the read_geo_file method
    :param save_file_name: File name to save the expression data as a csv file. Must end in '.csv'
    :return:  A dataframe of the expression data
    """
    if save_file_name is None:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.expression_data(geo_data))
    else:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.expression_data(geo_data, save_file_name))


def feature_data_from_geo_data(geo_data, save_file_name=None):
    """
    Retrieve the feature data from a GEO dataframe. Access the feature of the dataset seeing information about each gene recorded in the dataset e.g. 'Gene Symbols' and ''

    :param geo_data: The GEO dataframe most likely accessed using the read_geo_file method
    :param save_file_name: File name to save the feature data as a csv file. Must end in '.csv'
    :return: A dataframe of feature data
    """
    if save_file_name is None:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.feature_data(geo_data))
    else:
        with (ro.default_converter + pandas2ri.converter).context():
            return ro.conversion.get_conversion().rpy2py(r_func.feature_data(geo_data, save_file_name))


def print_all_columns(df):
    """
    Print all columns present in a GEO dataframe. Useful for finding the label col in the Phenotype DF

    :param df: A GEO dataframe e.g. from the read_geo_file method or expression_data_from_geo_data method
    :return: Prints all columns
    """
    print('Columns:')
    if isinstance(df, pd.DataFrame):
        cols = df.columns
        for i in cols:
            print(i)
    else:
        print(r_func.see_columns(df))


def print_unique_values_in_geo_column(df, col_name):
    """
    See the  unique values in a GEO column. Useful for seeing all phenotypes recorded

    :param df: A dataframe of GEO data e.g. result of get_data_from_geo method
    :param col_name: Name of the column we wish to see the values of
    :return: Prints all unique values
    """
    if isinstance(df, pd.DataFrame):
        print(f"Unique values in: {col_name}")
        print(df[col_name].unique())

    else:
        r_func.unique_values_in_col(df, col_name)


def differential_expression_analysis(geo_data, phenotype_col_name,
                                     list_of_phenotype_groups, excel_file_name,
                                     e_bayes_proportion=0.01, decide_test_adj_method='BH',
                                     decide_test_p_val=0.05):
    """
    Run differential expression analysis

    :param geo_data: GEO data can be obtained using the read_geo_file method
    :param phenotype_col_name: String name of column that contains the phenotypes
    :param list_of_phenotype_groups: List of lists separating phenotypic labels into group e.g. [['Cancer Tissue'], ['Non-Cancer Tissue'], ['Control', 'Healthy']]
    :param excel_file_name: String file name to save the results should end with the extension '.xlsx'
    :param e_bayes_proportion: Float representing assumed percentage of differentially expressed genes. Defaults to 0.01
    :param decide_test_adj_method: String specifying the P-Value adjustment method
    :param decide_test_p_val: Float of the desired P-Value of the test
    :return: A dictionary containing the results of the differential expression analysis and if selected the saved Excel file
    :key 'contrast_fit': result of contrast.fit method on a linear model and the different phenotypic groups
    :key 'empirical_bayes_fit': result of the eBayes method on the result of contrast_fit
    :key 'table_empirical_bayes': table of the eBayes method for all genes
    :key 'differentially_expressed_genes': dataframe showing which genes showing whether genes are differentially expressed
    :key 'calculation_info': record of the user supplied parameters for this calculation
    :key 'phenotype_data': phenotype data used in the calculation
    :key 'expression_data': expression data used in the calculation
    :key 'feature_data': feature data used in the calculation
    """
    if decide_test_adj_method not in decide_adj_test_algorithms:
        raise Exception(f'Invalid "decide_test_adj_method" arg ({decide_test_adj_method})\nUse one of the following: {decide_adj_test_algorithms}')

    output = r_func.differential_expression_analysis(geo_data, phenotype_col_name, list_of_phenotype_groups,
                                                     e_bayes_proportion, decide_test_adj_method,
                                                     decide_test_p_val, excel_file_name)

    return helper_convert_r_to_python(output)


def read_single_rna_seq_file(file_name, file_sep_delimiter='\t', header=True, excel_file_name=None):
    """
    Read a single RNA Seq file

    :param file_name: RNA seq file name
    :param file_sep_delimiter: String delimiter used to separate data in the RNA seq file. Defaults to tab
    :param header: Boolean value whether the file has a header. Defaults to True
    :param excel_file_name: String Excel file name to save the RNA seq file must end in xlsx. Defaults to None and does not save the file
    :return: A dictionary containing 3 items:
    :key 'phenotype': dataframe of the phenotype information.
    :key 'expression': dataframe of the expression information.
    :key 'feature': dataframe of the feature information.
    """
    if excel_file_name is None:
        output = r_func.read_single_rna_seq_file(file_name, file_sep_delimiter, header)
    else:
        output = r_func.read_single_rna_seq_file(file_name, file_sep_delimiter, header, excel_file_name)

    return helper_convert_r_to_python(output)


def read_multiple_rna_seq_files(path_to_files, files=None, file_sep_delimiter='\t', header=True,
                                names_and_count_cols=[1, 2], excel_file_name=None):
    """
    Read multiple RNA Seq files

    :param path_to_files: String the path to directory containing RNA Seq files
    :param files: List of RNA seq files useful if RNA seq files are in directory with other files. Defaults to None and uses all files in the path_to_files directory
    :param file_sep_delimiter: String delimiter used to separate data in the RNA seq file. Defaults to tab
    :param header: Boolean value whether the file has a header. Defaults to True
    :param names_and_count_cols: A list of the columns that are the gene names and counts respectively. Defaults to [1, 2]. Note that in R columns being counting from 1 not 0
    :param excel_file_name: String Excel file name to save the RNA seq file must end in 'xlsx'.  Defaults to None and does not save the file
    :return: A dictionary containing 3 items:
    :key 'phenotype': dataframe of the phenotype information
    :key 'expression': dataframe of the expression information
    :key 'feature': dataframe of the feature information
    """
    if files is None:
        files = os.listdir(path_to_files)

    if excel_file_name is None:
        output = r_func.read_multiple_rna_seq_files(files, path_to_files, names_and_count_cols, file_sep_delimiter,
                                                    header)
    else:
        output = r_func.read_multiple_rna_seq_files(files, path_to_files, names_and_count_cols, file_sep_delimiter,
                                                    header,
                                                    excel_file_name)

    return helper_convert_r_to_python(output)


def rna_seq_differential_expression_analysis(phenotype_data, expression_data, feature_data,
                                             phenotype_col_name, list_of_phenotype_groups, normalising_method='TMM',
                                             e_bayes_proportion=0.01, decide_test_adj_method='BH',
                                             decide_test_p_val=0.05, file_name=None):
    """
    Differential expression analysis on RNA Seq data

    :param phenotype_data: Dataframe of the phenotype data
    :param expression_data: Dataframe of the expression data
    :param feature_data: Dataframe of the feature data
    :param phenotype_col_name: String name of column that contains the phenotypes
    :param list_of_phenotype_groups: List of lists separating phenotypic labels into group e.g. [['Cancer Tissue'], ['Non-Cancer Tissue'], ['Control', 'Healthy']]
    :param normalising_method: Normalisation method used. Defaults to 'TMM'. Alternate choices are 'RLE', 'upperquartile' or 'none'
    :param e_bayes_proportion: Float representing assumed percentage of differentially expressed genes. Defaults to 0.01
    :param decide_test_adj_method: String specifying the P-Value adjustment method
    :param decide_test_p_val: Float of the desired P-Value of the test
    :param file_name:
    :return: A dictionary containing the results of the differential expression analysis and if selected the saved Excel file
    :key 'contrast_fit': result of contrast.fit method on a linear model and the different phenotypic groups
    :key 'empirical_bayes_fit': result of the eBayes method on the result of contrast_fit
    :key 'table_empirical_bayes': table of the eBayes method for all genes
    :key 'differentially_expressed_genes': dataframe showing which genes showing whether genes are differentially expressed
    :key 'calculation_info': record of the user supplied parameters for this calculation
    :key 'phenotype_data': phenotype data used in the calculation
    :key 'expression_data': expression data used in the calculation
    :key 'feature_data': feature data used in the calculation
    :key 'normalised_expression_data': the normalised expression data
    :key 'normalising_expression_weights': the normalising weights
    """
    if normalising_method not in normalising_method_algorithms:
        raise Exception(f'Invalid "normalising_method" arg ({normalising_method})\nUse one of the following: {normalising_method_algorithms}')

    if decide_test_adj_method not in decide_adj_test_algorithms:
        raise Exception(f'Invalid "decide_test_adj_method" arg ({decide_test_adj_method})\nUse one of the following: {decide_adj_test_algorithms}')

    with (ro.default_converter + pandas2ri.converter).context():
        phenotype_data = ro.conversion.get_conversion().py2rpy(phenotype_data)
        feature_data = ro.conversion.get_conversion().py2rpy(feature_data)
        expression_data = ro.conversion.get_conversion().py2rpy(expression_data)

    if file_name is None:
        output = r_func.rna_seq_differential_expression_analysis(phenotype_data, expression_data, feature_data,
                                                                 phenotype_col_name, list_of_phenotype_groups,
                                                                 normalising_method,
                                                                 e_bayes_proportion, decide_test_adj_method,
                                                                 decide_test_p_val)
    else:
        output = r_func.rna_seq_differential_expression_analysis(phenotype_data, expression_data, feature_data,
                                                                 phenotype_col_name, list_of_phenotype_groups,
                                                                 normalising_method,
                                                                 e_bayes_proportion, decide_test_adj_method,
                                                                 decide_test_p_val, file_name)

    return helper_convert_r_to_python(output)


def list_of_differentially_expressed_genes(df, gene_name_col, differential_analysis_col, sep_rows_string=None,
                                           save_df_csv_file_name=None, save_list_txt_file_name=None):
    """
    Create a list of differentially expressed genes

    :param df: Differentially expressed dataframe the result of differential_expression_analysis method
    :param gene_name_col: Column in dataframe that contains the names of genes i.e. values returned
    :param differential_analysis_col: Column in dataframe that has the result of the differential analysis
    :param sep_rows_string: String used to separate multiple genes within the same cell of the gene_name_col commonly '///'. Explodes the row out giving each gene its own row. If each gene already has its own row leave as 'None'
    :param save_df_csv_file_name: File name to save the data frame of differentially expressed genes should end with the extension '.csv'
    :param save_list_txt_file_name: File name to save the differentially expresssed genes as a list should end with the extensioh '.csv'
    :return: A dictionary with two items:
    :key 'list': list of the differentially expressed genes
    :key 'dataframe': dataframe consisting of two columns, gene column and differential expression analysis
    """

    df = df[[gene_name_col, differential_analysis_col]]
    df = df.dropna()

    df = df.loc[df[differential_analysis_col] != 0]

    if sep_rows_string is not None:
        df[gene_name_col] = df[gene_name_col].apply(lambda x: x.split(sep_rows_string))
        df = df.explode(gene_name_col)
        df[gene_name_col] = df[gene_name_col].str.strip()

    if save_df_csv_file_name is not None:

        if not save_df_csv_file_name.endswith('.csv'):
            raise Exception(
                f"save_df_csv_file_name needs to end with the extension '.csv' \nnot: {save_df_csv_file_name}")

        df.to_csv(save_df_csv_file_name, index=False)

    if save_list_txt_file_name is not None:
        lstring = '\n'.join(list(df[gene_name_col]))

        if not save_list_txt_file_name.endswith('.txt'):
            raise Exception(
                f"save_list_txt_file_name needs to end with the extension '.txt' \nnot: {save_list_txt_file_name}")

        with open(save_list_txt_file_name, 'w') as f:
            f.write(lstring)

    return {'dataframe': df, 'list': list(df[gene_name_col])}


