import pandas as pd
import numpy as np
from scipy import stats
from sklearn import linear_model
import warnings
from sklearn.exceptions import ConvergenceWarning
from drugstone import deep_search
import os

warnings.filterwarnings("ignore", category=ConvergenceWarning)


# -------------- Drug Search Methods --------------#
drug_search_algorithms = ["adjacentDrugs", "trustrank", "multisteiner", "keypathwayminer",
                          "closeness", "degree",  "proximity",  "betweenness"]


def drug_search(gene_list, algorithm='trustrank', output_csv_file_name=None):
    """
    Search for drugs related a list of supplied genes

    :param gene_list: List of genes to search for. If String treat as a path to file where genes are listed with a new line between them. Otherwise pass a list of genes
    :param algorithm: Algorithm to use to search for drugs. Defaults to 'trustrank' other possible choices - "adjacentDrugs", "trustrank", "multisteiner", "keypathwayminer", "closeness", "degree", "proximity", "betweenness"
    :param output_csv_file_name: String name of csv file to save results. Defaults to None and does not save the dataframe
    :return: Dataframe showing the related drugs and if selected a csv file of the results
    """

    if algorithm not in drug_search_algorithms:
        raise Exception(f'Invalid "algorithm" arg ({algorithm})\nUse one of the following: {drug_search_algorithms}')

    if isinstance(gene_list, str):
        with open(gene_list, 'r') as f:
            gene_list = f.read().splitlines()

    parameters = {'algorithm': algorithm}
    search = deep_search(gene_list, parameters)
    results = search.get_result()

    if output_csv_file_name is not None:
        results.download_drugs_csv(name=output_csv_file_name)

    return results


def drug_search_polobag_analysis_files(polobag_analysis_directory, output_directory=None, algorithm='trustrank'):
    """
    Apply the drug search method to the analysed PoLoBag files

    :param polobag_analysis_directory: String path to the directory containing the PoLoBag analysis files
    :param output_directory: String path to directory to save the drug search outputs. Defaults to None and saves the outputs to the current directory
    :param algorithm: String algorithm to use when searching for related drugs. Defaults to 'trustrank' other possible choices - "adjacentDrugs", "trustrank", "multisteiner", "keypathwayminer", "closeness", "degree", "proximity", "betweenness"
    :return: Performs the drug search and saves the files
    """
    if algorithm not in drug_search_algorithms:
        raise Exception(f'Invalid "algorithm" arg ({algorithm})\nUse one of the following: {drug_search_algorithms}')

    if output_directory is None:
        output_directory = polobag_analysis_directory
    else:
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

    for filename in os.listdir(polobag_analysis_directory):
        if not filename.endswith(".xlsx"):
            continue

        out_file = filename.split("_polobag_analysis.xlsx")[0]
        out_file = out_file.split('module_')[-1]
        out_file = "".join(['module_', out_file, '_active_gene_drug_search.csv'])

        out_file = os.path.join(output_directory, out_file)

        filename = os.path.join(polobag_analysis_directory, filename)

        df = pd.read_excel(filename, sheet_name='active_genes_found')
        gl = df['active_genes_found'].tolist()

        drug_search(gl, algorithm, out_file)


# -------------- PoloBag Methods --------------#
def PoLoBag(n12, n22, nM, nB, alpha, infilename, outfilename):
    """
    PoLoBag Algorithm

    :param n12: Float number of linear features in each bootstrap sample
    :param n22: Float number of nonlinear features in each bootstrap sample
    :param nM: Float bootstrap sample size
    :param nB: Int total number of bootstrap samples in the ensemble
    :param alpha: Float Lasso regularization parameter
    :param infilename: Name of input file
    :param outfilename: Name of output file
    :return: Results of the PoLoBag algorithm saved to the outfilename
    """
    # Read input expression file
    D = {}
    genes = []
    print('Reading File')
    with open(infilename, 'r') as f:
        geneCount = -1
        for line in f:
            geneCount += 1
            if geneCount == 0:  # header
                continue
            lineData = line.rstrip().split("\t")  # tab separated input file
            genes.append(lineData[0])
            D[lineData[0]] = stats.zscore(np.array(lineData[1:]).astype(float))
    genes = np.unique(genes)

    regs = genes  # Potential regulators

    edges = []
    w = []
    print("Regression for targets genes")
    print(f'Total Number of Genes: {len(genes)}')
    gene_num = 1
    for t in genes:  # Regression problem for each target gene
        yt = np.transpose(D[t][np.newaxis])
        Xt = yt[:, []]
        for reg in regs:
            if not reg == t:  # no auto-regulation
                Xt = np.hstack((Xt, np.transpose(D[reg][np.newaxis])))
        ntR = Xt.shape[1]
        if ntR == 0:
            continue
        valm = np.median(np.hstack((yt, Xt)))
        tau = 3
        yt = yt - tau * valm
        Xt = Xt - tau * valm
        nlin = int(n12 * np.sqrt(ntR))
        nnlin = int(n22 * np.sqrt(ntR))
        wt = np.zeros((ntR, 1))  # Weight
        wtM = np.zeros((ntR, 1))  # Weight magnitude
        wtS = np.zeros((ntR, 1))  # Weight sign
        stw = np.zeros((ntR, 1)) + 1e-10  # How many times a feature was chosen in a sample
        m = Xt.shape[0]
        for b in range(nB):  # For each sample
            rowindex = np.random.choice(m, int(nM * m))
            ytb = yt[rowindex, :]
            XtbF = Xt[rowindex, :]
            # Separate idtb for linear and nonlinear features
            idtblin = []
            idtbnlin = []
            if nlin > 0:
                idtblin = np.random.choice(ntR, nlin, replace=False)
            Xtb = XtbF[:, idtblin]  # linear features
            if nnlin > 0:
                allindex = [x for x in range(ntR) if x not in idtblin]
                index = np.random.choice(allindex, 2 * nnlin, replace=False)
                Xtb21F = XtbF[:, index[:nnlin]]
                Xtb22F = XtbF[:, index[nnlin:]]
                vala = Xtb21F * Xtb22F
                vals1 = np.sign(Xtb21F)
                vals12 = np.sign(vala)
                atbk = 0.5 * np.sqrt(np.abs(vala))
                ftb21 = vals1 * (1 + vals12) * atbk
                ftb22 = vals1 * (1 - vals12) * atbk
                Xtb = np.hstack((Xtb, ftb21[:, :int(0.5 * nnlin)], ftb22[:, int(0.5 * nnlin):]))  # nonlinear features
                for l in range(int(0.5 * nnlin)):
                    idtbnlin.append(str(index[l]) + '#' + str(index[l + nnlin]))
                for l in range(int(0.5 * nnlin), nnlin):
                    idtbnlin.append('-' + str(index[l]) + '#' + str(index[l + nnlin]))  # an indicator for category

            clf = linear_model.Lasso(alpha=alpha, fit_intercept=False)  # Lasso
            clf.fit(Xtb, ytb)
            wtb = clf.coef_
            if len(idtblin) > 0:  # for linear features
                wtM[idtblin, 0] += np.abs(wtb[0:len(idtblin)])
                wtS[idtblin, 0] += wtb[0:len(idtblin)]
                stw[idtblin, 0] += 1
            for l in range(len(idtbnlin)):  # for nonlinear features
                indi = idtbnlin[l].split("#")
                valai = np.sqrt(np.abs(wtb[l + len(idtblin)]))
                valsi = np.sign(wtb[l + len(idtblin)]) * valai
                if indi[0][0] == '-':
                    wtM[int(indi[0][1:]), 0] += valai
                    wtS[int(indi[0][1:]), 0] += valsi
                    stw[int(indi[0][1:]), 0] += 1
                    wtM[int(indi[1]), 0] += valai
                    wtS[int(indi[1]), 0] += -valsi
                    stw[int(indi[1]), 0] += 1
                else:
                    wtM[int(indi[0]), 0] += valai
                    wtS[int(indi[0]), 0] += valsi
                    stw[int(indi[0]), 0] += 1
                    wtM[int(indi[1]), 0] += valai
                    wtS[int(indi[1]), 0] += valsi
                    stw[int(indi[1]), 0] += 1

        wt = np.sign(wtS) * wtM / stw  # Compute weights as average
        j = 0
        for reg in regs:
            if not reg == t:
                val = wt[j, 0]
                if (abs(val) > 0.0):
                    edges.append(reg + "\t" + t)
                    w.append(val)
                j += 1
    gene_num += 1

    sortindex = np.flip(np.argsort(np.abs(w))[np.newaxis], 1).astype(int)  # Sort by absolute value of edge weights
    sedgevals = [w[s] for s in sortindex[0, :]]
    sedgevals = sedgevals / abs(sedgevals[0])
    sedges = [edges[s] for s in sortindex[0, :]]
    ofile = open(outfilename, 'w')
    print('regulator' + '\t' + 'target' + '\t' + 'score', file=ofile)
    for s in range(len(sedges)):
        print(sedges[s] + "\t" + str(sedgevals[s]), file=ofile)  # write to output edge list file

    ofile.close()


def save_expression_data_for_polobag(expression_data, save_txt_file_name, module_genes=None,
                                     samples_to_ignore=None, genes_to_ignore=None):
    """
    Format expression data for PoLoBag. Saves expression data to the PoLoBag format where genes are the index of a dataframe and columns are the samples

    :param expression_data: Dataframe of expression data. Index should correspond to genes and columns are the samples
    :param save_txt_file_name: String file name to save the expresison data should end with the extension '.txt'
    :param module_genes: List of genes (index values) to include in the data. If left as None defaults to all genes present
    :param samples_to_ignore: List of samples (columns) to ignore. If left as None defaults to keeping all columns
    :param genes_to_ignore: List of genes (index values) to ignore. If left as Non defaults to keeping all indices
    :return: A txt file of the saves expression data
    """
    if module_genes is not None:
        expression_data = expression_data.loc[expression_data.index.isin(module_genes)]

    if samples_to_ignore is not None:
        expression_data = expression_data.drop(samples_to_ignore, axis=1)

    if genes_to_ignore is not None:
        expression_data = expression_data.drop(genes_to_ignore)

    expression_data.to_csv(save_txt_file_name, sep='\t')


def prepare_expression_data_polobag(expression_data, feature_data=None, feature_gene_id_col=None,
                                    feature_gene_sep_string=None, phenotype_data=None, phenotype_col_name=None,
                                    expression_data_output_txt_file_name=None):
    """
    Helper method to prepare input for PoLoBag. This method is used to alternate the expression data index (genes) with the symbol names that appear in the feature data. The phenotype data is used to exclude those samples that were not used in the differential expression analysis

    :param expression_data: Dataframe of expression data
    :param feature_data: Dataframe of feature data. If left as None no transformations based on feature data are applied
    :param feature_gene_id_col: String column name of gene symbols in the feature data. If feature data is not None this parameter is required
    :param feature_gene_sep_string: String character used if multiple genes are kept in the same row commonly '///'. This separates them creating an individual row per gene. If left as None no separation is carried out
    :param phenotype_data: Dataframe of the phenotype data
    :param phenotype_col_name: String column name that shows phenotypic assignments all Nan phenotypes are removed from the expression data. If phenotype data is not None this value is required commonly his value will be 'phenotype'
    :param expression_data_output_txt_file_name: String file name to save the expression data if supplied must end with the extension '.txt'
    :return: The modified expression data
    """

    if feature_data is not None:
        if feature_gene_id_col is None:
            raise Exception('Need an ID column to replace current expression index')

        feature_data = feature_data[[feature_gene_id_col]]
        feature_data = feature_data.dropna()

        if feature_gene_sep_string is not None:
            feature_data[feature_gene_id_col] = feature_data[feature_gene_id_col].apply(lambda x: x.split(feature_gene_sep_string))

            feature_data = feature_data.explode(feature_gene_id_col)
            feature_data[feature_gene_id_col] = feature_data[feature_gene_id_col].str.strip()

        expression_data = pd.merge(expression_data, feature_data, left_index=True, right_index=True)
        expression_data.index = expression_data[feature_gene_id_col]
        expression_data.index.name = None
        expression_data = expression_data.drop(feature_gene_id_col, axis=1)

    if phenotype_data is not None:
        if phenotype_col_name is None:
            raise Exception('Need a phenotype column to drop samples omitted from differential analysis')

        phenotype_data = phenotype_data.loc[phenotype_data[phenotype_col_name].isna()]
        samples_to_remove = list(set(phenotype_data.index))
        expression_data = expression_data.drop(samples_to_remove, axis=1)

    if expression_data_output_txt_file_name is not None:
        save_expression_data_for_polobag(expression_data, expression_data_output_txt_file_name)

    return expression_data


def run_polobag_on_modules(modules, expression_data, output_directory, feature_data=None, feature_gene_id_col=None,
                feature_gene_sep_string=None, phenotype_data=None, phenotype_col_name=None, n12=0.5,
                n22=2.5, nM=0.5, nB=500, alpha=0.1):
    """
    Run PoLoBag on modules found from the DOMINO algorithm

    :param modules: Either a list of modules or Str of the file name of the modules
    :param expression_data: Expression data if no phenotype or feature data is supplied this expression data should have been altered using the prepare_expression_data_for_polobag method
    :param output_directory: String of the directory where to save output. If the directory doesn't exist it is created
    :param feature_data: Feature data if altering the expression data using the prepare_expression_data_for_polobag method. Deafault to None if the expression data has already been formatted
    :param feature_gene_id_col: String gene id col in the feature dataframe. Default to None if no feature data is provided
    :param feature_gene_sep_string: String character used if multiple genes are kept in the same row commonly '///'. This separates them creating an individual row per gene. If left as None no separation is carried out
    :param phenotype_data: Dataframe of the phenotype data. Defaults to None
    :param phenotype_col_name: String column name that shows phenotypic assignments all Nan phenotypes are removed from the expression data. If phenotype data is not None this value is required commonly his value will be 'phenotype'
    :param n12: Float number of linear features in each bootstrap sample. Defaults to 0.5
    :param n22: Float number of nonlinear features in each bootstrap sample. Defaults to 2.5
    :param nM: Float bootstrap sample size. Defaults to 0.5
    :param nB: Int total number of bootstrap samples in the ensemble. Defaults to 500
    :param alpha: Float Lasso regularization parameter. Defaults to 0.1
    :return: PoLoBag output per file and a txt file denoting which runs were successful
    """

    if isinstance(modules, str):
        def splitter(li):
            li = li.split('[')[-1]
            li = li.split(']')[0]
            li = li.split(', ')
            return li

        with open(modules, 'r') as f:
            mods = f.read().splitlines()

        modules = list(map(splitter, mods))

    if isinstance(expression_data, str):
        expression_data = pd.read_csv(expression_data, sep='\t', index_col=0)

    if feature_data is not None or phenotype_data is not None:
        expression_data = prepare_expression_data_polobag(expression_data, feature_data, feature_gene_id_col,
                                                          feature_gene_sep_string, phenotype_data, phenotype_col_name)

    if not os.path.exists(output_directory):
        os.mkdir(output_directory)

    exp_files = os.path.join(output_directory, 'polobag_expression_files')
    os.mkdir(exp_files)

    mod_dicts = []
    count = 0
    for mod in modules:
        tmp_string = "".join(['module_', str(count), '_expression_data.txt'])
        file_name = os.path.join(exp_files, tmp_string)
        save_expression_data_for_polobag(expression_data, file_name, mod)

        tmp_dic = {'module_number': count, 'module_genes': mod, 'expression_data_file_name': tmp_string}
        mod_dicts.append(tmp_dic)

        count += 1

    log_file_name = "polobag_expression_modules_log.csv"
    log_file_name = os.path.join(exp_files, log_file_name)

    log_df = pd.DataFrame(mod_dicts)
    log_df = log_df.sort_values(by=['module_number'], ignore_index=True)
    log_df.to_csv(log_file_name, index=False)

    pb_files = os.path.join(output_directory, 'polobag_output_files')
    os.mkdir(pb_files)

    error_files = []
    successful_files = []

    for filename in os.listdir(exp_files):
        if not filename.endswith(".txt"):
            continue

        out_file_name = filename.split('_expression_data.txt')[0]
        out_file_name = out_file_name.split('module_')[-1]
        out_file_name = "".join(["module_", out_file_name, "_polobag_results.tsv"])
        out_file_name = os.path.join(pb_files, out_file_name)

        filename = os.path.join(exp_files, filename)

        try:
            PoLoBag(n12, n22, nM, nB, alpha, filename, out_file_name)
        except Exception as e:
            print(str(e))
            print('\n')
            error_files.append(filename)
            print(f'{filename} : Unsuccessful')
            print('\n' * 2)
            continue
        else:
            print(f'{filename} : Successful')
            print('\n' * 2)
            successful_files.append(filename)

    pb_run_info = ['PoLoBag ran', '\n', 'Parameters:', ": ".join(["n12", str(n12)]), ": ".join(["n22", str(n22)]),
                   ": ".join(["nM", str(nM)]), ": ".join(["nB", str(nB)]), ": ".join(["alpha", str(alpha)]), '\n']

    if len(successful_files) > 0:
        pb_run_info.append('Successful runs:')
        pb_run_info.extend(successful_files)
        pb_run_info.append('\n')

    if len(error_files) > 0:
        pb_run_info.append('Unsuccessful runs:')
        pb_run_info.extend(error_files)

    pb_info_file = os.path.join(pb_files, 'polobag_run_information.txt')
    pb_info_string = "\n".join(pb_run_info)

    with open(pb_info_file, 'w') as f:
        f.write(pb_info_string)


def polobag_analysis_single_module(polobag_df, network_df, active_genes, excel_file_name=None,
                                    polobag_target_col='target', polobag_regulator_col='regulator',
                                    network_df_col1='interactor_a', network_df_col2='interactor_b',
                                    network_df_merge_exists=False, active_gene_drug_search_csv_file=None,
                                    drug_search_algorithm='trustrank'):
    """
    Analyse PoloBag results of a single module to see which predictions appear in the supplied network

    :param polobag_df:  The PoLoBag dataframe that is generated as part of hte run_polobag_method. The dataframe can be supplied directly or a String file path to the dataframe
    :param network_df: Dataframe of the ground truth network, most likely the one used in DOMINO to find the modules
    :param active_genes: List of genes of interest most commonly the active genes supplied as an input for DOMINO
    :param excel_file_name: File name to save the analysis should end with the extension '.xlsx'. If left as None no file is saved
    :param polobag_target_col: String name of the target column of the PoLoBag dataframe. Default is 'target'
    :param polobag_regulator_col: String name of the regulator column of the PoLoBag dataframe. Default is 'regulator'
    :param network_df_col1: String name of the column that depicts the first node in an edge of the network dataframe. Defaults to 'interactor_a' as found in the APID dataframe
    :param network_df_col2: String name of the column that depicts the second node in an edge of the network dataframe. Defaults to 'interactor_b' as found in the APID dataframe
    :param network_df_merge_exists: Boolean value if the 'merged' column exists in the network. This is a helper parameter for the polobag_analysis_multiple_files method. Defaults as False and created the column as part of the method
    :param active_gene_drug_search_csv_file: String file name to save the results of the drug search. If left as None no drug search is carried out
    :param drug_search_algorithm: Algorithm to use to search for drugs. Defaults to 'trustrank' other possible choices - "adjacentDrugs", "trustrank", "multisteiner", "keypathwayminer", "closeness", "degree", "proximity", "betweenness"
    :return: A dictionary
    :key 'polobag_network': Edges predicted by PoLoBag that are found in the supplied network
    :key 'active_gene_network': Edges in the PoLoBag network that contain at least one active gene
    :key 'active_genes': A list of the active genes supplied to the method
    """

    if isinstance(polobag_df, str):
        polobag_df = pd.read_csv(polobag_df, sep='\t')

    polobag_df[polobag_target_col] = polobag_df[polobag_target_col].astype(str)
    polobag_df[polobag_regulator_col] = polobag_df[polobag_regulator_col].astype(str)

    polobag_df['merged'] = [''.join(sorted(filter(None, x)))
                            for x in polobag_df[[polobag_target_col, polobag_regulator_col]].to_numpy()]

    if not network_df_merge_exists:
        network_df[network_df_col1] = network_df[network_df_col1].astype(str)
        network_df[network_df_col2] = network_df[network_df_col2].astype(str)

        network_df['merged'] = [''.join(sorted(filter(None, x)))
                                for x in network_df[[network_df_col1, network_df_col2]].to_numpy()]

    polobag_df = polobag_df.loc[polobag_df['merged'].isin(network_df['merged'])]

    a_gene_network = polobag_df.loc[(polobag_df[polobag_regulator_col].isin(active_genes)) |
                                    (polobag_df[polobag_target_col].isin(active_genes))].copy()

    a_gene_network['active_gene_1'] = np.where(a_gene_network[polobag_regulator_col].isin(active_genes),
                                               a_gene_network[polobag_regulator_col], "")

    a_gene_network['active_gene_2'] = np.where(a_gene_network[polobag_target_col].isin(active_genes),
                                               a_gene_network[polobag_target_col], "")

    a_gene_network['active_genes'] = a_gene_network['active_gene_1'].astype(str) + ' ' + \
                                     a_gene_network['active_gene_2'].astype(str)

    a_gene_network['active_genes'] = a_gene_network['active_genes'].str.strip()

    active_genes_found = list(set(a_gene_network['active_gene_1']).union(set(a_gene_network['active_gene_2'])))
    try:
        active_genes_found.remove("")
    except ValueError:
        pass

    active_genes_found_df = {'active_genes_found': active_genes_found}
    active_genes_found_df = pd.DataFrame(active_genes_found_df)

    a_gene_network = a_gene_network.drop(['active_gene_1', 'active_gene_2'], axis=1)

    output = {'polobag_network': polobag_df, 'active_genes_network': a_gene_network,
              'active_genes_found': active_genes_found_df}

    if excel_file_name is not None:
        with pd.ExcelWriter(excel_file_name) as writer:
            for name, df in output.items():
                df.to_excel(writer, sheet_name=name, index=False)

    if active_gene_drug_search_csv_file is not None:
        drug_search(active_genes_found, drug_search_algorithm, active_gene_drug_search_csv_file)

    return output


def polobag_analysis_multiple_files(path_to_polobag_output_directory, network_df, active_genes, analysis_output_directory,
                                    polobag_target_col='target', polobag_regulator_col='regulator',
                                    network_df_col1='interactor_a', network_df_col2='interactor_b',
                                    active_gene_drug_search_output_directory=None, drug_search_algorithm='trustrank'):
    """
    Analyse multiple PoloBag results to see which predictions appear in the supplied network

    :param path_to_polobag_output_directory: String path to the directory containing PoLoBag results. This is the output to the run_polobag method
    :param network_df: Dataframe of the ground truth network, most likely the one used in DOMINO to find the modules
    :param active_genes: List of genes of interest most commonly the active genes supplied as an input for DOMINO
    :param analysis_output_directory: String path to the output directory where to store the results
    :param polobag_target_col: String name of the target column of the PoLoBag dataframe. Default is 'target'
    :param polobag_regulator_col: String name of the regulator column of the PoLoBag dataframe. Default is 'regulator'
    :param network_df_col1: String name of the column that depicts the first node in an edge of the network dataframe. Defaults to 'interactor_a' as found in the APID dataframe
    :param network_df_col2: String name of the column that depicts the second node in an edge of the network dataframe. Defaults to 'interactor_b' as found in the APID dataframe
    :param active_gene_drug_search_output_directory: String path to directory of where to store the results of drug search of active genes found in the PoLoBag network. If left as None no drug search is carried out
    :param drug_search_algorithm: Algorithm to use to search for drugs. Defaults to 'trustrank' other possible choices - "adjacentDrugs", "trustrank", "multisteiner", "keypathwayminer", "closeness", "degree", "proximity", "betweenness"
    :return: Analysed PoLoBag files and drug search files to the specified directories
    """

    if not os.path.exists(analysis_output_directory):
        os.mkdir(analysis_output_directory)

    network_df[network_df_col1] = network_df[network_df_col1].astype(str)
    network_df[network_df_col2] = network_df[network_df_col2].astype(str)

    network_df['merged'] = [''.join(sorted(filter(None, x)))
                            for x in network_df[[network_df_col1, network_df_col2]].to_numpy()]

    if isinstance(active_genes, str):
        with open(active_genes, 'r') as f:
            active_genes = f.read().splitlines()

    for filename in os.listdir(path_to_polobag_output_directory):
        if not filename.endswith(".tsv"):
            continue

        out_file = filename.split("_polobag_results.tsv")[0]
        out_file = out_file.split('module_')[-1]

        if active_gene_drug_search_output_directory is None:
            drug_file = active_gene_drug_search_output_directory
        else:
            if not os.path.exists(active_gene_drug_search_output_directory):
                os.mkdir(active_gene_drug_search_output_directory)

            drug_file = "".join(['module_', out_file, '_active_gene_drug_search'])
            drug_file = os.path.join(active_gene_drug_search_output_directory, drug_file)

        out_file = "".join(['module_', out_file, '_polobag_analysis.xlsx'])
        out_file = os.path.join(analysis_output_directory, out_file)

        filename = os.path.join(path_to_polobag_output_directory, filename)

        polobag_analysis_single_module(filename, network_df, active_genes, out_file, polobag_target_col,
                                        polobag_regulator_col, network_df_col1, network_df_col2,
                                        network_df_merge_exists=True, active_gene_drug_search_csv_file=drug_file,
                                        drug_search_algorithm=drug_search_algorithm)

        print(f'{filename}: Successful')




