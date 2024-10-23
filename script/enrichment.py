import os
import glob
import pandas as pd
import numpy as np
import ast
import scipy.special
from scipy.stats import hypergeom
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

# Function to process each CSV file in the folder structure as dam_file
def process_csv_in_folder(psat_file, folder_path, compounds_list_N):
    # Read the psat.csv as the main dataframe
    df = pd.read_csv(psat_file)
    df['Compounds'] = df['Compounds'].apply(convert_string_to_list)

    # Traverse through all subdirectories and find all CSV files (dam_files)
    for dam_file in glob.glob(os.path.join(folder_path, '**/*.csv'), recursive=True):
        print(f"Processing dam_file: {dam_file}")
        
        # Read the dam_file CSV
        dam = pd.read_csv(dam_file)

        # Perform operations on the dam_file to extract kegg_id and compounds_list_n
        kegg_id = pd.concat([dam['External Identifier'].dropna(), dam['External Identifier'].dropna()]).unique().tolist()
        compounds_list_n = [value for value in kegg_id if value.startswith('C')]
        no_comma = []
        for item in compounds_list_n:
            if ', ' in item:
                no_comma.extend(item.split(', '))
            else:
                no_comma.append(item)
        compounds_list_n = sorted(set(no_comma))

        # Update df with the hits based on the dam_file compounds
        df_updated = update_dataframe_with_hits(df.copy(), compounds_list_N, compounds_list_n)
        result = calculate_p_values(df_updated)
        result = result.sort_values(by='P value', ascending=True)

        # Save enrichment.csv in the same subfolder as dam_file
        output_enrichment_csv = os.path.join(os.path.dirname(dam_file), 'enrichment.csv')
        result.to_csv(output_enrichment_csv, index=False)

        # Filter for visualization and save enrichment.pdf
        result_filtered = result[result['Num_n'] != 0].head(15)
        create_visualization(result_filtered, os.path.dirname(dam_file))

        # Save degree.csv in the same subfolder as dam_file
        generate_node_degree_csv(result_filtered, os.path.dirname(dam_file))

def create_visualization(df, output_dir):
    # Visualization logic here (same as original plotting code)
    df['Color'] = -np.log10(df['P value'])
    df['Size'] = df['Num_n'] * 50
    
    # Initialize Graph
    G = nx.Graph()

    formatted_labels = {}
    for index, row in df.iterrows():
        formatted_label = format_label(row['Description'], 30)
        formatted_labels[row['Description']] = formatted_label
        G.add_node(formatted_label, hits_n=row['Hits_n'], size=row['Size'] * 10, color=row['Color'])

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            intersection = set(df.iloc[i]['Hits_n']).intersection(df.iloc[j]['Hits_n'])
            if intersection:
                G.add_edge(formatted_labels[df.iloc[i]['Description']], formatted_labels[df.iloc[j]['Description']], weight=len(intersection))

    pos = nx.kamada_kawai_layout(G)

    fig, ax = plt.subplots(figsize=(14, 8))
    node_sizes = [G.nodes[node]['size'] for node in G.nodes]
    node_colors = [G.nodes[node]['color'] for node in G.nodes]
    edge_weights = [G[u][v]['weight'] for u, v in G.edges]

    nodes = nx.draw_networkx_nodes(G, pos, node_size=node_sizes, node_color=node_colors, cmap=plt.cm.Wistia, edgecolors='darkorange', linewidths=1.5, ax=ax)
    labels = nx.draw_networkx_labels(G, pos, ax=ax, font_size=6, font_color='#404040')

    edge_colors = [get_edge_color(weight) for weight in edge_weights]
    nx.draw_networkx_edges(G, pos, width=edge_weights, edge_color=edge_colors, alpha=1, ax=ax)

    # Add legends and colorbars (same as original code)
    output_pdf = os.path.join(output_dir, 'enrichment.pdf')
    plt.savefig(output_pdf, bbox_inches='tight')
    plt.close()

def generate_node_degree_csv(df, output_dir):
    # Generate node degree and connected nodes information (same as original)
    G = nx.Graph()
    formatted_labels = {}
    for index, row in df.iterrows():
        formatted_label = format_label(row['Description'], 30)
        formatted_labels[row['Description']] = formatted_label
        G.add_node(formatted_label, hits_n=row['Hits_n'], size=row['Size'] * 10, color=row['Color'])

    for i in range(len(df)):
        for j in range(i + 1, len(df)):
            intersection = set(df.iloc[i]['Hits_n']).intersection(df.iloc[j]['Hits_n'])
            if intersection:
                G.add_edge(formatted_labels[df.iloc[i]['Description']], formatted_labels[df.iloc[j]['Description']], weight=len(intersection))

    node_info = []
    for node in G.nodes:
        neighbors = list(G[node])
        node_info.append({
            'Description': node,
            'Degree': G.degree(node),
            'Connected Nodes': ', '.join(neighbors)
        })

    node_df = pd.DataFrame(node_info)
    node_df.sort_values(by='Degree', ascending=False, inplace=True)
    node_df['Connected Nodes'] = node_df['Connected Nodes'].replace('\n', ' ', regex=True)
    node_df['Description'] = node_df['Description'].replace('\n', ' ', regex=True)

    output_degree_csv = os.path.join(output_dir, 'degree.csv')
    node_df.to_csv(output_degree_csv, index=False)
    return node_df

# Helper Functions (same as your original logic)
def update_dataframe_with_hits(df, compounds_list_N, compounds_list_n):
    df['Hits_N'] = df['Compounds'].apply(lambda x: list(set(x).intersection(compounds_list_N)))
    df['Num_N'] = df['Hits_N'].apply(len)
    
    df['Hits_n'] = df['Compounds'].apply(lambda x: list(set(x).intersection(compounds_list_n)))
    df['Num_n'] = df['Hits_n'].apply(len)
    
    return df

def convert_string_to_list(string):
    try:
        return ast.literal_eval(string)
    except ValueError:
        return []

def calculate_p_values(df):
    p_values = []
    N = df['Num'].sum()
    n = df['Num_n'].sum()

    for index, row in df.iterrows():
        K = row['Num']
        k = row['Num_n']

        p = hypergeom.sf(k-1, N, K, n)
        p_values.append(p)

    df['P value'] = p_values
    return df

def format_label(description, lim):
    words = description.split()
    current_line = ""
    formatted_description = ""
    
    for word in words:
        if len(current_line) + len(word) + 1 > lim:
            if formatted_description:
                formatted_description += "\n"
            formatted_description += current_line
            current_line = word
        else:
            if current_line:
                current_line += " "
            current_line += word

    if current_line:
        if formatted_description:
            formatted_description += "\n"
        formatted_description += current_line

    return formatted_description

def get_edge_color(weight):
    if weight >= 5:
        return '#C40C0C'
    elif weight == 4:
        return '#FF6500'
    elif weight == 3:
        return '#FF8A08'
    elif weight == 2:
        return '#ebb134'
    else:
        return '#f2de9d'

# Main entry point
psat_file = 'psat.csv'  # Path to psat.cs
folder_path = '/Users/bowen/Desktop/multi_omics_analysis/2408_pea 2/scripts/240924_Enrichment_Result/TS_sME'
# two tissues
compounds_list_N = ['C00018', 'C00021', 'C00021', 'C00025', 'C00042', 'C00042', 'C00072', 'C00072', 'C00073', 'C00078', 'C00099', 'C00122', 'C00122', 'C00134', 'C00134', 'C00144', 'C00147', 'C00149', 'C00152', 'C00163', 'C00170', 'C00179', 'C00183', 'C00263', 'C00263', 'C00299', 'C00307', 'C00311', 'C00330', 'C00349', 'C00362', 'C00408', 'C00408', 'C00431', 'C00436', 'C00449', 'C00493', 'C00555', 'C00559', 'C00624', 'C00628', 'C00643', 'C00777', 'C00805', 'C00830', 'C00942', 'C00971', 'C00990', 'C01186', 'C01353', 'C01384', 'C01494', 'C01601', 'C01799', 'C01877', 'C01877', 'C01924', 'C01933', 'C02155', 'C02298', 'C02372', 'C02427', 'C02504', 'C02614, C02612', 'C02647', 'C02700', 'C02714', 'C02728', 'C02926', 'C03090', 'C03194', 'C03340', 'C03343', 'C03401', 'C03440', 'C03564', 'C03618', 'C03646', 'C03884', 'C03943', 'C04076', 'C04092', 'C04366', 'C05578', 'C05635', 'C05635', 'C05711', 'C05715', 'C05942', 'C05949', 'C05951', 'C05958', 'C06029', 'C06032', 'C06186', 'C06231', 'C06231', 'C06424', 'C06442', 'C06469', 'C06469', 'C07354', 'C08306', 'C09163', 'C09315', 'C10833', 'C11221', 'C11918', 'C12269', 'C12468', 'C12633', 'C14100', 'C14422', 'C14766', 'C14828', 'C14828', 'C14836', 'C14837', 'C15996', 'C16196', 'C16196', 'C16196', 'C16196', 'C16196', 'C16272', 'C16319', 'C16326', 'C16344', 'C16403', 'C16414', 'C16461', 'C16537', 'C16594', 'C16673', 'C17242', 'C18202', 'C19617', 'C19618', 'C19618', 'C19911', 'C20253', 'C20781', 'C21248', 'C22137', 'C22137', 'C22137', 'C22137', 'C22137', 'C22137', 'C22138', 'C00352', 'C02305', 'C01102', 'C00606', 'C00350', 'C00506', 'C00051', 'C00051', 'C03296', 'C00035', 'C04895', 'C03771', 'C00062', 'C02218', 'C20492', 'C00020', 'C20911', 'C20639', 'C00064', 'C00327', 'C00002', 'C00065', 'C05938', 'C00049', 'C01015', 'C00188', 'C03824', 'C00956', 'C02512', 'C00157', 'C21646', 'C01047', 'C06114', 'C03684', 'C06772', 'C00361', 'C06193', 'C16445', 'C05519', 'C03955', 'C00334', 'C01165', 'C05651', 'C15700', 'C05730', 'C05636', 'C00041', 'C05516', 'C05945', 'C00441', 'C16432', 'C01046', 'C06735', 'C20248', 'C04020', 'C00979', 'C20920', 'C03194', 'C01234', 'C05665', 'C03656', 'C00242', 'C17754', 'C20470', 'C00570', 'C02587', 'C00431', 'C01205', 'C22140', 'C00262', 'C03800', 'C00555', 'C16701', 'C12115', 'C00835', 'C00148', 'C18377', 'C12455', 'C05526', 'C06185', 'C05714', 'C00141', 'C20913', 'C00106', 'C12137', 'C05911', 'C00328', 'C02518', 'C04137', 'C01682', 'C00647', 'C05670', 'C01092', 'C00018', 'C17755', 'C15563', 'C06052', 'C01888', 'C17756', 'C14516', 'C16827', 'C00123', 'C00491', 'C01211', 'C00055', 'C16186', 'C01279', 'C17235', 'C12134', 'C00974', 'C15699', 'C00077', 'C02442', 'C02323', 'C08731', 'C00047', 'C02105', 'C00156', 'C00482', 'C02298', 'C02379', 'C11332', 'C21634', 'C01987', 'C03758', 'C10434', 'C10254', 'C21762', 'C00755', 'C02657', 'C09789', 'C00388', 'C10945', 'C00082', 'C10305', 'C05627', 'C00315', 'C01563', 'C02718', 'C00048', 'C03290', 'C00847', 'C00542', 'C01077', 'C00822', 'C01444', 'C00160', 'C19631', 'C02341', 'C03508', 'C05411', 'C01042', 'C00417', 'C00383', 'C02225', 'C01013', 'C01020', 'C02220', 'C00186', 'C01596', 'C00940', 'C12108', 'C05568', 'C01118', 'C00168', 'C01109', 'C17951', 'C00254', 'C00188', 'C20904', 'C20313', 'C02226', 'C02170', 'C00232', 'C00263', 'C00322', 'C20914', 'C05984', 'C20333', 'C00022', 'C00033', 'C03340', 'C20310', 'C17757', 'C05789', 'C03273', 'C01732', 'C01817', 'C05786', 'C01053', 'C22140', 'C20387', 'C05791', 'C03005', 'C00164', 'C00975', 'C01197', 'C00642', 'C00511', 'C05593', 'C05595', 'C00568', 'C00587', 'C05852', 'C03590', 'C01771', 'C14096', 'C03217', 'C00931', 'C00486', 'C05840', 'C04075', 'C00108', 'C00135', 'C00036', 'C00108', 'C02774', 'C17513', 'C16322', 'C00490', 'C00026', 'C00122', 'C16311', 'C05582', 'C17227', 'C01839', 'C08491', 'C20397', 'C00902', 'C11950', 'C08317', 'C08278', 'C17230', 'C16309', 'C02265', 'C11924', 'C11854', 'C02035', 'C02165', 'C11940', 'C11854', 'C05957', 'C11863', 'C02367', 'C04785', 'C14762', 'C16346', 'C16308', 'C16316', 'C19615', 'C05954', 'C11905', 'C16320', 'C18218', 'C00249', 'C21818', 'C04577', 'C07289', 'C19678', 'C14765', 'C19616', 'C04742', 'C01226', 'C14777', 'C11875', 'C16325', 'C02367', 'C16300', 'C06427', 'C08362', 'C01595', 'C00219', 'C03242', 'C00249', 'C00869']
process_csv_in_folder(psat_file, folder_path, compounds_list_N)
