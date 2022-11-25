import pandas as pd
import os
import taxopy
import warnings

warnings.filterwarnings("ignore")

## Variables
# Change paths to NCBI database for taxopy (used for the taXaminer analysis) and taXaminer results
tax_db = taxopy.TaxDb(nodes_dmp='/share/project2/sebastian/taXaminer-main/nodes.dmp', names_dmp='/share/project2/sebastian/taXaminer-main/names.dmp', keep_files=True)
input_path = '/share/project2/sebastian/Data/'


# Group name of considered data
groups = [['best_hit', '']]


# exclude genes in query phylum
exclude_query_phylum = True
include_cLCA = False
exclude_inLineage = True

# Output path
output_path = '/home/sebastian/projects/heatmap/Subset/best_hit/genus_level/'


# Level of visualization x axis and grouping of y axis (must be found by taxopy)
vis_level = 'genus'
query_vis_level = 'phylum'

# wanted families
families = ['Enterococcaceae', 'Staphylococcaceae', 'Enterobacteriaceae', 'Moraxellaceae', 'Pseudomonadaceae', 'Actinomycetia']

# Errors
missing_files = []
not_found_id = []
not_found_phylum = []
double_species = []
finished = []

# Left out events
inLineage = []

# Create dataframe for heatmap
for group, path_name in groups:
    if path_name == '':
        input_data = input_path
    else:
        input_data = os.path.join(input_path, path_name)
    query_ID = []
    query_phylum = []
    query_name = []
    event_phyla = ['species', 'query_vis_level', 'has_actinomycetia', 'has_pathogens']
    heatmap = pd.DataFrame(0, index=query_ID, columns=event_phyla)
    for root, subdirectories, files in os.walk(input_data):
        for subdirectory in subdirectories:
                file = os.path.join(root, subdirectory)
                file = os.path.join(file, 'taxonomic_assignment/')
                if not os.path.exists(file):
                    missing_files.append(f'Missing taxonomic_assignment directory: {subdirectory,file}')
                    continue
                path = os.path.join(file, 'summary.txt')
                if os.path.exists(path):
                    with open(file + 'summary.txt') as f:
                        content = f.readlines()
                        id = content[0].split(' ')[-1]
                        id = int(id[:-2])
                        query = taxopy.Taxon(id, tax_db)
                        inLineage.append(query.rank_name_dictionary)
                        name = query.name
                        if name in heatmap.index:
                            double_species.append(query.name)
                            break
                        new_query = pd.DataFrame(0, index=[name], columns=event_phyla)
                        heatmap = pd.concat([heatmap, new_query], ignore_index=False, axis=0)
                        heatmap.at[name, 'species'] = query.rank_name_dictionary['species']
                        heatmap.at[name, 'has_actinomycetia'] = 'F'
                        heatmap.at[name, 'has_pathogens'] = 'F'
                        try:
                            heatmap.at[name, 'query_vis_level'] = query.rank_name_dictionary [query_vis_level]
                        except:
                            heatmap.at[name, 'query_vis_level'] = 'not found'
                        query_name.append(query.rank_name_dictionary['species'])
                    taxaminer_result = pd.read_csv(file + 'gene_table_taxon_assignment.csv', index_col=0)
                    assigned_taxaminer_result = taxaminer_result.loc[(taxaminer_result['best_hitID'] != 'Unassigned')]                    
                    results_LCA = assigned_taxaminer_result.loc['best_hit' != assigned_taxaminer_result['best_hitID']]['best_hitID'].value_counts()
                    results_cLCA = assigned_taxaminer_result.loc[assigned_taxaminer_result['corrected_lca'] == assigned_taxaminer_result['taxon_assignment']]['taxon_assignment'].value_counts()
                    for n, count in results_LCA.items():
                        # exclude hits which lay inside lineage (when exclude_inLineage=True)
                        if (n in query.name_lineage) and (exclude_inLineage):
                            inLineage.append (n + ': ' + str (count) + ' found in LCA ' + query.name)
                            continue
                        tax_ids = [n] 
                        if len(tax_ids) > 0:     
                            event_id = tax_ids[0]
                            try:
                                event_tax = taxopy.Taxon(event_id, tax_db)
                            except taxopy.exceptions.TaxidError:
                                 try:
                                 	event_id2 = tax_db.oldtaxid2newtaxid[event_id]
                                 	event_tax = taxopy.Taxon(event_id2, tax_db)
                                 except TypeError:
                                 	#print(query)
                                 	event_tax = query # event_phyl wird zu query, damit es rausfliegt.
                            if (event_tax.rank_name_dictionary.get('class') in families) or (event_tax.rank_name_dictionary.get('family') in families):  #Das hier ist der Filter, hier kommt nur durch was in Famlies abgedeckt ist !
                                if event_tax.rank_name_dictionary.get('class') in families:
                                    heatmap.at[name, 'has_actinomycetia'] = 'T'
                                else:
                                    heatmap.at[name, 'has_pathogens'] = 'T'
#                                print(event_tax)                         
                                event_phylum = event_tax.rank_name_dictionary.get('genus')
                                print(event_phylum)                         
                                if event_phylum not in event_phyla:
                                    if event_phylum is not None:
                                    	event_phyla.append(event_phylum)
                                    	heatmap.insert(len(heatmap.columns), event_phylum, 0)
                                    	heatmap.at[name, event_phylum] = count
                                    else:
                                    	continue
                                else:
                                    if 'na' in set(heatmap[event_phylum]):
                                        heatmap.at[name, event_phylum] = count
                                    else:
                                        heatmap.at[name, event_phylum] = heatmap.loc[name, event_phylum] + count
                            else:
                                if event_tax.name not in not_found_phylum:
                                    not_found_phylum.append(event_tax.name)
                                else:
                                    continue
                            continue
                        else:
                            if n not in not_found_id:
                                not_found_id.append(n)
                    finished.append(f'Finished {query.name,path}')
                else:
                    missing_files.append(f'Missing summary file: {subdirectory,path}')
    heatmap = heatmap.fillna(0)
    heatmap.to_csv(output_path+'heatmap_taXaminer_'+group+'.csv')
    event_phyla.remove('species')
    event_phyla.remove('query_vis_level')
    event_phyla.remove('has_pathogens')
    event_phyla.remove('has_actinomycetia')
    with open(output_path + 'query_phyla_' + group + '.txt', 'w') as f:
        f.write("\n".join(map(str, query_name)))
    f.close()
    with open(output_path + 'event_phyla_' + group + '.txt', 'w') as f:
        f.write("\n".join(map(str, event_phyla)))
    f.close()
with open(output_path + 'not_found_with_ID.txt', 'w') as f:
    f.write("\n".join(map(str, not_found_id)))
f.close()
with open(output_path + 'not_found_phylum.txt', 'w') as f:
    f.write("\n".join(map(str, not_found_phylum)))
f.close()
with open(output_path + 'double_species.txt', 'w') as f:
    f.write("\n".join(map(str, double_species)))
f.close()
with open(output_path + 'missing_files.txt', 'w') as f:
    f.write("\n".join(map(str, missing_files)))
f.close()
with open(output_path + 'inLineage.txt', 'w') as f:
    f.write("\n".join(map(str, inLineage)))
f.close()
with open(output_path + 'finished.txt', 'w') as f:
    f.write("\n".join(map(str, finished)))
f.close()
