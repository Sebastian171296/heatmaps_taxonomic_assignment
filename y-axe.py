import taxopy
import pandas as pd

tax_db = taxopy.TaxDb(nodes_dmp='/share/project2/sebastian/taXaminer-main/nodes.dmp', names_dmp='/share/project2/sebastian/taXaminer-main/names.dmp', keep_files=True)

heatmap = pd.read_csv('/home/sebastian/projects/heatmap/Subset/taxo_assignment/species_of_interest/heatmap_taXaminer_subset_taxo_assignment.csv',index_col=0)


sampling = {
"Isopoda": "hand sorting",

"Chilopoda": "hand sorting",

"Diplopoda": "hand sorting",

"Symphyla": "soil sample",

"Diplura": "soil sample",

"Pauropoda": "soil sample",

"Collembola": "soil sample or exhaustor",

"Astigmata": "culture",

"Gamasin": "soil sample or culture",

"Oribatida": "soil sample",

"Enchytraeidae": "soil sample or culture",

"Lumbricina": "hand sorting",

"Nematoda": "soil sample or culture",

"Tardigrada": "culture"
}


list_to_write = []
orders = []
for n in heatmap['species']:
	id = taxopy.taxid_from_name(n, tax_db)
	id = id[0]
	taxon = taxopy.Taxon(id, tax_db)
	linage = taxon.name_lineage
	
	for key in sampling:
		result = linage.count(key)
		if result > 0:
			heatmap['species'] = heatmap['species'].replace([n], f"[{sampling[key]}] {n}")
		else: 
			continue

heatmap.to_csv('/home/sebastian/projects/heatmap/Subset/taxo_assignment/species_of_interest/heatmap_taXaminer_subset_taxo_assignment.csv')





