from ete3 import Tree, NCBITaxa
import sys, os
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database() ###Zum Updaten der Datenbank

datei = open("event_ids.txt", "r")
Data = datei.read()
IDs = Data.split("\n")
int_list = [int(x) for x in IDs]
tree = ncbi.get_topology(int_list)


sorted_ids=[]
for leaf in tree.iter_leaves():
  sorted_ids = sorted_ids + [int(leaf.name)]



sorted_names=[]
for e in sorted_ids:
  name = ncbi.get_taxid_translator([e])
  sorted_names.append(name[e])


with open('event_phyla_subset_taxo_assignment.txt', 'w') as f:
    f.write("\n".join(map(str, sorted_names)))
f.close()



	

	




