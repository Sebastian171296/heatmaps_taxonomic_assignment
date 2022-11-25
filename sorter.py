from ete3 import Tree, NCBITaxa
import sys, os
ncbi = NCBITaxa()
#ncbi.update_taxonomy_database() ###Zum Updaten der Datenbank

datei = open("taxidlist.txt", "r")
Data = datei.read()
IDs = Data.split("\n")
IDs = IDs[:-1]
int_list = [int(x) for x in IDs]
tree = ncbi.get_topology(int_list)


sorted_ids=[]
for leaf in tree.iter_leaves():
  sorted_ids = sorted_ids + [int(leaf.name)]



sorted_names=[]
for e in sorted_ids:
  name = ncbi.get_taxid_translator([e])
  sorted_names.append(name[e])





with open('event_phyla_best_hit.txt', 'w') as f:
  for item in sorted_names:
    f.write(item + '\n')



	

	




