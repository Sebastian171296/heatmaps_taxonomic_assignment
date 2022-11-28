import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from numpy.ma import masked_array
import taxopy


tax_DB = taxopy.TaxDb(nodes_dmp='/share/project2/sebastian/taXaminer-main/nodes.dmp', names_dmp='/share/project2/sebastian/taXaminer-main/names.dmp', keep_files=True)

groups= ['best_hit']

# wanted families
families = ['Enterococcaceae', 'Staphylococcaceae', 'Enterobacteriaceae', 'Moraxellaceae', 'Pseudomonadaceae']

# Define in and ouput path
input_path = '/home/sebastian/projects/heatmap/Subset/best_hit/genus_level/'
output_path = '/home/sebastian/projects/heatmap/Subset/best_hit/genus_level/'

vis_level = 'taxon'
query_vis_level = 'class'
reduced = True

actinomycetia = []
pathogens = []

# colours, if classes are more, just add more colours to this list
colours = ['black', 'cyan', "magenta", "red", "green"]


#inlcude new option 'ispathogen' for efficient sorting
for group in groups:
    heatmap = pd.read_csv(input_path+'heatmap_taXaminer_'+group+'.csv', index_col=0)
    query_sorted = []
    event_sorted = []
    try:
        with open(input_path + 'sorted_query_' + group + '.txt') as f:
            lines = f.readlines()
            for l in lines:
                l = l.replace('\n', "")
                query_sorted.append(l)
    except:
        print('No sorted y axis')
        with open(input_path + 'query_phyla_' + group + '.txt') as f:
            lines = f.readlines()
            for l in lines:
                l = l.replace('\n', "")
                query_sorted.append(l)
    f.close()
    try:
        with open(input_path + 'sorted_event_' + group + '.txt') as f:
            lines = f.readlines()
            for l in lines:
                l = l.replace('\n', "")
                event_sorted.append(l)
    except:
        print('No sorted x axis')
        with open(input_path + 'event_phyla_' + group + '.txt') as f:
            lines = f.readlines()
            for l in lines:
                l = l.replace('\n', "")
                event_sorted.append(l)
    f.close()
    if 'root' in heatmap.index:
        event_sorted.append('root')
    heatmap['species_ordered'] = pd.Categorical(heatmap['species'],categories=query_sorted,ordered=True)
    heatmap.sort_values('species_ordered', inplace=True)
    heatmap = heatmap.drop('species', axis=1)
    heatmap = heatmap.drop('species_ordered', axis=1)
    heatmap = heatmap.drop('has_actinomycetia', axis=1)
    heatmap = heatmap.drop('has_pathogens', axis=1)
    heatmap_short = heatmap.copy()
    heatmap_short = heatmap_short.drop('query_vis_level', axis=1)
    heatmap_short = heatmap_short[event_sorted]
 
    if reduced:
        heatmap_short = heatmap_short.loc[(heatmap_short != 0).any(1)]
    # Define colour for different classes of query
    classes = heatmap['query_vis_level'].unique()
    tick_colour_def = pd.DataFrame(0, index=[], columns=['black'])
    i = 0

    for c in classes:
        new_query = pd.DataFrame ('', index=[c], columns=['black'])
        tick_colour_def = pd.concat([tick_colour_def, new_query], ignore_index=False, axis=0)
        tick_colour_def.at[c, 'black'] = colours[i]
        i= i+1
    tick_colour = []
    for l, row in heatmap_short.iterrows():
        labels = heatmap.at[l,'query_vis_level']
        tick_colour.append(tick_colour_def.at[labels, 'black'])
    xtick_colour = []
    for i in heatmap_short.columns:
        tax_ids = taxopy.taxid_from_name (i, tax_DB)
        tax_id = tax_ids[0]
        tax = taxopy.Taxon(tax_id, tax_DB)
        if tax.rank_name_dictionary.get('family') in families:
            pathogens.append(tax.name)
            xtick_colour.append('red')
        elif tax.rank_name_dictionary.get('class') == 'Actinomycetia':
            actinomycetia.append (tax.name)
            xtick_colour.append ('blue')
 
    heatmap_short_actinomycetia = heatmap_short.copy()
    heatmap_short_pathogens = heatmap_short.copy()
    heatmap_short_actinomycetia[pathogens] = np.nan
    heatmap_short_pathogens[actinomycetia] = np.nan
    no_taxassignment = masked_array(heatmap_short, heatmap_short > 0)
    taxassign_actinomycetia = masked_array(heatmap_short_actinomycetia, heatmap_short_actinomycetia <= 0)
    taxassign_pathogen = masked_array(heatmap_short_pathogens, heatmap_short_pathogens <= 0)
    # make plot
    label_size = min(np.round(len(heatmap_short.index), 0), np.round(len(heatmap_short.columns)/100000,0))
    fig, ax = plt.subplots(figsize=(np.round(len(heatmap_short.index)+1,0), np.round(len(heatmap_short.columns)/6.6,0))) #plot größe
    plt.rc('font', size=label_size+1)
    shw1 = ax.imshow(no_taxassignment, cmap='binary', aspect=30)
    shw2 = ax.imshow(taxassign_actinomycetia, cmap='Blues',aspect=30, norm='log', vmax=200)
    shw3 = ax.imshow(taxassign_pathogen, cmap='Reds',aspect=30, norm='log', vmax=200)
    
    bar = plt.colorbar(shw2,ax = [ax],shrink=0.2, location='top', anchor=(0.8, -1.395))
    bar.outline.set_linewidth(0.2)
    bar.ax.tick_params(labelsize=label_size/2, width=0.2, length=0.5, pad=0.1)
    bar.ax.minorticks_off()
    
    bar2 = plt.colorbar(shw3,ax = [ax],shrink=0.2, location='top', anchor=(0.2, -0.318))
    bar2.outline.set_linewidth(0.2)
    bar2.ax.tick_params(labelsize=label_size/2, width=0.2, length=0.5, pad=0.1)
    bar2.ax.minorticks_off()
    
    plt.yticks(np.arange(0, len(heatmap_short.index), 1), heatmap_short.index, style='italic' )
    plt.xticks(np.arange(0, len(heatmap_short.columns), 1), heatmap_short.columns, rotation=90, fontsize=0.00001) ### Datenpunkte unter der X-Achse 
    plt.tick_params(axis='x', which='major', labelsize = label_size/4000, width=0.3, length=0.5, pad=0.1) # pad= um die beschriftung näher ranzuholen
    plt.tick_params(axis='y', which='major', labelsize = label_size/4000, width=0.3, length=0.5, pad=0.1)
    plt.xlabel('Number of genes with taxonomic assignment belongig to the particular ' + vis_level, size = label_size)
    plt.ylabel('Query organism', size = label_size)
    plt.title('\n\n\n\n\n$\it{taXaminer}$ results for Soil_invertebrates', fontweight="bold")
  
    for axis in ['top','bottom','left','right']:
    	ax.spines[axis].set_linewidth(0.2)
    #bar.set_label('Actinomycetia genes', size=label_size+1)
    #bar2.set_label('Pathogen genes', size=label_size+1)
    # print colouring query organisms
    #for ticklabel, tickcolor in zip(plt.gca().get_yticklabels(), tick_colour):
     #   ticklabel.set_color(tickcolor)
    for ticklabel, tickcolor in zip(plt.gca().get_xticklabels(), xtick_colour):
        ticklabel.set_color(tickcolor)
    #for i in range(heatmap_short.shape[0]):
        #for j in range(heatmap_short.shape[1]):
            #if heatmap_short.iat[i, j]>0:
                #ax.text(j, i, "{:.0f}".format(heatmap_short.iat[i, j]), ha="center", va="center", size=label_size/2) ## Kachelbeschriftung
    ax.set_aspect(1)   
    #plt.text(1+(label_size/200), 1-(label_size/400), 'Colours:\n', horizontalalignment='left', verticalalignment='top', transform=ax.transAxes, color='black', size=label_size)
    i = 1-2*(label_size/400)
    plt.savefig(output_path +'heatmap_taXaminer_'+group+'.png', bbox_inches='tight', dpi=1000)
