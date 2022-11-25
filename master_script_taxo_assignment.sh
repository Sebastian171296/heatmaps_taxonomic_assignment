#!/bin/bash

eval "$(conda shell.bash hook)"

find_in_conda_env(){
    conda env list | grep "${@}" >/dev/null 2>/dev/null
}

if find_in_conda_env ".*ete3.*" ; then
	echo "ete3 conda environment is present."
else 
	echo "ete3 conda environment is missing."
	conda create -n ete3 python=3
	conda activate ete3
	conda install -c etetoolkit ete3 ete_toolchain
	ete3 build check
	conda deactivate
fi

if find_in_conda_env ".*heatmap.*" ; then
	echo "heatmap conda environment is present."
	conda activate heatmap
else 
	echo "heatmap conda environment is missing."
	conda create --name heatmap
	conda activate heatmap
	conda install -c anaconda pandas
	conda install -c bioconda taxopy
	conda install -c conda-forge matplotlib
	conda install -c anaconda numpy
fi

#running sub-scripts 
echo "creating dataframe"
python make_heatmap_taxo_assignment.py
echo "dataframe was created"
echo "sorting x-axis"
conda activate ete3 
python sorter.py
echo "x-axes was sorted" 
conda deactivate
conda activate heatmap

echo "sorting y-axis"
python y-axe.py 
echo "y-axis was sorted" 
echo "creating heatmap" 
python Visualization_taxo_assignment.py
