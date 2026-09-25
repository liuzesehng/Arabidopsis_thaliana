awk 'NR==FNR {a[$1]; next} $1 in a' coex_experiment/yellow_module_genes.txt coex2/total/yellow_module_genes.txt | cut -f1 > coex_experiment/yellow_module_genes_intersect.txt
