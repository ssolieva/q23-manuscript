# Files for analysis of sequence search results.

## To parse through output log file of search hits:
	---> parse_all_log.py

Extract information from the search hits, such as donor names, genes, cdrh3 sequences, etc.

	---> logfile_toolbox.py

Functions needed for parse_all_log.py
	
## To analyze search hits: 
	---> calculate_frequencies.py

Calculate frequency of hits per donor. 

	---> frequency_toolbox.py

Functions needed for calculate_frequencies.py

	---> example-remove-redundant-seqs.ipynb

Notebook to show an example of how redundant sequences were removed. 

	---> dgene-enrichment-v033.ipynb  

IGHD gene enrichment analysis and plot for V033. 

	---> CH01_frequency_per_donor.ipynb, PG9_frequency_per_donor.ipynb, RM_frequency_per_donor.ipynb 

Graphing summary frequencies per donor per search.
