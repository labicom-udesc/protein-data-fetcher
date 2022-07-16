# protein-data-fetcher
Fetches fasta and pbd files, runs [psi/pss]pred and rosetta fragpicker on a protein

The protein-data-fetcher project is used to generate the protein fragments.

Use the frag_picker.sh file to configure the protocol used to generate the fragments.

Use the protein_fetcher_config.yaml file to configure the directory of secondary structure predictors used.

To run the project use:
"python3 protein_fetcher.py protein_name"
ex: python3 protein_fetcher.py 1crn

-----------------------------

