1) Run: augment_catalogue_data.py:
		- Reads "Merged_Catalogue.csv" and appends data from:
			- "earthquake_query_2017-2022.csv"
			- "combined_au_mw.dat"
		- sets Australian ML region
		- exports "merged_cat.pkl"
		
2) Run: revise_mls_2023.py:
		- applies ML correction according to [Allen (2021)](https://link.springer.com/article/10.1007/s10950-021-10004-5), but updated according to [add ref to catalogue record here]

3) Run: regress_ml_legacy_target.py:
		- develops relationship between ML_legacy and ML_target for pre-instrumental earthquakes
		- outputs 'ml_revision_reg.csv' 
		
4) Rerun: revise_mls_2023.py:
		- applies new regression coefficients determined in regress_ml_legacy_target.py
		- exports "merged_cat_revised_ml.pkl" and "merged_cat_revised_ml.csv"

5) Run Mw conversion codes and output coefficients: 
		- regress_mw_vs_ml.py;
		- regress_mw_vs_ms.py;
		- regress_mw_vs_mb.py
		
6) Run: select_pref_mags.py:
		- applies logic in catalogue documentation to select preferred mags
		- exports "merged_cat_pref_mags.pkl" and "merged_cat_pref_mags.csv"
		
7) Run: ../decluster_2023.py