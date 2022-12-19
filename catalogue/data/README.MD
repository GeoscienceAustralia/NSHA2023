List of catalogues used to determine MFDs for the 2023 NSHA:
	
* NSHA23CAT.MW.V0.1.csv: Raw Australian earthquake catalogue developed for the 2023 NSHA:
	- Preserves all magnitude types and converted magnitudes
	- Uses hierarchical structure as outlined in Allen et al (in prep) to determine preferred magnitude type
	- Uses Ghasemi & Allen (2017) as preferred ML-MW conversions
	- Uses updated mb-MW and MS-MW conversions as outlined in Allen et al (in prep)
	- Preferred magnitude column is MX_REVML
	- Further details on the catalogue's development is [here](https://github.com/GeoscienceAustralia/NSHA2023/blob/master/info/NSHA23_description_of_files_and_metadata.XLS) (to be updated)

* NSHA23CAT_V0.1_hmtk_mx_orig.csv: full catalogue developed for the 2023 NSHA in [Hazard Modellers Toolkit (HMTK)](https://github.com/GEMScienceTools/hmtk) format:
	- Based on NSHA23CAT.MW.V0.1.csv
	- Preferred magnitude is [mx_origML](https://github.com/GeoscienceAustralia/NSHA2023/tree/master/catalogue/data/2012_catalogue_readme.txt) as required for the Leonard (2008) declustering algorithm
	- Catalogue adds several additional magnitude columns to standard HMTK format.  Details are [here](https://github.com/GeoscienceAustralia/NSHA2023/blob/master/info/NSHA23_description_of_files_and_metadata.XLS) and preferred Mw as per 2023 conversions
	- Not declustered
	
* NSHA23CAT_V0.1_hmtk_declustered.csv: declustered catalogue developed for the 2023 NSHA in HMTK format:
	- Based on NSHA23CAT.MW.V0.1.csv
	- Declustered using Leonard ([2008](https://pubs.geoscienceworld.org/ssa/bssa/article/98/3/1458/341978)) algorithm using [mx_origML](https://github.com/GeoscienceAustralia/NSHA2023/blob/master/info/NSHA23_description_of_files_and_metadata.XLS)
	- Preferred magnitude is in MW as per 2023 conversions
	- Catalogue adds several additional magnitude columns to standard HMTK format.
	
* Catalogues used to generate magnitude conversion equations (Allen et al, in prep):
	- NSHA23CAT.ML-MW.csv
	- NSHA23CAT.mb-MW.csv
	- NSHA23CAT.MS-MW.csv
	
* ISC-GEM_V9.1_hmtk_full.csv: full [International Seismological Centre-Global Earthquake Model (ISC-GEM) earthquake catalogue](http://www.isc.ac.uk/iscgem/overview.php) in HMTK format:
	- Catalogue to be used for adjacent offshore sources
	- Refomated from original catalogue format as described in Storchak et al. ([2015](https://www.sciencedirect.com/science/article/pii/S003192011400154X))
	- Version 9.1 released on 2022-06-27 for the period of 1904-2018
	- Appends GA catalogue from 2019-2022 for magnitude MW >= 5.5 (run io_catalogues.ga_query2htmk)
		- mb(ISC) converted to Mw using conversion equations developed by Ghasemi et al. ([2020](https://link.springer.com/article/10.1007/s10518-020-00966-1)) for Papua New Guinea region

* ISC-GEM_V9.1_hmtk_clip.csv: full [International Seismological Centre-Global Earthquake Model (ISC-GEM) earthquake catalogue](http://www.isc.ac.uk/iscgem/overview.php) in HMTK format:
	- Catalogue to be used for adjacent offshore sources
	- Refomated from original catalogue format as described in Storchak et al. ([2015](https://www.sciencedirect.com/science/article/pii/S003192011400154X))
	- Version 9.1 released on 2022-06-27 for the period of 1904-2018
	- Appends GA catalogue from 2019-2022 for magnitude MW >= 5.5
	
* GIS files for catalogues are provided [here](https://github.com/GeoscienceAustralia/NSHA2023/tree/master/catalogue/shapefiles)