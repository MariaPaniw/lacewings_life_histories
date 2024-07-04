# Comparative life-history responses of lacewings to changes in temperature

This readme.txt file was generated on 2024-07-04 by Maria Paniw


GENERAL INFORMATION

1. Title of Dataset/Repisitory: Data and Analyses from: Comparative life-history responses of lacewings to changes in temperature.

2. Author Information
	A. Principal Investigator Contact Information
		Name: Maria Paniw
		Institution: Doñana Biological Station (EBD-CSIC)
		Address: Seville, 41001 Spain
		Email: m.paniw@gmail.com

	B. Associate or Co-investigator Contact Information
		Name: Hanna Serediuk
		Institution: Doñana Biological Station (EBD-CSIC)
		Address: Seville, 41001 Spain
		Email: anna.serediuk@gmail.com

C. Associate or Co-investigator Contact Information
		Name: Sanne M Evers
		Institution: Doñana Biological Station (EBD-CSIC)
		Address: Seville, 41001 Spain
		Email: sanne.m.evers@gmail.com 

D. Associate or Co-investigator Contact Information
		Name: John Jackson
		Institution: Doñana Biological Station (EBD-CSIC)
		Address: Seville, 41001 Spain
		Email: jjackson0308@gmail.com 


4. Date of data collection: 2023-01-15 until 2023-05-31 for lieterature review. All studies published before 2023-03-31 were considered  

5. Geographic location of data collection: global

6. Information about funding sources that supported the collection of the data: 

Data collection was funded by the Consejo Superior de Investigaciones Científicas (CSIC) with the grant UCRAN20052. MP was additionally funded by the grant RYC2021-033192-I by MCIN/AEI/10.13039/501100011033 by “European Union NextGenerationEU/PRTR”.
 

SHARING/ACCESS INFORMATION

1. Licenses/restrictions placed on the data: none

2. Links to publications that cite or use the data: NA

3. Links to other publicly accessible locations of the data: [here](https://drive.google.com/drive/folders/1jWptakO8ea5g_97oe0XaUMgc4MzH4rsK?usp=sharing)

4. Links/relationships to ancillary data sets: NA

5. Was data derived from another source? No

6. Recommended citation for this dataset: 
Serediuk, H., Jackson, J., Evers, S.M., Paniw, M. Data and Analyses from: Comparative life-history responses of lacewings to changes in temperature. Zenodo. xxxx

DATA & FILE OVERVIEW & WORKFLOW OF R SCRIPTS: 

Scripts are meant to be run in the following order: 

1.	Get coefficients for meta-analysis: **metaanalysis_aquisition_exploration_may24.R**.

    - a.	*Input*: LH_Neuroptera_red.csv; data extracted from artciles through literature review
    - b.	*Output*: neuroptera.RData
 
2.	Run the meta-analysis: **metaanalysis_study_regressions_may24.R**

    - a.	*Input*: neuroptera.RData
    - b.	*Output*: coef_study.RData: Figure 2 main text
 
3.	Perform multivariate analyses on co-variation in 6 life-history processes: **mcmc_pca_Neuroptera_main_text.R**

    - a.	 *Input*: LH_Neuroptera_red.csv
    - b.	*Output*: Table 1; Figures 3 & 4 main text; Figure S3.1; Table S3.1 
 
4.	Perform multivariate analyses on co-variation in 4 life-history processes: **mcmc_pca_Neuroptera_dev_only.R**

    - a.	 *Input*: LH_Neuroptera_red.csv
    - b.	*Output*: Table 2 main text; Figures S3.2-S3.4; Table S3.2 

5.	Perform multivariate analyses on co-variation in 2 life-history processes: **mcmc_pca_Neuroptera_surv_repro_only.R**

    - a.	 *Input*: LH_Neuroptera_red.csv
    - b.	*Output*: Table 3 main text; Figures S3.5 & S3.6; Table S3.2
      

METHODOLOGICAL INFORMATION

1. Description of methods used for collection/generation of data: 

We searched Web of Science and Scopus for literature for studies (published before March 31, 2023) that quantified the effects of temperature on different life-history processes in Neuroptera as described in the previous section. Additionally, we used Google Scholar to search for articles, including grey literature, that may not have been indexed in the aforementioned databases. The latter included book chapters, monographs, or regional journals. For search terms, see the main text of the article. 

To compare the variation in life-history processes (i.e., development times, survival, and reproduction) under different environmental regimes among species, we performed a meta-analysis.

To explore general patterns of among-species differences in life-history strategies independent of species temperature sensitivities, we analysed the covariation in life-history processes (i.e., development times, survival, and reproduction), after accounting for the effect of temperature and other confounding factors on this covariation.


2. Methods for processing the data: 

All data were processed in R

3. Instrument- or software-specific information needed to interpret the data: 

R statistical software, version 4.1.2 (and packages as described in the R scripts)

4. Standards and calibration information, if appropriate: NA

5. Environmental/experimental conditions: as described in primary literature

6. Describe any quality-assurance procedures performed on the data: All R scripts are fully commented and have been checked; all reviewed information was quality-checked by all co-authors.

7. People involved with sample collection, processing, analysis and/or submission: NA


