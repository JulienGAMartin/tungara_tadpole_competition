This readme file was generated on [21-10-2022] by [Cammy Beyts]

GENERAL INFORMATION

Title of Dataset: Trinidad_2019_ACT_EXP_PRED

Author/Principal Investigator Information
Name: Cammy Beyts
ORCID: 0000-0002-4729-2982
Institution: Edinburgh University
Address: The Roslin Institute and R(D)SVS, University of Edinburgh, Easter Bush, UK
Email: cammy.beyts@ed.ac.uk

Author/Associate or Co-investigator Information
Name: Maddalena Cella
Institution: Digital Futures
Address: Warnford Court, 29 Throgmorton St, London EC2N 2AT

Name: Nick Colegrave
Institution: University of Edinburgh
Address: Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, UK

Name: Roger Downie
Institution: University of Glasgow
Address: Institute of Biodiversity Animal Health & Comparative Medicine, University of Glasgow, UK

Name: Julien G.A. Martin
ORCID: 0000-0001-7726-6809
Institution: University of Ottawa
Address: Department of Biology, University of Ottawa, Canada

Name: Patrick Walsh
Institution: University of Edinburgh
Address: Institute of Ecology and Evolution, School of Biological Sciences, University of Edinburgh, UK


Date of data collection: 2019

Geographic location of data collection: Trinidad, Trinidad and Tobago

Information about funding sources that supported the collection of the data: 
This project was funded by a NERC doctoral training partnership grant (NE/L002558/1) and a Davis Expedition Fund grant (E08668) awarded to CB. Funding was also provided by The School of Biology, The University of Edinburgh through funding received by PW.


SHARING/ACCESS INFORMATION

Licenses/restrictions placed on the data: None

Links to publications that cite or use the data: None

Links to other publicly accessible locations of the data: 
https://github.com/cammybeytstungara_tadpole_competition

Links/relationships to ancillary data sets: None

Was data derived from another source? No
If yes, list source(s):

Recommended citation for this dataset: 

Beyts C, Cella M, Colegrave N, Martin, GA, Walsh P. 2022. Data from: The effect of heterospecific and conspecific competition on inter-individual differences in tungara frog tadpole (Engystomops pustulosus) behavior. Behav Ecol. 

DATA & FILE OVERVIEW

File List:

1) Trinidad_2019.R (R script for data analysis)
2) Trinidad_2019_ACT_EXP_PRED.txt (txt file with data to be read into R)
3) SVL_brms_2022c.rda (brms model of tadpole body size)
4) ACT_EXP_PRED_brms_2022a.rds (brms model of tadpole behavior)

Relationship between files, if important: 

Additional related data collected that was not included in the current data package: 

Are there multiple versions of the dataset?
If yes, name of file(s) that was updated: No
Why was the file updated? 
When was the file updated? 


METHODOLOGICAL INFORMATION

Description of methods used for collection/generation of data: 
We reared tungara frog tadpoles (Engystomops pustulosus) either in isolation, with a conspecific tadpole or with a heterospecific tadpole. We looked at the effect of treatment on body size and behavior by measuring tadpole growth and the distance focal E. pustulosus tadpoles swam in familiar, novel and predator risk contexts six times during development.  We used univariate and multivariate hierarchical mixed effect models to investigate the effect of treatment on body size, mean behavior, variance among individuals, the variance within individuals and the repeatability of behaviour in the familiar, novel and predator risk contexts as well as how tadpole behavior was correlated between contexts. 

Methods for processing the data: Data was processed using raw data in RStudio 2021.09.1+372 

Instrument- or software-specific information needed to interpret the data: RStudio 2021.09.1+372

Standards and calibration information, if appropriate: 

Environmental/experimental conditions: 

Describe any quality-assurance procedures performed on the data: 

People involved with sample collection, processing, analysis and/or submission: 
Cammy Beyts (data collection, processing, analysis and submission)
Maddelena Cella (data collection, processing)
Roger Dowinie (data collection)
Julien Martin (analysis)
Nick Colegrave (analysis)
Patrick Walsh (analysis)

DATA-SPECIFIC INFORMATION FOR: Trinidad_2019_ACT_EXP_PRED.txt

Number of variables: 10

Number of cases/rows: 967

Variable List:

Set: Batch number animals tested in
TestingOrder: Order animals were tested in 	
TadpoleID: Tadpole identification number
Rep: Trial number (1-6)
Treatment: Treatment name (no competition, conspecific competition, heterospecific competition)
FocalNest: Nest identification number of focal tadpole	
SVL: Tadpole snout vent length = body size (mm)
ACT.Dist.Pixels: Distance travelled in pixels in the activity assay (familar context)
EXP.Dist.Pixels: Distance travelled in pixels in the exploration assay (novel context)
PRED.Dist.Pixels: Distance travelled in pixels in the predation-risk assay (predatory risk context)

Missing data codes: None

Specialized formats or other abbreviations used: 
ACT = activity assay (familar context)
EXP = exploration assay (novel context)
PRED = predation-risk assay (predatory risk context)
nocomp = no competiton treatment
consp = conspecific treatment 
hetero = heterospeciifc treatment