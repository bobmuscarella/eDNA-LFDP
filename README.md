*Data and scripts for the paper:*

# Soil eDNA reflects regionally dominant species rather than local composition of tropical tree communities

Francesc Borràs Sayas<sup>1</sup>, Ottavia Iacovino<sup>2</sup>, María Uriarte<sup>3</sup>, Jess K. Zimmerman<sup>4</sup>, Christopher J. Nytch<sup>4</sup>, Glenn Dunshea<sup>2*</sup>, and Robert Muscarella<sup>1*</sup>

<sup>1</sup> Department of Ecology and Genetics, Uppsala University, Uppsala, Sweden

<sup>2</sup> Department of Natural History, Norwegian University of Science and Technology, Trondheim 7491, Norway

<sup>3</sup> Department of Ecology, Evolution and Environmental Biology, Columbia University, New York, USA

<sup>4</sup> Department of Environmental Sciences, University of Puerto Rico, Rio Piedras, Puerto Rico

<sup>*</sup> Corresponding authors: robert.muscarella@ebc.uu.se, glenn.dunshea@ntnu.no

Please contact Robert Muscarella for enquires about:
The LFDP, LFDP census data, confusion matrix analyses & simulations

Please contact Glenn Dunshea for enquires about:
The eDNA experimental design, laboratory workflows and bioinformatics

# Instructions
1. Clone this repository

2. Download data files (Zenodo link)

3. Tranfer the datafiles from the unzipped Zenodo file to the directory(ies) of the same name in the cloned eDNA-LFDP repo

4. Open the R project in RStudio and go to the 'Files' tab

5. The scripts are sequentially numbered from 00 onwards to reproduce all analyses and figures presented in this study

6. NOTE: In the 00.dada2+taxonomy+phyloseq.R script, there are instructions at thet start for setting working directory paths, software packages to install and shell scripts run from the terminal - not all of this script is an R script. After DADA2 denoising of sequence data, a portion of the LULU curation of ASV tables is also run from the terminal, which is clearly flagged in the 00 script. All scripts after this point are run from R.

7. Scripts can be:
	a. Either run from first principles from the raw sequencing data received from the sequencing service provider, or
	b. At the start of some scripts after 00, there is the option to open an R image file where the scripts prior to the current script have already been run, so the relevant R objects are already available.



