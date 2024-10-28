# Quinlan_et_al_2024_Molec_Ecol
# Leveraging transcriptional signatures of diverse stressors for bumble bee conservation

## Gabriela M. Quinlan1*, Heather M. Hines1,2, Christina M. Grozinger1 
1.	Penn State University, Department of Entomology, Center for Pollinator Research, Huck Institutes of the Life Sciences University Park, PA, US 16802
2.	Penn State University, Department of Biology University Park, PA, US 16802

Corresponding author: gmq5021@psu.edu

## Abstract 
Organisms in nature are subjected to a variety of stressors, often simultaneously. Foremost among stressors of key pollinators are pathogens, poor nutrition, and climate change. Landscape transcriptomics can be used to decipher the relative role of stressors, provided there are unique signatures of stress that can be reliably detected in field specimens. In this study, we identify biomarkers of bumble bee (Bombus impatiens) responses to key stressors by first subjecting bees to various short-term stressors (cold, heat, nutrition, and pathogen challenge) in a laboratory setting and assessing their transcriptome responses. Using random forest classification on this whole transcriptome data, we were able to discriminate each stressor. Our best model (tissue-specific model trained on a subset of important genes) correctly predicted known stressors with 92% accuracy. We then applied this random forest model to wild-caught bumble bees sampled across a heatwave event at two sites in central Pennsylvania, US, expected to differ in baseline temperature and floral resource availability. Transcriptomes of bees sampled during the heat-wave’s peak showed signatures of heat stress, while bees collected in the relatively cooler morning periods showed signatures of starvation and cold stress. We failed to pick up on signals of heat stress shortly after the heatwave, suggesting this set of biomarkers are more useful for identifying acute stressors than long-term monitoring of chronic, landscape-level stressors. We highlight future directions to fine-tune landscape transcriptomics toward the development of better stress biomarkers that can be used both for conservation and improving understanding of stressor impacts on bees.

## Data
All data required for running these models is provided as a supplement Excel file in the manuscript <doi: XXX>
custom.txt -- Custom GO database, can be created from Table S11 with code in fullScript.R, but this provides a shortcut.

## Code 
preChecks.R -- R code for creating preliminary PCA plot and heatmap for lab stressor pannel transcriptome data.
fullScript.R -- R code for data upload, formatting, running all models, and creating all figures.

