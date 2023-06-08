# Rabadan_Tomer2022
Code associated with Rabadan et al 2022
https://www.nature.com/articles/s41467-022-31073-1


Python code for Calcium imaging data analysis
===================================

0. Module file containing function definitions
- monnet_utils.py

1. For segmentation of Spheroids in MoNNets and aggregated activity calculations
- Max Projections (along time-axis) of Ca2+ imaging data
   GenMaxProj.ipynb

- Registration of before and after treatment Ca2+ imaging datasets
   Registration_Before_and_After_Treatment.ipynb

- Semi-automated segmentation and aggregated calcium signal calculations
   SegmentSpheroids.ipynb
   
2. Neuronal activity source extraction
- to extract spatial footprints of neuronal activity sources
   ActivitySourceExtractions.ipynb

- DF/F calculations
   DFF_Neuronal_Sources.ipynb

3. Analysis of MoNNet network activity
- quantifying network dynamics over time
   NetworkActivityAnalysis.ipynb

- Pharmacological Treatment Comparative analysis
   PharmacologicalCharacterization.ipynb

- SCZ-associated phenotype rescue data analysis
   SCZ_Phenotype_RescueAnalysis.ipynb


MATLAB script for hierarchical consensus analysis
=================================================
using the package hierarchicalConsensus: https://github.com/LJeub/HierarchicalConsensus
consensus_clust.m

R scripts for RNAseq data analysis
===============================
RNA_analysis.R
