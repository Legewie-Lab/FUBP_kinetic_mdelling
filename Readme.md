

This repository provides MATLAB code for two different models of gene splicing: exon definition and intron definition.

The exon_definition_ode.m and intron_definition_ode.m files contain the differential equations for the respective models, while the solv_ode_exondef.m and solv_ode_introndef.m files solve the differential equations for each model.

To simulate a population of 10,000 in silico exons (Fig. Xxx in <citation of the paper>), you can use the sampling_model_exondef.m and sampling_model_introndef.m scripts. These scripts provide the solution to the sampling and output Excel tables with different conditions, as well as plots for delta PSI.

Additionally, this repository contains Excel data files (exon_definition_data1.7.xlsx, exon_definition_data3.5.xlsx, intron_definition_data1.7.xlsx, and intron_definition_data3.5.xlsx) containing stimulation values for exon and intron definition models with perturbations of 1.7 or 3.5.

If you're interested in learning more about the specifics of gene splicing and the differences between the exon definition and intron definition models, please see the included documentation , the pdf files named 'Topology of exon definition model' and 'Topology of intron definition model'.
 
