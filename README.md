# PancreaticBetaCellMetabolism
Code Depository for the manuscript “Kinetic and data-driven modeling of pancreatic -cell central carbon metabolism and insulin secretion”. Will update upon revision and publication. 

To run all analyses described in manuscript, use “AllAnalysesCollectedTimeCourse”. 
That file uses “AllPredictions_CollectedTimeCourse”, “GetAveRxnRate”, and “getRXN_Velocity_INS1E” files,  which all use the “INS1E_Pathway” ODE file. Analyses with the Metformin perturbations were done by calling the aforementioned files with the “INS1E_Pathway_Metformin” ODE model file. 

To perform model sensitivity analysis, we used the “model_efast” file
To perform parameter fitting, we used the “ODE_driver_CollectedTimeCourse”, “ODE_ObjFunc_CollectedTimeCourse”, “Run_simulation_fitting”, and “simulation_protocol_fitting_CollectedTimeCourse” files. Those files require the PESTO parameter estimation toolbox to be downloaded and accessible on the matlab path. 

