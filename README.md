# ManuscriptScripts_DCISRecurrenceVsProgression
Scripts to perform the statistical analysis and generate the plots for the manuscript:

***Evolutionary Measures Show that Recurrence of DCIS is Distinct from Progression to Breast Cancer***. *Fortunato, Mallo, Cisneros, et al. (2024)* 

# Instructions
To replicate these analyses, follow these steps:
1. Modify the *configFile.in* file to configure it for your system and rename it as *configFile*. An example configuration file, *configFile.ex*, can be found in the repository.
2. Export an environmental variable with the absolute path of this repository `export ManuscriptScripts_DCISRecurrenceVsProgression="/path/to/the/repository/ManuscriptScripts_DCISRecurrenceVsProgression/"`
3. Obtain or generate the data files indicated in the configFile
4. Execute the *allPaperPlots.R* script

# Input data
With the current version of these scripts, you need to obtain the following datafiles from the authors. Instructions on how to generate each of them will be added soon.

- aim1GeneticDataFile
- aim1IntScoreDataFile
- aim1IntRNDScoreDataFile
- aim1DivergenceDataFile
- aim1EMDistDataFile
- aim1CDIScoreDataFile
- aim4SNVDataFile
- aim4CNADataByPatientFile
- aim4CNADataBySampleFile
- aim4ERDataFile
- aim4GLUTDataFile
- aim4ERPosERDataFile
- aim4ERPosGLUTDataFile
- aim4PhenotypicHeterDataFile
- theClinicalDataFile
