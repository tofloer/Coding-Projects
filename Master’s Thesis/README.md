# Computational Part of the Master Thesis: Predicting Interdependent Diseases in Individuals Using Survey Data and Machine Learning Algorithms

This README is intended to give an overview for the provided files and to explain how to reproduce the results of the thesis.


## Files 

- 'MA_Computations_Docu.html' is a HTML file that documents the code and the results. It does not require any software installation and opens in the browser.
- 'GEDA14.sav' is the used SPSS version of the GEDA dataset.
- 'MA_Computational_Part.ipynb' is a Jupyter Notebook that entails runable code and documentation.
- 'MA_Imputation.R' is a R script that is used for imputation.
- 'Build_Tools.png' is a picture that is included in the HTML file and the Jupyter Notebook  

## Software 

In order to run the Jupyter Notebook and the R script such that the MA results are preserved, the following software and versions have to be installed in advance:

- Anaconda software distribution version 2021.11 from (https://repo.anaconda.com/archive/), i.e. Anaconda3-2021.11-Windows-x86_64.exe
- 'Desktopentwicklung mit C++' in the 'Microsoft Build Tools für C++' from (https://visualstudio.microsoft.com/de/visual-cpp-build-tools/) (default package selection)
- R version 4.1.2 from (https://cran.r-project.org/bin/windows/base/old/)
- RStudio version 2021.9.2.382 from (https://www.rstudio.com/products/rstudio/older-versions/#2021092)

Additional information are provided in HTML file 'MA_Computations_Docu.html'


## Usage 

It is recommended to keep all files in one directory. This is advantageous since script and Notebook automatically identify the storage location as working direcory. 

Jupyter Notebook becomes available by the Anaconda software distribution installation. It is web-based and opens in the brower. One can navigate to the file location and run 'MA_Computational_Part.ipynb'.

The imputation is conducted via 'MA_Imputation.R'. The Notebook writes the dataset 'X.csv' to the directory from which it is loaded via RStudio. The imputed data is exported to the directory as 'X_imp.csv' and then loaded into the notebook.
