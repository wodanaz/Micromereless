MS ID#: DEVELOP/2024/203152
MS TITLE: Reprogramming of cells during embryonic transfating: overcoming a reprogramming block
AUTHORS: Alejandro Berrio, Esther Miranda, Abdull J Massri, Anton Afanassiev,
         Geoffrey Schiebinger, Gregory A Wray, and David R. McClay

This README contains the instructions to replicate each figure in the above manuscript.
Created by: Alejandro Berrio
Date: 11/06/2024
Usage:

Ensure all required packages are installed.
Follow the instructions provided in the following file for running the scripts.

License:
This work is licensed under the Creative Commons Attribution-NonCommercial 4.0 International License.
To view a copy of this license, visit http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to
Creative Commons, PO Box 1866, Mountain View, CA 94042, USA.

Contact:
For any questions or issues, please contact alejo.berrio@duke.edu


## Description
This repository contains data and R scripts used for the analysis of Waddington data. The main script, `Figures_and_analysis.R`, performs data processing and generates figures for the study.

## Files 
- `precomputed_control.Rda`: Precomputed control data. 
- `precomputed_micromereless.Rda`: Precomputed micromereless data.
- `precomputed_integrated.Rda`: Precomputed integrated data.
- `mmless.zip`: Compressed file containing micromereless data.
- `control.zip`: Compressed file containing control data.
- `Figures_and_analysis.R`: R script for data analysis and figure generation.

## Instructions

### Download and Extract Files
1. **Download the data**: Download the `*.Rda` ,  `mmless.zip` and `control.zip` files from Zenodo (10.5281/zenodo.14051567)
2. **Extract the data**: Extract the contents of `mmless.zip` and `control.zip` to your working directory.

### Running the R Script
1. **Open R or RStudio**: Ensure you have R or RStudio installed on your computer.
2. **Set the working directory**: Set your working directory to the location where you extracted the data. You can do this in R using the `setwd()` function. For example:
    ```r
    setwd("path/to/your/directory")
    ```
3. **Install required packages**: Make sure you have the necessary R packages installed. You can install them using:
    ```r
    install.packages(c("tidyverse", "data.table", "Seurat", "viridis", "xlsx", "Matrix", "scCustomize", "cowplot, "pheatmap"))

    # To Install CIDER, follow the instructions from https://github.com/zhiyhu/CIDER
    
    ```
4. **Run the script**: Source the `Figures_and_analysis.R` script to perform the analysis and generate figures. You can do this by running each step in:
    ```r
    "Figures_and_analysis.R"
    ```

## Notes
- Ensure that the paths in the `Figures_and_analysis.R` script are correctly set to the location of your data files.
- If you encounter any issues, please check the script for any hardcoded paths and update them accordingly.

## Contact
For any questions or issues, please contact Alejandro Berrio at alejo.berrio@duke.edu
