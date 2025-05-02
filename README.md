# Impact of protein conformational diversity on AlphaFold3 predictions

Main data preparation script can be found here: [Google Colab](https://colab.research.google.com/drive/16k_oZTqws-6uCptCGIRy0ykzQ-ltVHtL?usp=sharing)
This repository contains code and scripts to perform root-mean-square deviation (RMSD) analysis and visualize protein structural alignments. It focuses on comparing APO and HOLO forms of proteins, analyzing differences, and generating insightful visualizations using R.

## Features

- **Data Preparation**: Preprocess RMSD datasets.
- **Statistical Tests**: Perform Wilcoxon signed-rank tests to determine significant differences between APO and HOLO conformers.
- **Visualization**:
  - Density plots of RMSD values.
  - Scatter plots comparing APO and HOLO alignments.
  - Box plots grouped by protein families and flexibility.
  - Conformational diversity impact on plDDT score.

## Requirements

- R (≥ 4.0)
- R packages:
  - `ggplot2`
  - `dplyr`
  - `tidyr`
  - `readr`
  - `ggpubr`

Install the necessary packages using:
```r
install.packages(c("ggplot2", "dplyr", "tidyr", "readr", "ggpubr", "gridExtra"))
```

## Usage

1. Clone the repository:
   ```bash
   git clone https://github.com/Komceks/Impact-of-protein-conformational-diversity-on-AlphaFold3-predictions.git
   ```

2. Open the R script `Plots.R` in RStudio or run it in an R environment.

3. Ensure to set working directory.

4. Execute the script to:
   - Preprocess and analyze the data.
   - Generate and save plots to the `Stats/` directory.

5. View saved figures for insights.

## Analysis Sections

### **Figure 0**: Density Plot
- Visualize RMSD distributions between APO and HOLO states.

### **Figure 1A & 1B**: Scatter Plots
- Scatter plots of predicted APO structure similarity to experimental APO vs HOLO (full and zoomed).

### **Figure 2A & 2B**: Scatter Plots
- Scatter plots of predicted HOLO similarity to experimental APO vs HOLO (full and zoomed).

### **Figure 4A & 4B**: Scatter Plot of plDDT vs. RMSD
- Correlate predicted plDDT scores with RMSD values.

### **Figure 5A-D**: Scatter Plot of Maximum and minimum RMSD trends in Protein Clusters

### **Figure 6A & 6B**: Scatter Plot of RMSD by protein families
- Analyze RMSD by protein families (homogeneous vs. heterogeneous).

### **Figure 7A & 7B**: Scatter Plot of RMSD by protein flexibility
- Analyze RMSD by protein flexibility (flexible vs. rigid).

### **Figure 8A & 8B**: Box plot of RMSF vs pLDDT

## Outputs

Generated figures are saved in the `Stats/` directory.
