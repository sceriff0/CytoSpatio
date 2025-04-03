# CytoSpatio: Learning Cell Type Spatial Relationships
Haoran Chen and Robert F. Murphy\
Carnegie Mellon University\
V0.2.0 April 3, 2025

CytoSpatio is designed to decipher the complex spatial relationships between different cell types. Using generative, multirange, multitype point process models, it captures intricate spatial interactions between cell types across various ranges simultaneously. By supplying a simple three-column input of cell coordinates and types, CytoSpatio produces interaction coefficients that delineate both the inherent and apparent spatial relationships among cell types. In addition, it generates synthetic tissue images which preserve the spatial relationships observed in training images.

Reference: Chen, Haoran, and Robert F. Murphy. "CytoSpatio: Learning cell type spatial relationships using multirange, multitype point process models." bioRxiv (2024): 2024-10.


## Requirements

- **R Version**: 3.6.3 or later
- **Required R Packages**:
   `spatstat, spatstat.utils, spatstat.data, ggplot2, dplyr, permute, data.table, igraph, deldir, sf`

  All packages will be automatically installed.

## Usage

### Input Format

Your data should be formatted in a CSV file with the following columns:

- The x-coordinate of the cell.
- The y-coordinate of the cell.
- The cell type of the cell.

An example input file is contained in the "example" folder.

### Parameters

- **TR (Total Range)**: The maximum range within which two cells can interact. Default is 500 pixels.
- **IR (Interval Range)**: The distance of each range with same interaction between two cell types. Default is 100 pixels.
- **HR (Hardcore Range)**: The smallest allowable distance between two cells to maintain spatial integrity. Cells cannot be closer than this distance. Default: 1 pixel.

### Command

From your R terminal:

```R
source("cytospatio.R")
cytospatio(input_file = "path_to_input_file", output_dir = "path_to_output_dir", TR, IR, HR)
```
## Examples

### Example input file

There is an example input file "cell_data.csv" in the "example" folder.  It contains data for 43,900 cells.

### Example R script

There is also an example R script to run CytoSpatio using this example file.  It writes its output to an "output" folder.

From your terminal:

```Rscript exampleCytoSpatioscript.R
```

### Example Jupyter Notebook

The CytoSpatioNotebook.ipynb notebook allows specifying input and output files as well as ranges for model construction.  It runs CytoSpatio and then displays some of the outputs.

From your terminal, go to the folder containing the notebook and enter
```jupyter notebook
```

Default values are contained in the notebook, so you can run it by selected "Run all" from the "Cell" menu

### Contact

Robert F. Murphy - murphy@cmu.edu\
Haoran Chen - haoranch@cs.cmu.edu

