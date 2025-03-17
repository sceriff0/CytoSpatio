# CytoSpatio: Learning Cell Type Spatial Relationships
Haoran Chen and Robert F. Murphy\
Carnegie Mellon University\
V0.1.0 Feb 11, 2025

CytoSpatio is designed to decipher the complex spatial relationships between different cell types. Using generative, multirange, multitype point process models, it captures intricate spatial interactions between cell types across various ranges simultaneously. By supplying a simple three-column input of cell coordinates and types, CytoSpatio produces interaction coefficients that delineate both the inherent and apparent spatial relationships among cell types. In addition, it generates synthetic tissue images which preserve the spatial relationships observed in training images.

Reference: Chen, Haoran, and Robert F. Murphy. "CytoSpatio: Learning cell type spatial relationships using multirange, multitype point process models." bioRxiv (2024): 2024-10.


## Requirements

- **R Version**: 3.6.3
- **Required R Packages**:
   `spatstat, spatstat.utils, spatstat.data, ggplot2, dplyr, permute, data.table, igraph, deldir, sf`

  All packages will be automatically installed.

## Usage

### Input Format

Your data should be formatted in a CSV file with the following columns:

- The x-coordinate of the cell.
- The y-coordinate of the cell.
- The cell type of the cell.

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

### Contact

Robert F. Murphy - murphy@cmu.edu\
Haoran Chen - hrchen@cmu.edu

