# Raw Data

This folder is intended to contain the original, unprocessed datasets used in this project.

## Data Source

The data used in this analysis originates from the Environment Agency (EA) open data on freshwater macroinvertebrate monitoring in England. These datasets include long-term ecological observations collected across a national network of river sites.

Key datasets used:

* Site metadata (locations, catchments, coordinates)
* Sampling records (dates, methods, analysis types)
* WHPT family-level abundance data
* TAXA-level abundance data

## Data Availability

The full raw datasets are **not included in this repository** due to:

* Large file sizes
* External data ownership and licensing considerations

Instead, this project is designed to be **fully reproducible** using the original data sources.

## How to Reproduce

To run this project:

1. Download the required datasets from the official Environment Agency or relevant data portals.

2. Place the `.parquet` files in this folder (`data/raw/`).

3. Ensure filenames match those expected in the scripts:

   * `INV_OPEN_DATA_SITE.parquet`
   * `INV_OPEN_DATA_METRICS.parquet`
   * `R_INV_WHPT_METRICS_B.parquet`
   * `INV_OPEN_DATA_TAXA.parquet`

4. Run the scripts in the `scripts/` folder in order.

## Notes

* Data is loaded using the `arrow` package for efficient handling of large datasets.
* The pipeline is designed to work with **lazy loading**, meaning data is only read into memory when required.
* No modifications are made to raw data files; all transformations occur in downstream scripts.

## Ethical and Practical Considerations

* This project works with environmental monitoring data intended for public good and research purposes.
* Care has been taken to preserve the integrity of the original data during processing.
* Users should ensure compliance with the original data provider’s terms of use when accessing and using the datasets.

