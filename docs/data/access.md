# Private data access and setup

The repository intentionally does not distribute study data, source slides, annotations, or model
weights. Obtain access through the project maintainer, then create a data-root directory matching
[the data contract](../../data/README.md).

Before starting a full run, verify:

1. The four mouse H5AD matrices and pathway JSON libraries are present.
2. The four mouse segmentation files are exactly the v4 GeoJSONs named in `data/catalog.csv`.
3. The requested result root is writable and outside the Git checkout when possible.
4. The Harmony R packages are installed in the analysis environment.
5. Segmentation slide paths and annotation paths are supplied through a private manifest derived
   from `segmentation/data/manifest_v4.example.csv`.

Never add any of these artifacts to Git. The hygiene check rejects them by extension and size.
