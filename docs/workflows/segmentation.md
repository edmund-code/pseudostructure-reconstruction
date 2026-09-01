# Panoptic kidney segmentation

The `segmentation/` project predicts semantic classes plus local interior/boundary maps, then
uses marker-controlled watershed and tiled stitching to make instances. It intentionally has no
centroid-HV output; see the architecture decisions in `docs/architecture/`.

Use a dedicated environment:

```bash
conda env create -f segmentation/environment.yml
conda activate kidney-panoptic
cd segmentation
pytest tests
```

For data-backed work, copy `data/manifest_v4.example.csv` outside the repository, replace each
placeholder with an authorized absolute or mounted path, then invoke the scripts with a config
that points to that private manifest. Whole-slide inference also requires a separately obtained
v4 `best_pq.pt` checkpoint.

Default tests are CPU/unit-level checks. Training, Hugging Face encoder download, and whole-slide
inference are deliberately not CI jobs.
