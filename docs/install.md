---
label: Install
icon: desktop-download
order: 100
---

# Install HARPY

Until this pipeline gets completed and hosted on Bioconda, it will be available by cloning/downloading [the repository](https://github.com/pdimens/harpy). The dependencies can be installed into a conda environment using the provided `harpyenv.yaml`:
```bash
conda env create --name harpy --file misc/harpyenv.yaml
```

The version of [EMA](https://github.com/arshajii/ema) bundled in this repository (`ema-h`) is a [fork](https://github.com/EdHarry/ema/tree/haplotag) of the orignal EMA modified to work with Generation 1 haplotag beadtags (AxxCxxBxxDxx).