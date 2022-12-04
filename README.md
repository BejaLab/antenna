Carotenoid antenna workflow
===========================

This is the [snakemake](https://snakemake.readthedocs.io/) workflow behind the paper [Chazan et al (2022)](https://www.biorxiv.org/content/10.1101/2022.08.24.505090v1) "Energy transfer in ubiquitous rhodopsin pumps with xanthophyll antennas". For supplementary data, see the [Figshare repository](https://doi.org/10.6084/m9.figshare.20502384).

Run as: `snakemake -c{number of cores} --use-conda` (the dependencies are under the conda control). Check the file `workflow/Snakefile` for specific rules.
