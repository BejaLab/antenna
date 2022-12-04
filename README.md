Carotenoid antenna workflow
===========================

This is the [snakemake](https://snakemake.readthedocs.io/) workflow behind the paper [Chazan et al (2022)](https://www.biorxiv.org/content/10.1101/2022.08.24.505090v1) "Energy transfer in ubiquitous rhodopsin pumps with xanthophyll antennas". For supplementary data, see the [Figshare repository](https://doi.org/10.6084/m9.figshare.20502384).

Run as: `snakemake -c{number of cores} --use-conda` (the dependencies are under the conda control). Check the file `workflow/Snakefile` for specific rules.

To run the pipeline you have to download the heavy input data:

```shell
# Create direcotries
mkdir -p data/gem/split data/jgi data/om-rgc
#
# Download GEM data
wget https://portal.nersc.gov/GEM/genomes/faa.tar
tar xf faa.tar
find faa -name '*.faa.gz' | xargs zcat | seqkit split -s 1000000 -O data/gem/split/
#
# Download OM-RGC data
wget -P data/om-rgc/ https://www.ebi.ac.uk/biostudies/files/S-BSST297/Files/OM-RGC_v2.tsv.gz
wget -P data/om-rgc/ https://www.ebi.ac.uk/biostudies/files/S-BSST297/Files/OM-RGC_v2_gene_profile_metaG.tsv.gz
gunzip data/om-rgc/*.tsv.gz
#
# Download JGI rhodopsins
wget -O Data_jgi.tar.gz https://figshare.com/ndownloader/files/38419106
tar xfz Data_jgi.tar.gz
```
