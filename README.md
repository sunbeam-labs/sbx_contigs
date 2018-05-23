# Sunbeam Contigs extension

This is an extension to select [contigs](https://github.com/sunbeam-labs/sunbeam/blob/dev/rules/assembly/assembly.rules) based on the [blast summary](https://github.com/sunbeam-labs/sunbeam/blob/dev/rules/annotation/annotation.rules) for a give taxon of interest, e.g. Escherichia coli, calculate the coverage by mapping reads back and generate coverage plot.

## Installing

With you sunbeam conda environment activated, simply clone the extension directory into the sunbeam/extensions/ folder, installing requirements, and add new options to you existing configuration file:

  ```bash
   git clone https://github.com/zhaoc1/sbx_contigs.git extensions/sbx_contigs
   conda install --file extensions/sbx_contigs/requirements.txt
   cat extensions/sbx_contigs/config.yml >> sunbeam_config.yml
   ```
 
## Running

  ```bash
  sunbeam run --configfile=sunbeam_config.yml give_my_contigs
  ```


 
