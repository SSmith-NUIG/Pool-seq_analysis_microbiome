# Pool-seq_analysis_microbiome
Microbiome pipeline to analyse Pool-seq data from honey bees

## Step one
Run ```Alignment.sh``` to align your raw fastqs to the honey bee reference genome and exracts 
reads which did not align. These are the microbiome reads we will use for the analysis

## Step two

Run ```Kaiju.sh``` to use Kaiju for taxonomic classification of the reads.

## Step three

Run ```kaiju2phyloseq.py``` to transform the individual kaiju output files into one phyloseq object.

## Step four

Run ```Pyloseq.R``` to analyse the pyloseq object and create all plots using in microbiome chapter of thesis.
