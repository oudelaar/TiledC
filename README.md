# TiledC

Overview

The most straightforward way to analyse Tiled-C data is to use the HiCPro pipeline (Servant et al. Genome Biol 2015. https://github.com/nservant/HiC-Pro) with the options for Capture-Hi-C analysis. 

It is also possible to analyse Tiled-C data using a combination of the CCseqBasic pipeline, custom scripts and ICE normalisation. Instructions for this approach are described below. 

If you have enriched for multiple regions, you can analyse these tiles simultaneously in Step 1, but will need to do every region separate in the following steps. 

Step 1: CCseqBasic CM5

This step uses the CCseqBasiq pipeline generated by Jelena Telenius to analyse the fastq files and report bam files containing reads in which at least one interaction with a restriction fragment within the tile(s) is detected. These reads are filtered for ploidy and PCR duplication related artefacts.  

Information about the CCseqBasiq pipeline can be found here:
https://github.com/Hughes-Genome-Group/CCseqBasicF/releases
http://userweb.molbiol.ox.ac.uk/public/telenius/CCseqBasicManual/index

Tiled-C analysis is supported in the CM5 version (https://github.com/Hughes-Genome-Group/CCseqBasicM), of which scripts are available upon request. 

Brief instructions:
- Use the flag --tiled
- Use an oligo file reporting the targeted region(s) of interest with no proximity exclusion; e.g. TiledC_oligofile_alpha_mouse.txt 
- Skip blat filtering using a dummy folder

Step 2: Generate raw matrix

A) Convert bam file generated by pipe to sam file:

samtools view -h /path_to_F6_folder/name_COMBINED_reported_capture_reads_CM5.bam > /path/yourname.sam &

B) Convert sam file to raw sparse contact matrix:

Use script TiledC_sam2rawmatrix.pl. Instructions can be found in the comment lines in the script. This script also outputs a bed file that specifies the coordinates of each bin number. 

Example:

nohup perl /path/TiledC_sam2rawmatrix.pl -sam /path/name_COMBINED_reported_capture_reads_CM5.sam -c 11 -s 29902951 -e 33226736 -r /path/mm9_dpnII_coordinates.txt -n name -b 2000 &

Step 3: ICE normalisation (Imakaev et al. Nature Methods, 2012)

The TiledC_sam2rawmatrix.pl creates a matrix in a specific directory hierarchy that is compatible with the HiC-Pro pipeline. You can therefore use HiC-Pro to perform ICE normalisation.

HiC-Pro documentation:
https://github.com/nservant/HiC-Pro

Briefly:
- Edit the HiC-Pro configuration file. As we’re just using the ICE normalisation, you only need to modify the settings for #Contact Maps and #ICE Normalization at the bottom of the file. Make sure the bin_size matches the bin sizes in your matrix folder. You can specify multiple bin sizes (if you ran the TiledC_sam2rawmatrix.pl multiple times with different bin sizes) separated by a space.
- Run command with -s ice_norm flag to specify you only want to run ICE normalisation and point the -i flag to the path of the ‘matrix’ directory generated by the Tiled_sam2rawmatrix_MO.pl script.

Step 4: Visualisation

The normalised matrices can be visualised with a simple python script, e.g. TiledC_matrix_visualisation.py 