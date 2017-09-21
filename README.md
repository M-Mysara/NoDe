# NoDe
A fast error-correction algorithm for pyrosequencing amplicon reads
## Mysara, M., Leys, N., et al., 2015. NoDe: a fast error-correction algorithm for pyrosequencing amplicon reads. BMC Bioinformatics, 16(1), p.88.
## Abstract
The popularity of new sequencing technologies has led to an explosion of possible applications, including biodiversity studies. However each of these sequencing technologies suffer from sequencing errors originating from different factors.
For 16S rRNA metagenomics studies the 454 pyrosequencing technology is already for many years one of the most frequently used platforms, but sequencing errors still leads to important data analysis issues (e.g. in OTU clustering and biodiversity estimation). The new error correction algorithm proposed in this work is trained to identify those positions in 454 sequencing reads that are likely to have an error, and subsequently clusters those error-prone reads with correct reads resulting in error-free consensus reads. The NoDe software can be downloaded here:
http://science.sckcen.be/en/Institutes/EHS/MCB/MIC/Bioinformatics/NoDe

## Installation Requirement
Only Perl needs to be installed to be able to run NoDe. All other software packages that are required to run NoDe are included in the downloaded file (NoDe.run).

## Included Software
Following software is used by the NoDe algorithm, however you do NOT need to install it separately as these software modules are included in the NoDe software.

- Mothur:

    Available at http://www.mothur.org/wiki/Download_mothur. Note about mothur: the command called "pre.cluster" is modified to be compatible with NoDe algorithm.
- WEKA 3.6: 
    Available online at http://www.cs.waikato.ac.nz/ml/weka/.

Both WEKA and Mothur are distributed under the GNU licence.
## Syntax
./CATCh.run {options}

### Use the Following Mandatory Options:

!!! make sure you use an underscore (and not a hyphen) to specificy the option you want to set.

- _s The sff.txt file
                Can be produced via the mothur command [sff.info (sff = file.sff, sfftxt = T)]

- _n Name file with the redundancy
Tab separated file with the read ID on the first column and the IDs of identical reads on the second column [separated by comma]).
Can be produced after Dereplication of the reads via Mothur command [unique.seqs(fasta = file.fasta)]

- _f Aligned sequences file
                Aligned multi-fasta files without any ambiguity, only "AGTC" bases.

- _o Output Path             
 Use the Following Non-Mandatory Options:

- _p number of processors, [default 1]
- _g group file (in case of having different sample-groups)
 Tab separated file, with read ID in the first column and read sample (group) name on the second.

- _d differences tolerated
[Default: will be automatically calculated with a cut-off of 98% distance. 
(it is recommended to increase the differences by "1" each additional 100 bp in the read average length).

- To see the program help type NoDe.run (with no parameters).

## Output Files
The NoDe program generates different text output files distributed over two folders "NoDe_Final" and "NoDe_Temp"

Inside each of them another folder titled the same random number representing the process ID

### NoDe_Final:
Contain the final alignment and name file after denoising and preclustering named:

Results.NoDe.names
Results.NoDe.fasta
### NoDe_Temp:
- Splited alignment files (0.fasta, 1.fasta ,etc) and name files (0.names ,1.names, etc) according to the number of processors used.
- NoDe model (classifier) input and output named (0.Test.arff, 0.Test.Final, etc) respectively.
- Splited alignment files with the marked positions (0.Results, 1.Results, etc).
- Merged aligment file (Results.fasta) and names (Results.names)
- SLP-modified output: Results.precluster.fasta, Result.precluster.names
-Final output after removing of the marked positions (Results.precluster.pick.fasta, Results.precluster.pick.names)
logfile
The program produces a log file of the steps being run, where it is possible to monitor the different steps, and trackdown possible errors:
NoDe: OUTPUT_PATH/NoDe_Temp/****/NoDe_****.logfile 
(**** represents a random process ID).

## Testing

Type:

NoDe_Main.run _f ./sample.fasta _n ./sample.names _s ./sample.sfftxt _o OUTPUT_PATH

It will produce in the output path a file containing the results i.e. two files named Sample.NoDe.fasta and Sample.NoDe.names
The Test dataset is a part of the data published in (Schloss et al, 2012)  and is part of Project SRP002397 in NCBI short archive.

-------

## Citing
If you are going to use NoDe, please cite it with the included software (Mothur, WEKA):

    Mysara, M., Leys, N., et al., 2015. NoDe: a fast error-correction algorithm for pyrosequencing amplicon reads. BMC Bioinformatics, 16(1), p.88.
    Schloss PD, Westcott SL, Ryabin T, Hall JR, Hartmann M, Hollister EB, et al. (2009). Introducing mothur: open-source, platform-independent, community-supported software for describing and comparing microbial communities. Applied and environmental microbiology 75:7537–41.
    Hall M, National H, Frank E, Holmes G, Pfahringer B, Reutemann P, et al. (2009). The WEKA Data Mining Software?: An Update. SIGKDD Explorations 11:10–18.
