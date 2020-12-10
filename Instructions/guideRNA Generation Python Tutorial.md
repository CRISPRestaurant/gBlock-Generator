# Determining gRNA Candidates using CHOPCHOP and Python 2.7

#### WARNING: This tutorial is mainly for the Saccharomyces Cerevisiae species. Nuances shall apply for other species and their respective configuration files.

1. If you were like me and panicked a little pondering why Python 2.7 was still relevant, [here's a few lines](https://mothergeo-py.readthedocs.io/en/latest/development/how-to/venv-win.html) to setup a Python 2 virtual environment.

```
pip install virtualenv
cd my-project
virtualenv --python C:\Path\To\Python2\python.exe python2-venv
cd python2-venv\Scripts
activate
```
2. Download our more structured [CHOPCHOP-containing folder](https://github.com/CRISPRestaurant/guideRNA).

3. Obtain the ```.gene_table``` file for your organism by visiting [this site](http://genome.ucsc.edu/cgi-bin/hgTables?hgsid=859909031_IKpLA1AwyQ5WqaEjCJXPycmLjRkz&clade=other&org=0&db=0&hgta_group=genes&hgta_track=refGene&hgta_table=refFlat&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=), selecting the CHOPCHOP API's recommended inputs below, and clicking "get output".

Option | Recommended Input
--- | ---
group | Gene and Gene Predictions
track | SGD Genes
table | sgdGene
region | genome
output format | all fields from selected table

Option | Remember Input for Later
--- | ---
assembly | Month and Year

4. Right click the new page, click "Save Page As", navigate to wherever your ```Gene_Table``` [folder](https://github.com/CRISPRestaurant/guideRNA/tree/master/Gene_Tables) is,and type ```<your organism name>.gene_table``` into the Filename input, and save.

5. Obtain the ```.2bit``` file for your organism by visiting [this site](https://hgdownload.soe.ucsc.edu/downloads.html) and looking for your species and the exact same assembly month and year as what you had inputed for obtaining the ```.gene_table```. Then navigate to ```Annotations``` and then ```Fileserver (bigBed, maf, fa, etc) annotations Also see Data Access```.

6. Locate the ```.2bit``` file and save it into the ```2Bit_Genes``` [folder](https://github.com/CRISPRestaurant/guideRNA/tree/master/2Bit_Genes).

7. Obtain your EBWT files in two manners depending on if you are running Windows or Linux.

For Linux run these set of commands.

```
cd chopchop
./twoBitToFa ../2Bit_Genes/<your organism name>.2bit ../FASTA_Genes/<your organism name>.fa
./bowtie/bowtie-build ../FASTA_Genes/<your organism name>.fa <your organism name>
```
Move all your .ebwt files into ```EBWT_Genes```

For Windows go [here](https://chopchop.cbu.uib.no/genomes/) and download all .ebwt files concerning your organism into the ```EBWT_Genes``` [folder](https://github.com/CRISPRestaurant/guideRNA/tree/master/EBWT_Genes)

8. Run the CHOPCHOP program to design the gRNAs by running the command below.

```
./chopchop/chopchop.py -G <your organism name> -o <your gRNA output directory> -Target <your gene name>
```

```<your organism name>``` must match the name of the ```.gene_table``` file

```<your gene name>``` must belong to the ```name``` column inside the ```.gene_table``` file

9. Import any modules that are required and keep importing until the program spits out a list of gRNAs.

For biopython, since version 1.75 is the latest version Python 2 will accept, you must install this version.

```pip install biopython==1.75```
