# guide_rna_library Library
[Source code](https://github.com/CRISPRestaurant/gBlock-Generator/blob/master/cooked_up_libraries/guide_rna_library.py) of library

## Methods

__Obtain guideRNA sequences using [CHOPCHOP's API folder](https://github.com/CRISPRestaurant/gBlock-Generator/tree/master/chopchop)__
**WARNING: This folder requires [extensive setup](https://github.com/CRISPRestaurant/gBlock-Generator/blob/master/Instructions/guideRNA%20Generation%20Python%20Tutorial.md) and can only be done on Linux**

```python
guideRNAs = obtainGuideRNAsOfflineCHOPCHOP(organism_name, target_gene, guide_rna_results_storage_folder, chopchop_dir)
# organism_name is the name of the organism consistent with the file name for the organism's genetic data files
# target_gene is the name of the gene present in the organism's gene table file
# guide_rna_results_storage_folder is this library's specific folder that archives results in a manner explained below
# chopchop_dir is the path to the chopchop directory that is used to run the chopchop/chopchop.py
```

__Obtain guideRNA sequences using [CHOPCHOP's website](https://chopchop.cbu.uib.no/)__

```python
guideRNAs = obtainGuideRNAsOnlineCHOPCHOP(organism_name, target_gene, guide_rna_results_storage_folder, time_allotment = 120, selenium_driver = "chromedriver/chromedriver")
# organism_name is the name of the organism consistent with the file name for the organism's genetic data files
# target_gene is the name of the gene present in the organism's gene table file
# guide_rna_results_storage_folder is this library's specific folder that archives results in a manner explained below
# time_allotment is the time in seconds the method will wait before it quits if a web element does not respond
# selenium_driver is the location of the chromedriver
```

__Obtain guideRNA sequences using archived run (example of consistent archived folder)__

```python
guideRNAs = obtainGuideRNAsAncientRunCHOPCHOP(organism, gene, guide_rna_results_storage_folder)
# organism is the name of the organism consistent with the file name for the organism's genetic data files
# gene is the name of the gene present in the organism's gene table file
# guide_rna_results_storage_folder is this library's specific folder that this method will retrieve guideRNA data from

# If the archived data cannot be found in the guide_rna_results_storage_folder, False will be returned
```

__Obtain age (in days) of guideRNA archived results__

```python
guideRNAs = obtainAgeOfSavedGuideRNAsinDays(organism, gene, guide_rna_results_storage_folder)
# organism is the name of the organism consistent with the file name for the organism's genetic data files
# gene is the name of the gene present in the organism's gene table file
# guide_rna_results_storage_folder is this library's specific folder that this method will retrieve guideRNA data from

# If the archived data cannot be found in the guide_rna_results_storage_folder, False will be returned
```

__Obtain archived data if it exists and is young enough or else obtain guideRNA data using [CHOPCHOP's website](https://chopchop.cbu.uib.no/)__

```python
guideRNAs = obtainGuideRNAsSavedOrOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, time_allotment = 120, selenium_driver = "chromedriver/chromedriver", acceptable_age = sys.maxsize)
# organism is the name of the organism consistent with the file name for the organism's genetic data files
# gene is the name of the gene present in the organism's gene table file
# guide_rna_results_storage_folder is this library's specific folder that this method will retrieve guideRNA data from
# time_allotment is the time in seconds the method will wait before it quits if a web element does not respond
# selenium_driver is the location of the chromedriver
# acceptable_age is the maximum age (in days) the archived guideRNA data can be. If the archived data is older, the method will fetch the guideRNA data online, resetting the age of the archived data to 0 days
```

## guide_rna_results_storage_folder structure (best not to touch this)

```
guide_rna_results_storage_folder
|
|___Organism #1
|   |
|   |___Gene #1
|   |   |   guideRNAs.csv
|   |   |   Knockout Output Summary.txt
|   |   |___Off Targets
|   |          1.offtargets
|   |          2.offtargets
|   |          ...
|   |
|   |___Gene #2
|
|___Organism #2
    |
    |___Gene #3
    |
    |___Gene #4
```
