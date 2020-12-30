import sys
sys.path.append("../../cooked_up_libraries")

import guide_rna_library
# When running this file change the following attributes.

# chopchop/config.json

# "PRIMER3": "../../chopchop/primer3_core",
# "BOWTIE": "../../chopchop/bowtie/bowtie",
# "TWOBITTOFA": "../../chopchop/twoBitToFa",
# "TWOBIT_INDEX_DIR": "../../CHOPCHOP_Necessary_Genetic_Data/2Bit_Genes",
# "BOWTIE_INDEX_DIR": "../../CHOPCHOP_Necessary_Genetic_Data/EBWT_Genes",
# "GENE_TABLE_INDEX_DIR": "../../CHOPCHOP_Necessary_Genetic_Data/Gene_Tables"

# When running it back under the direct parent folder of gBlock Generation change it back
# to this configuration

# chopchop/config.json

# "PRIMER3": "chopchop/primer3_core",
# "BOWTIE": "chopchop/bowtie/bowtie",
# "TWOBITTOFA": "chopchop/twoBitToFa",
# "TWOBIT_INDEX_DIR": "CHOPCHOP_Necessary_Genetic_Data/2Bit_Genes",
# "BOWTIE_INDEX_DIR": "CHOPCHOP_Necessary_Genetic_Data/EBWT_Genes",
# "GENE_TABLE_INDEX_DIR": "CHOPCHOP_Necessary_Genetic_Data/Gene_Tables"

organism = "sacCer3"
gene = "YDR210W"
guide_rna_results_storage_folder = "../../chopchop/temp"
selenium_driver = "../../chromedriver/chromedriver"

# obtainGuideRNAsOnlineCHOPCHOP has shown to work!
# guide_rna_library.obtainGuideRNAsOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, selenium_driver = selenium_driver)

# obtainGuideRNAsAncientRunCHOPCHOP has shown to work!
# print(guide_rna_library.obtainGuideRNAsAncientRunCHOPCHOP(organism, gene, guide_rna_results_storage_folder))

print(guide_rna_library.obtainAgeOfSavedGuideRNAsinDays(organism, gene, guide_rna_results_storage_folder))
