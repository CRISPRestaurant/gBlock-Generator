import os
from math import *
from random import *
import primer3
import openpyxl as pyxl
import pandas as pd

# libraries cooked up for greater organization and specialized potential in the future
import cooked_up_libraries.dna_properties as dna_properties
import cooked_up_libraries.gene_table as gene_table
import cooked_up_libraries.restriction_enzymes as restriction_enzymes
import cooked_up_libraries.guide_rna_library as guide_rna_library
import cooked_up_libraries.fasta_library as fasta_library

# Ensures all the previous three filters were passed for a given sequence
def firstgBlockFilterOutcome(sequence):
    return dna_properties.noBsaIRestrictionSites(sequence) and dna_properties.noHomopolymerMoreThanFourTs(sequence) and dna_properties.noHomopolymerMoreThanFiveAsAndGs(sequence)

#### Assign Priority System for Guide RNAs ####

# Assigns a higher priority for guides that are earlier
# in the ORF
def obtainEarlierORFPriority(position):
    return (1) / log(position)

# Assigns a higher priority to guides that do not have off-target
# sites
def obtainOffTargetSitePriority(off_target_output_folder, output_number):
    file = open("%s/%s.offtargets" % (off_target_output_folder, output_number), "r")

    first_line = True

    for line in file:
        if first_line:
            first_line = False
        else:
            if "There are no predicted off-targets." in line:
                return 1
            else:
                split_list = line.split(";")
                off_targets = len(split_list)

                return e ** (-1 * off_targets)

# Prioritizes a list of guideRNAs based on the last four
# priority metrics

# sequences is a list in the form of [[first sequence number, first sequence nucleotide ordering, first sequence full_chromosome
# and position, first sequence sense or anti-sense, . . . more metadata . . .]]

# output_folder is where temp files regarding offtarget sites for the guideRNA will find its home
def prioritize_guide_RNAs(sequences, off_target_output_folder):
    # priority_dict is contained as {guideRNA: priority_value}
    priority_dict = {}
    sequence_properties = {}
    sequenced_list = []

    for sequence in sequences:
        chromosome_and_position = sequence[2].split(":")
        chromosome = chromosome_and_position[0]
        position = int(chromosome_and_position[1])

        gc_priority = dna_properties.obtainGCPriority(sequence[1])
        last_purine_priority = dna_properties.obtainLastPurinePriority(sequence[1])
        earlier_orf_priority = obtainEarlierORFPriority(position)
        off_target_site_priority = obtainOffTargetSitePriority(off_target_output_folder, sequence[0])

        priority_value = gc_priority * last_purine_priority * earlier_orf_priority * off_target_site_priority

        priority_dict[sequence[1]] = priority_value
        sequence_properties[sequence[1]] = sequence

    # sequenced_list datatype is [(guideRNA nucleotide sequence, priority number)] and is sorted from highest priority
    # to lowest priority
    sequenced_list = sorted(priority_dict.items(), key = lambda x: x[1], reverse = True)

    # We don't need priority number in the output so we're extracting the nucleotide sequence from sequenced_list
    output_sequenced_list = []
    for sequence_and_priority in sequenced_list:
        output_sequenced_list.append(sequence_properties[sequence_and_priority[0]])

    return output_sequenced_list

# END OF METHODS THAT HELP OBTAIN THE PRIORITY OF GUIDERNAS

# Obtain a series of primers from a target strand
def obtainPrimers(target, left_primer_length, right_primer_start, right_primer_length, location):
    raw_results = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': target,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[0, left_primer_length, right_primer_start, right_primer_length]]
        },
        {
            'PRIMER_OPT_SIZE': 20,
            # 'PRIMER_PICK_INTERNAL_OLIGO': 1,
            # 'PRIMER_INTERNAL_MAX_SELF_END': 8,
            # 'PRIMER_MIN_SIZE': 18,
            # 'PRIMER_MAX_SIZE': 25,
            # 'PRIMER_OPT_TM': 60.0,
            'PRIMER_MIN_TM': 55.0,
            'PRIMER_MAX_TM': 63.0,
            'PRIMER_MIN_GC': 20.0,
            'PRIMER_MAX_GC': 80.0,
            # 'PRIMER_MAX_POLY_X': 100,
            # 'PRIMER_INTERNAL_MAX_POLY_X': 100,
            # 'PRIMER_SALT_MONOVALENT': 50.0,
            # 'PRIMER_DNA_CONC': 50.0,
            # 'PRIMER_MAX_NS_ACCEPTED': 0,
            # 'PRIMER_MAX_SELF_ANY': 12,
            # 'PRIMER_MAX_SELF_END': 8,
            # 'PRIMER_PAIR_MAX_COMPL_ANY': 12,
            # 'PRIMER_PAIR_MAX_COMPL_END': 8,
            # 'PRIMER_PRODUCT_SIZE_RANGE': [[75,100],[100,125],[125,150],
            #                               [150,175],[175,200],[200,225]],
        }
    )

    # Saves explicity primer3 output data to analyze errors
    df = pd.DataFrame(list(raw_results.items()), columns = ["Key", "Value"])
    df.to_csv(location, index = False)

    num_left_primers = raw_results["PRIMER_LEFT_NUM_RETURNED"]
    num_right_primers = raw_results["PRIMER_RIGHT_NUM_RETURNED"]

    min_primers = min(num_left_primers, num_right_primers)

    output = []

    for i in range(min_primers):
        # Obtain left primer sequence
        left_primer = raw_results["PRIMER_LEFT_" + str(i) + "_SEQUENCE"]

        # Obtain left primer melting temperature
        left_primer_tm = raw_results["PRIMER_LEFT_" + str(i) + "_TM"]

        # Obtain left primer position
        left_primer_position = raw_results["PRIMER_LEFT_" + str(i)][0]

        # Obtain right primer sequence
        right_primer = raw_results["PRIMER_RIGHT_" + str(i) + "_SEQUENCE"]

        # Obtain right primer melting temperature
        right_primer_tm = raw_results["PRIMER_RIGHT_" + str(i) + "_TM"]

        # Obtain right primer position
        right_primer_position = raw_results["PRIMER_RIGHT_" + str(i)][0]

        output.append([[left_primer, left_primer_tm, left_primer_position], [right_primer, right_primer_tm, right_primer_position]])

    return output

# Obtains HI-CRISPR gBlock that does not contain the BsaI endpoints and only contains the
# left homology arm, stop codon, SbfI, right homology arm, and guideRNA
# Also obtains
def obtaingBlockCandidatesWithoutBsaISite(organism, gene, gene_table_for_organism, guide_rna_results_storage_folder, chopchop_dir, online):
    output = []

    # Very crude processing of primer information
    primer_target_round_one = {}
    # Very refined version with primer sequence and melting temperature
    primer_target_final_round = {}

    SbfI = "CCTGCAGG"

    guide_rnas = []
    if online:
        guide_rnas = guide_rna_library.obtainGuideRNAsSavedOrOnlineCHOPCHOP(organism, gene, guide_rna_results_storage_folder, time_allotment = 6000)
    else:
        guide_rnas = guide_rna_library.obtainGuideRNAsOfflineCHOPCHOP(organism, gene, guide_rna_results_storage_folder)

    organism_fasta = fasta_library.processOrganismFasta(organism)

    first_filter_passed_guideRNAs = []
    guide_rna_to_assembled_sequence_dict = {}

    ## This section processes the guideRNAs and prioritizes
    ## them based on a variety of metrics
    for guide_rna in guide_rnas:
        if guide_rna[3] == "+" and float(guide_rna[10]) > 0:
            full_chromosome = guide_rna[2].split(":")

            chromosome = full_chromosome[0]

            # CHOPCHOP's output position does not follow computer science conventions of 0
            # being the first character in a string
            position = int(full_chromosome[1]) - 1

            left_hr_start_position = position - 50
            left_hr_end_position = position - 1

            left_hr = organism_fasta[chromosome][left_hr_start_position: left_hr_end_position]
            guide_plus_pam = guide_rna[1]
            guide = guide_plus_pam[0: (len(guide_plus_pam) - 3)]
            one_extra = organism_fasta[chromosome][position + len(guide_plus_pam)]

            right_hr_start_position = position + len(guide_plus_pam) + 1
            right_hr_end_position = position + len(guide_plus_pam) + 51

            right_hr = organism_fasta[chromosome][right_hr_start_position: right_hr_end_position]

            # Obtaining string to obtain primers
            start_position_first_for_primer = position - 50
            end_position_first_for_primer = position + len(guide_plus_pam) + 50

            length_of_arm = len(guide_plus_pam) + 101
            length_for_ends = 1100 - length_of_arm

            potential_left_end_length = int(length_for_ends / 4)
            potential_right_end_length = length_for_ends - potential_left_end_length

            length_of_chromosome = len(organism_fasta[chromosome])

            availability_at_the_beginning = start_position_first_for_primer
            availability_at_the_end = length_of_chromosome - 1 - end_position_first_for_primer

            if availability_at_the_beginning < potential_left_end_length:
                add_to_right_end = potential_left_end_length - availability_at_the_beginning
                potential_left_end_length = availability_at_the_beginning
                potential_right_end_length += add_to_right_end

            if availability_at_the_end < potential_right_end_length:
                add_to_left_end = potential_right_end_length - availability_at_the_end
                potential_right_end_length = availability_at_the_end
                potential_left_end_length += add_to_left_end

            # Information for primer deduction formatted in
            # [[start of amplicon, end of left homology arm], [start of right homology arm, end of amplicon], chromosome number, length of available place to find left primers, length of available place to find right primers]
            primer_target = [[position - 50 - potential_left_end_length, position], [position + len(guide_plus_pam) + 1, position + len(guide_plus_pam) + 51 + potential_right_end_length], chromosome, potential_left_end_length, potential_right_end_length]

            # End Obtaining string to obtain primers

            first_filter_block = left_hr + SbfI + right_hr + guide

            if firstgBlockFilterOutcome(first_filter_block):
                first_filter_passed_guideRNAs.append(guide_rna)
                primer_target_round_one[guide_rna[1]] = primer_target
                guide_rna_to_assembled_sequence_dict[guide_rna[1]] = [left_hr, right_hr]


    ranked_guideRNAs = prioritize_guide_RNAs(first_filter_passed_guideRNAs, "%s/%s/%s/Off Targets" % (guide_rna_results_storage_folder, organism, gene))
    ## END SECTION

    # Chooses four highest ranked guideRNAs
    for i in range(0, 4):
        if i < len(ranked_guideRNAs):
            guide_rna = ranked_guideRNAs[i]
            position = int(guide_rna[2].split(":")[1])
            exon_ranges = gene_table_for_organism[gene]["exonRanges"]

            # Simple procedure to find exon start position that the guideRNA is residing in
            exon_start_position = 0

            for exon_range in exon_ranges:
                if position >= int(exon_range[0]) and position <= int(exon_range[1]):
                    exon_start_position = int(exon_range[0])
                    break

            difference = position - exon_start_position
            remainder = difference % 3

            prefix = ""

            if remainder != 0:
                prefix = dna_properties.generateRandomSequence(3 - remainder)

            homology_arms = guide_rna_to_assembled_sequence_dict[guide_rna[1]]

            # build variable is easy recognition and analysis by human observers
            build = homology_arms[0] + "|" + prefix + "|" + "TAA" + "|" + SbfI + "|" + homology_arms[1] + "|||" + guide_rna[1]
            # build real is the real sequence that should be placed into order forms and et cetera
            build_real = homology_arms[0] + prefix + SbfI + homology_arms[1]
            output.append(build)

            primer_build_info = primer_target_round_one[guide_rna[1]]

            amplicon = organism_fasta[primer_build_info[2]][primer_build_info[0][0]: primer_build_info[0][1]] + prefix + "TAA" + "CCTGCAGG" + organism_fasta[chromosome][primer_build_info[1][0]: primer_build_info[1][1]]
            left_primer_search_space_length = primer_build_info[3]

            right_primer_start_position_search = primer_build_info[3] + len(prefix) + 8 + 100 + 100
            right_primer_search_space_length = primer_build_info[4] - 100

            location_of_detailed_primer_file = "Test Output Garbage Dump/Primer Results/Detailed Primer Results for Gene " + gene + " and guideRNA " + guide_rna[1] + ".csv"
            primers = obtainPrimers(amplicon, left_primer_search_space_length, right_primer_start_position_search, right_primer_search_space_length, location_of_detailed_primer_file)
            primer_target_final_round[build] = [amplicon, primers]

    return [output, primer_target_final_round]

def obtaingBlocks(organism, gene_list, plasmid, gene_table_for_organism, guide_rna_results_storage_folder, chopchop_dir, online):
    # Obtains two sets of four nucleotides needed to glue the ends of the gBlock onto the plasmid
    start_end_BsaI_sites = restriction_enzymes.findBsaISites(plasmid)
    # Pair of start and end nucleotides for BsaI ligation in order of gBlock
    BsaI_four_nucleotides_ordered = []
    # length of nucleotides before and after special 4 nucleotide defining feature
    BsaI_arm_lengths = []
    # List of four nucleotides for BsaI that have already been guessed and should not be generated again for later gBlocks
    BsaI_four_nucleotides_chaos = []

    BsaI_four_nucleotides_chaos.append(start_end_BsaI_sites[0])
    BsaI_four_nucleotides_chaos.append(start_end_BsaI_sites[1])
    if len(gene_list) == 1:
        BsaI_four_nucleotides_ordered = [start_end_BsaI_sites]
        BsaI_arm_lengths = [[24, 26]]
    else:
        for i in range(len(gene_list)):
            if i == 0:
                start_four_nucleotide = start_end_BsaI_sites[0]
                end_four_nucleotide = dna_properties.generateRandomSequenceWithBlackList(4, BsaI_four_nucleotides_chaos)

                BsaI_four_nucleotides_chaos.append(end_four_nucleotide)

                BsaI_four_nucleotides_ordered.append([start_four_nucleotide, end_four_nucleotide])
                BsaI_arm_lengths.append([24, 29])
            elif i == len(gene_list) - 1:
                start_four_nucleotide = BsaI_four_nucleotides_ordered[-1][1]
                end_four_nucleotide = start_end_BsaI_sites[1]

                BsaI_four_nucleotides_ordered.append([start_four_nucleotide, end_four_nucleotide])
                BsaI_arm_lengths.append([34, 26])
            else:
                start_four_nucleotide = BsaI_four_nucleotides_ordered[-1][1]
                end_four_nucleotide = dna_properties.generateRandomSequenceWithBlackList(4, BsaI_four_nucleotides_chaos)

                BsaI_four_nucleotides_chaos.append(end_four_nucleotide)

                BsaI_four_nucleotides_ordered.append([start_four_nucleotide, end_four_nucleotide])
                BsaI_arm_lengths.append([34, 26])

    gBlocks = []
    associated_primers = {}

    counter = 1
    for i in range(len(gene_list)):
        gBlocks_raw = obtaingBlockCandidatesWithoutBsaISite(organism, gene_list[i], gene_table_for_organism, guide_rna_results_storage_folder, chopchop_dir, online)
        gBlocks_without_BsaI = gBlocks_raw[0]
        primer_output = gBlocks_raw[1]

        BsaI_nucleotides = BsaI_four_nucleotides_ordered[i]
        BsaI_arm_length = BsaI_arm_lengths[i]
        gBlocks_for_gene = []

        for j in range(len(gBlocks_without_BsaI)):
            build_complete = False
            total_build = ""

            while not build_complete:
                BsaI_left_arm = BsaI_nucleotides[0] + "|" + dna_properties.generateRandomSequence(1) + "|GAGACC"
                BsaI_left_arm = dna_properties.generateRandomSequence(5) + "|" + BsaI_left_arm
                BsaI_left_arm = BsaI_left_arm + "|" + dna_properties.generateRandomSequence(BsaI_arm_length[0] - len("".join(BsaI_left_arm.split("|"))))

                BsaI_right_arm = "GGTCTC|" + dna_properties.generateRandomSequence(1) + "|" + BsaI_nucleotides[1]
                BsaI_right_arm = BsaI_right_arm + "|" + dna_properties.generateRandomSequence(5)
                BsaI_right_arm = dna_properties.generateRandomSequence(BsaI_arm_length[1] - len("".join(BsaI_right_arm.split("|")))) + "|" + BsaI_right_arm

                total_build = BsaI_left_arm + "|" + gBlocks_without_BsaI[j] + "|" + BsaI_right_arm
                total_build_without_pipes = "".join(total_build.split("|"))

                if dna_properties.noHomopolymerMoreThanFourTs(total_build_without_pipes) and dna_properties.noHomopolymerMoreThanFiveAsAndGs(total_build_without_pipes):
                    build_complete = True

            gBlocks_for_gene.append(total_build)
            associated_primers[total_build] = primer_output[gBlocks_without_BsaI[j]]
            counter = counter + 1


        gBlocks_for_gene.insert(0, gene_list[i])
        gBlocks.append(gBlocks_for_gene)
    return [gBlocks, associated_primers]

organism = "sacCer3"
gene = "YBL011W"
gene_table_for_organism = gene_table.obtainGeneTableInfo(organism)
plasmid = "pCRCT"

gene_list = ["YDR210W", "YOL031C", "YPR030W"]
gBlocks_and_primers = obtaingBlocks(organism, gene_list, plasmid, gene_table_for_organism, "GuideRNA Archives", "chopchop", True)
file_name = "_".join(gene_list)

file = open("Test Output Garbage Dump/HI-CRISPR Knockout/" + file_name + ".txt", "w")

primer_wb = pyxl.Workbook()
ws = primer_wb[primer_wb.sheetnames[0]]

# Assigns A1 to "Name"
wcell_name = ws.cell(1, 1)
wcell_name.value = "Name"

# Assigns B1 to "Sequence"
wcell_sequence = ws.cell(1, 2)
wcell_sequence.value = "Sequence"

# Assigns C1 to "Scale"
wcell_scale = ws.cell(1, 3)
wcell_scale.value = "Scale"

# Assigns D1 to "Purification"
wcell_purification = ws.cell(1, 4)
wcell_purification.value = "Purification"

gene_primers_added_to_excel = []
# Outputs all the gBlock generation and primer verifiers into different
# text files and excel documents
#
# Look at Test Output Garbage Dump/HI-CRISPR Knockout for gBlocks and corresponding primers and amplicons.
for gene_position in range(len(gBlocks_and_primers[0])):
    gene = gBlocks_and_primers[0][gene_position]
    file.write(gene[0] + "\n\n")
    for i in range(1, len(gene)):
        file.write("gBlock #" + str(i) + "\n\n")
        file.write(gene[i] + "\n\n")

        primer_information = gBlocks_and_primers[1][gene[i]]
        print(primer_information)
        file.write(primer_information[0] + "\n\n")

        primer_list = primer_information[1]

        for j in range(len(primer_list)):
            file.write("Left Primer #" + str(j + 1) + ": " + primer_list[j][0][0] + "\n")
            file.write("Left Primer #" + str(j + 1) + " Melting Temperature: " + str(primer_list[j][0][1]) + "\n\n")

            file.write("Right Primer #" + str(j + 1) + ": " + primer_list[j][1][0] + " (search for " + dna_properties.getAntiSenseStrand(primer_list[j][1][0]) + ")" + "\n")
            file.write("Right Primer #" + str(j + 1) + ": " + str(primer_list[j][1][1]) + "\n\n")

            if j == 0 and not gene[0] in gene_primers_added_to_excel:
                left_primer_row_num = 2 * gene_position + 2
                right_primer_row_num = 2 * gene_position + 3

                # Assigns left primer row values

                wcell_lp_name = ws.cell(left_primer_row_num, 1)
                wcell_lp_name.value = gene[0] + "_LEFT_PRIMER"

                wcell_lp_sequence = ws.cell(left_primer_row_num, 2)
                wcell_lp_sequence.value = primer_list[j][0][0]

                wcell_lp_scale = ws.cell(left_primer_row_num, 3)
                wcell_lp_scale.value = "25nm"

                wcell_lp_purification = ws.cell(left_primer_row_num, 4)
                wcell_lp_purification.value = "STD"

                # Assigns right primer row values

                wcell_rp_name = ws.cell(right_primer_row_num, 1)
                wcell_rp_name.value = gene[0] + "_RIGHT_PRIMER"

                wcell_rp_sequence = ws.cell(right_primer_row_num, 2)
                wcell_rp_sequence.value = primer_list[j][1][0]

                wcell_rp_scale = ws.cell(right_primer_row_num, 3)
                wcell_rp_scale.value = "25nm"

                wcell_rp_purification = ws.cell(right_primer_row_num, 4)
                wcell_rp_purification.value = "STD"

                gene_primers_added_to_excel.append(gene[0])

    file.write("\n\n")

primer_wb.save("Test Output Garbage Dump/HI-CRISPR Knockout/" + file_name + ".xlsx")
file.close()
