# Parses the genes, their information, and their exon positions from an organism
# gene file
def obtainGeneTableInfo(organism_name):
    output_dict = {}
    property_list = []

    first_line = True

    file = open("CHOPCHOP_Necessary_Genetic_Data/Gene_Tables/" + organism_name + ".gene_table")

    for line in file:
        line = line[0: (len(line) - 1)]

        if first_line:
            property_list = line.split("\t")
            first_line = False
        else:
            line_split = line.split("\t")

            output_dict[line_split[1]] = {}
            for i in range(0, 9):
                if i != 1:
                    output_dict[line_split[1]][property_list[i]] = line_split[i]

            exon_start_positions = line_split[9].split(",")
            exon_start_positions = exon_start_positions[0: (len(exon_start_positions) - 1)]

            exon_end_positions = line_split[10].split(",")
            exon_end_positions = exon_end_positions[0: (len(exon_end_positions) - 1)]

            exon_ranges = []

            for i in range(len(exon_start_positions)):
                exon_ranges.append([exon_start_positions[i], exon_end_positions[i]])

            output_dict[line_split[1]]["exonRanges"] = exon_ranges

            output_dict[line_split[1]]["proteinID"] = line_split[11]

    return output_dict
