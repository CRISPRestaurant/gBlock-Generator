# Finds BsaI sites and their four nucleotide lego extension
def findBsaISites(plasmid_name):
    file = open("Plasmids/" + plasmid_name + ".fasta", "r")
    total_genome = ""

    first_line = True

    start_bsai_positions = []
    start_bsai_four_nucleotides = []

    end_bsai_positions = []
    end_bsai_four_nucleotides = []

    for line in file:
        line = line[0: (len(line) - 1)]

        if first_line:
            first_line = False
        else:
            line = line.upper()

            start_bsai_position = line.find("GAGACC")
            end_bsai_position = line.find("GGTCTC")

            if start_bsai_position != -1:
                start_bsai_positions.append(len(total_genome) + start_bsai_position)

                start_four_nucleotide_position = start_bsai_position - 5
                start_four_nucleotide = line[start_four_nucleotide_position: start_four_nucleotide_position + 4]

                start_bsai_four_nucleotides.append(start_four_nucleotide)

            if end_bsai_position != -1:
                end_bsai_positions.append(len(total_genome) + end_bsai_position)

                end_four_nucleotide_position = end_bsai_position + 7
                end_four_nucleotide = line[end_four_nucleotide_position: end_four_nucleotide_position + 4]

                end_bsai_four_nucleotides.append(end_four_nucleotide)

            total_genome += line

    return [start_bsai_four_nucleotides[0], end_bsai_four_nucleotides[0]]
