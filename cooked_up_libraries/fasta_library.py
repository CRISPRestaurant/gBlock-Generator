# Obtain the genome information for an organism
def processOrganismFasta(organism_name):
    organism_fasta_file = open("CHOPCHOP_Necessary_Genetic_Data/Fasta_Genes/" + organism_name + ".fa", "r")

    output = {}

    chromosome = ""
    chromosome_build = ""

    total_genome = ""

    for line in organism_fasta_file:
        line = line[0: (len(line) - 1)]

        if ">" in line:
            if chromosome_build != "":
                output[chromosome] = chromosome_build
                total_genome += chromosome_build

            chromosome = line[1:]
            chromosome_build = ""
        else:
            chromosome_build += line

    return output
