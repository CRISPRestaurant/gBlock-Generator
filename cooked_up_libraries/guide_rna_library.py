import os
# Uses the CHOPCHOP API to obtain guideRNA candidates for an
# organism and a gene
def obtainGuideRNAs(organism_name, target_gene):
    command = "chopchop/chopchop.py -G " + organism_name + " -o chopchop/temp -Target " + target_gene
    command_stream = os.popen(command)
    command_output = command_stream.read()

    command_output_split = command_output.split("\n")

    for i in range(len(command_output_split)):
        command_output_split[i] = command_output_split[i].split("\t")

    command_output_split = command_output_split[0: len(command_output_split) - 1]
    return command_output_split
