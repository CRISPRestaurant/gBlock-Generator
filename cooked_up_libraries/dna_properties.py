from random import *
import cooked_up_libraries.bio_string_parser as bsp
import primer3

# Generates a random DNA sequence of desired length
def generateRandomSequence(length):
    output = ""
    for i in range(length):
        output += choice(["A", "T", "G", "C"])
    return output

# Generates a random DNA sequence of desired length after ensuring it is not in the
# blacklisted_sequences
# Will output "ATC" if length = 3 and blacklisted_sequences is ["GTA", "TAC"]
# Will not output "ATC" if length = 3 and blacklisted_sequences is ["ATC", "ACC"]
def generateRandomSequenceWithBlackList(length, blacklisted_sequences):
    output = generateRandomSequence(length)

    while output in blacklisted_sequences:
        output = generateRandomSequence(length)

    return output

def generateRandomSequenceWithoutRepeatNucleotidesOfMinimalLength(length, max_repeat):
    output = ""

    for i in range(length):
        repeat_at_maximum = False
        character_if_at_repeat_maximum = "haha you don't know me yet"

        if i >= max_repeat:
            last_max_repeat = bsp.substring(output, len(output) - 1, max_repeat, left_endpoint = False)
            unique_list = []

            for character in last_max_repeat:
                if character not in unique_list:
                    unique_list.append(character)

            if len(unique_list) == 1:
                repeat_at_maximum = True
                character_if_at_repeat_maximum = unique_list[0]

        if repeat_at_maximum:
            output += generateRandomSequenceWithBlackList(1, [character_if_at_repeat_maximum])
        else:
            output += generateRandomSequence(1)

    return output

def generateRandomSequenceWithoutRepeatNucleotidesOfMinimalLengthAndWithBlackList(length, max_repeat, blacklisted_sequences):
    output = generateRandomSequenceWithoutRepeatNucleotidesOfMinimalLength(length, max_repeat)

    while output in blacklisted_sequences:
        output = generateRandomSequenceWithoutRepeatNucleotidesOfMinimalLength(length, max_repeat)

    return output


# Return the complementary strand for a sequence
def getComplementaryStrand(sequence):
    complement = {"A": "T", "T": "A", "G": "C", "C": "G"}

    output = ""
    for i in sequence:
        output += complement[i]

    return output

# Return the anti-sense strand for a given sequence
def getAntiSenseStrand(gRNA):
    complement = getComplementaryStrand(gRNA)

    return complement[::-1]

##################### MAJOR SPECIALIZED PROPERTIES REGARDING CRISPR CAS-9 #############################

# Ensures no BsaI restriction sites
def noBsaIRestrictionSites(sequence):
    if "GGTCTC" in sequence or "GAGACC" in sequence:
        return False
    else:
        return True

# Ensures no sub-sequence more than four Ts
def noHomopolymerMoreThanFourTs(sequence):
    loader = ""

    for i in sequence:
        if i == "T":
            loader += "T"

            if len(loader) > 4:
                return False
        else:
            loader = ""

    return True

# Ensures no sub-sequence more than 5 As and 5 Gs
def noHomopolymerMoreThanFiveAsAndGs(sequence):
    loader = ""

    for i in sequence:
        if i == "A":
            loader += "A"

            if len(loader) > 5:
                return False
        else:
            loader = ""

    for i in sequence:
        if i == "G":
            loader += "G"

            if len(loader) > 5:
                return False
        else:
            loader = ""

    return True

################## IMPORTANT DNA PROPERTIES FOR RANKING GUIDE RNAS ######################

# Returns the GC content for a sequence
def obtainGCContent(sequence):
    g_c_count = 0

    for i in sequence:
        if i == "C" or i == "G":
            g_c_count += 1

    return float(g_c_count) / len(sequence)

# Assigns a GC priority for a given sequence, giving
# higher priority to sequences with GC content between
# 35% and 75%
def obtainGCPriority(sequence):
    g_c_content = obtainGCContent(sequence)

    if g_c_content >= 0.35 and g_c_content <= 0.75:
        return 2
    else:
        return 1

# Returns whether GC content of primer candidate is within
# optimal range between 40% and 60%
def isOptimalGCPrimer(sequence):
    gc_content = obtainGCContent(sequence)

    if gc_content >= 0.40 and gc_content <= 0.60:
        return True
    else:
        return False

def isMeltingTemperatureWithinRange(sequence, lower_end, upper_end):
    melting_temperature = primer3.calcTm(sequence)

    return melting_temperature >= lower_end and melting_temperature <= upper_end

def hasHairpinStructure(sequence):
    return primer3.calcHairpin(sequence).structure_found

# Assigns a priority based on the purine content in
# the last four nucleotides of the sequence
def obtainLastPurinePriority(sequence):
    last_four = sequence[(len(sequence) - 4): ]

    purine_count = 0

    for i in sequence:
        if i == "A" or i == "G":
            purine_count += 1

    return 2 ** purine_count
