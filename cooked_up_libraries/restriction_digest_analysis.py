import cooked_up_libraries.bio_string_parser as bsp
import cooked_up_libraries.dna_properties as dnap
import matplotlib.pyplot as plt
import primer3

# Obtain a series of primers from a target strand
def obtainPrimers(target, left_primer_start, left_primer_length, right_primer_start, right_primer_length):
    raw_results = primer3.bindings.designPrimers(
        {
            'SEQUENCE_ID': 'MH1000',
            'SEQUENCE_TEMPLATE': target,
            'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': [[left_primer_start, left_primer_length, right_primer_start, right_primer_length]]
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
            'PRIMER_MAX_GC': 80.0
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

    #print(raw_results)
    # Saves explicity primer3 output data to analyze errors
    #df = pd.DataFrame(list(raw_results.items()), columns = ["Key", "Value"])
    #df.to_csv(location, index = False)

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
        left_primer_position_and_length = raw_results["PRIMER_LEFT_" + str(i)]
        left_primer_position = left_primer_position_and_length[0]
        left_primer_length = left_primer_position_and_length[1]

        # Obtain right primer sequence
        right_primer = raw_results["PRIMER_RIGHT_" + str(i) + "_SEQUENCE"]

        # Obtain right primer melting temperature
        right_primer_tm = raw_results["PRIMER_RIGHT_" + str(i) + "_TM"]

        # Obtain right primer position
        right_primer_position_and_length = raw_results["PRIMER_RIGHT_" + str(i)]
        right_primer_position = right_primer_position_and_length[0]
        right_primer_length = right_primer_position_and_length[1]

        output.append([[left_primer, left_primer_tm, left_primer_position, left_primer_length], [right_primer, right_primer_tm, right_primer_position, right_primer_length]])

    return output

def obtainPrimersForMeaningfulDigestAnalysis(chromosome_sequence, restriction_site_region_start_position, restriction_site_region_end_position, base_pair_position_before_digest_cut, degree_of_asymmetry, minimum_spacing_length, primer_search_space_length, primer_length, maximum_amplicon_length):
    space_to_the_left = restriction_site_region_start_position # Correct
    space_to_the_right = len(chromosome_sequence) - restriction_site_region_end_position - 1 # Correct

    is_left_minimum = min(space_to_the_left, space_to_the_right) == space_to_the_left

    spacing_from_left_primer_search_space_to_restriction_region = 0
    spacing_from_right_primer_search_space_to_restriction_region = 0

    left_primer_search_space_shift_rate = 1
    right_primer_search_space_shift_rate = 1

    if is_left_minimum:
        spacing_from_left_primer_search_space_to_restriction_region = minimum_spacing_length
        spacing_from_right_primer_search_space_to_restriction_region = minimum_spacing_length * degree_of_asymmetry

        right_primer_search_space_shift_rate *= degree_of_asymmetry
    else:
        spacing_from_left_primer_search_space_to_restriction_region = minimum_spacing_length * degree_of_asymmetry
        spacing_from_right_primer_search_space_to_restriction_region = minimum_spacing_length

        left_primer_search_space_shift_rate *= degree_of_asymmetry

    end_position_of_left_primer_search_space = restriction_site_region_start_position - spacing_from_left_primer_search_space_to_restriction_region - 1 # Correct
    start_position_of_right_primer_search_space = restriction_site_region_end_position + spacing_from_right_primer_search_space_to_restriction_region + 1 # Correct

    start_position_of_left_primer_search_space = end_position_of_left_primer_search_space - primer_search_space_length + 1
    end_position_of_right_primer_search_space = start_position_of_right_primer_search_space + primer_search_space_length - 1

    primers = []
    while (restriction_site_region_start_position - spacing_from_left_primer_search_space_to_restriction_region - primer_search_space_length >= 0) and (restriction_site_region_end_position + spacing_from_right_primer_search_space_to_restriction_region + primer_search_space_length < len(chromosome_sequence)) and (restriction_site_region_end_position - restriction_site_region_start_position + spacing_from_right_primer_search_space_to_restriction_region + spacing_from_left_primer_search_space_to_restriction_region + 2 * primer_search_space_length + 1 <= maximum_amplicon_length):

        amplicon_region = chromosome_sequence[start_position_of_left_primer_search_space: end_position_of_right_primer_search_space + 1]

        left_primer_start_search_space_region = 0
        left_primer_search_space_length = primer_search_space_length

        right_primer_start_search_space_region = start_position_of_right_primer_search_space - start_position_of_left_primer_search_space
        right_primer_search_space_length = primer_search_space_length

        primers = obtainPrimers(amplicon_region, left_primer_start_search_space_region, left_primer_search_space_length, right_primer_start_search_space_region, right_primer_search_space_length)

        if primers != []:
            break

        spacing_from_left_primer_search_space_to_restriction_region += left_primer_search_space_shift_rate
        spacing_from_right_primer_search_space_to_restriction_region += right_primer_search_space_shift_rate

        end_position_of_left_primer_search_space -= left_primer_search_space_shift_rate
        start_position_of_left_primer_search_space -= left_primer_search_space_shift_rate

        start_position_of_right_primer_search_space += right_primer_search_space_shift_rate
        end_position_of_right_primer_search_space += right_primer_search_space_shift_rate

    primer_information = []
    for left_and_right_primer in primers:
        left_primer_results = left_and_right_primer[0]
        right_primer_results = left_and_right_primer[1]

        left_primer_start_position = start_position_of_left_primer_search_space + left_primer_results[2]
        right_primer_end_position = start_position_of_left_primer_search_space + right_primer_results[2]

        left_dna_cut_length = base_pair_position_before_digest_cut - left_primer_start_position + 1
        right_dna_cut_length = right_primer_end_position - base_pair_position_before_digest_cut

        primer_information.append([left_primer_results, right_primer_results, [left_dna_cut_length, right_dna_cut_length]])

    return primer_information
