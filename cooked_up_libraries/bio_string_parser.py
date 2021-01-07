def substring(string, end_position, length, left_endpoint = True):
    if left_endpoint:
        return string[end_position: end_position + length]
    else:
        return string[end_position - length + 1: end_position + 1]

def all_possible_substrings_by_length(string, length):
    output = []
    for i in range(0, len(string) - length + 1):
        output.append(string[i: i + length])

    return output
