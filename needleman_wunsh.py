import argparse
import numpy as np
import time


def read_single_fasta_sequence(path):
    file = open(path, 'r')
    lines = file.readlines()
    file.close()
    sequence_id = lines[0]
    sequence_content = [a for a in "".join(lines[1:]) if a!= '\n']
    return sequence_id, sequence_content


def read_config(path):
    file = open(path, 'r')
    lines = file.readlines()
    file.close()
    num_lines = len(lines)
    if num_lines != 4:
        raise Exception(f"Config should consist of 4 lines. {num_lines} were found.")
    parameters = {'GAP_PENALTY': None, 'SAME_AWARD': None, 'DIFFERENCE_PENALTY': None, 'MAX_SEQ_LENGTH': None}
    for line in lines:
        for key in parameters.keys():
            if line[:len(key)] == key:
                parameters[key] = int(line[len(key)+1:])
    for k, v in parameters.items():
        if v is None:
            raise Exception(f"No configuration for parameter {k} was found.")
    gap_p = parameters['GAP_PENALTY']
    same_a = parameters['SAME_AWARD']
    diff_p = parameters['DIFFERENCE_PENALTY']
    max_seq_len = parameters['MAX_SEQ_LENGTH']
    return gap_p, same_a, diff_p, max_seq_len


VERTICAL = 1
HORIZONTAL = 2
CORNER = 4
VERTICAL_OR_CORNER = VERTICAL + CORNER
HORIZONTAL_OR_CORNER = HORIZONTAL + CORNER


def get_score_and_array(sequence_a, sequence_b, gap_p, diff_p, same_a):
    sequence_a = sequence_a[::-1]
    sequence_b = sequence_b[::-1]
    na = len(sequence_a)
    nb = len(sequence_b)

    incoming_directions = np.zeros((na+1, nb+1), dtype=np.int32)
    arr_scores = np.arange(0, nb+1)*(-gap_p)
    incoming_directions[0,1:] = np.ones(nb, dtype=np.int32)*HORIZONTAL
    incoming_directions[1:,0] = np.ones(na, dtype=np.int32)*VERTICAL
    for i in range(1, na+1):
        memorized_corner_score = arr_scores[0]
        arr_scores[0] -= gap_p
        for j in range(1, nb+1):
            horizontal_score = arr_scores[j-1]-gap_p
            vertical_score = arr_scores[j]-gap_p
            corner_score = memorized_corner_score + \
                (-diff_p if sequence_a[i-1] != sequence_b[j-1] else same_a)
            memorized_corner_score = arr_scores[j]
            arr_scores[j] = max(horizontal_score, vertical_score, corner_score)
            directions = 0
            if horizontal_score == arr_scores[j]:
                directions += HORIZONTAL
            if vertical_score == arr_scores[j]:
                directions += VERTICAL
            if corner_score == arr_scores[j]:
                directions += CORNER
            incoming_directions[i, j] = directions

    score = arr_scores[nb]
    return score, incoming_directions


def get_alignments(sequence_a, sequence_b, incoming_directions):
    sequence_a = sequence_a[::-1]
    sequence_b = sequence_b[::-1]
    na = len(sequence_a)
    nb = len(sequence_b)

    stack_gathered_alignment_length = [0]
    stack_point_i = [na]
    stack_point_j = [nb]
    stack_outcoming_direction = [0]

    alignment_a = ""
    alignment_b = ""

    while(len(stack_gathered_alignment_length) > 0):
        alignment_length = stack_gathered_alignment_length.pop()
        i = stack_point_i.pop()
        j = stack_point_j.pop()
        outcoming_direction = stack_outcoming_direction.pop()

        alignment_a = alignment_a[:alignment_length]
        alignment_b = alignment_b[:alignment_length]

        alignment_a += sequence_a[i] if outcoming_direction & VERTICAL_OR_CORNER else "-"
        alignment_b += sequence_b[j] if outcoming_direction & HORIZONTAL_OR_CORNER else "-"

        if i==0 and j==0:
            yield (alignment_a[1:], alignment_b[1:])
        else:
            if incoming_directions[i, j] & VERTICAL:
                stack_gathered_alignment_length.append(alignment_length+1)
                stack_point_i.append(i-1)
                stack_point_j.append(j)
                stack_outcoming_direction.append(VERTICAL)
            if incoming_directions[i, j] & HORIZONTAL:
                stack_gathered_alignment_length.append(alignment_length+1)
                stack_point_i.append(i)
                stack_point_j.append(j-1)
                stack_outcoming_direction.append(HORIZONTAL)
            if incoming_directions[i, j] & CORNER:
                stack_gathered_alignment_length.append(alignment_length+1)
                stack_point_i.append(i-1)
                stack_point_j.append(j-1)
                stack_outcoming_direction.append(CORNER)


if __name__ == '__main__':
    start_timestamp = time.time()
    description = 'Global alignment of DNA sequences using Needleman-Wunsh algorithm'

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('-a', type=str,
                        help='path to the first sequence in FASTA format')
    parser.add_argument('-b', type=str,
                        help='path to the second sequence in FASTA format')
    parser.add_argument('-c', type=str,
                        help='path to configuration file')
    parser.add_argument('-o', type=str,
                        help='output file path')
    args = parser.parse_args()

    sequence_a_id, sequence_a = read_single_fasta_sequence(args.a)
    sequence_b_id, sequence_b = read_single_fasta_sequence(args.b)
    gap_p, same_a, diff_p, max_seq_len = read_config(args.c)

    if len(sequence_a) > max_seq_len:
        raise Exception(f"First sequence was longer than it was allowed in config.")

    if len(sequence_b) > max_seq_len:
        raise Exception(f"Second sequence was longer than it was allowed in config.")

    output_file = open(args.o, 'w')

    score, incoming_directions = get_score_and_array(
                                    sequence_a=sequence_a,
                                    sequence_b=sequence_b,
                                    gap_p=gap_p,
                                    diff_p=diff_p,
                                    same_a=same_a)
    
    output_file.writelines([f'SCORE = {score}'])
    print(f'Score: {score}')
    num_alignments = 0
    for alignment_a, alignment_b in get_alignments(
                                    sequence_a=sequence_a,
                                    sequence_b=sequence_b,
                                    incoming_directions=incoming_directions):
        output_file.writelines(['\n', f'\n{alignment_a}', f'\n{alignment_b}'])
        num_alignments += 1
    print(f"Time of execution: {time.time() - start_timestamp :.4f} s")
    print(f"{num_alignments} alignments found")
    output_file.close()
