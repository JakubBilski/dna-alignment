from needleman_wunsh import get_alignments, get_score_and_array


if __name__ == '__main__':
    sequence_a = "AAAAA"
    sequence_b = "AAAAA"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=2, diff_p=1, same_a=1)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 1
    assert score == 5
    assert alignments[0] == ("AAAAA", "AAAAA")

    sequence_a = "AA"
    sequence_b = "AA"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=0, diff_p=0, same_a=0)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 13
    assert score == 0

    sequence_a = "TT"
    sequence_b = "AA"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=0, diff_p=0, same_a=0)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 13
    assert score == 0

    sequence_a = "AAAATTTT"
    sequence_b = "CCCCGGGG"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=2, diff_p=0, same_a=1)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 1
    assert score == 0
    assert alignments[0] == ("AAAATTTT", "CCCCGGGG")

    sequence_a = "ACGTT"
    sequence_b = "CGTTA"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=2, diff_p=5, same_a=1)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 1
    assert score == -2+4-2
    assert alignments[0] == ("ACGTT-", "-CGTTA")

    sequence_a = "GCATGCG"
    sequence_b = "GATTACA"
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=1, diff_p=1, same_a=1)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 3
    assert score == 0
    assert ("GCATG-CG", "G-ATTACA") in alignments
    assert ("GCAT-GCG", "G-ATTACA") in alignments
    assert ("GCA-TGCG", "G-ATTACA") in alignments

    sequence_a = ""
    sequence_b = ""
    score, array = get_score_and_array(sequence_a, sequence_b, gap_p=200, diff_p=100, same_a=100)
    alignments = list(get_alignments(sequence_a, sequence_b, array))
    assert len(alignments) == 1
    assert score == 0
    assert ("", "") in alignments
