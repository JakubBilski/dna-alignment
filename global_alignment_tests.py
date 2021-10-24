def get_alignments(sequence_a, sequence_b, gap_p, miss_p, good_p, max_paths):
    return 0, []


if __name__ == '__main__':
    score, alignments = get_alignments("AAAAA", "AAAAA", gap_p=-2, miss_p=-1, good_p=1, max_paths=5)
    assert len(alignments) == 1
    assert score == 5

    score, alignments = get_alignments("AAAAA", "AAAAA", gap_p=-2, miss_p=-1, good_p=1, max_paths=5)
    assert alignments[0] == ("AAAAA", "AAAAA")

    score, alignments = get_alignments("AA", "AA", gap_p=0, miss_p=0, good_p=0, max_paths=13)
    assert len(alignments) == 13
    assert score == 0
    score, alignments = get_alignments("TT", "AA", gap_p=0, miss_p=0, good_p=0, max_paths=13)
    assert len(alignments) == 13
    assert score == 0

    score, alignments = get_alignments("AAAATTTT", "CCCCGGGG", gap_p=-2, miss_p=0, good_p=1, max_paths=100)
    assert len(alignments) == 1
    assert score == 0
    assert alignments[0] == ("AAAATTTT", "CCCCGGGG")

    score, alignments = get_alignments("ACGTT", "CGTTA", gap_p=-2, miss_p=-5, good_p=1, max_paths=100)
    assert len(alignments) == 1
    assert score == -2+4-2
    assert alignments[0] == ("ACGTT-", "-CGTTA")

    score, alignments = get_alignments("GCATGCG", "GATTACA", gap_p=-1, miss_p=-1, good_p=1, max_paths=100)
    assert len(alignments) == 3
    assert score == 0
    assert ("GCATG-CG", "G-ATTACA") in alignments
    assert ("GCAT-GCG", "G-ATTACA") in alignments
    assert ("GCA-TGCG", "G-ATTACA") in alignments

    score, alignments = get_alignments("GCATGCG", "GATTACA", gap_p=-1, miss_p=-1, good_p=1, max_paths=1)
    assert len(alignments) == 1
    assert score == 0

    score, alignments = get_alignments("GCATGCG", "GATTACA", gap_p=-1, miss_p=-1, good_p=1, max_paths=2)
    assert len(alignments) == 2
    assert score == 0

    score, alignments = get_alignments("GCATGCG", "GATTACA", gap_p=-1, miss_p=-1, good_p=1, max_paths=3)
    assert len(alignments) == 3
    assert score == 0
