import numpy as np
import re


def seq_valid(seq):
    return bool(re.match("^[AUGC]+$", seq))


def is_base_pair(i, j):
    """
    returns 1 if i,j are base pairs, 0 otherwise
    """
    pair = [seq[i], seq[j]]
    if pair == ['A', 'U'] or pair == ['U', 'A'] or pair == ['C', 'G'] or pair == ['G', 'C']:
        return 1
    else:
        return 0

def score(seq, N):
    """
    score the DP matrix according to the OPT formula
    """
    L = len(seq)
    # fill matrix with zeros
    matrix = np.zeros((L, L))
    # np.fill_diagonal(matrix, 0) # necessary?
    # fil the DP matrix
    for L in range(1, N):
        for i in range(0, N - L):
            j = i + L
            if j - i >= 1:
                # case 1: i,j pair -> add i,j pair onto best structure found for subsequence i+1, j-1
                case1 = matrix[i + 1, j - 1] + is_base_pair(i, j)
                # case 2: i unpaired -> add unpaired position i onto best structure found for subsequence i+1, j
                case2 = matrix[i + 1, j]
                # case 3: j unpaired -> add unpaired position j onto best structure found for subsequence i, j-1
                case3 = matrix[i, j - 1]
                # case 4: bifurcation: combine tow optimal substructures i,k and k+1,j
                case4 = 0
                tmp = []
                for k in range(i + 1, j):
                    tmp.append(matrix[i, k] + matrix[k + 1, j])
                    case4 = max(tmp)
                matrix[i, j] = max(case1, case2, case3, case4)
            else:
                matrix[i, j] = 0
    return matrix


def traceback(matrix, seq, i, j, pairs):
    if i < j:  # main diagonal not reached
        if matrix[i, j] == matrix[i + 1, j - 1] + is_base_pair(i, j):
            # add the pair (i,j) to our list of pairs
            pairs.append([i, j, str(seq[i]), str(seq[j])])
            # move the recursion to the button-left diagonal
            traceback(matrix, seq, i + 1, j - 1, pairs)
        # if i is unpaired, there will be no change in score when we take it out, so we just recurse to the next index
        elif matrix[i, j] == matrix[i + 1, j]:
            traceback(matrix, seq, i + 1, j, pairs)
        # if j is unpaired, there will be no change in score when we take it out, so we just recurse to the next index
        elif matrix[i, j] == matrix[i, j - 1]:
            traceback(matrix, seq, i, j - 1, pairs)
        else:
            # try pairing j with a matching index k.
            for k in range(i + 1, j):
                if matrix[i, j] == matrix[i, k] + matrix[k + 1, j]:
                    traceback(matrix, seq, k + 1, j, pairs)
                    traceback(matrix, seq, i, k, pairs)
                    break
    return pairs


def print_results(pairs, matrix):
    print("Results for sequence: \n" + seq)
    parenthesis_list = ['.'] * len(seq)
    # print Matrix
    print("Scored Matrix:")
    print(matrix)
    print()
    # print list of pairs
    print("List of pairs:")
    for pair in pairs:
        print(pair)
        parenthesis_list[pair[0]] = '('
        parenthesis_list[pair[1]] = ')'
    print()
    # print Parenthesis representation of matching pairs
    print("Parenthesis representation of matching pairs:")
    parenthesis_str = ""
    print(seq)
    print(parenthesis_str.join(parenthesis_list))


def run(seq):
    matrix = score(seq, len(seq))
    pairs = traceback(matrix, seq, 0, len(seq) - 1, [])
    print_results(pairs, matrix)


if __name__ == '__main__':
    seq = "GGCAGUACCAAGUCGCGAAAGCGAUGGCCUUGCAAAGGGUAUGGUAAUAAGCUGCC"
    run(seq)
    # seq = input("Enter an RNA sequence consisting of A, U, G, C bases only & Press ENTER:\n")
    # if seq_valid(seq):
    #     run(seq)
    # else:
    #     print("The RNA sequence must contain only A, U, G, C bases")
