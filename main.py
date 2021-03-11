import numpy as np


def base_pair(i, j):
    pair = [seq[i], seq[j]]

    if pair == ['A', 'U'] or pair == ['U', 'A']:
        return (1)
    elif pair == ['C', 'G'] or pair == ['G', 'C']:
        return (1)
    else:
        return (0)


def score(seq, N):
    L = len(seq)
    s = np.zeros((L, L))
    for L in range(1, N):
        for i in range(0, N - L):
            j = i + L
            if j - i >= 1:
                case1 = s[i + 1, j - 1] + base_pair(i, j)
                case2 = s[i + 1, j]
                case3 = s[i, j - 1]
                case4 = 0
                tmp = []
                for k in range(i + 1, j):
                    tmp.append(s[i, k] + s[k + 1, j])
                    case4 = max(tmp)
                    s[i, j] = max(case1, case2, case3, case4)
            else:
                s[i, j] = 0
    return s


def traceback(s, seq, i, j, pair):
    if i < j:
        if s[i, j] == s[i + 1, j - 1] + base_pair(i, j):
            pair.append([i, j, str(seq[i]), str(seq[j])])
            traceback(s, seq, i + 1, j - 1, pair)
        elif s[i, j] == s[i + 1, j]:
            traceback(s, seq, i + 1, j, pair)
        elif s[i, j] == s[i, j - 1]:
            traceback(s, seq, i, j - 1, pair)
        else:
            for k in range(i + 1, j):
                if s[i, j] == s[i, k] + s[k + 1, j]:
                    traceback(s, seq, k + 1, j, pair)
                    traceback(s, seq, i, k, pair)
                    break
    return pair


def run(seq):
    N = len(seq)
    s = score(seq, N)
    pair = traceback(s, seq, 0, len(seq) - 1, [])
    print(seq)
    print(s)
    print(pair)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    seq = "GGGAAAUCC"
    s = np.zeros([len(seq), len(seq)])
    np.fill_diagonal(s, 0)
    run(seq)
