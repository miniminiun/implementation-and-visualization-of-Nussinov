import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap

pairs = [['A', 'U'], ['G', 'C'], ['U', 'A'], ['C', 'G']]

cdict = {'red': [[0.0, 0.0, 1.0],
                 [0.25, 1.0, 1.0],
                 [0.75, 1.0, 1.0],
                 [1, 1.0, 1.0]],
         'green':[[0.0, 0.0, 1.0],
                  [0.25, 1.0, 0.5],
                  [0.75, 0.5, 0.0],
                  [1.0, 0.0, 0.0]],
         'blue':[[0.0, 1.0, 1.0],
                 [0.25, 1.0, 0.5],
                 [0.75, 0.5, 0.0],
                 [1.0, 0.0, 0.0]]}

def findmaxk(i, j, dp):
    max = -1
    pointk = -4
    for k in range(i + 1, j):
        if dp[i][k] + dp[k + 1][j] > max:
            max = dp[i][k] + dp[k + 1][j]
            pointk = k
    return max, pointk


def max(a, b):
    if a > b:
        return a
    else:
        return b

def paintgraph(paren, dp, path, i, j):
    if i >= j:
        return
    else:
        dp[i][j] = 1.0
        if path[i][j] == -1:
            paintgraph(paren, dp, path, i + 1, j)
        elif path[i][j] == -2:
            paintgraph(paren, dp, path, i, j - 1)
        elif path[i][j] == -3:
            paren[i] = '('
            paren[j] = ')'
            paintgraph(paren, dp, path, i + 1, j - 1)
        elif path[i][j] == -4:
            paintgraph(paren, dp, path, i + 1, j - 1)
        else:
            for m in range(path[i][j] + 1, j):
                dp[i][m] = 0.5
            for m in range(i + 1, path[i][j] + 1):
                dp[m][j] = 0.5
            paintgraph(paren, dp, path, i, path[i][j])
            paintgraph(paren, dp, path, path[i][j] + 1, j)
        return



def main():
    f = open("RNA.txt", "r")
    line = f.readline()
    line = line.replace('\n', '')
    dp = [[0 for i in range(len(line))] for j in range(len(line))]
    path = [[0 for i in range(len(line))] for j in range(len(line))]
    paren = ['-' for i in range(len(line))]

    for i in range(len(line)):
        for j in range(len(line) - 1, -1, -1):
            maxvalue = 0;
            if j >= i:
                continue;
            else:
                maxk, k = findmaxk(j, i, dp)
                maxdia = 0
                pair = [line[j], line[i]]
                if not pair in pairs:
                    maxdia = dp[j + 1][i - 1]
                else:
                    maxdia = dp[j + 1][i - 1] + 1
                maxleft = dp[j][i - 1]
                maxdown = dp[j + 1][i]
                """dp[j][i] = max(maxk,max(maxdia,max(maxleft,maxdown)))"""
                if maxdown > maxleft:
                    if maxdown > maxdia:
                        if maxdown > maxk:
                            path[j][i] = -1
                            dp[j][i] = maxdown
                        else:
                            path[j][i] = k
                            dp[j][i] = maxk
                    else:
                        if maxdia > maxk:
                            if maxdia == dp[j + 1][i - 1] + 1:
                                path[j][i] = -3
                            else:
                                path[j][i] = -4
                            dp[j][i] = maxdia
                        else:
                            path[j][i] = k
                            dp[j][i] = maxk
                else:
                    if maxleft > maxdia:
                        if maxleft > maxk:
                            path[j][i] = -2
                            dp[j][i] = maxleft
                        else:
                            path[j][i] = k
                            dp[j][i] = maxk
                    else:
                        if maxdia > maxk:
                            if maxdia == dp[j + 1][i - 1] + 1:
                                path[j][i] = -3
                            else:
                                path[j][i] = -4
                            dp[j][i] = maxdia
                        else:
                            path[j][i] = k
                            dp[j][i] = maxk

    for i in range(len(line)):
        print(i, " ", end='')
    print("")
    for i in range(len(line)):
        print(i, " ", end='')
        for j in range(len(line)):
            print(dp[i][j], " ", end='')
        print("")
    print("")
    for i in range(len(line)):
        print(i, " ", end='')
    print("")
    for i in range(len(line)):
        print(i, " ", end='')
        for j in range(len(line)):
            print(path[i][j], " ", end='')
        print("")
    print("")
    bicolormap = LinearSegmentedColormap('bicolormap', segmentdata=cdict)

    dp2 = [[0 for i in range(len(line))] for j in range(len(line))]
    paintgraph(paren, dp2, path, 0, len(line) - 1)
    """for i in range(len(line)):
        print(dp2[i])"""

    col = []
    for i in range(len(line)):
        col.append(line[i])

    fig = plt.figure(figsize=(16, 6))
    ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
    ax2 = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])

    table1 = plt.table(
        cellText=[paren],
        colLabels=col,
        cellLoc='center'
    )
    the_table = plt.table(
        cellText=dp2,
        rowLabels=col,
        rowLoc='center',
        colLabels=col,
        colLoc='center',
        cellColours=bicolormap(dp2),
        cellLoc='center',
        loc='center'
    )
    """print(bicolormap(dp2))"""

    plt.show()

if __name__ == "__main__":
    main()
