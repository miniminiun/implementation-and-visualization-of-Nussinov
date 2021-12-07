import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.animation import FuncAnimation
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
class RNA:
    line = 0
    dp = None
    path = None
    paren = None
    dp2 = None
    fig = None
    ax1 = None
    ax2 = None
    col = []
    the_table = None
    table1 = None
    leftpath = []
    downpath = []
    node = []

    def __init__(self):
        f = open("RNA.txt", "r")
        self.line = f.readline()
        self.line = self.line.replace('\n', '')
        length = len(self.line)
        self.dp = [[0 for i in range(length)] for j in range(length)]
        self.path = [[0 for i in range(length)] for j in range(length)]
        self.paren = ['-' for i in range(length)]
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, figsize=(10,6), gridspec_kw={'height_ratios': [1, len(self.line)]})

    def findmaxk(self, i, j):
        max = -1
        pointk = -4
        for k in range(i + 1, j):
            if self.dp[i][k] + self.dp[k + 1][j] > max:
                max = self.dp[i][k] + self.dp[k + 1][j]
                pointk = k
        return max, pointk


    def max(self, a, b):
        if a > b:
            return a
        else:
            return b

    def paintgraph(self, i, j):
        if i >= j:
            return
        else:
            self.dp2[i][j] = 1.0
            if self.path[i][j] == -1:
                self.paintgraph(i + 1, j)
            elif self.path[i][j] == -2:
                self.paintgraph(i, j - 1)
            elif self.path[i][j] == -3:
                #self.paren[i] = '('
                #self.paren[j] = ')'
                self.paintgraph(i + 1, j - 1)
            elif self.path[i][j] == -4:
                self.paintgraph(i + 1, j - 1)
            else:
                for m in range(self.path[i][j] + 1, j):
                    self.dp2[i][m] = 0.5
                for m in range(i + 1, self.path[i][j] + 1):
                    self.dp2[m][j] = 0.5
                self.paintgraph(i, self.path[i][j])
                self.paintgraph(self.path[i][j] + 1, j)
            return

    def update(self, num):
        """for (i, j) in self.node:
            print(i, j)
            self.the_table._cells[(i, j)].set_facecolor('coral')"""
        if self.downpath or self.leftpath:
            tempdownpath = []
            for (i, j) in self.downpath:
                self.the_table._cells[(i + 1, j)].set_facecolor('coral')
                if self.dp2[i + 1][j] == 1.0:
                    self.node.append((i + 1, j))
                else:
                    tempdownpath.append((i + 1, j))
            self.downpath = tempdownpath
            templeftpath = []
            for (i, j) in self.leftpath:
                #print(i, j, self.dp2[i][j - 1])
                self.the_table._cells[(i + 1, j)].set_facecolor('coral')
                if self.dp2[i][j - 1] == 1.0:
                    self.node.append((i, j - 1))
                else:
                    templeftpath.append((i, j - 1))
            self.leftpath = templeftpath
        else:
            #print("here2")
            tempnode = []
            for (i, j) in self.node:
                self.the_table._cells[(i + 1, j)].set_facecolor('red')
                if self.path[i][j] == -3:
                    self.table1._cells[(1, i)]._text.set_text('(')
                    self.table1._cells[(1, j)]._text.set_text(')')
                    if i + 1 < j - 1:
                        tempnode.append((i + 1, j - 1))
                elif self.path[i][j] == -4:
                    if i + 1 < j - 1:
                        tempnode.append((i + 1, j - 1))
                elif self.path[i][j] == -1:
                    if i + 1 < j:
                        tempnode.append((i + 1, j))
                elif self.path[i][j] == -2:
                    if i < j - 1:
                        tempnode.append((i, j - 1))
                else:
                    self.downpath.append((i + 1, j))
                    self.leftpath.append((i, j - 1))
            self.node = tempnode

    def start(self):
        #path == -3 means diagnal and add 1, -4 means diagnal without add 1, -1 means down, -2 means left
        for i in range(len(self.line)):
            for j in range(len(self.line) - 1, -1, -1):
                maxvalue = 0;
                if j >= i:
                    continue;
                else:
                    maxk, k = self.findmaxk(j, i)
                    maxdia = 0
                    pair = [self.line[j], self.line[i]]
                    if not pair in pairs:
                        maxdia = self.dp[j + 1][i - 1]
                    else:
                        maxdia = self.dp[j + 1][i - 1] + 1
                    maxleft = self.dp[j][i - 1]
                    maxdown = self.dp[j + 1][i]
                    """dp[j][i] = max(maxk,max(maxdia,max(maxleft,maxdown)))"""
                    if maxdown > maxk:
                        if maxdown > maxleft:
                            if maxdown > maxdia:
                                self.path[j][i] = -1
                                self.dp[j][i] = maxdown
                            else:
                                if maxdia == self.dp[j + 1][i - 1] + 1:
                                    self.path[j][i] = -3
                                else:
                                    self.path[j][i] = -4
                                self.dp[j][i] = maxdia
                        else:
                            if maxleft > maxdia:
                                self.path[j][i] = -2
                                self.dp[j][i] = maxleft
                            else:
                                if maxdia == self.dp[j + 1][i - 1] + 1:
                                    self.path[j][i] = -3
                                else:
                                    self.path[j][i] = -4
                                self.dp[j][i] = maxdia
                    else:
                        if maxk > maxleft:
                            if maxk > maxdia:
                                self.path[j][i] = k
                                self.dp[j][i] = maxk
                            else:
                                if maxdia == self.dp[j + 1][i - 1] + 1:
                                    self.path[j][i] = -3
                                else:
                                    self.path[j][i] = -4
                                self.dp[j][i] = maxdia
                        else:
                            if maxleft > maxdia:
                                self.path[j][i] = -2
                                self.dp[j][i] = maxleft
                            else:
                                if maxdia == self.dp[j + 1][i - 1] + 1:
                                    self.path[j][i] = -3
                                else:
                                    self.path[j][i] = -4
                                self.dp[j][i] = maxdia

        """for i in range(len(line)):
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
        print("")"""
        bicolormap = LinearSegmentedColormap('bicolormap', segmentdata=cdict)

        self.dp2 = [[0 for i in range(len(self.line))] for j in range(len(self.line))]
        self.paintgraph(0, len(self.line) - 1)
        """for i in range(len(line)):
            print(dp2[i])"""


        for i in range(len(self.line)):
            self.col.append(self.line[i])

        #fig = plt.figure(figsize=(16, 6))
        #ax = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])
        #ax2 = fig.add_subplot(111, frameon=True, xticks=[], yticks=[])


        #fig, (ax1, ax2) = plt.subplots(2, figsize=(16, 6))
        #self.fig.set_tight_layout(True)
        self.ax1.set_axis_off()
        self.ax2.set_axis_off()

        self.table1 = self.ax1.table(
            cellText=[self.paren],
            colLabels=self.col,
            cellLoc='center',
            loc='upper left'
        )

        self.the_table = self.ax2.table(
            cellText=self.dp,
            rowLabels=self.col,
            rowLoc='center',
            colLabels=self.col,
            colLoc='center',
            cellLoc='center',
            #cellColours=bicolormap(self.dp2),
            loc='upper left'
        )

        """print(bicolormap(dp2))"""
        self.node.append((0, len(self.line) - 1))
        ani = FuncAnimation(self.fig, self.update, repeat=False)

        plt.show()

if __name__ == "__main__":
    rna = RNA()
    rna.start()
