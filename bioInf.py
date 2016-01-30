#author: Martin Mihal
#about:
    #python program - homework in course bioinformatic
    #we have given phylogenetic tree in file and we have to find
    #most possible evolution case


import math
E = math.exp(1)
BASE = ["A", "C", "G", "T"] #all posiblle

def jukes(changeFrom, changeTo, t, alpha):
    '''
    :param changeFrom:
    :param changeTo:
    :param t:
    :param alpha:
    :return: probability that given change occured
    '''
    if changeFrom == changeTo:         # no change in evolution
        return (1 + 3*E**(-4/3*alpha*t)) / 4
    else:                              # was change in evolution
        return (1 - E**(-4/3*alpha*t))/4

def felsenstein(tree, alpha = 0.1):
    '''
    felsenstein algorithm for computing the likelihood of an evolutionary tree from nucleic acid sequence data
    :param tree:
    :param alpha:
    :return: prob that given tree occured
    '''
    A = {} #2D dict[v,a], probality data in subtree with root v, if in node v is base a
    for node in tree.getLeafs(): #initializatio of leaf probality from given sequence
        for b in BASE:
            if node.name not in A:
                A[node.name] = {}
            if node.base == "-" or  node.base == "N":
                A[node.name][b] = 1
            else:
                A[node.name][b] =int(node.base == b)

    for node in tree.postorder():
        if not node.isLeaf(): # leafs was handle above
            for b in BASE:
                sl= 0
                if node.hasLeftChildren(): #firstly, handle left children
                    leftNode = node.childrens[0]
                    for L in BASE:
                        jukesProb = jukes(L, b, leftNode.lengthToParent, alpha)
                        sl += jukesProb * A[leftNode.name][L]
                sr = 0
                if node.hasRightChildren():  # secondly,handle right children
                    rightNode = node.childrens[1]
                    for R in BASE:
                        jukesProb = jukes(R, b, rightNode.lengthToParent, alpha)
                        sr += jukesProb * A[rightNode.name][R]

                if node.name not in A:
                    A[node.name] = {}
                A[node.name][b] = sl * sr

    result = 0
    for b in BASE:
        #print('{0:f}'.format(0.25 * A["Root"][b]))
        result += 0.25 * A["Root"][b]

    return result


class Tree():
    '''
    'clasic' tree class with a few method needed for felsenstein
    '''
    root = None

    def __init__(self, root):
        self.root = root

    def postorder(self, node=None, result = None):
        '''
        :param node:
        :param result:
        :return: list of nodes in tree in postorder order
        '''
        if node is None:
            node = self.root
            result = []
        if len(node.childrens) > 0:
            left = node.childrens[0]
            self.postorder(left, result)
        if len(node.childrens) > 1:
            right = node.childrens[1]
            self.postorder(right, result)
        result.append(node)
        return result

    def getLeafs(self):
        '''
        :return: all leaf nodes  in tree
        '''
        result = []
        for node in tree.postorder():
            if node.isLeaf():
                result.append(node)
        return result

    def setBaseOfLeafs(self, dictOfBases):
        for node in tree.postorder():
            if node.isLeaf():
                node.base = dictOfBases[node.name]

    def clearBases(self):
        for node in tree.postorder():
            node.base = None

class Node():
    '''
    class for representing node in tree
    '''
    name = None
    lengthToParent = None
    parentNode = None
    childrens = []
    base = None

    def __init__(self, name, lengthToParent, parentNode, base = None):
        self.lengthToParent = lengthToParent
        self.name = name
        self.parentNode = parentNode
        if parentNode is not None:
            self.parentNode.childrens.append(self)
        self.childrens = []
        self.base = None

    def addChildren(self, node):
        self.childrens.append(node)
        node.parent = self

    def isLeaf(self):
        return len(self.childrens) == 0

    def hasLeftChildren(self):
        return len(self.childrens) > 0

    def hasRightChildren(self):
        return len(self.childrens) > 1

    def __str__(self):
        return self.name

#creating tree from data
rootNode = Node("Root", 0, None)
tree = Tree(rootNode)
nodes = {} #name to node object
nodes["Root"] = rootNode
with open("data.txt") as f:
    for row in f.readlines():
        row = row.rstrip('\r\n')
        name, parent, length = row.split(" ")[0], row.split(" ")[1], float(row.split(" ")[2])
        parentNode = nodes.get(parent)
        if parentNode is None:
            print("None", name, parent, length)
        node = Node(name, length, parentNode)
        nodes[name] = node

#first try with all base as A
for node in tree.postorder():
    if node.isLeaf():
        node.base = "A"
# try it with different alpha factor and compare results
result1 = felsenstein(tree, 1)
result2 = felsenstein(tree, 0.2)
print("s alpha 1", result1) # vysledok bol 0.02607126198072631
print("s alpha 0.2", result2) # vysledok bol 0.1571909061005545


def frange(x, y, jump):
  '''
  iteating with float number
  :param x:
  :param y:
  :param jump:
  :return:
  '''
  while x < y:
    yield x
    x += jump

#task3
def findBestAlpha(sequences, tree, W):
    '''
    find best value for alpha factor
    :param sequences: dict pet -> alligment
    :param tree:
    :param W:
    :return: most possible alpha value
    '''
    bestAlpha, maxProbalityYet = None, float("-inf")
    for alpha in frange(0.1, 2.1, 0.1):#every coulmn send to tree, calculate probality
        probs = []
        for c in range(W):
            dictToTree = {}
            for nameOfPet in sequences:
                dictToTree[nameOfPet] = sequences.get(nameOfPet)[c]
            tree.clearBases()
            tree.setBaseOfLeafs(dictToTree)
            probOfColumn = felsenstein(tree, alpha)
            probs.append(probOfColumn)
        resultForAlphaLog = 0
        resultForAlpha = 1
        for p in probs:
            resultForAlpha = resultForAlpha * p
            resultForAlphaLog += math.log(p, 10)
        if resultForAlphaLog > maxProbalityYet:
            maxProbalityYet = resultForAlphaLog
            bestAlpha = alpha
    return bestAlpha

#task4
dictOfSequences = {} # to every pets his sequence from file
with open("cftr.txt") as f:
    currentPet = None
    for row in f.readlines():
        if  ">" in row:
            currentPet = row.split(">")[1].rstrip('\r\n')
            dictOfSequences[currentPet] = ""
            print(currentPet)
        else:
            if currentPet is None:
                continue
            dictOfSequences[currentPet] += row.rstrip('\r\n')
W = 100 #divide into sequences of length of 100 and send to findBestAlpha()
nOfWindows = 0
startW, endW = 0, W
while endW < 1100:
    argumentDict = {}
    for pet in dictOfSequences:
        seq = dictOfSequences.get(pet)
        argumentDict[pet] = seq[startW:endW]
    tree.clearBases()
    bestAlpha = findBestAlpha(argumentDict, tree, W)
    nOfWindows +=1
    startW += W
    endW += W
