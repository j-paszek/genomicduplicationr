from treeop import Tree, str2tree


def readgs(filename, verbose=""):
    t = [l.strip() for l in open(filename, "r").readlines() if len(l.strip()) and l.strip()[0] != '#']
    if verbose.find("3") != -1:
        print(t)
    return Tree(str2tree(t[0])), Tree(str2tree(t[1]))


