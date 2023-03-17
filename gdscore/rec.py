import itertools
import math


def rec(gtrees, st):
    dup = []  # duplication nodes == nodes with defined interval
    stnodespostorder = st.root.get_nodes()
    gtnodespostorder = list(itertools.chain.from_iterable(gt.root.get_nodes() for gt in gtrees))

    # limit gene tree nodes only to duplications
    for g in gtnodespostorder:
        if g.interval:
            dup.append(g)

    for s in stnodespostorder:
        s.topinterval = []
        s.botinterval = []
        s.allintervals = []
        s.active = False
        s.min_len_int = math.inf

    for g in dup:
        g.interval[0].botinterval.append(g)
        g.interval[1].topinterval.append(g)

    # set min_len_int value for every node s
    # which is the distance to the nearest top of an interval that pass thru s
    for s in stnodespostorder:
        if s.topinterval:
            s.min_len_int = 0   # node s is an ending of some interval
            s.active = True     # this is a place for an episode, all dup nodes with passing intervals are set to here
            if s.allintervals:
                for i in s.allintervals:
                    i.interval[1].topinterval.remove(i)  # passing interval is removed from interval top list
            if s.botinterval:
                for i in s.botinterval:
                    i.interval[1].topinterval.remove(i)  # starting interval is removed from interval top list
        else:
            if s.botinterval:       # add to parent all intervals that start here
                for i in s.botinterval:
                    s.parent.allintervals.append(i)
            if s.allintervals:      # add to parent all intervals that pass thru here
                for i in s.allintervals:
                    s.parent.allintervals.append(i)

    # Counts selected episode node
    ec_score = 0
    for s in stnodespostorder:
        if s.active:
            ec_score += 1
            print(s)

    return ec_score
