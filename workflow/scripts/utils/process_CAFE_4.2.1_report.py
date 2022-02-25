import sys
import re
from collections import defaultdict


def readCAFEreport(CAFEreportInputFile):

    try:
        if 'input_data/' in CAFEreportInputFile:
            CAFEreport = open(CAFEreportInputFile, 'r')
            CAFEreportInputFile = CAFEreportInputFile.replace('input_data/', '')
        else:
            CAFEreport = open('input_data/' + CAFEreportInputFile, 'r')
    except OSError:
        print('\nInvalid argument to CAFE report input filename, cannot open path.\nPlease specify the input filename present in the input_data directory.\n\nUsage examples:\n\n\
    \tpython3 main.py <input_filename> -cafe <CAFEreport_input_filename>\n')
        sys.exit()

    lines = CAFEreport.readlines()

    CAFEreport.seek(0)
    CAFEreport.close()

    return(lines)


def getCAFEpvalues(CAFEreportInputFile):

    lines = readCAFEreport(CAFEreportInputFile)
    CAFEpvaluesDict = defaultdict()

    i = 0
    for line in lines:
        if line[0:4] == "\'ID\'":
            og_first_line = i + 1
            break
        else:
            i += 1

    if og_first_line:
        pass
    else:
        print('ERROR: could not find first OG line in CAFE report with processCAFEreport.py script. Check report formatting.')

    for line in lines[og_first_line:]:
        splitline = line.strip().split('\t')
        orthoGroup = splitline[0]
        CAFEpvalue = splitline[2]
        CAFEpvaluesDict[orthoGroup] = CAFEpvalue

    return(CAFEpvaluesDict)


# Function used later downstream
def getEventCounts(leaves, node):
    if(':' in node):
        counts_node = float(node[node.index('_') + 1:node.index(':')])
    else:
        counts_node = float(node[node.index('_') + 1:])
    if ',' in leaves:
        leaf1 = leaves.split(',')[0]
        counts_leaf1 = float(leaf1[leaf1.index('_') + 1:leaf1.index(':')])
        leaf2 = leaves.split(',')[1]
        counts_leaf2 = float(leaf2[leaf2.index('_') + 1:leaf2.index(':')])
        counts_leaves = [counts_leaf1, counts_leaf2]
    else:
        counts_leaves = [float(leaves[leaves.index('_') + 1:leaves.index(':')])]
    countsDict = {'leaves': counts_leaves, 'node': counts_node}
    return(countsDict)


# Function used later downstream
def addEventCounts(countsDict, eventCounts_Dict):
    for counts_leaf in countsDict['leaves']:
        if counts_leaf > countsDict['node']:
            eventCounts_Dict['nExpansions'] += 1
        elif counts_leaf < countsDict['node']:
            eventCounts_Dict['nContractions'] += 1
        else:
            eventCounts_Dict['nStable'] += 1
    return(eventCounts_Dict)


def getCAFEevents(CAFEreportInputFile):

    lines = readCAFEreport(CAFEreportInputFile)
    CAFEeventCountsDict = defaultdict(dict)

    for line in lines[10:]:
        splitline = line.split('\t')
        orthoGroup = splitline[0]
        tree = splitline[1]

        eventCountsDict = {'nExpansions': 0, 'nContractions': 0, 'nStable': 0}

        # Regular expression, don't touch
        match = re.search(r'\([^\(]*?\).*?[\,)]', tree)

        while match is not None:
            leaves = match.group().split(')')[0].strip('(')
            node = match.group().split(')')[1].strip(',')
            countsDict = getEventCounts(leaves, node)
            eventCountsDict = addEventCounts(countsDict, eventCountsDict)
            # Removing the computed leaves from tree, but keeps the node
            # to compute the events related to that node with its parent
            # in a recursive way
            tree = tree[:match.span()[0]] + tree[tree.find(')') + 1:]
            match = re.search(r'\([^\(]*?\).*?[\,)]', tree)

        # The very last iteration has a different '(leaves)node' syntax,
        # so it must be dealt apart, with its own RegEx.
        last_match = re.search(r'\([^\(]*?\).*?$', tree)
        if last_match is not None:
            leaves = last_match.group().split(')')[0].strip('(')
            node = last_match.group().split(')')[1].strip(',')
            countsDict = getEventCounts(leaves, node)
            eventCountsDict = addEventCounts(countsDict, eventCountsDict)

        CAFEeventCountsDict[orthoGroup] = eventCountsDict

    return(CAFEeventCountsDict)
