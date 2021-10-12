import sys
import re
from collections import defaultdict

input_file = snakemake.input.cafe_results
output_file = open(snakemake.output.copy_number_variation_table, 'w')


def readCAFEreport(CAFEreportInputFile):

    with open(input_file) as f:
        lines = f.readlines()

    return lines


# String pattern search function to detect copy number counts in the the tree elements (node or leaf1 or leaf2)
def getEventCounts(leaves, node):
    if ';' in node:
        node = node.strip(';')


    def get_match(tree_element):
        count_pattern = '>_'
        significant_count_pattern = '>*_'
        count_index = 2
        if significant_count_pattern in tree_element:
            count_pattern = significant_count_pattern
            count_index = 3

        return [count_pattern, count_index]


    count_pattern, count_index = get_match(node)

    if ':' in node:
        counts_node = float(node[node.index(count_pattern)+count_index : node.index(':')])
    else:
        counts_node = float(node[node.index(count_pattern)+count_index : ])

    count_pattern, count_index = get_match(leaves)

    if ',' in leaves:
        leaf1 = leaves.split(',')[0]
        count_pattern, count_index = get_match(leaf1)
        counts_leaf1 = float(leaf1[leaf1.index(count_pattern)+count_index : leaf1.index(':')])

        leaf2 = leaves.split(',')[1]
        count_pattern, count_index = get_match(leaf2)
        counts_leaf2 = float(leaf2[leaf2.index(count_pattern)+count_index : leaf2.index(':')])

        counts_leaves = [counts_leaf1, counts_leaf2]
    else:
        counts_leaves = [float(leaves[leaves.index(count_pattern)+count_index : leaves.index(':')])]
    countsDict = {'leaves': counts_leaves, 'node': counts_node}

    return countsDict


# Function to add event counts according to the event's nature (EXP, CON or STA).
# Here the Stabilities are filtered as to not include the null stabilities (i.e. children nodes that maintain 0 gene copy numbers like the parent node).
def addEventCounts(countsDict, eventCounts_Dict):

    for counts_leaf in countsDict['leaves']:
        if counts_leaf > countsDict['node']:
            eventCounts_Dict['n_expansions'] += 1
        elif counts_leaf == countsDict['node'] and counts_leaf != 0:
            eventCounts_Dict['n_stabilities'] += 1
        elif counts_leaf < countsDict['node']:
            eventCounts_Dict['n_contractions'] += 1

    return eventCounts_Dict


# Regular expression function to extract a signle node-children subtree from each tree-line in the NEXUS cafe5 output file.
# The goal of this function is to parse the whole tree by iteratevily removing each single node-children subtree so
# to extract gene copy number counts at each bifurcation. 
def getCAFEevents(CAFEreportInputFile):

    lines = readCAFEreport(CAFEreportInputFile)
    CAFEeventCountsDict = defaultdict(dict)

    for line in lines:

        if 'TREE ' not in line:  ## Checking that we are dealing with a TREE line of the NEXUS file, and not a formatting line
            continue

        splitline = line.strip().split(' = ')
        orthogroup = splitline[0].split(' ')[1]
        tree = splitline[1]

        eventCountsDict = {'n_expansions': 0, 'n_stabilities': 0, 'n_contractions': 0}

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

        CAFEeventCountsDict[orthogroup] = eventCountsDict

    return CAFEeventCountsDict


copy_number_variations = getCAFEevents(input_file)


# Process output files
output_file.write('orthogroup\tEXP\tSTA\tCON\n')

for orthogroup in sorted(copy_number_variations.keys()):
    output_file.write(f"{orthogroup}\t{copy_number_variations[orthogroup]['n_expansions']}\t{copy_number_variations[orthogroup]['n_stabilities']}\t{copy_number_variations[orthogroup]['n_contractions']}\n")

# Close files
output_file.close()
