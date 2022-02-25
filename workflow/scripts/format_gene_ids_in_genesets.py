with open('genelists.txt') as inf:
    lines = inf.readlines()
with open('genelists_2.txt', 'w') as out:
    out.write(lines[0])
    for line in lines[1:]:
        new_line = []
        for b in line.split('\t'):
            if b != '' and b != '\n':
                c = '7460:' + b
                new_line.append(c)
            else:
                new_line.append(b)
        out.write('\t'.join(new_line))
