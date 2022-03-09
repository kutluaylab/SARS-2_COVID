import sys

inFile = open(sys.argv[1],'r')
sense_counts=['bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous']
antisense_counts=['bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous']
combined_counts=['bp\tA\tG\tC\tT\tdel\tins\tinserted\tambiguous']

for line in inFile:
        data = line.strip().split('\t')
        bp = data[1]
        bases = data[4].upper()
        ref = data[2].upper()

        types = {'sense':{'A':0,'G':0,'C':0,'T':0},'antisense':{'A':0,'G':0,'C':0,'T':0},'combined':{'A':0,'G':0,'C':0,'T':0},'-':0,'+':[],'X':[]}

        i = 0
        while i < len(bases):
                base = bases[i]
                if base == '^' or base == '$':
                        i += 1
                elif base == '-':
                        i += 1
                elif base == '*':
                        types['-'] += 1
                elif base == '+':
                        i += 1
                        addNum = int(bases[i])
                        addSeq = ''
                        for a in range(addNum):
                                i += 1
                                addSeq += bases[i]

                        types['+'].append(addSeq)
                elif base == '.':
                        types['sense'][ref] += 1
                        types['combined'][ref] += 1
                elif base == ',':
                        types['antisense'][ref] += 1
                        types['combined'][ref] += 1
                else:
                        if types['combined'].has_key(base):
                                types['combined'][base] += 1
                        else:
                                types['X'].append(base)

                i += 1

        adds = '.'
        if len(types['+']) > 0:
                adds = ','.join(types['+'])
                print adds

        amb = '.'
        if len(types['X']) > 0:
                amb = ','.join(types['X'])
 
        sense_out = [bp,types['sense']['A'],types['sense']['G'],types['sense']['C'],types['sense']['T'],types['-'],len(types['+']),adds,amb]
        antisense_out = [bp,types['antisense']['A'],types['antisense']['G'],types['antisense']['C'],types['antisense']['T'],types['-'],len(types['+']),adds,amb]
        combined_out = [bp,types['combined']['A'],types['combined']['G'],types['combined']['C'],types['combined']['T'],types['-'],len(types['+']),adds,amb]

        sense_counts.append('\t'.join([str(x) for x in sense_out]))
        antisense_counts.append('\t'.join([str(x) for x in antisense_out]))
        combined_counts.append('\t'.join([str(x) for x in combined_out]))
        #print'\t'.join([str(x) for x in out])

with open("sense_counts", "w") as outfile:
    outfile.write("\n".join(sense_counts))

with open("antisense_counts", "w") as outfile:
    outfile.write("\n".join(antisense_counts))

with open("combined_counts", "w") as outfile:
    outfile.write("\n".join(combined_counts))
