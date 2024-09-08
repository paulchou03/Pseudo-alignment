from collections import defaultdict
from Bio import SeqIO
import pandas as pd

#kmer length
kmer = 21

#referce complement function
def reverse_complement(reads):
    reversed = reads[::-1]
    rc = ""
    for i in reversed:
        if i == 'A':
            rc += 'T'
        elif i == 'T':
            rc += 'A'
        elif i == 'C':
            rc += 'G'
        elif i == 'G':
            rc += 'C'
    return rc


#creating the index (hashmap) here of kmers of the transcripts
map_transcript = defaultdict(set)
for line in SeqIO.parse("reference.fasta", "fasta"):
        transcript = line.id
        t_sequence = str(line.seq)

        for i in range(len(t_sequence) - kmer + 1):
            t_read = t_sequence[i:i + kmer]
            map_transcript[t_read].add(transcript)
            rc_read = reverse_complement(t_read)
            map_transcript[rc_read].add(transcript)

readlist = []
for row in SeqIO.parse("reads.fasta", "fasta"):
    read = str(row.seq)
    if 'N' in read:
        continue #skipping the entire read when we see 'N'. this is one of the assumption we make
    else:
        readlist.append(read)

#setting up the tsv
column_names = ['counts', 'number of items in equivalence class', 'isoforms in equivalence class']
df = pd.DataFrame(columns=column_names)

#for every read in the readlist
for i in readlist:
    list1 = []
    #for every kmer in the read
    for k in range(len(i) - kmer + 1):
        list2 = []
        k_read = i[k:k + kmer]
        if k_read in map_transcript:
            list2.extend(map_transcript[k_read])
            list1.append(list2)
        else:
            list1 = [] #there is no isoform so the intersection will just be blank, so we can just break
            break
    if list1: #if not empty
        sets = map(set, list1)
        intersection = set.intersection(*sets)
        intersection_list = list(intersection)
        item_nums = len(intersection_list) #finds the intersection of the kmers of the read
        if item_nums > 0:
            isoform = ','.join(intersection_list)
            if not df.empty and (df['isoforms in equivalence class'] == isoform).any(): #if the eq class already exists in the output, add 1 to counter
                index = df.index[df['isoforms in equivalence class'] == isoform]
                df.loc[index, 'counts'] += 1
            else: #add the row to the result
                df = pd.concat([df, pd.DataFrame([{'counts': 1, 'number of items in equivalence class': item_nums,
                                                   'isoforms in equivalence class': isoform}])], ignore_index=True)

    else: #if empty
        if not df.empty and (df['isoforms in equivalence class'] == 'NA').any():
            index = df.index[df['isoforms in equivalence class'] == 'NA']
            df.loc[index, 'counts'] += 1
        else:
            df = pd.concat([df, pd.DataFrame([{'counts': 1, 'number of items in equivalence class': 0, 'isoforms in equivalence class': 'NA'}])], ignore_index = True)

    df = df.sort_values(by='counts', ascending=False).reset_index(drop=True)
    df.to_csv('output.tsv', sep='\t', index=False)


