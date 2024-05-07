from typing import List, Union
import numpy.typing as npt
import numpy as np


def base_count(fastafile: str) -> List[int]:
    # 課題 1-1
    # read and skip first line
    i = 0
    s = ''
    with open('1_sequence_processing/'+fastafile) as f:
        s = f.read()
    for c in s:
        if c == '\n':
            break
        i+=1
    # count
    res = [0, 0, 0, 0]
    offset = 0
    for i in range(i,len(s)):
        c = s[i]
        if c == 'A':
            offset = 0
        elif c == 'T':
            offset = 1
        elif c == 'G':
            offset = 2
        elif c == 'C':
            offset = 3
        else:
            continue
        res[offset] = res[offset]+1


    return res # A, T, G, C

def gen_rev_comp_seq(fastafile: str) -> str:
    # 課題 1-2
    # read and skip first line
    i0 = 0
    s = ''
    with open('1_sequence_processing/'+fastafile) as f:
        s = f.read()
    for c in s:
        if c == '\n':
            break
        i0+=1
    # rev_comp
    res = ""
    for i in range(i0,len(s)):
        c = s[i]
        if c == 'A':
            res = 'T'+res
        elif c == 'T':
            res = 'A'+res
        elif c == 'G':
            res = 'C'+res
        elif c == 'C':
            res = 'G'+res
        else:
            continue
    return res

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    # read and skip first line
    i0 = 0
    s = ''
    with open('1_sequence_processing/'+fastafile) as f:
        s = f.read()
    for c in s:
        if c == '\n':
            break
        i0+=1
    temp = s[i0+1:]
    s = ''
    for c in temp:
        if c != '\n':
            s+=c
    cnt = 0
    w_cnt = int((len(s)-window)/step)
    if w_cnt < 0:
        return []
    res = []

    # do double loop
    for i in range(0,w_cnt+1):
        for j in range(0,window):
            c = s[i*step+j]
            if c == 'G' or c == 'C':
                cnt+=1
        res.append(cnt/window*100)
        cnt = 0
    return res

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4

    # read and skip first line
    i0 = 0
    s = ''
    with open('1_sequence_processing/'+fastafile) as f:
        s = f.read()
    for c in s:
        if c == '\n':
            break
        i0+=1
    temp = s[i0+1:]
    s = ''
    for c in temp:
        if c != '\n':
            s+=c
    # get reverse complemetary 
    rev_com = ""
    for c in motif:
        if c == 'A':
            rev_com = 'T'+rev_com
        elif c == 'T':
            rev_com = 'A'+rev_com
        elif c == 'G':
            rev_com = 'C'+rev_com
        elif c == 'C':
            rev_com = 'G'+rev_com
        else:
            continue

    res = []
    # search from front
    start = s.find(motif)
    while start>=0:
        res.append('F'+str(start+1))
        start = s.find(motif,start+len(motif))
    start = 0
    # search from back with reverse sequence
    start = s.rfind(rev_com)
    while start>=0:
        res.append('R'+str(start+len(motif)))
        start = s.rfind(rev_com,0,start)
    return res

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    # read and skip first line
    i0 = 0
    s = ''
    with open('1_sequence_processing/'+fastafile) as f:
        s = f.read()
    for c in s:
        if c == '\n':
            break
        i0+=1
    temp = s[i0+1:]
    s = ''
    for c in temp:
        if c != '\n':
            s+=c

    # get reverse complemetary 
    rev_com = ""
    for c in s:
        if c == 'A':
            rev_com = 'T'+rev_com
        elif c == 'T':
            rev_com = 'A'+rev_com
        elif c == 'G':
            rev_com = 'C'+rev_com
        elif c == 'C':
            rev_com = 'G'+rev_com
        else:
            continue
    
    # T:0 C:1 A:2 G:3
    dna_code=[
        [
        ['F','F','L','L'],['S','S','S','S'],['Y','Y','_','_'],['C','C','_','W']
        ],
        [
        ['L','L','L','L'],['P','P','P','P'],['H','H','Q','Q'],['R','R','R','R']
        ],
        [
        ['I','I','I','M'],['T','T','T','T'],['N','N','K','K'],['S','S','R','R']
        ],
        [
        ['V','V','V','V'],['A','A','A','A'],['D','D','E','E'],['G','G','G','G']
        ],

    ]
    res = []
    def to_number(c):
        if c == 'T':
            return 0
        elif c == 'C':
            return 1
        elif c == 'A':
            return 2
        elif c == 'G':
            return 3
        return -1
    def to_amino_acid(str):
        return dna_code[to_number(str[0])][to_number(str[1])][to_number(str[2])]
    # search original sequence
    start = s.find('ATG')
    while start>=0:
        res_s = 'M'
        start = start+3
        while start+3<len(s):
            aa=to_amino_acid(s[start:start+3])
            res_s+=aa
            start = start+3
            if aa == '_':
                break
        res.append(res_s)
        start = s.find('ATG',start)
    # search reverse complemetary sequence
    start = rev_com.find('ATG')
    while start>=0:
        res_s = 'M'
        start = start+3
        while start+3<len(rev_com):
            aa=to_amino_acid(rev_com[start:start+3])
            res_s+=aa
            start = start+3
            if aa == '_':
                break
        res.append(res_s)
        start = rev_com.find('ATG',start)
    return res

if __name__ == "__main__":
    filepath = "data/NT_113952.1.fasta"
    # 課題 1-1
    print(base_count(filepath))
    # 課題 1-2
    print(gen_rev_comp_seq(filepath))
    # 課題 1-3
    print(calc_gc_content(filepath))
    # 課題 1-4
    print(search_motif(filepath, "ATG"))
    # 課題 1-5
    print(translate(filepath))
