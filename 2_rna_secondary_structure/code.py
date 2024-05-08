from typing import List, Tuple, Union
import numpy.typing as npt
import numpy as np

def isPair(a,b):
    def to_number(a):
        if a=='A':
            return 0
        elif a=='U' or a=='T':
            return 4
        elif a=='G':
            return 1
        elif a=='C':
            return 3
        return -1
    return (to_number(a)+to_number(b)) == 4
def readFile(fastafile):
    # read and skip first line
    i0 = 0
    s = ''
    with open('2_rna_secondary_structure/'+fastafile) as f:
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
    return s
def enumerate_pairs(fastafile: str) -> List[Tuple[int, int]]:
    # 課題 2-1
    s = readFile(fastafile)
    res = []
    for i in range(0,len(s)):
        for j in range(i+1,len(s)):
            if isPair(s[i],s[j]):
                res.append((i+1,j+1))
    return res

def enumerate_possible_pairs(fastafile: str, min_distance: int=4) -> List[Tuple[int, int]]:
    # 課題 2-2
    s = readFile(fastafile)
    # 
    res = []
    for i in range(0,len(s)):
        for j in range(i+min_distance,len(s)):
            if isPair(s[i],s[j]):
                res.append((i+1,j+1))
    return res

def enumerate_continuous_pairs(fastafile: str, min_distance: int=4, min_length: int=2) -> List[Tuple[int, int, int]]:
    # 課題 2-3
    s = readFile(fastafile)
    res = []
    for i in range(0,len(s)):
        for j in range(i+min_distance,len(s)):
            l = int((j-i+1)/2)
            con = 0
            # check outter
            if i-1>=0 and j+1<len(s):
                if isPair(s[i-1],s[j+1]):
                    # already contained in previous result
                    continue
            # check inner
            for k in range(0,l):
                if isPair(s[i+k],s[j-k]):
                    con += 1
                else:
                    break
            if con >= min_length:
                res.append((i+1,j+1,con))
    return res

def create_dotbracket_notation(fastafile: str, min_distance: int=4, min_length: int=2) -> str:
    # 課題 2-4
    pairs = enumerate_continuous_pairs(fastafile,min_distance,min_length)
    s = readFile(fastafile)
    res = ['.']*len(s)
    # 
    def isInRange(r,p):
        return r[0]<=p[0]-1 and r[1]>p[1]-1
    range_stack = [(int(0),len(s))]
    def onRightOfRange(r,p):
        return r[1]<=p[0]-1
    p_index = 0
    while len(range_stack) > 0 and p_index < len(pairs):
        r = range_stack.pop()
        while p_index < len(pairs):
            p = pairs[p_index]
            p_index += 1
            if isInRange(r,p):
                # make pairs
                for i in range(0,p[2]):
                    res[p[0]-1+i] = '('
                    res[p[1]-1-i] = ')'
                # push new range
                # right
                if r[1]-p[1] >=2:
                    range_stack.append((p[1],r[1]))
                # middle
                if p[1]-p[0]+1-p[2]*2 >= 2:
                    range_stack.append((p[0]-1+p[2],p[1]-1-p[2]+1))
                break
            elif onRightOfRange(r,p):
                p_index -= 1
                break

    return ''.join(res)

if __name__ == "__main__":
    filepath = "data/AUCGCCAU.fasta"
    # 課題 2-1
    print(enumerate_pairs(filepath))
    # 課題 2-2
    print(enumerate_possible_pairs(filepath))
    # 課題 2-3
    print(enumerate_continuous_pairs(filepath, 2))
    # 課題 2-4
    print(create_dotbracket_notation(filepath, 2))


