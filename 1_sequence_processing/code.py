from typing import List, Union
import numpy.typing as npt
import numpy as np


def base_count(fastafile: str) -> List[int]:
    # 課題 1-1
    # skip first line
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
    return ""

def calc_gc_content(fastafile: str, window: int=1000, step: int=300) -> Union[npt.NDArray[np.float_], List[float]]:
    # 課題 1-3
    # 値を出力するところまで。matplotlibを使う部分は別途実装してください。
    return []

def search_motif(fastafile: str, motif: str) -> List[str]:
    # 課題 1-4
    return []

def translate(fastafile: str) -> List[str]:
    # 課題 1-5
    return []

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
