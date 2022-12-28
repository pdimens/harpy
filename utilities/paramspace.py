#! /usr/bin/env python3

with open('params.tsv', "w") as file:
        file.write(
'''
model   useBX   k   s   nGen
pseudoHaploid   TRUE    10  5   50
pseudoHaploid   TRUE    10  1   50
pseudoHaploid   TRUE    15  10  100
'''
    )