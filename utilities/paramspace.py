#! /usr/bin/env python3

with open('stitch.params', "w") as file:
        file.write(
'''
model\tuseBX\tk\ts\tnGen
pseudoHaploid\tTRUE\t10\t5\t50
pseudoHaploid\tTRUE\t10\t1\t50
pseudoHaploid\tTRUE\t15\t10\t100
'''
    )