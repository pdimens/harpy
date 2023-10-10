#! /usr/bin/env python3

with open('stitch.params', "w") as file:
        file.write(
'''
model\tuseBX\tbxlimit\tk\ts\tnGen
diploid\tTRUE\t10\t100000\t5\t50
diploid\tTRUE\t10\t100000\t1\t50
diploid\tTRUE\t15\t100000\t10\t100
'''
    )