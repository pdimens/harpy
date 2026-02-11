"""
Functions for GIH-haplotagging preprocessing
"""

def needs_stagger(filename):
    counts = {}
    total_count = 0
    in_col3 = False

    with open(filename, 'r') as file:
        for line in file:
            line = line.strip()
            if line == 'col3=startpost':
                in_col3 = True
                continue
            elif line == 'col4=endpos':
                in_col3 = False
                continue

            if in_col3:
                try:
                    count, position = map(int, line.split())
                except ValueError:
                    continue  # Skip lines that don't have two integers
                counts[position] = counts.get(position, 0) + count
                total_count += count

    if 51 in counts and counts[51] > total_count / 2:
        return False
    else:
        return True
