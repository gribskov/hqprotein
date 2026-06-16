"""=====================================================================================================================
clusters.py
Read in clusters for each nid. Check each isoform group to see what cluster the members belong to

2026-06-16 gribskov
====================================================================================================================="""
import sys

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    cl_file = sys.argv[1]
    cl = open(cl_file, 'r')
    print(f'cluster file: {cl_file}')

    # skip header
    for line in cl:
        if line.startswith('begin'): break

    nid2og = {}
    nid_n = 0
    og_member = {}
    og_n = 0
    for line in cl:
        # data ends opens with (mclmatrix\nbegin and ends with )
        if line.startswith(')'): break

        tag = line[:7].rstrip()
        if tag:
            og = int(tag)
            og_member[og] = []
            og_n += 1
        print(f'{og:<7d} {line[7:-1]}')

        id = line[7:].rstrip().split()
        for this_id in id:
            if this_id.startswith('$'): continue

            # print(f'id: |{this_id}|\tog: {og}')
            nid2og[this_id] = og
            og_member[og].append(this_id)
            nid_n += 1

    cl.close()
    print(f'genes: {nid_n}')
    print(f'orthogroups: {og_n}')

    exit(0)
