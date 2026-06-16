"""=====================================================================================================================
clusters.py
Read in clusters for each nid. Check each isoform group to see what cluster the members belong to

2026-06-16 gribskov
====================================================================================================================="""
import sys


def cluster_read(cluster_fname):
    """-----------------------------------------------------------------------------------------------------------------
    Read the cluster information from Orthofinder output (workingDirectory/clusters_Orthofinder_I1.5.txt)
    Format
    # cline: mcl /scratch/negishi/nngubane/Orthofinder/OF_39249197/annotated/OrthoFinder/Results_Jun14/WorkingDirectory/OrthoFinder_graph.txt "-I" "1.5" "-o" "/scratch/negishi/nngubane/Orthofinder/OF_39249197/annotated/OrthoFinder/Results_Jun14/WorkingDirectory/clusters_OrthoFinder_I1.5.txt" "-te" "12" "-V" "all"
    (mclheader
    mcltype matrix
    dimensions 251971x16582
    )
    (mclmatrix
    begin
    0      6641 14275 18284 18285 25672 25673 25674 36351 36352 41597
             52809 57742 61006 61007 62424 64419 65117 65466 66472 66473
             ...
             201639 201640 202525 204837 204838 206386 208251 214831
             226400 226402 226403 234409 242334 247141 $
    1      20036 27395 36379 91902 92538 93489 99739 100492 100493 100494
            ...

    :param cluster_fname:
    :return: nid2og     dict with nid as key and og as value
             og_member  dict with og as key and list of nids of group emebers as values
    -----------------------------------------------------------------------------------------------------------------"""
    cl = open(cluster_fname, 'r')

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
        # print(f'{og:<7d} {line[7:-1]}')

        id = line[7:].rstrip().split()
        for this_id in id:
            if this_id.startswith('$'): continue

            # print(f'id: |{this_id}|\tog: {og}')
            nid2og[this_id] = og
            og_member[og].append(this_id)
            nid_n += 1

    cl.close()
    return nid2og, og_member


# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    cl_file = sys.argv[1]

    nid2og, og_member = cluster_read(cl_file)
    print(f'cluster file: {cl_file}')

    print(f'genes: {len(nid2og)}')
    print(f'orthogroups: {len(og_member)}')

    exit(0)
