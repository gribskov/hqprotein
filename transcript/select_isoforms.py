"""=====================================================================================================================
select_isoforms.py

Construct a table of multiple predicted isofroms based on Orthofinder WorkingDirectory/SequenceIDs
Format:
gene_gid: transcript_id
0_0: g1.t1
0_1: g2.t1
...
0_13: g14.t1
0_14: g14.t2
0_15: g15.t1

row number gives the ID in sequential numbering (used in some orthofinder outputs)
gene_id: species_gid - sequential numbering withing species
transcript.id: gene.transcript given in original GFF input (in this case braker3)

Usage:
select_isoforms.py SequenceIDs.txt                      # writes yaml to stdout
select_isoforms.py SequenceIDs.txt isoform.list.yaml    # writes yaml to isoform.list.yaml

2026-06-16 gribskov
====================================================================================================================="""
import sys
import yaml

# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    seqname = sys.argv[1]
    seq = open(seqname, 'r')
    print(f'sequence IDs: {seqname}')

    if len(sys.argv) < 3:
        # no output file given, write to stdout
        outfile = sys.stdout
    else:
        outfile = open(sys.argv[2], 'w')

    nid = 0
    gene_n = 0
    tgroup = 0
    ilist = []
    isoformlist = {}
    prev = ''
    for line in seq:
        nid += 1
        gid, tid = line.rstrip().split()
        gid = gid.replace(':', '')
        gene, isoform = tid.split('.')

        if gene == prev:
            # gene is the same so this is an isoform
            ilist.append({'nid': nid, 'gid': gid, 'tid': tid})
        else:
            # new gene is a different isoform than the one in ilist
            gene_n += 1
            if len(ilist) > 1:
                # isoformlist.append(ilist)
                isoformlist[f'tgroup_{tgroup}'] = ilist
                tgroup += 1
            ilist = [{'nid': nid, 'gid': gid, 'tid': tid}]

        prev = gene
        print(f'{nid}\t{gid}\t{tid}')

    # if len(ilist): isoformlist.append(ilist)
    if len(ilist): isoformlist[f'tgroup_{tgroup}'] = ilist
    seq.close()

    print(f'\nisoforms read: {nid}')
    print(f'genes: {gene_n}')
    print(f'multiple isoforms: {len(isoformlist)} ({len(isoformlist) / gene_n * 100:.2f}%)')

    yaml.dump(isoformlist, outfile, explicit_start=True, default_flow_style=None, sort_keys=False)
    outfile.close()

    exit(0)
