"""=====================================================================================================================
splice_scan.py
use splice donor and acceptor pssms to predict splice sites in gene regions

2026-03-21 gribskov
====================================================================================================================="""
from collections import defaultdict
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from include.gff.gff2 import GxfSet
from include.sequence.fasta import Fasta
from pssm import PSSM


def pssm_scan(site, sequence):
    """-----------------------------------------------------------------------------------------------------------------
    calculate the match score for a pssm at each position along the sequence. score is sum of pssm values
    TODO should it be sum of logs?

    :param site: PSSM object    a site defined as a position specific scoring matrix
    :param sequence: string     sequence of a DNA region to analyze
    :return: list               positional scores
    -----------------------------------------------------------------------------------------------------------------"""
    pssm = site.matrix
    pssmlen = site.cols

    alphabet = 'ACGT'
    a2i = {alphabet[i]: i for i in range(len(alphabet))}

    # ap bases to indices, produces an integer list where each position is the index in a2i
    seq_indices = [a2i[base] for base in sequence]
    seq_len = len(sequence)

    # one-hot encoding for DNA sequence
    one_hot = np.zeros((seq_len, 4))
    one_hot[np.arange(seq_len), seq_indices] = 1

    # 2D convolution to apply PSSM efficiently, PSSM needs to be flipped for convolution (correlation)
    scores = np.zeros(seq_len - pssmlen + 1)
    for i in range(len(site.rows)):
        # scores += np.convolve(one_hot[:, i], pssm[i, :][::-1], mode='valid')
        scores += np.correlate(one_hot[i, :], pssm[i, :], mode='valid')

    return scores


def score_vectorized(dna_seq, pssm):
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    seq_indices = np.array([mapping[base] for base in dna_seq])

    # Create windows of the sequence indices
    # (e.g., if motif is len 3, windows are [0,1,2], [1,2,3], etc.)
    motif_len = pssm.shape[1]
    shape = (seq_indices.size - motif_len + 1, motif_len)
    strides = (seq_indices.strides[0], seq_indices.strides[0])
    windows = np.lib.stride_tricks.as_strided(seq_indices, shape=shape, strides=strides)

    # Advanced indexing to pull PSSM values directly
    # This effectively "looks up" the score for each base in each window
    pos_idx = np.arange(motif_len)
    scores = np.sum(pssm[windows, pos_idx], axis=1)

    return scores


def score_sequence_with_correlate(dna_seq, pssm):
    """
    Computes PSSM scores using np.correlate (no flipping required).
    pssm: 4xL matrix
    dna_seq: String of DNA characters
    """
    mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    seq_indices = [mapping[base] for base in dna_seq]

    # Create one-hot matrix (4, seq_len) to match PSSM rows
    one_hot = np.zeros((4, len(dna_seq)))
    one_hot[seq_indices, np.arange(len(dna_seq))] = 1

    # We sum the correlations of each row
    # mode='valid' ensures we only get scores where the PSSM fits fully
    motif_len = pssm.shape[1]
    scores = np.zeros(len(dna_seq) - motif_len + 1)

    for i in range(4):
        # np.correlate does NOT flip the kernel
        scores += np.correlate(one_hot[i, :], pssm[i, :], mode='valid')

    return scores


# ======================================================================================================================
# Main
# ======================================================================================================================
if __name__ == '__main__':
    gff_file = 'data/z.tritici.IP0323.reannot.gff3'
    gff = GxfSet(file=gff_file, fmt='gff')
    genome = Fasta(filename='data/z.tritici.IP0323.fasta')
    donorfile = 'donor_zt.dat'
    acceptorfile = 'acceptor_zt.dat'

    padding = 500
    donor = PSSM()
    donor.read(donorfile)
    donor.logtransform()
    acceptor = PSSM()
    acceptor.read(acceptorfile)
    acceptor.logtransform()

    ####################################################################################################################
    # read exon features from GFF file. Store in transcript_set indexed by transcript ID
    ####################################################################################################################
    feature_n = gff.feature_get(['exon'])
    junction = defaultdict(list)
    feature_count = 0
    current = ''
    transcript_set = defaultdict(list)
    for feature in gff.features:
        # find exons and store grouped by transcript
        feature_count += 1
        # limit input for debug
        if feature_count > 50: break

        parent = feature.attribute['Parent']
        transcript_set[parent].append(feature)

    print(f'\nexons found: {feature_count}')
    print(f'transcript sets: {len(transcript_set)}')

    ####################################################################################################################
    # read one sequence at a time from the reference genome file, and process transcripts in that sequence
    # for initial testing just find a single + strand transcript with multiple exons
    ####################################################################################################################
    # for sequence in genome:
    #     print(f'\nProcessing {sequence.id}')
    #     for tid in transcript_set:
    #         transcript = transcript_set[tid]
    #         if len(transcript) < 2 or transcript[0].strand == '-':
    #             # TODO remove - for testing/debugging only
    #             continue
    #
    #         print(f'tid: {tid} \t strand: {transcript[0].strand}\t exons: {len(transcript)}')
    #         begin = transcript[0].start
    #         end = transcript[0].end
    #         print(f'\t{transcript[0].start}\t{transcript[0].end}\t{transcript[0].strand}')
    #         for exon in transcript[1:]:
    #             print(f'\t{exon.start}\t{exon.end}\t{exon.strand}')
    #             begin = min(begin, exon.start)
    #             end = max(end, exon.end)
    #
    #         print(f'\tbegin: {begin} \t end: {end}')
    #         seq = sequence.seq[begin-padding:end+padding]
    #         print(seq)
    #         pssm_scan(donor, seq)
    #     # TODO remove - for testing only, limit to one transcript
    #     break
    #
    # # TODO remove - for testing only, limit to one chromosome
    # break

    with open('data/transcript.seq', 'r') as f:
        seq = f.read().rstrip()
    # pssm_scan(donor, seq)
    pdonor = score_sequence_with_correlate(seq, donor.matrix)
    pacceptor = score_sequence_with_correlate(seq, acceptor.matrix)

    # Plot probabilities
    pdonor = np.power(2.0, pdonor)
    pacceptor = np.power(2.0, pacceptor)
    pssmlen = donor.cols
    xlen = len(seq) - pssmlen + 1
    x = np.linspace(0, xlen, xlen)
    # y = np.sin(x)

    # figure and axes
    fig, ax = plt.subplots()
    line, = ax.plot(x, pdonor, label='donor')
    line, = ax.plot(x, pacceptor, color='red', label='acceptor')

    # Relimit the axes based on the current data in the artists
    ax.relim()
    ax.autoscale_view()

    # chr_1	ingenannot	mRNA	116892	117149	.	+	.	ID=ZtIPO323_000030.1;Name=ZtIPO323_000030.1;Parent=ZtIPO323_000030;ev_tr=None;aed_ev_tr=1.0000;ev_tr_penalty=undef;ev_pr=None;aed_ev_pr=1.0000;ev_lg=None;aed_ev_lg=1.0000;ev_lg_penalty=undef;product=unnamed protein product;Dbxref=InterPro:-,MobiDBLite:mobidb-lite;locus_tag=ZtIPO323_000030;
    #
    # chr_1	ingenannot	exon	116892	116904	.	+	.	ID=exon:ZtIPO323_000030.1;Parent=ZtIPO323_000030.1;locus_tag=ZtIPO323_000030;
    #
    # chr_1	ingenannot	exon	116959	117149	.	+	.	ID=exon:ZtIPO323_000030.2;Parent=ZtIPO323_000030.1;locus_tag=ZtIPO323_000030;
    # begin 500 donor 512 acceptor 567 begin=500, end =757
    # plt.vlines(x=[495, 507, 757], ymin=0, ymax=1e-7, colors='black', lw=2)
    # fig.patches.extend([plt.Rectangle((500,0.5e-7), 12, 1e-7,
    #                                   fill=True, color='black', alpha=0.4, figure=fig)])

    rectangle = Rectangle((495,5.74e-18), 12, -1e-8,fill=True, color='black', alpha=0.2,
        linewidth=0,
        edgecolor='none',
        facecolor='none'
    )
    ax.add_patch(rectangle)
    rectangle = Rectangle((557, 5.74e-18), 257, -1e-8, fill=True, color='black', alpha=0.2,
                          linewidth=0,
                          edgecolor='none',
                          facecolor='none'
                          )
    ax.add_patch(rectangle)
    plt.title("Splice Donor/Acceptor sites")
    plt.xlabel("Sequence Position")
    plt.ylabel("P(site|sequence)")
    plt.legend()
    plt.show()

    exit(0)
