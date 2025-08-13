#!/usr/bin/env python
#
# (c) 2017 by Joffrey Fitz (joffrey.fitz@tuebingen.mpg.de),
# Max Planck Institute for Developmental Biology,
# Tuebingen, Germany
#

import h5py
import numpy


f = h5py.File('/datapool/home/2023102768/lico_share_dir/life-gongl/zesheng/Arabidopsis_thaliana/list/snp/1001_SNP_MATRIX/imputed_snps_binary.hdf5','r')

# Coordinates for At1g01070.1
start_pos = 16568746
end_pos = 16575692

# Get all SNP positions for all chromosomes (len=10709949)
positions = f['positions'][:]

# Array of tupels with start/stop indices for each chromosome
chr_regions = f['positions'].attrs['chr_regions']
indices_for_chr2 = chr_regions[1]

# Subset positions for SNPs on Chr1
positions_on_chr2 = positions[indices_for_chr2[0]:indices_for_chr2[1]]


# Index for start_pos/end_pos
idx_start_pos = numpy.where( numpy.logical_and(positions_on_chr2>=start_pos, positions_on_chr2<=end_pos) )[0][0]
idx_end_pos = numpy.where( numpy.logical_and(positions_on_chr2>=start_pos, positions_on_chr2<=end_pos) )[0][-1]

# Get SNPs
snps_in_region = f['snps'][idx_start_pos:idx_end_pos]

# Get index of specific accessions
accs = f['accessions'][:]

# Open a file to save the results
with open('snp_results.txt', 'w') as file:
    # Write the header line
    file.write("position\tsample\tSNP\n")

    # Loop through positions from start_pos to end_pos
    for position in range(start_pos, end_pos + 1):
        try:
            # Find index of a specific position 
            ix = numpy.where(positions_on_chr2 == position)[0][0]
            # Get the corresponding SNP for that position
            snp = f['snps'][ix]
            for acc in accs:
                ix_of_acc = numpy.where(accs == acc)[0][0]
                # Get SNP for pos and acc
                snp_for_specific_pos_and_acc = snp[ix_of_acc]
                acc_str=acc.decode('utf-8')
                result = f"{position}\t{acc_str}\t{snp_for_specific_pos_and_acc}"
                print(result)
                file.write(result + '\n')
        except IndexError:
            result = f"Position {position} not found in chromosome 2."
            print(result)


