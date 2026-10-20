#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import itertools
import os

from dipy.io.streamline import load_tractogram
import nibabel as nib
import numpy as np
from scipy.ndimage import map_coordinates

from scilpy.io.utils import (add_overwrite_arg,
                             add_verbose_arg,
                             add_reference_arg,
                             assert_headers_compatible,
                             assert_inputs_exist,
                             assert_outputs_exist)
from scilpy.tractanalysis.streamlines_metrics import compute_tract_counts_map


def _build_arg_parser():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument('in_tractogram',
                   help='Input tractogram file to decompose')
    p.add_argument('in_labels',
                   help='Input labels file containing the connectivity signatures.')
    p.add_argument('out_dir',
                   help='Output directory to save the decomposed TDI files.')
    add_reference_arg(p)
    add_overwrite_arg(p)
    add_verbose_arg(p)

    return p


def main():
    parser = _build_arg_parser()
    args = parser.parse_args()
    assert_inputs_exist(parser, [args.in_tractogram, args.in_labels])

    sft = load_tractogram(args.in_tractogram, 'same')
    labels_img = nib.load(args.in_labels)
    assert_headers_compatible(sft, [labels_img])
    labels = labels_img.get_fdata().astype(np.uint16)

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]
    comb_list = list(itertools.combinations(unique_labels, 2))
    comb_list.extend([(i, i) for i in unique_labels])

    all_out_files = [os.path.join(args.out_dir, f'{i}_{j}.nii.gz')
                     for i, j in comb_list]
    assert_outputs_exist(parser, args, all_out_files)

    sft.to_vox()
    sft.to_corner()

    streamlines = sft.streamlines
    n_streamlines = len(streamlines)

    if n_streamlines > 0:
        seq = streamlines
        start_vox = seq._data[seq._offsets].astype(np.int16)
        end_vox = seq._data[seq._offsets + seq._lengths - 1].astype(np.int16)

        labels_start = map_coordinates(labels, start_vox.T, order=0, mode='nearest')
        labels_end = map_coordinates(labels, end_vox.T, order=0, mode='nearest')

        min_labels = np.minimum(labels_start, labels_end)
        max_labels = np.maximum(labels_start, labels_end)
    else:
        min_labels = np.array([], dtype=np.uint16)
        max_labels = np.array([], dtype=np.uint16)

    for i, j in comb_list:
        indices = np.where((min_labels == i) & (max_labels == j))[0]
        if len(indices) > 0:
            filtered_sft = sft[indices]
            density = compute_tract_counts_map(filtered_sft.streamlines,
                                               sft.dimensions)
        else:
            density = np.zeros(sft.dimensions, dtype=np.int32)

        nib.save(nib.Nifti1Image(density.astype(np.int32), sft.affine),
                 os.path.join(args.out_dir, f'{i}_{j}.nii.gz'))


if __name__ == '__main__':
    main()
