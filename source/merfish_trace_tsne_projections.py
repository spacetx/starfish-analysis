#!/usr/bin/env python

import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from starfish import data
from starfish.types import Features, Indices
from starfish.image import Filter
from starfish.spots import SpotFinder
from copy import deepcopy
from MulticoreTSNE import MulticoreTSNE as TSNE

use_test_data = os.getenv("USE_TEST_DATA") is not None
experiment = data.MERFISH(use_test_data=use_test_data)
primary_image = experiment.fov().primary_image
primary_image.show_stack({Indices.CH: 0})

ghp = Filter.GaussianHighPass(sigma=3)
high_passed = ghp.run(primary_image, verbose=True, in_place=False)

dpsf = Filter.DeconvolvePSF(num_iter=15, sigma=2, clip=True)
deconvolved = dpsf.run(high_passed, verbose=True, in_place=False)

glp = Filter.GaussianLowPass(sigma=1)
low_passed = glp.run(deconvolved, in_place=False, verbose=True)

if use_test_data:
    scale_factors = {
        (t[Indices.ROUND], t[Indices.CH]): t['scale_factor']
        for t in experiment.extras['scale_factors']
    }
else:
    scale_factors = {
        (t[Indices.ROUND], t[Indices.CH]): t['scale_factor']
        for index, t in primary_image.tile_metadata.iterrows()
    }

scaled_image = deepcopy(low_passed)

for indices in primary_image._iter_indices():
    data = scaled_image.get_slice(indices)[0]
    scaled = data / scale_factors[indices[Indices.ROUND.value],
                                  indices[Indices.CH.value]]
    scaled_image.set_slice(indices, scaled)

psd = SpotFinder.PixelSpotDetector(
    codebook=experiment.codebook,
    metric='euclidean',
    distance_threshold=0.5176,
    magnitude_threshold=1.77e-5,
    min_area=2,
    max_area=np.inf,
    norm_order=2,
    crop_size=(0, 40, 40)
)

initial_spot_intensities, prop_results = psd.run(scaled_image)

spot_intensities = initial_spot_intensities.loc[
    initial_spot_intensities[Features.PASSES_THRESHOLDS]
]

traces = initial_spot_intensities.stack(traces=('r', 'c'))
tsne = TSNE(n_jobs=8)
tsne_xy = tsne.fit_transform((traces / traces.sum(axis=1)).values)

f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(
    tsne_xy[:, 0][traces.passes_thresholds],
    tsne_xy[:, 1][traces.passes_thresholds],
    s=50, alpha=0.7,
    c=pd.Categorical(traces[traces.passes_thresholds].target).codes,
    cmap=plt.cm.CMRmap
)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('merfish_tsne_passing_codes.png', dpi=150)

f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(
    tsne_xy[:, 0][~traces.passes_thresholds],
    tsne_xy[:, 1][~traces.passes_thresholds],
    s=50, alpha=0.5, c='k'
)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('merfish_tsne_failing_codes.png', dpi=150)