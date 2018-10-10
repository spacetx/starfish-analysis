#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from MulticoreTSNE import MulticoreTSNE as TSNE

from starfish import data

from starfish.image import Filter
from starfish.spots import SpotFinder

sns.set_context('talk')
sns.set_style('ticks')


# Note that the data here correspond to DARTFISHv1 2017.
# The group is actively working on improving the protocol.
use_test_data = os.getenv("USE_TEST_DATA") is not None
exp = data.DARTFISH(use_test_data=use_test_data)

stack = exp.fov().primary_image

sc_filt = Filter.ScaleByPercentile(p=100)
z_filt = Filter.ZeroByChannelMagnitude(thresh=.05, normalize=False)

norm_stack = sc_filt.run(stack)
zero_norm_stack = z_filt.run(norm_stack)

magnitude_threshold = 0.5
area_threshold = (5, 30)
distance_threshold = 0.6

psd = SpotFinder.PixelSpotDetector(
    codebook=exp.codebook,
    metric='euclidean',
    distance_threshold=distance_threshold,
    magnitude_threshold=magnitude_threshold,
    min_area=area_threshold[0],
    max_area=area_threshold[1]
)

initial_spot_intensities, results = psd.run(zero_norm_stack)
spot_intensities = initial_spot_intensities[initial_spot_intensities.passes_thresholds]

# note that tSNE is not deterministic and I forgot to set a random seed for
# the slides. :(
traces = spot_intensities.stack(traces=('r', 'c'))
tsne = TSNE(n_jobs=8)
tsne_xy = tsne.fit_transform(traces)

f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(tsne_xy[:, 0], tsne_xy[:, 1], s=50, alpha=0.7, c=pd.Categorical(traces.target).codes, cmap=plt.cm.nipy_spectral)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('dartfish_tsne_passing_codes.png', dpi=150)
