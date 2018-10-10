#!/usr/bin/env python

import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from showit import image

from starfish import data
from starfish.types import Features, Indices
from starfish.image import Filter
from starfish.image import Registration
from starfish.spots import SpotFinder
import warnings
from skimage.color import rgb2gray

use_test_data = os.getenv("USE_TEST_DATA") is not None
experiment = data.ISS(use_test_data=use_test_data)

fov = experiment.fov()
primary_image = fov.primary_image
dots = fov['dots']
nuclei = fov['nuclei']
images = [primary_image, nuclei, dots]


# filter raw data
masking_radius = 15
filt = Filter.WhiteTophat(masking_radius, is_volume=False)
for img in images:
    filt.run(img, verbose=True, in_place=True)

registration = Registration.FourierShiftRegistration(
    upsampling=1000,
    reference_stack=dots,
    verbose=True)
registered_image = registration.run(primary_image, in_place=False)

# parameters to define the allowable gaussian sizes (parameter space)
min_sigma = 1
max_sigma = 10
num_sigma = 30
threshold = 0.01

p = SpotFinder.GaussianSpotDetector(
    min_sigma=min_sigma,
    max_sigma=max_sigma,
    num_sigma=num_sigma,
    threshold=threshold,
    measurement_type='mean',
)

# detect triggers some numpy warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")

    # blobs = dots; define the spots in the dots image,
    # but then find them again in the stack.
    blobs_image = dots.max_proj(Indices.ROUND, Indices.Z)
    intensities = p.run(registered_image, blobs_image=blobs_image)

decoded = experiment.codebook.decode_per_round_max(intensities)

genes, counts = np.unique(
    decoded.loc[decoded[Features.PASSES_THRESHOLDS]][Features.TARGET],
    return_counts=True)
table = pd.Series(counts, index=genes).sort_values(ascending=False)


GENE1 = 'HER2'
GENE2 = 'VIM'

# blobs figure
rgb = np.zeros(registered_image.tile_shape + (3,))
rgb[:,:,0] = nuclei.max_proj(Indices.ROUND, Indices.CH, Indices.Z)
rgb[:,:,1] = dots.max_proj(Indices.ROUND, Indices.CH, Indices.Z)
do = rgb2gray(rgb)
do = do/(do.max())

image(do,size=10)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', FutureWarning)
    is_gene1 = decoded.where(decoded[Features.AXIS][Features.TARGET] == GENE1, drop=True)
    is_gene2 = decoded.where(decoded[Features.AXIS][Features.TARGET] == GENE2, drop=True)

plt.plot(is_gene1.x, is_gene1.y, 'or')
plt.plot(is_gene2.x, is_gene2.y, 'ob')
plt.title(f'Red: {GENE1}, Blue: {GENE2}')
plt.savefig('iss_spots_blob_detection.png', dpi=150)

# pixels figure
psf = SpotFinder.PixelSpotDetector(
    experiment.codebook, metric='euclidean',
    magnitude_threshold=0.001, distance_threshold=0.517,
    min_area=15, max_area=300
)
intensities_psf, ccdr = psf.run(registered_image)

rgb = np.zeros(registered_image.tile_shape + (3,))
rgb[:,:,0] = nuclei.max_proj(Indices.ROUND, Indices.CH, Indices.Z)
rgb[:,:,1] = dots.max_proj(Indices.ROUND, Indices.CH, Indices.Z)
do = rgb2gray(rgb)
do = do/(do.max())

image(do,size=10)
with warnings.catch_warnings():
    warnings.simplefilter('ignore', FutureWarning)
    is_gene1 = intensities_psf.where(intensities_psf[Features.AXIS][Features.TARGET] == GENE1, drop=True)
    is_gene2 = intensities_psf.where(intensities_psf[Features.AXIS][Features.TARGET] == GENE2, drop=True)

plt.plot(is_gene1.x, is_gene1.y, 'or')
plt.plot(is_gene2.x, is_gene2.y, 'ob')
plt.title(f'Red: {GENE1}, Blue: {GENE2}')
plt.savefig('iss_spots_pixel_detection.png', dpi=150)
