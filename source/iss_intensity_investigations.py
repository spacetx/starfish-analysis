import requests
import starfish
from MulticoreTSNE import MulticoreTSNE as TSNE
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns

# get the file
url = 'http://czi.starfish.data.public.s3.amazonaws.com/browse/formatted/20180924/iss_breast/results/decoded_intensities_fov_001.nc'
r = requests.get(url, allow_redirects=True)
open('intensities.nc', 'wb').write(r.content)

# now load the intensities
intensities = starfish.IntensityTable.load('intensities.nc')

traces = intensities.stack(traces=('r', 'c'))
tsne = TSNE(n_jobs=8)
tsne_xy = tsne.fit_transform(traces / traces.sum(axis=1))

# plot both failing and passing spots
f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(tsne_xy[:, 0][~traces.passes_thresholds], tsne_xy[:, 1][~traces.passes_thresholds], s=75, alpha=0.5, c='k')
plt.scatter(tsne_xy[:, 0][traces.passes_thresholds], tsne_xy[:, 1][traces.passes_thresholds], s=50, alpha=0.7, c=pd.Categorical(traces.target).codes[traces.passes_thresholds], cmap=plt.cm.tab20)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('tsne_with_failing_codes.png', dpi=150)

# plot passing spots
f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(tsne_xy[:, 0][traces.passes_thresholds], tsne_xy[:, 1][traces.passes_thresholds], s=50, alpha=0.7, c=pd.Categorical(traces.target).codes[traces.passes_thresholds], cmap=plt.cm.tab20)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('tsne_without_failing_codes.png', dpi=150)

# plot subsampled codes
subsample = []
for gene, group in traces.groupby('target'):
    inds = np.random.randint(0, group.shape[0], 50)
    subsample.append(group.isel(features=inds))

subsample = xr.concat(subsample, dim='features')

tsne = TSNE(n_jobs=8)
subsample_xy = tsne.fit_transform(subsample / subsample.sum(axis=1))

f, ax = plt.subplots(figsize=(10, 10))
plt.scatter(subsample_xy[:, 0], subsample_xy[:, 1], s=75, alpha=0.7, c=pd.Categorical(subsample.target).codes, cmap=plt.cm.tab20)
ax.axes.set_axis_off()
f.tight_layout()
plt.savefig('tsne_subsampled_no_failing_codes.png', dpi=150)

# plot nearest centroid for each failing spot
valid_average_traces = traces[traces.passes_thresholds].groupby('target').mean(axis=0)

nn = NearestNeighbors(n_neighbors=1)
nn.fit(valid_average_traces)
dists, inds = nn.kneighbors(traces[~traces.passes_thresholds])

closest_target_failing_traces = pd.Series(*np.unique(inds, return_counts=True)[::-1])
closest_target_failing_traces.index = valid_average_traces.target.values[closest_target_failing_traces.index]

f, ax = plt.subplots(figsize=(10, 4))
closest_target_failing_traces.plot.bar()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.despine()
f.tight_layout()
plt.savefig('failed_closest_codes.png', dpi=150)

# plot number of spots per gene
f, ax = plt.subplots(figsize=(10, 4))
target_counts = traces[traces.passes_thresholds].groupby('target').count()
pd.Series(target_counts, index=target_counts.target).plot.bar()
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
sns.despine()
f.tight_layout()
plt.savefig('spots.png', dpi=150)

# plot regression between number of genes detected and number of failing codes
# closest to each spot
f, ax = plt.subplots()
sns.regplot(
    np.log(closest_target_failing_traces),
    np.log(target_counts.loc[closest_target_failing_traces.index])
)
sns.despine()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('$log_{10}$ closest target', fontsize=18)
plt.ylabel('$log_{10}$ gene counts', fontsize=18)
plt.tight_layout()
plt.savefig('iss_regression_plot.png', dpi=150)


# plot highest intensity traces for each code
f, ax = plt.subplots()
pd.Series(
    *np.unique(traces[traces.passes_thresholds].groupby('target').mean(axis=0).argmax(axis=1), return_counts=True)[::-1]
).plot.bar(facecolor='royalblue')
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.xlabel('max-intensity bit', fontsize=18)
plt.ylabel('num genes', fontsize=18)
plt.tight_layout()
sns.despine()
plt.savefig('num_max_channels.png', dpi=150)