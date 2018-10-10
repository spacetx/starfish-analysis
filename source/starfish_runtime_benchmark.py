"""
to generate these hard-coded runtime results, the python scripts that
correspond to the notebooks in the starfish repository, found in
`starfish/notebooks/py` for ISS, MERFISH, DARTFISH, and 3D smFISH (allen
smFISH) were run using:

time python3 starfish/notebooks/py/<notebook_name>

And the following plot was generated based on the outputs.
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

processing_time = pd.Series(
        [54.5, 60 * 2 + 59.3, 46.2, 11 * 60 + 51],
        index=['ISS', 'MERFISH', 'DARTFISH', '3D smFISH']
)

(processing_time / 60).plot.bar(facecolor='royalblue', edgecolor='k')
sns.despine()
plt.ylabel('time (m)', fontsize=16)
plt.title('Single FOV processing time', fontsize=18)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.tight_layout()
plt.savefig('processing_times.png', dpi=150)
