# Starfish-Analysis

Contains analyses executed using the starfish package

Because starfish is an evolving library, we the
authors cannot guarantee that all analyses will be
updated as the library evolves. However, we are committed
to reproducible science. Because of this, We version both
the libraries and the datasets used in starfish, which
means that if the correct version of starfish is installed,
the scripts in this library should work at any point in the
future.

To enable this, we assume that the user has cloned the
starfish repository from source and installed it in develop
mode using a virtualenv in `.venv` as instructed in the
starfish documentation

```bash
git clone https://github.com/spacetx/starfish.git
cd starfish
python -m venv .venv
source .venv/bin/activate
pip install -e .
```

Currently, all scripts in this repository depend on a single
starfish version. To run the scripts, install additional
dependencies for these scripts:

```
pip install -r MAYBE_EXTRA_REQUIREMENTS.TXT
```

Then, install the specific version of starfish that was used
to generate these plots.

```
git checkout 7f643c4
pip install -e .
```

Now all the scripts should run!
