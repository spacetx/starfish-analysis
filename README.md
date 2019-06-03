# starfish-analysis
This repo contains recipes for processing spaceTX data using starfish on a single FOV


## Format Specification
A recipe should be a python .py file containing a single method 
```python
def process_fov(fov_num: int, experiment_str: str):
    """Process a single field of view of starfish data
    Parameters
    ----------
    fov_num : int
        the field of view to process
    experiment_str : int
        path of experiment json file

    Returns
    -------
    DecodedSpots :
        tabular object containing the locations of detected spots.
    """
```

At the top of the file the recipe should also include the url location of the full experiment data. 

An example `recipe.py` file can be found [here](recipes/iss/recipe.py)


## Special processing 
If an experiment requires processing data on a per round/ch/zplane basis within a fov the single method within
the provided recipe file should include the round/ch/zplane number as an extra parameter to the `process_fov` method. 

An example of a `recipe.py` that processes per round within a fov can be found [here](recipes/smFish/recipe.py)
