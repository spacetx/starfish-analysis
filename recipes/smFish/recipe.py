import starfish
from starfish.types import Axes
from starfish.core.imagestack import indexing_utils

"https://d2nhj9g34unfro.cloudfront.net/browse/20190111/allen_mouse_panel_1/experiment.json"


def process_fov(fov_num: int, round_num: int,  experiment_str: str):
    """Process a single field of view of ISS data
    Parameters
    ----------
    fov_num : int
        The field of view to process
    experiment_str : int
        path of experiment json file
    round_num : int
        The round number within the given fov to process

    Returns
    -------
    DecodedSpots :
        tabular object containing the locations of detected spots.
    """

    fov_str: str = f"fov_{int(fov_num):03d}"

    # load experiment
    experiment = starfish.Experiment.from_json(experiment_str)

    fov = experiment[fov_str]
    img = fov.get_images(starfish.FieldOfView.PRIMARY_IMAGES, rounds=[round_num])

    codebook = indexing_utils.index_keep_dimensions(experiment.codebook, {Axes.ROUND: round_num})

    # filter
    clip1 = starfish.image.Filter.Clip(p_min=50, p_max=100)
    clip1.run(img)

    bandpass = starfish.image.Filter.Bandpass(lshort=.5, llong=7, threshold=0.0)
    bandpass.run(img)

    glp = starfish.image.Filter.GaussianLowPass(sigma=(1, 0, 0), is_volume=True)
    glp.run(img)

    clip2 = starfish.image.Filter.Clip(p_min=99, p_max=100, is_volume=True)
    clip2.run(img)

    # find spots
    tlmpf = starfish.spots.DetectSpots.TrackpyLocalMaxPeakFinder(
        spot_diameter=5,  # must be odd integer
        min_mass=0.02,
        max_size=2,  # this is max radius
        separation=7,
        noise_size=0.65,  # this is not used because preprocess is False
        preprocess=False,
        percentile=10,  # this is irrelevant when min_mass, spot_diameter, and max_size are set properly
        verbose=True,
        is_volume=True,
    )

    intensities = tlmpf.run(img)

    # decode
    decoded = codebook.decode_per_round_max(intensities)

    # save results
    df = decoded.to_decoded_spots()
    return df