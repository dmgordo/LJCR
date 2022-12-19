# The La Jolla Combinatorics Repository

This repository is a port of the datasets on [my
website](https://dmgordon.org).  It contains datasets and existence
results on four combinatorial objects I've studied over the years:

- Covering Designs
- Difference Sets
- Circulant Weighing Matrices
- Signed Difference Sets

There is a large literature on all of these but the last, which were
defined by me in a recent paper.
For mathematical background on some of these objects, see [my papers
on arXiv](https://arxiv.org/search/?query=Gordon%2C+Daniel+M&searchtype=author&abstracts=show&order=-announced_date_first&size=50). 


For each one there is a Jupyter Notebook, Sage code, and one or more
json files containing data about the objects.  Each notebook describes
the corresponding object, and runs through how to access and
manipulate the data.

The same data is available on [my website](https://dmgordon.org) as a
MySQL database.  Both will be available for the near term, but
eventually the website will go away, while this is
a permanent location for the data.

This is also an experiment in making a FAIR (Findable, Accessible,
Interoperable, Reusable) mathematical database (see [The FAIR Guiding
Principles for scientific data management and
stewardship](https://doi.org/10.1038/sdata.2016.18) for details.
After researching different approaches used by researchers in
many different scientific areas, implementing it as a Jupyter Notebook
seemed like the best current option.

The repository files may be downloaded and run locally, viewed
statically in [nbviewer](nbviewer.org), or run in your browser using
[binder](mybinder.org).  The last option may be slow, especially for
covering designs, which has a far larger database (2.7GB) than the
other datasets.

## Covering Designs

| Dataset | Nbviewer | binder |
| --               | ---         | ---  |
| Covering Designs | [covering_designs.ipynb](https://nbviewer.jupyter.org/github/dmgordo/LJCR/coverings/covering_designs.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dmgordo/LJCR/master?filepath=coverings%2Fcovering_designs.ipynb) |
| Difference Sets | [difference_sets.ipynb](https://nbviewer.jupyter.org/github/dmgordo/LJCR/diffsets/difference_sets.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dmgordo/LJCR/master?filepath=diffsets%2Fdifference_sets.ipynb) |
| Circulant Weighing Matricess | [circulant_weighing_matrices.ipynb](https://nbviewer.jupyter.org/github/dmgordo/LJCR/cwm/circulant_weighing_matrices.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dmgordo/LJCR/master?filepath=cwm%2Fcirculant_weighing_matrices.ipynb) |
| Signed Difference Sets | [signed_difference_sets.ipynb](https://nbviewer.jupyter.org/github/dmgordo/LJCR/sds/signed_difference_sets.ipynb) | [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dmgordo/LJCR/master?filepath=sds%2Fsigned_difference_sets.ipynb) |



## Contact 
If you encounter any problems with this repository, please report them
to <dmgordo@gmail.com>.
