# 2023 OpenFF Workshop Vignettes

A collection of notebooks highlighting functionality of OpenFF tools.

## Installation

### Local install

Each notebooks can be run locally in an environment created from the included YAML file:

```shell
$ mamba env create --file environment.yaml
Preparing transaction: ...working... done
Verifying transaction: ...working... done
Executing transaction: ...working... done
$ mamba activate vignettes-env
```

Note that, if running locally, you can skip the first cell in each notebook that looks like this:

```python
!wget -q https://raw.githubusercontent.com/openforcefield/2023-workshop-vignettes/main/colab_setup.ipynb
%run colab_setup.ipynb
%env LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib/
```

## Web-based install-free usage

Each notebooks is also mirrored on a Google Colab instance linked below. This is a browser-based service that bypasses the need for local installations, but requires a bit of setup before the code can be run.

* [PDB file to OpenMM simulation](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-PDB-to-simulation.ipynb)
* [SMILES to parameters](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-SMILES-to-parameters.ipynb)
* [Vectorized representations](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-vectorized-representations.ipynb)
* [Ligand in water](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-ligand-in-water.ipynb)
* [Protein-ligand complex](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-protein-ligand.ipynb)
* [Retrieve a torsion drive from QCArchive](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-vectorized-representations.ipynb)
* [Parametrize a non-canonical amino acid](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/R-custom_substructures_and_nagl.ipynb)
* [Modify a ligand with RDKit reactions](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-rdkit-ligand-modification.ipynb)
