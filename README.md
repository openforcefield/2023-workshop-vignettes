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

Note that, if running locally, you can skip the first cells in each notebook that look something like this:

```python
!wget -q https://raw.githubusercontent.com/openforcefield/2023-workshop-vignettes/main/colab_setup.ipynb
%run colab_setup.ipynb
```

## Web-based install-free usage

Each notebooks is also mirrored on a Google Colab instance linked below. This is a browser-based service that bypasses the need for local installations, but requires a few minutes of setup before the code can be run. (No action should be needed outside of running the cells.)

* [🟢 PDB file to OpenMM simulation](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-PDB-to-simulation.ipynb)
* [🟢 PDB file to OpenMM simulation with bespoke parameters](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-PDB-bespokefit-to-simulation.ipynb)
* [🟢 SMILES to parameters](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-SMILES-to-parameters.ipynb)
* [🟢 Ligand in water](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-ligand-in-water.ipynb)
* [🟢 Parametrize a non-canonical amino acid](https://gist.github.com/Yoshanuikabundi/66007cb9966b1455a259baaf7cd7e7c3)
* [🟢 Modify a ligand with RDKit reactions](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-rdkit-ligand-modification.ipynb)
* [🟢 Retrieve a torsion drive from QCArchive](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-retrieve-qcarchive-torsiondrive.ipynb)
* [🟢 Vectorized representations](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/G-vectorized-representations.ipynb)
* [🟡 Protein-ligand complex via Interchange combination](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/Y-interchange-combination-export.ipynb)
* [🟡 Importing a prepared OpenMM system via Interchange](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/Y-from-openmm-xml.ipynb)
* [🟡 Preparing mixed solvents and exporting to GROMACS](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/Y-interchange-gromacs-export.ipynb)
* [🟡 Micelle self-assembly](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/Y-micelle-self-assembly.ipynb)
* [🔴 Custom substructure parametrization and charge assignment with OpenFF NAGL](https://colab.research.google.com/github/openforcefield/2023-workshop-vignettes/blob/main/R-custom_substructures_and_nagl.ipynb)
