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

Note that, if running locally, you can skip or comment out the first cell in each notebook that looks like this:

```python
!wget https://raw.githubusercontent.com/openforcefield/2023-workshop-vignettes/update-install-instructions/colab_setup.ipynb
%run colab_setup.ipynb
```

## Web-based Install-free usage

Each notebooks is also mirrored on a Google Colab instance. This is a browser-based service that bypasses the need for local installations, but requires a bit of setup before the code can be run.
