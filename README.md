# fealden

![gpl3.0](https://img.shields.io/github/license/Paradoxdruid/academia-admin-automation.svg?color=success "Licensed under GPL 3.0")  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/fealden/badge)](https://www.codefactor.io/repository/github/paradoxdruid/fealden)  [![Codacy Badge](https://app.codacy.com/project/badge/Grade/520de6c0a1aa463b8b12f4ddc746b4d3)](https://www.codacy.com/gh/Paradoxdruid/fealden/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Paradoxdruid/fealden&amp;utm_campaign=Badge_Grade) [![CodeQL](https://github.com/Paradoxdruid/fealden/actions/workflows/codeql.yml/badge.svg)](https://github.com/Paradoxdruid/fealden/actions/workflows/codeql.yml) [![CI](https://github.com/Paradoxdruid/fealden/actions/workflows/CI.yml/badge.svg)](https://github.com/Paradoxdruid/fealden/actions/workflows/CI.yml) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black) [![mypy](https://img.shields.io/badge/%20type_checker-mypy-%231674b1?style=flat)](https://github.com/python/mypy)

## Introduction

**Fealden** is a command line tool to generate optimized structure switching DNA biosensors for fluorescent or electrochemical detection of trace amounts of a large number of biomolecule targets.  An abbreviated [bibliography](#bibliography) of works utilizing such nucleic acid-based biosensors is available below.

-------------------------

## Dependencies

**Fealden** is written for [Python 3.9+](https://www.python.org/), and depends upon external secondary structure prediction routines.

It can use either the excellent [RNAstructure package](https://rna.urmc.rochester.edu/RNAstructure.html) from the Mathews Lab, or [UNAfold v3.8 and mfold v3.6](http://www.unafold.org/) by Markham and Zuker.

<details>
  <summary>Dependency Installation Tips</summary>

### RNAstructure

* available from the Mathews lab at <https://rna.urmc.rochester.edu/RNAstructure.html>
* fealden requires the "text" interface version, for your operating system
* has a python_interface available, but requires C++ compilation step
* RNAstructure version 6.4 needs corrections in `rna_sources.h` in the `python_interface` folder:
  * all references to `TurboFold` directory need to be replaced with `src` directory
  * can use command `sed -i 's/TurboFold\/src\//g' rna_sources.h` to correct
* requires installation of [swig](https://github.com/swig/swig)
  * such as `pip install swig` or `conda install swig`
* once `swig` installed and `rna_sources.h` corrected, enter the `python_interface` directory and run:
  * `make swig`
  * `make interface-from-distutils`
  * update `.env` file (see below) with path to RNAstructure

### UNAfold and mfold

* version 3.8 of UNAfold is available from sourceforge: <https://rnaspace.sourceforge.net/software/unafold-3.8.tar.gz>
* once unzipped, enter the `unafold-3.8` directory and run:
  * `./configure --prefix=/A/GOOD/PATH/FOR/USER` (for instance, `/home/user/unafold-final`)
  * `make`
  * `make install`
  * update `.env` file (see below) with path to `HYBRID_SS_MIN` in Unafold
* You will also need the program `sir_graph` included in version 3.6 of mfold, available at <http://www.unafold.org/download/mfold-3.6.tar.gz>
* once unzipped, enter the `mfold-3.6` directory and run:
  * `./configure --prefix=/A/GOOD/PATH/FOR/USER` (for instance, `/home/user/mfold-final`)
  * `make`
  * `make install`
  * update `.env` file (see below) with path to `SIR_GRAPH` in Unafold

</details>

Before use, you will need to create a `.env` file following the format in [structure.py](fealden/structure.py).

<details>
  <summary>Example .env file</summary>

```env
FEALDEN_BACKEND=mfold   # either 'mfold' or 'rnastructure'
HYBRID_SS_MIN=/home/username/unafold-new/bin/hybrid-ss-min
SIR_GRAPH=/home/username/mfold/bin/sir_graph
RNASTRUCTURE=/home/username/RNAstructure
```

</details>

-------------------------

## Usage

To use **Fealden**:
`python -m fealden "TATATAA" 1`

(Where `"TATATAA"` is the input binding/recognition element (such as an aptamer), and `1` indicates whether the binding element is predominantly double-stranded (0) or single-stranded (1) in the binding-active state.)

**Fealden** generates a .csv file of optimized biosensor sequences along with scoring metrics.

-------------------------

## Contributors

**Fealden** is developed as academic software by the **[Bonham Lab](http://www.bonhamlab.com)** and Dr. Andrew J. Bonham at the [Metropolitan State University of Denver](http://www.msudenver.edu).  It is licensed under the GPL v3.0.  

Contributors include: Dr. Andrew J. Bonham / [@Paradoxdruid]( https://github.com/Paradoxdruid ) (initial and ongoing development), Jody Stephens / [@23jodys]( https://github.com/23jodys ) (early implementation), Becky Addison (early implementation), Aviva Bulow / [@aviva-bulow]( https://github.com/aviva-bulow ) (code rewrite and development of current approach), and Austin Haider / [@WallFacerGibbs]( https://github.com/wallfacergibbs ) (further development).

-------------------------

## Bibliography

* [Transcription Factor Beacons for the Quantitative Detection of DNA Binding Activity](http://dx.doi.org/10.1021/ja204775k)
* [Quantification of Transcription Factor Binding in Cell Extracts Using an Electrochemical, Structure-Switching Biosensor](http://dx.doi.org/10.1021/ja2115663)
* [Electrochemical Aptamer Scaffold Biosensors for Detection of Botulism and Ricin Proteins](http://dx.doi.org/10.1007/978-1-4939-6958-6_2)
* [Structure-switching biosensors: inspired by Nature](https://doi.org/10.1016/j.sbi.2010.05.001)
