# fealden

![gpl3.0](https://img.shields.io/github/license/Paradoxdruid/academia-admin-automation.svg "Licensed under GPL 3.0")  [![Language grade: Python](https://img.shields.io/lgtm/grade/python/g/Paradoxdruid/fealden.svg?logo=lgtm&logoWidth=18)](https://lgtm.com/projects/g/Paradoxdruid/fealden/context:python)  [![CodeFactor](https://www.codefactor.io/repository/github/paradoxdruid/fealden/badge)](https://www.codefactor.io/repository/github/paradoxdruid/fealden)  [![Codacy Badge](https://app.codacy.com/project/badge/Grade/520de6c0a1aa463b8b12f4ddc746b4d3)](https://www.codacy.com/gh/Paradoxdruid/fealden/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=Paradoxdruid/fealden&amp;utm_campaign=Badge_Grade) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/ambv/black)

### Introduction 
**Fealden** is a command line tool to generate optimized structure switching DNA biosensors for fluorescent or electrochemical detection of trace amounts of a large number of biomolecule targets.  An abbreviated [bibliography](#bibliography) of works utilizing such nucleic acid-based biosensors is available below.

-------------------------

### Dependencies

**Fealden** depends upon the excellent [RNAstructure package](https://rna.urmc.rochester.edu/RNAstructure.html) from the Mathews Lab, and is written for [Python 3.6+](https://www.python.org/).  

Before use, put the path to the installed RNAstructure package in `config.ini`.

------------------------

### Usage

To use **Fealden**:
`fealden "TATATAA" 1`

(Where `"TATATAA"` is the input binding/recognition element (such as an aptamer), and `1` indicates whether the binding element is predominantly single-stranded (0) or double-stranded (1) in the binding-active state.)

**Fealden** generates a .csv file of optimized biosensor sequences along with scoring metrics.

----------------------

### Contributors

**Fealden** is developed as academic software by the **[Bonham Lab](http://www.bonhamlab.com)** and Dr. Andrew J. Bonham at the [Metropolitan State University of Denver](http://www.msudenver.edu).  It is licensed under the GPL v3.0.  

Contributors include: Dr. Andrew J. Bonham / [@Paradoxdruid]( https://github.com/Paradoxdruid ) (initial and ongoing development), Jody Stephens / [@23jodys]( https://github.com/23jodys ) (early implementation), Becky Addison (early implementation), Aviva Bulow / [@aviva-bulow]( https://github.com/aviva-bulow ) (code rewrite and development of current approach), and Austin Haider / [@WallFacerGibbs]( https://github.com/wallfacergibbs ) (further development).

---------------------

### Bibliography

* [Transcription Factor Beacons for the Quantitative Detection of DNA Binding Activity](http://dx.doi.org/10.1021/ja204775k)
* [Quantification of Transcription Factor Binding in Cell Extracts Using an Electrochemical, Structure-Switching Biosensor](http://dx.doi.org/10.1021/ja2115663)
* [Electrochemical Aptamer Scaffold Biosensors for Detection of Botulism and Ricin Proteins](http://dx.doi.org/10.1007/978-1-4939-6958-6_2)
* [Structure-switching biosensors: inspired by Nature](https://doi.org/10.1016/j.sbi.2010.05.001)
