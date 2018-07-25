# fealden

### Introduction 
**Fealden** is a command line tool to generate optimized structure switching DNA biosensors for fluorescent or electrochemical detection of trace amounts of a large number of biomolecule targets.  An abbreviated [bibliography](#Bibliography) of works utilizing such nucleic acid-based biosensors is available below.

-------------------------

### Dependencies

**Fealden** depends upon the excellent [RNAstructure package](https://rna.urmc.rochester.edu/RNAstructure.html) from the Mathews Lab, and is written for [Python 3](https://www.python.org/).  

------------------------

### Usage

To use **Fealden**:
`python ./fealden.py "TATATAA" 1`

(Where `"TATATAA"` is the input binding/recognition element (such as an aptamer), and `1` indicates whether the binding element is predominantly single-stranded (0) or double-stranded (1) in the binding-active state.)

**Fealden** generates a .csv file of optimized biosensor sequences along with scoring metrics.

----------------------

### Authors

**Fealden** is developed as academic software by the **[Bonham Lab](http://www.bonhamlab.com)** and Dr. Andrew J. Bonham at the [Metropolitan State University of Denver](http://www.msudenver.edu).  It is licensed under the LGPL v3.0.  

Contributors include: Dr. Andrew J. Bonham (initial and ongoing development), Jody Stephens (early implementation), Becky Addison (early implementation), Aviva Bulow (code rewrite and development of current approach), and Austin Haider (further development).

---------------------

### Bibliography

* [Transcription Factor Beacons for the Quantitative Detection of DNA Binding Activity](http://dx.doi.org/10.1021/ja204775k)
* [Quantification of Transcription Factor Binding in Cell Extracts Using an Electrochemical, Structure-Switching Biosensor](http://dx.doi.org/10.1021/ja2115663)
* [Electrochemical Aptamer Scaffold Biosensors for Detection of Botulism and Ricin Proteins](http://dx.doi.org/10.1007/978-1-4939-6958-6_2)
* [Structure-switching biosensors: inspired by Nature](https://doi.org/10.1016/j.sbi.2010.05.001)
