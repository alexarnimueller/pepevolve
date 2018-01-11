# pepevolve
Evolutionary algorithm to generate offspring of one or morph two (peptide) sequences into each other.

## pepevolve.py
Generate offspring for a given sequence. This script is designed to work with peptide
sequences. However, any sequence can be provided, given a suitable distance matrix is referenced.

### Usage:
``` bash
python pepevolve.py <parent> --lambd <int> --sigma <float> --matrixfile <file> --skip <str> --seed <int>
```
#### Parameters:
- `parent`: {str} parent string from which you wish to generate children
- `lambd`: {int} number of offsprings to generate for the given parent
- `sigma`: {float} width of the gaussian distribution for tuning the distance to the parent
- `matrixfile`: {str} filename of the distance matrix to use
- `skip`: {str} letters (AA) to skip when sampling
- `seed`: {int} random seed to use when sampling (makes runs reproducible)

#### Example:
``` bash
python pepevolve.py GLFDIVKKVVGALGSL --lambd 10 --sigma 0.1 --matrixfile grantham.txt --skip CM --seed 42
```

#### Output:
generated sequences with corresponding distances and sigma values written to the file `restult.txt`

## pepmorph.py
Evolutionary algorithm to morph two given sequences into each other with the help of an evolutionary algorithm. This
script is designed to work with peptide sequences. However, any sequence can be provided, given a suitable distance
matrix is referenced.

### Usage:
``` bash
python pepmorph.py <start> <target> --lambd <int> --sigma <float> --matrixfile <file> --skip <str> --seed <int>
```

#### Parameters:
- `parent`: {str} parent string from which you wish to generate children
- `target`: {str} target string to which the parent should be morphed to (same length as parent!)
- `lambd`: {int} number of offsprings to generate in between the two sequences
- `sigma`: {float} width of the gaussian distribution for tuning the distance to the previous sequence
- `matrixfile`: {str} filename of the distance matrix to use
- `skip`: {str} letters (AA) to skip when sampling
- `seed`: {int} random seed to use when sampling (makes runs reproducible)

#### Example:
``` bash
python pepmorph.py KLLKLLKKLLKLLK GLFDIVKKVVGALGSL --lambd 10 --sigma 0.1 --matrixfile grantham.txt --skip CM --seed 42
```

#### Output:
generated sequences with corresponding distances and sigma values written to the file `restult.txt`

