# CS2BPS
Chemical Shifts guided Base-Pairing Status prediction

## Installation

set up a conda environment for running CS2BPS

```bash
conda create -n r-tf pip python=3.6
source activate r-tf
```

install python tensorflow (CPU version) and keras

```bash
conda install -c anaconda tensorflow
conda install -c conda-forge keras
```

install r-tensorflow using conda

```
conda install -c r r-tensorflow --name r-tf
```

## Usage

```bash
Rscript cs2bps.r [options] path_to_chemical_shift_file
Options:
	-i ID, --id=ID
		ID tag of the test RNA
	-s SPEED, --speed=SPEED
		whether to perform fast or slow imputation
	-p PROGRAM, --program=PROGRAM
		impute or predict
	-o OUTPUT, --output=OUTPUT
		name of imputed chemical shifts file or names of predicted base pairing probability file
	-w WHOLE_SET_PREDICTION, --whole_set_prediction=WHOLE_SET_PREDICTION
		whether to output predictions from all classifiers
	-h, --help
		Show this help message and exit
```

## License
