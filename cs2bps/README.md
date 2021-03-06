# CS2BPS
Chemical Shift to Base Pairing Status predictions

## Installation

set up a conda environment for running CS2BPS

```bash
conda create -n r-environment r-essentials r-base
conda activate r-environment
```

install MICE (for imputation)

```bash
conda install  -y -c conda-forge r-mice 
```

## Usage

```bash
# load environment
conda activate r-environment

# Print help message
Rscript cs2bps.r --help # this will automatically detect and install necessary libraries

# Impute only
Rscript cs2bps.r test_miR20b/measured_shifts_2N7X.dat -i 2N7X -p impute -o test_miR20b/

# Impute and then run cs2bps
Rscript cs2bps.r test_miR20b/measured_shifts_2N7X.dat -i 2N7X -p cs2pbs -o test_miR20b/
```

## License
MIT License

Copyright (c) 2019 Aaron T. Frank

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
