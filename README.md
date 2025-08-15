# hdx-rates-mixtures

Code for calculation of Hydrogen-Deuterium forward and back intrinsic exchange rates in mixtures.
The method implemented, described in ..., is an extension of Englander's method (https://hx2.med.upenn.edu/download.html).

To calculate rates from terminal, run the Python script `rates.py` from the folder containing the input sequence, specifying the required inputs.
The results shown in the "example" folder have been generated running the command

`python ../rates.py --seq example.seq --pH 7 --temp 293 --eta 0.9 --out rates_090D.csv`.

As an alternative, the scripts (or functions thereof) can be imported as modules, see the examples in the "notebooks" folder.

Required arguments:
- `--seq` protein sequence, must be a string or a path-like object.
- `--pH` pH read value
- `--temp` temperature (in K)

Optional arguments:
- `--deut` deuteration level, must be a float in [0,1], with
  - 0 = 100% H2O
  - 1 = 100% D2O
  Default to 1.
- `--ref = 'PDLA' OR '3Ala'` reference for intrinsic exchange rates calculation.
  Default to `'3Ala'`.
- `--time = 'h' OR 'm' OR 's'` inverse units of the computed rates.
  Default to `'s'`.
- `--shift` shift in the index of the first amino acid, must be an integer.
  Default to 0.
- `--out` output file name

The output is a CSV file with three items per row, separated by a comma.
Each line corresponds to a residue of the protein sequence. The index of the dataframe contains the residue numbers.
The header has the following keywords:
- `res` residue symbol (one-letter code)
- `kint` intrinsic exchange rate of the forward reaction, i.e. H to D
- `kback` intrinsic exchange rate of the reverse reaction, i.e. D to H
