# FlexAID

Compile FlexAID using CMake:
```
git clone --branch flexaid-cpp https://github.com/NRGlab/FlexAID
cd FlexAID
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --target FlexAID -j 4
```

FlexAID requires a config and ga file to run. These can be generated using `ProcessLigand` (installed with pypi: `pip install processligand-py`).
When using ProcessLigand make sure `atom_index=90000` on the ligand.

# Required Config file codes


| Code     | Description              | Value                                                             | 
|:---------|:-------------------------|:------------------------------------------------------------------|
| `INPLIG` | Ligand input file        | Absolute path to ligand .inp file                                 |
| `METOPT` | Optimization method      | `GA`                                                              |
| `OPTIMZ` | Ligand Flexible residues | One line for each flexible bond in the ligand                     |
| `PDBNAM` | Target input file        | Absolute path to target .inp.pdb file                             |
| `RNGOPT` | Binding site file        | `GLOBAL` or `LOCCLF` + Absolute path to binding site `_sph_` file |

## More details for OPTMIZ:
This line appears at least once for the rigid docking search. Each line contains the ID of the residue to be optimized (AAA – NNNN), followed by an integer number.
This number is the number of the rotatable bond to be optimized or a zero for the ligand to be docked. For example,
`OPTIMZ 132 – 0` defines that residue 132 chain “ “ is the ligand to be docked.

Adding the following lines, you would be setting flexible the first rotatable bond of the ligand and the second flexible bond of the residue whose number is 76, chain A:

`OPTIMZ 132 – 1`

`OPTIMZ 76 A 2`

When using ProcessLigand the residue number is typically `9999` and at least 2 lines of `OPTIMZ` are required:

`OPTIMZ 9999 – -1`

`OPTIMZ 9999 – 0`


Additionally, one line is required for each line with `FLEDIH` in the ProcessLigand output.

---

## Optional Config file codes


| Code     | Description                                                              | Value                                                                                       | 
|:---------|:-------------------------------------------------------------------------|:--------------------------------------------------------------------------------------------|
| `ACSWEI` |                                                                          |                                                                                             |
| `BPKENM` |                                                                          |                                                                                             |
| `CLRMSD` |                                                                          |                                                                                             |
| `CLUSTA` |                                                                          | `FO` or `DP` or `CF` (Typically not set)                                                    |
| `COMPLF` | Complementarity function to use                                          | `SPH` or `VCT`                                                                              |
| `CONSTR` |                                                                          |                                                                                             |
| `DEECLA` |                                                                          |                                                                                             |
| `DEEFLX` |                                                                          |                                                                                             |
| `DEFTYP` |                                                                          |                                                                                             |
| `DEPSPA` |                                                                          |                                                                                             |
| `EXCHET` |                                                                          |                                                                                             |
| `FLEXSC` | Target flexibility                                                       | One line per flexible residue (Residue number, chain, Residue name). Example: ` 196  A HIS` |
| `HTPMOD` | Makes printing and file writing minimal for use in a high throughput way | N/A                                                                                         |
| `IMATRX` | Matrix file to be loaded                                                 | Absolute path to matrix file                                                                |
| `INCHOH` |                                                                          |                                                                                             |
| `INTRAF` |                                                                          |                                                                                             |
| `MAXRES` | Maximum number of results to output                                      | 10                                                                                          |
| `NMAAMP` |                                                                          |                                                                                             |
| `NMAEIG` |                                                                          |                                                                                             |
| `NMAMOD` |                                                                          |                                                                                             |
| `NOINTR` |                                                                          |                                                                                             |
| `NORMAR` |                                                                          |                                                                                             |
| `NRGOUT` | Time FlexAID waits before aborting when `NRGSUI` option is specified     | 60 (seconds)                                                                                |
| `NRGSUI` | Writes a .update file and waits for it to be deleted before continuing   | N/A                                                                                         |
| `OMITBU` |                                                                          |                                                                                             |
| `OUTRNG` |                                                                          |                                                                                             |
| `PERMEA` | Permeability                                                             | 0.9                                                                                         |
| `RMSDST` | Reference for calculating RMSD                                           | Absolute path to ligand _ref.pdb file                                                       |
| `ROTOBS` |                                                                          |                                                                                             |
| `ROTOUT` |                                                                          |                                                                                             |
| `ROTPER` |                                                                          |                                                                                             |
| `SCOLIG` |                                                                          |                                                                                             |
| `SCOOUT` |                                                                          |                                                                                             |
| `SLVTYP` | User specified atom type for solvent                                     | 40                                                                                          |
| `SLVPEN` |                                                                          |
| `SPACER` | Spacer length                                                            | 0.375                                                                                       |
| `STATEP` | Path to folder where Pause and Abort files can be written.               | Absolute path                                                                               |
| `TEMPER` |                                                                          |                                                                                             |
| `TEMPOP` | Temp folder path                                                         | Absolute path to temp folder (typically inside the `STATEP` folder)                         |
| `USEACS` |                                                                          |                                                                                             |
| `USEDEE` |                                                                          |                                                                                             |
| `VARANG` | Delta angle                                                              | 5.0                                                                                         |
| `VARDIS` |                                                                          |                                                                                             |
| `VARDIH` | Delta dihedral                                                           | 5.0                                                                                         |
| `VARFLX` | Delta flexibility                                                        | 10.0                                                                                        |
| `VCTPLA` |                                                                          |                                                                                             |
| `VCTSCO` |                                                                          |                                                                                             |
| `VINDEX` |                                                                          |                                                                                             |

---

## GA Codes

| Code       | Description                                                   | Value                | 
|:-----------|:--------------------------------------------------------------|:---------------------|
| `NUMCHROM` | Number of chromosomes                                         | (int)                |
| `NUMGENER` | Number of generations                                         | (int)                |
| `ADAPTVGA` |                                                               |                      |
| `ADAPTKCO` |                                                               | (list) with 4 floats |
| `CROSRATE` |                                                               |                      |
| `MUTARATE` | Mutation rate                                                 |                      |
| `POPINIMT` | Population initialization method                              | `RANDOM` or `IPFILE` |
| `FITMODEL` | Fitness model                                                 | `PSHARE` or `LINEAR` |
| `SHAREALF` | Sharing Parameter α                                           |                      |
| `SHAREPEK` | Sharing parameter peaks                                       |                      |
| `SHARESCL` | Sharing scaling                                               |                      |
| `REPMODEL` | Reproduction technique code                                   | `STEADY`, `BOOM`     |
| `BOOMFRAC` | Population boom size  (fraction of the number of chromosomes) | 0 to 1 (float)       |
| `PRINTCHR` | Number of best chromosome to print each generation            | (int)                |
| `PRINTINT` | Print generation progress as well as current best cf          | 0 or 1               |
| `OUTGENER` |                                                               |                      |