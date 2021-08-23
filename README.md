# FlexAID

FlexAID, implements a docking algorithm that can use small-molecules and peptides as ligands and proteins/nucleic acids as targets. It permits full ligand flexibility as well as target side-chain flexibility. FlexAID utilizes a soft scoring function, i.e. one that is not highly dependent on specific geometric criteria, based on surface complementarity. The energy parameters of the scoring function were derived from the classification of a largedataset of native and near native (less than 2Å RMSD) conformations for nearly 1500 complexes from the PDBbind database as true positive examples. These were countered over successive rounds of Monte Carlo optimization over an ever increasing and successively more difficult sets (increasingly lower energy decoys) with RMSD above 2Å.

## References 

[Gaudreault F, Najmanovich RJ. FlexAID: Revisiting Docking on Non-Native-Complex Structures. J Chem Inf Model. 2015;55: 1323–1336.](https://doi.org/10.1021/acs.jcim.5b00078)

[Morency LP, Gaudreault F, Najmanovich R. Applications of the NRGsuite and the Molecular Docking Software FlexAID in Computational Drug Discovery and Design. Methods in Molecular Biology (Clifton, N.J.). 2018;1762: 367-388.](https://dx.doi.org/10.1007/978-1-4939-7756-7_18)
