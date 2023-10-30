# seq_to_seq_and_dock_AMU
## Title
```Neural networks in the design of molecules with affinity to selected protein domains.```

# Workflow
Open ```Order_of_usage.md``` to see detalied description of each folder.

## Training data preparation
<p align="center">
  <img src="/scheme/data_preparation.png">
</p>


## Model training
<p align="center">
  <img src="/scheme/scheme_model_v2.png">
</p>

# Acknowledgments

Many thanks to Ph. D. Esben Jannik Bjerrum ``` https://github.com/EBjerrum ``` who is the author of Cheminformania.com blog. The encoder-decoder code (model) was inspired from code at 
```
https://www.cheminformania.com/master-your-molecule-generator-seq2seq-rnn-models-with-smiles-in-keras/
```

Special thanks to Ph. D. Rafał Bachorz ``` https://github.com/rafalbachorz ``` for invaluable help during code preparation. 

## Prediction

<p align="center">
  <img src="/scheme/scheme_prediction_v2.png">
</p>

## SMILES filtration (initials and generated)

<p align="center">
  <img src="/scheme/scheme_filtration.png">
</p>

## Docking

<p align="center">
  <img src="/scheme/scheme_docking.png">
</p>


# Dependencies:

```
https://github.com/chembl/ChEMBL_Structure_Pipeline
```
```
https://github.com/XDamianX-coder/pyscreener_old
```
```
https://github.com/lich-uct/syba
```



# Installation of environment
1. Install conda from [Conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html), it can be Anaconda or Miniconda - for purpose of this work Linux Miniconda has been used.
2. Git is necessary to dowload repository. [Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git)
3. When installed ```git clone https://github.com/XDamianX-coder/seq_to_seq_and_dock_AMU``` can be used to download this repo and use ```cd seq_to_seq_and_dock_AMU```.
4. Create usable environment by using ```conda env create -f environment.yml```
5. After environment is created use ```conda activate cheminf``` to use currect domain.
6. Run ```conda install -c rdkit -c lich syba``` to install SYBA classifier. (SYnthetic BAyesian classifier (SYBA))
7. Installation of [ChEMBL structure pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline) can be done by using ```git clone https://github.com/chembl/ChEMBL_Structure_Pipeline```
8. Then use cd ChEMBL_Structure_Pipeline and ```pip install .``` to make ChEMBL Structure Pipeline ready. 
9. Installation of [Pyscreener](https://github.com/XDamianX-coder/pyscreener_old) can be done by using ```git clone https://github.com/XDamianX-coder/pyscreener_old```
10. Little modification in name is necessary - from pyscreener_old to pyscreener.
11. One also must be sure that all dependencies from [Pyscreener](https://github.com/XDamianX-coder/pyscreener_old) are installed, such like [ADFR Suite](https://ccsb.scripps.edu/adfr/downloads/) and [Vina](http://vina.scripps.edu/) in case of this environment to work properly.
12. To run pythonic molecular docking user should be sure that lines to pyscreener are in the docking script ```Python_molecular_docking.py```. 
These are in 2nd line and should be prepared in manner given below:
```
import sys
sys.path.append('PATH_TO_PYSCREENER')
```
Two additional folders should be present one called ```ChEMBL_Structure_Pipeline``` and second ```pyscreener```.


# Citation
```
@article{nowak_neural_2023,
	title = {Neural Networks in the Design of Molecules with Affinity to Selected Protein Domains},
	volume = {24},
	rights = {All rights reserved},
	issn = {1422-0067},
	url = {https://www.mdpi.com/1422-0067/24/2/1762},
	doi = {10.3390/ijms24021762},
	abstract = {Drug design with machine learning support can speed up new drug discoveries. While current databases of known compounds are smaller in magnitude (approximately 108), the number of small drug-like molecules is estimated to be between 1023 and 1060. The use of molecular docking algorithms can help in new drug development by sieving out the worst drug-receptor complexes. New chemical spaces can be efficiently searched with the application of artificial intelligence. From that, new structures can be proposed. The research proposed aims to create new chemical structures supported by a deep neural network that will possess an affinity to the selected protein domains. Transferring chemical structures into {SELFIES} codes helped us pass chemical information to a neural network. On the basis of vectorized {SELFIES}, new chemical structures can be created. With the use of the created neural network, novel compounds that are chemically sensible can be generated. Newly created chemical structures are sieved by the quantitative estimation of the drug-likeness descriptor, Lipinski’s rule of 5, and the synthetic Bayesian accessibility classifier score. The affinity to selected protein domains was verified with the use of the {AutoDock} tool. As per the results, we obtained the structures that possess an affinity to the selected protein domains, namely {PDB} {IDs} 7NPC, 7NP5, and 7KXD.},
	pages = {1762},
	number = {2},
	journal = {International Journal of Molecular Sciences},
	shortjournal = {{IJMS}},
	author = {Nowak, Damian and Bachorz, Rafał Adam and Hoffmann, Marcin},
	date = {2023-01-16},
        year = {2023},
	langid = {english},
}
```
