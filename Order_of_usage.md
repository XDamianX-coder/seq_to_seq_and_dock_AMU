# Usage order

## Model

First of all data that will feed the model should be collected. It can be done by running ```Data_preparation.ipynb``` which is present in ```model/data``` directory.
To create model use ```model``` folder and run ```Molecule_generator-Generative_neural_network.py```
All necessary files will be created.
To plot model's architecture use ```Model_plot.ipynb``` which is presnt inside ```model``` folder.

## Prediction and data selection

To make predictions and select SMILES to be docked enter ```prediction_and_selection``` and go to 
```ROR_gamma_active_QED_Lipinski.ipynb```, run it, next step is done inside ```SYBA_class_ROR_gamma_activ_AFTER_QED_LIPINSKI.ipynb```.
When initial molecules are selected new molecules generation can be launched in  ```Prediction_1_0.1_tensor_scaling.ipynb``` or / and ```Prediction_2_0.2_tensor_scaling.ipynb``` then enter  ```Selection_from_0_1_tensor_scaling.ipynb``` or ```Selection_from_0_2_tensor_scaling.ipynb``` respectively.
Then combine obtained excel files with ```Combine_Generated_and_Selected_structures.ipynb```.
To assing generation mode to each selected structure enter and run ```Assigning_prediction_mode_to_selected_SMILES.ipynb``` as a result excel sheet will containt 0-1 columns indicating which mode has been used to get the structure.
Use additional classifier SYBA with ```SYBA_classifier_additional_filter.ipynb```.
If one is interested in checking if training data somehow converge with generated ```Checking_if_ROR_gamma_activ_are_in_ZINC_db.ipynb``` can be entered and launched.
To check how SYBA classfier assign scores to Kinase Inhibitors enter ```SYBA_classifier_additional_filter-ROR_gamma_activ.ipynb```.
To visualize structures that will be docked use ```All_generated_SMILES_visualization.ipynb```.

## Tanimoto similarity

To check similarity between generated structures, after selection, and training structures along Kiniase inhibitors structures enter ```Tanimoto_similarity-SYBA_selection.ipynb``` it will result in comparison of molecules that fulfilled SYBA threshold (above 0 score) other option is to compare all structures which fulfilled QED and Lipinski's Rule of 5 in ```Tanimoto_similarity_All_generated_and_selected.ipynb```.
It can be rerun after molecular docking procedure if some structures will not be able to be docked, to do so enter molecular docking procedure add ```_blind_try``` to the end of the filename, just before ```.xlsx``` and move back to ```All_generated_SMILES_visualization.ipynb``` and set ```first_try_done``` value to ```True```. Then Tanimoto similarity can be done one more time, but without structures that from different reasons cannot be docked. 

## Molecular docking

To dock selected structures (which fulfilled QED, Lipinski's Rule of 5 and SYBA threshold) use ```molecular_docking``` folder.
First of all use ```Python_molecular_docking.py``` excel sheet will be created.
To clean results use ```Clean_results.ipynb```.
To select structures with lowest energy binding energies use ```Selection_of_most_prominent_structures.ipynb```.
In case of manual docking with AutoDock 1.5.6 one can use ```SMILES_to_3D_PDB.ipynb``` and this will create useful 3D PDB structure.
Generation mode 0.1 or 0.2 tensor rotation molecules can be distinguish via ```Assigning_prediction_mode_to_obtained_and_selected_molecule.ipynb``` notebook.

## Overall workflow - click to enlarge the picture

<p align="center">
  <img src="/scheme/Impartant_scheme.png">
</p>
