# IT IS AN APPROACH TO GET ENERGY BINDING DATA FROM MOLECULAR DOCKING PROCEDURE

import sys
sys.path.append('/home/dnowak/seq_to_seq_and_dock_AMU/pyscreener') # a path to pyscreener

from chembl_structure_pipeline import standardize_mol, get_parent_mol

import pickle
import pandas as pd

import pyscreener
import prody
from openbabel import pybel

from pdbfixer import PDBFixer
from simtk.openmm.app import PDBFile

from Bio.PDB import PDBParser, PDBIO, Select
from Bio.PDB.PDBIO import PDBIO

from rdkit import Chem

from rdkit.Chem import PandasTools


import numpy as np

from pyscreener.docking import vina


debug = False

# returns coordinates (MIN, MAX)
def getMacromoleculeBox(PDBfile):
    print("Processed receptor: "+PDBfile)
    parser = PDBParser(QUIET = True)
    macromolecule = parser.get_structure(PDBfile,PDBfile)
    for model in macromolecule.get_models():
        print("model ", model.id)
        extreemeCoords = {}
        coordMin = [999, 999, 999]
        coordMax = [-999, -999, -999]
        for chain in model.get_chains():
            print("chain ", chain.id, len(list(chain.get_residues())), "residues")
            residuesPerChain = []
            for residue in chain.get_residues():
                for atom in residue.get_atoms():
                #print(residue.get_resname(), atom.get_name(), atom.get_coord())
                    coords = atom.get_coord()

                    for iii, (cooMin, cooMax, coo) in enumerate(zip(coordMin, coordMax, coords)):
                        coordMin[iii] = min(cooMin, float(coo))
                        coordMax[iii] = max(cooMax, float(coo))

        extreemeCoords[model.id] = [coordMin, coordMax]
    return extreemeCoords
    
    
    
# Allows to analyze receptor and find unwanted residues
def analyzeReceptor(receptorFile):
    pdb = prody.parsePDB(receptorFile)
    proteinC = pdb.select('protein')
    
    chains = list(set(pdb.getChids()))
    ligandGeo = []
    for chain in chains:
        ligand = pdb[chain].select('not protein and not water and not ion')
        if ligand: ligandGeo.append((ligand.getCoords(), ligand.getElements()))
        
    residuesInReceptor = None
    if proteinC: residuesInReceptor = list(set(proteinC.getResnames()))
    return ligandGeo, residuesInReceptor
    
    
#removing of unwanted residues

def removeUnwantedResidues(pdbFile, pdbFileOut):
    pdb = PDBParser().get_structure(receptor.split('.')[0], pdbFile)
    residue_to_remove = []
    for model in pdb:
        if debug: print(model)
        for chain in model:
            if debug: print(chain)
            for residue in chain:
                if debug: print(residue.get_resname())
                if residue.id[0] != ' ':
                    residue_to_remove.append((chain.id, residue.id))


    for residue in residue_to_remove:
        model[residue[0]].detach_child(residue[1])

    io = PDBIO()
    io.set_structure(pdb)
    io.save(pdbFileOut)
    return len(residue_to_remove)
    
#returns a search space - a box
def getSearchSpace(receptorFile, margin, mode='fullProteinSearch'):
    pdb = prody.parsePDB(receptorFile)
    ligandC = pdb.select('not protein and not water and not ion')
    proteinC = pdb.select('protein')
    if ligandC: lCoords = ligandC.getCoords()
    pCoords = proteinC.getCoords()
    if mode=='fullProteinSearch':
        size = [abs(pCoords[:, i].min()) + abs(pCoords[:, i].max())+margin for i in range(3)]
        center = [(pCoords[:, i].min() + pCoords[:, i].max())/2 for i in range(3)]
    elif mode=='ligandSearch' and lCoords:
        center = [(lCoords[:, i].min() + lCoords[:, i].max())/2 for i in range(3)]
        size = [20.0, 20.0, 20.0]
    
    return center, size


#perform a molecular docking procedure with specified software and other parameters
def dock_smiles(ligandSmiles, receptorsCleaned, center, size, ncpu, software='vina'): #smima
    pyscreenerInstance = vina.Vina(software=software, receptors=receptorsCleaned, center = center, size=size)
    ligand = pyscreenerInstance.prepare_from_smi(smi=ligandSmiles, rdkitOptimized=True)
    extra = ['--exhaustiveness=8', '--seed=1234']
    dock_results = pyscreenerInstance.dock_ligand(ligand=ligand, software=software, receptors=receptorsCleaned, 
                   center = center, size=size, ncpu=ncpu, extra=extra)
    return dock_results


#geometrical center obtainment
def getLigandGeometricalCenter(molFileName):
    # get molecules
    mols = pybel.readfile('pdbqt', molFileName) #qt
    mols_list = [mol for mol in mols]
    coords = []
    for atom in mols_list[0].atoms:
        coords.append(atom.coords)
    
    coords = np.array(coords)
    geoCenter = [coords[:, idx].mean() for idx in range(3)]
    return geoCenter

#Get structure of selected protein
prody.fetchPDB('7npc', compressed=False)
prody.fetchPDB('7np5', compressed=False)
prody.fetchPDB('7kxd', compressed=False)

consideredReceptors = ['7npc', '7np5','7kxd']

# Potential drugs to be docked loading

to_be_docked = pd.read_excel(r'../prediction_and_selection/All_generated_SMILES_SYBA_filtration.xlsx')
to_be_docked = to_be_docked['SMILES']

ligand_present_in_raw_pdb_file = ['C1=CC(=C(C(=C1)Cl)C2=NOC(=C2COC3=CC=C(C=C3)C(=O)O)C4=CNC=C4)C(F)(F)F', 'CC(C)C1=C(C=CC(=C1)OC2=C(C=C(C=C2Cl)CO)Cl)OC']

to_be_docked = list(to_be_docked)

to_be_docked = to_be_docked.copy()
for i in range(len(ligand_present_in_raw_pdb_file)):
    to_be_docked.append(ligand_present_in_raw_pdb_file[i])
    

##################
##################
#WHOLE PROCEDURE OF MOLECULAR DOCKING 
#WITH CREATION OF RESULTS FILE

rawSufix = '.pdb'
fixSuffix = '_fixer.pdb'
cleSuffix = '_clean.pdb'

nCPU = 12

resultsDF = pd.DataFrame(columns = ['smiles']+consideredReceptors)

sampledSmiles = pd.DataFrame(to_be_docked) #[0:5]
sampledSmiles = sampledSmiles.values.flatten().tolist()

resultsDF['smiles'] = sampledSmiles

for receptor in consideredReceptors:
    print(receptor)
    lCoords, residuesInReceptor = analyzeReceptor(receptor+rawSufix)
    print('Residues (aminiacids) in receptor: ', end=' ')
    for residue in residuesInReceptor:
        print(residue, end=' ')
        
    fixer = PDBFixer(receptor+rawSufix)
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    PDBFile.writeFile(fixer.topology, fixer.positions, open(receptor+fixSuffix, 'w'))
    removeUnwantedResidues(receptor+fixSuffix, receptor+cleSuffix)
    print('\n')
    
    # full search space with margin
    center, size = getSearchSpace(receptor+cleSuffix, margin=10)
    print('size:, center ', size, center)
    

    receptors = [receptor+cleSuffix] #+'qt'
    receptors_docking = [receptor+cleSuffix+'qt']
    results = []
    
    ##Create pdbqt files...
    try:
        dock = dock_smiles(sampledSmiles[0], receptors, center, size, nCPU)
    except:
        pass
    for smiles in sampledSmiles: #[:1]
        
        smilesCleaned = Chem.MolToSmiles(get_parent_mol(Chem.MolFromSmiles(smiles), neutralize=True, check_exclusion=True, verbose=False)[0])
        try:
            dock_results = dock_smiles(smilesCleaned, receptors_docking, center, size, nCPU) #here should be list of all receptors in *pdbqt format
        except: #receptors
            dock_results = None
            
        if dock_results:
            docked_molecule = dock_results[0][0]['out']
            geoCenter = getLigandGeometricalCenter(docked_molecule.name)
            dockingScore = dock_results[0][0]['score']
            result = tuple([dockingScore, geoCenter])
            results.append(result)
        else:
            result = tuple([None, None])
            results.append(result)            
    
    pd.Index.size!=len(sampledSmiles)
    resultsDF[receptor] = results  
    
##################
##################

#Add column with image of each docked molecule
smilesCleaned = [Chem.MolToSmiles(get_parent_mol(Chem.MolFromSmiles(smiles), neutralize=True, check_exclusion=True, verbose=False)[0]) for smiles in sampledSmiles]
resultsDF['Mol Image'] = [Chem.MolFromSmiles(s) for s in smilesCleaned]
PandasTools.SaveXlsxFromFrame(resultsDF, 'dockingResults_ROR_gamma_SYBA_selected.xlsx', molCol='Mol Image')

#Save results to csv file
#resultsDF.to_csv('dockingResults_ROR_gamma.csv', index=True) #results without 2D representation of molecule

