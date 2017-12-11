# BioPhysics
Scripts and data for the homework in Biophysics course by Rubén Molina.

Usage:

    python3 main.py [-h] [--backonly] [--nowats] pdb_path msms_path

(msms_path is not needed if you use the file attached to here, since it's by default the one uploaded here)


Modules of the script:

## main.py
Main code for the detection and quantification of polar contacts

## plot.py
Module for making a plot with the type of interactions of each residue

## simpleinteractions.py
Module for checking a list of polar residues that interact

## energyinteractions.py
Module for creating a table with interactions and energies.

## ForceField.py
Module for vdw parameters management

## ResLib.py
Module for Residue Library management

## data/vdwprm
Simple vdw parameters (based on AMBER ff)

## msms.MacOSX.2.6.1
Program to compute the Accesible Surface Area. It is needed for the program in order to work.

## data/aaLib.lib
Library for obtaining amino acid atom types and partial charges (based on AMBER ff)

## WARNING THIS LIBRARY HAVE BEEN MODIFIED, SINCE IT WAS ERRONEUS AND DIFFERENT KIND OF ATOMS DID NOT APPEAR, FOR EXAMPLE HISTIDINE, WHICH HAS THE NAME HIE, HIP AND HID BUT NOT NAME HIS, WHICH APPEARS IN THE PDF FILE.
ALSO, HIDROGEN DID NOT APPEAR, SO THE HOMEWORK HAD TO BE DONE WITH AN STRUCTURE WITHOUT HIDROGENS 



External dependencies

    Bio.PDB.NeighborSearch (BioPython)
    Bio.PDB.PDBParser (Biopython)
    Bio.PDB.ResidueDepth (Biopython)



