from Bio.PDB.ResidueDepth import get_surface, residue_depth 
import math

def energies(hblist, aaLib, vdwParams,st, msms_path):
	print ()
	print ('This function is a little slow, please, have patience.')
	print ()
	print ("Residue interactions")
	print ('Res1','ID1','Res2','ID2','Eint','Vdw', 'Dielectric')
	print()

	surface = get_surface(st, MSMS=msms_path)
	
	respairs = []   #We create the list for pairs of residues
	for hb in hblist:
		r1 = hb[0].get_parent()  
		r2 = hb[1].get_parent()
		if [r1,r2] not in respairs:
			respairs.append([r1,r2])   #And we append to the list the pairs of residues
    
	for rpair in sorted(respairs, key=lambda i: i[0].id[1]): 

		eint=0.
		evdw=0.
		diel = 1.0

		for at1 in rpair[0].get_atoms():
			resid1 = rpair[0].get_resname()  #For each atom we take the name of residue and the id of the atom
			atid1 = at1.id
            
			atparam1 = aaLib.getParams(resid1,atid1)  #We take from library the parameters for each type of atom
			vdwprm1 = vdwParams.atTypes[atparam1.atType] 
			for at2 in rpair[1].get_atoms():

                ### AND HERE WE ADD THE CONDITION FOR THE CHANGE OF DIELECTRIC CONSTANT
				if residue_depth(rpair[0], surface) <= 5:
					if residue_depth(rpair[1],surface) <= 5:
						diel = 80.0

				resid2 = rpair[1].get_resname()
				atid2 = at2.id
				atparam2 = aaLib.getParams(resid2,atid2)  
				vdwprm2 = vdwParams.atTypes[atparam2.atType]
				eint = eint + 332.16 * atparam1.charg * atparam2.charg/diel/(at1-at2)
				eps = math.sqrt(vdwprm1.eps*vdwprm2.eps)
				sig = math.sqrt(vdwprm1.sig*vdwprm2.sig)
				evdw = evdw + 4 * eps *( (sig/(at1-at2))**12-(sig/(at1-at2))**6)

		print (resid1,rpair[0].id[1],resid2,rpair[1].id[1],eint,evdw, diel) 
         