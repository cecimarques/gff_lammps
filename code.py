# -*- coding: utf-8 -*-
"""
Created on Thu Jan 30 15:46:57 2025

@author: u2371853
"""

import numpy as np

opls = []
# INPUT : CHANGE ACCORDINGLY
# GROMACS "OPLS-AA atom types" for all atoms composing your molecule. 
# These should be informed in ascending order of the LAMMPS atom ID in the
# 1-molecule LAMMPS data file got from LigParGen. 
opls.append("opls_135")
for it_1 in range (1,4):
    opls.append("opls_136")
opls.append("opls_157")
opls.append("opls_154")
for it_1 in range(7,18):
    opls.append("opls_140")
opls.append("opls_155")

NA = len(opls)

# INPUT : CHANGE ACCORDINGLY
# Here are further things you need to change. The number of bonds, angles, dihedrals
# and improper found declared in the 1-molecule LAMMPS data file.
number_of_bonds = 17
number_of_angles = 31
number_of_dihedrals = 39
number_of_impropers = 10

# ---------------------------------------------------------------------------------
atom_name = []
charge = []
epsilon = []
sigma = []

# You should change this accordingly depends on how many OPLS-AA atomtypes you have
# saved in your ffnonbonded file.
lines_ffnonbonded = 813

# Now I am going to find which lines of the ffnonbonded.itp file I care about wi-
# thin the scope of this given compound AND save the "atom name", charges, epsilon
# and sigma that I find there for each of the atom IDs.
for it_1 in range(0, len(opls)):
    value = opls[it_1]
    # ----------------------------------------------------------------------------
    ofi = open("ffnonbonded", 'r')
    for it_2 in range(0, lines_ffnonbonded):
        dump = ofi.readline()
        # This variable will count how many non-zero string sequences (i.e. text) in 
        # the dump variable (i.e. line of the ffnonbonded file) I will find.
        count = 0
        # This variable will take value of 0 everytime I read a space bar or \t in the dump
        # file as well as will be 0 for the first time a read a non-space-bar or \t: 
        # in that latter case, the value of flag will be 1 for all non-space-bar or \t that
        # I read in the sequence (see below).
        flag = 0
        # This variable functions to me as an indicator on whether or not I should be saving
        # data belonging to a text in the dump. Naturally it should only be "activated" if I in-
        # deed met a line of the ffnonbonded.itp file that has the atom type corresponding to
        # the one that appears in the line it_1 of the list opls
        save_data = 0
        # This loops over all slots of the dump string
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                # This condition will occur in the first non-space bar occurance of a
                # text string of the line. I will make sure I define it as 1 so that
                # I dont fall here in the sequence of the loop it_3.
                if (flag == 0):
                    # Found a text in the line!
                    count = count + 1
                    # defining an empty "string" variable named "value".
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                # This concatenates the string in the slot it_3 of dump to the variable
                # "value".
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    # If the next slot of the dump is a separator, for sure I will be in the last
                    # slot that is not a separator and that belongs to this text.
                    # I will check if it is a text that matches what is the line it_1 of opls and
                    # possibly set the save_data value accordingly.
                    # The if condition above was necessary because in the last slot of
                    # dump I will have a problem in the if condition below (out of range):
                    # to account for the only possible scenario of having it_3 == (len(dump)-1)
                    # in which I MAY need to do something, I will write an if condition below. The
                    # latter will only fail to do what I want if for any reason the last character
                    # is *not* a space bar or a tab AND is a random character, such as ; for exam-
                    # ple that is glued to the 8th text (i.e. epsilon value) that appears in the 
                    # line: as far as I can tell, I dont see this happening when somewhat rapidly
                    # scrolling down the ffnonbonded file.
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if value == opls[it_1]:
                            save_data = 1
                        if (save_data == 1) & (count == 2):
                            atom_name.append(value)
                        if (save_data == 1) & (count == 5):
                            charge.append(value) 
                        if (save_data == 1) & (count == 7):
                            sigma.append(value)
                        if (save_data == 1) & (count == 8):
                            epsilon.append(value)
                            # lets reset the save_data value since I already saved all
                            # the data I needed to save and break the loop of it_3 as
                            # I no longer will need to run it for the atom correspondnig
                            # to opls[it_1]
                            save_data = 0
                            break
                if it_3 == (len(dump)-1):
                    if (save_data == 1) & (count == 8):
                        epsilon.append(value)
                        save_data = 0
                        break

# -----------------------------------------------------------------------------
# Now lets print the converted parameter values for each atom type that exists in
# the file.
epsilon = np.array(epsilon,float)
epsilon = epsilon.reshape(NA,1)
sigma = np.array(sigma,float)
sigma = sigma.reshape(NA,1)

# Unit and GROMACS->LAMMPS potential form conversion (the former is not needed gi-
# ven the form of the non-bonded potential indicated in the forcefield.itp file that
# is found in the same directory as the ffnonbonded.itp):
epsilon = epsilon/4.184
sigma = sigma*10

# Now I am printing these non-bonded parameters for me to directly copy and paste to
# the LAMMPS data file.
nb_output = np.concatenate((epsilon, sigma), axis = 1)
print("non-bonded parameters")
for it_1 in range(0, len(nb_output)):
    print(int(it_1 + 1),"\t", "{0:.6f}".format(nb_output[it_1,0]),"\t", "{0:.6f}".format(nb_output[it_1,1]))

# These I will save in a separate file manually for latter plugging into the third co-
# lumn of the pos section of the propagated LAMMPS data file.
charge = np.array(charge,float)
charge = charge.reshape(NA,1)

# Sanity check to see if the charges of a (neutral) molecule are summing to zero. 
# See result in the kernel.
charsum = 0
for it_1 in range(0, len(charge)):
    charsum = charsum + charge[it_1]
print("Molecule net charge:", charsum)

ofi = open("chargesff.out", 'w')   
for it_1 in range(len(charge)):
    for it_2 in range(0,1):
            ofi.write(str("{0:.6f}".format(charge[it_1,it_2])))
    ofi.write('\n')

# ------------------------------------------------------------------------------
# ------------------------------------------------------------------------------
# Now lets find the suitable bonded parameters, starting with the bonds
bonds = []
ofi = open("bonds", 'r')
for it_1 in range(0, number_of_bonds):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(4)):
            bonds.append(float(e))
bonds = np.array(bonds,float)
bonds = bonds.reshape(number_of_bonds,4)

# You should change this accordingly depends on how many OPLS-AA bondtypes you have
# saved in your ffbond file.
lines_ffbond = 300

print("\n")
print("Bond Potential Parameters")
for it_1 in range(0, len(bonds)):
    # These receive the line indexes of the information that concerns the bonded a-
    # toms in the atom_name array.
    ID1 = atom_name[int(bonds[it_1,2]-1)]
    ID2 = atom_name[int(bonds[it_1,3]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffbond", 'r')
    for it_2 in range(0, lines_ffbond):
        dump = ofi.readline()
        count = 0
        flag = 0
        # These variables will later take the value of the atom names that appear
        # in a given line of the ffbond file.
        ID1_fileff = 0
        ID2_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                            # No need to continue reading this line if I finally found
                            # thesecond ID declared in the bond in the ffbond file.
                            break
        # This condition will only occur if I happened to have broken the loop in the
        # context of the previous command line (no other possibility), which should o-
        # curr for all lines that I read in the ffbond file though. And the if conditi-
        # ons below that will only enter if I happen to have read the line of the ffbond 
        # file that contains the information on the bond between the two given "atom na-
        # mes"). In that latter case I will make sure I output the information already in 
        # the LAMMPS format.
        if count == 2:
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff):
                # These for now empty strings will later take the value of
                # the parameters.
                bit1 = str()
                bit2 = str()
                # This variable will count the text in the line.
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # I should only get to the command line below if I do not
                    # fall into the if condition above.
                    # --------------------------------------------------
                    # If the condition below is met, I am in the first charac-
                    # ter of the text.
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    # These lines will store all the characters of the text
                    # that appear in the 5th and 4th slot.
                    if (count2 == 5) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 4) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                # Converting the units and accounting for differences in poten-
                # tial form.
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*100*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2*10)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(bonds[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)
            if (ID2 == ID1_fileff) & (ID1 == ID2_fileff):
                # These for now empty strings will later take the value of
                # the parameters.
                bit1 = str()
                bit2 = str()
                # This variable will count the text in the line.
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # I should only get to the command line below if I do not
                    # fall into the if condition above.
                    # --------------------------------------------------
                    # If the condition below is met, I am in the first charac-
                    # ter of the text.
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    # These lines will store all the characters of the text
                    # that appear in the 5th and 4th slot.
                    if (count2 == 5) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 4) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                # Converting the units and accounting for differences in poten-
                # tial form.
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*100*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2*10)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(bonds[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)

# There will be repeated occurances if these occur in the .itp files OR if
# the atom names of the two IDs are the same. You will be able to spot re-
# petition by seeing that more than one set of parameter is printed for a sa-
# me bond type. If they are the same, I suppose it is no problem since it 
# implies that the info through the file for the given "atom names" is cohe-
# rent. If not, you will need to go check the ffbond file to see if there is 
# some comment on the lines that help you choose which set of parameters is
# more suitable for your system.
# Note also that obviously you will need to erase the repetition after copy-
# ing pasting everything to the LAMMPS file. This is great, because it will
# then force you to check also if something is missing, since no wildcard
# check is made.
# ------------------------------------------------------------------------------
angles = []
ofi = open("angles", 'r')
for it_1 in range(0, number_of_angles):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(5)):
            angles.append(float(e))
angles = np.array(angles,float)
angles = angles.reshape(number_of_angles,5)

# You should change this accordingly depends on how many OPLS-AA angletypes you have
# saved in your ffangle file.
lines_ffangle = 929

print("\n")
print("Angle Potential Parameters")
for it_1 in range(0, len(angles)):
    ID1 = atom_name[int(angles[it_1,2]-1)]
    ID2 = atom_name[int(angles[it_1,3]-1)]
    ID3 = atom_name[int(angles[it_1,4]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffangle", 'r')
    for it_2 in range(0, lines_ffangle):
        dump = ofi.readline()
        count = 0
        flag = 0
        ID1_fileff = 0
        ID2_fileff = 0
        ID3_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                        if count == 3:
                            ID3_fileff = value
                            break

        if count == 3:
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff) & (ID3 == ID3_fileff):
                bit1 = str()
                bit2 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 6) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 5) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(angles[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)
            if (ID3 == ID1_fileff) & (ID2 == ID2_fileff) & (ID1 == ID3_fileff):
                bit1 = str()
                bit2 = str()
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 6) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 5) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                # ---------------------------------------
                bit1 = float(bit1)
                bit1 = "{0:.6f}".format(bit1/(4.184*2))
                bit1 = str(bit1)
                bit2 = float(bit2)
                bit2 = "{0:.6f}".format(bit2)
                bit2 = str(bit2)
                # ---------------------------------------
                out = str(int(angles[it_1,1]))+str("\t")+bit1+str("\t")+bit2
                print(out)

# Same comment found in the end of the section for finding the bond potential para-
# meters hold here for the angles.
# ---------------------------------------------------------------------
dihedrals = []
ofi = open("dihedrals", 'r')
for it_1 in range(0, number_of_dihedrals):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(6)):
            dihedrals.append(float(e))
dihedrals = np.array(dihedrals,float)
dihedrals = dihedrals.reshape(number_of_dihedrals,6)

# You should change this accordingly depends on how many OPLS-AA atomtypes you have
# saved in your ffnonbonded file.
lines_ffdihedral = 1047

print("\n")
print("Dihedral Potential Parameters")
for it_1 in range(0, len(dihedrals)):
    ID1 = atom_name[int(dihedrals[it_1,2]-1)]
    ID2 = atom_name[int(dihedrals[it_1,3]-1)]
    ID3 = atom_name[int(dihedrals[it_1,4]-1)]
    ID4 = atom_name[int(dihedrals[it_1,5]-1)]
    # ----------------------------------------------------------------------------
    ofi = open("ffdihedral", 'r')
    for it_2 in range(0, lines_ffdihedral):
        dump = ofi.readline()
        count = 0
        flag = 0
        ID1_fileff = 0
        ID2_fileff = 0
        ID3_fileff = 0
        ID4_fileff = 0
        for it_3 in range(0,len(dump)):
            if (dump[it_3] == ' ') or (dump[it_3] == '\t'):
                flag = 0
            if (dump[it_3] != ' ') & (dump[it_3] != '\t'):
                # ------------------------------------------------------------
                if (flag == 0):
                    count = count + 1
                    value = str()
                    flag = 1
                # ------------------------------------------------------------
                value = value+''.join(dump[it_3])
                
                if it_3 != (len(dump)-1):
                    if (dump[it_3+1] == ' ') or (dump[it_3+1] == '\t'):
                        if count == 1:
                            ID1_fileff = value
                        if count == 2:
                            ID2_fileff = value
                        if count == 3:
                            ID3_fileff = value
                        if count == 4:
                            ID4_fileff = value
                            # No need to continue reading this line if I finally found
                            # thesecond ID declared in the bond in the ffbond file.
                            break
        if count == 4:
            if (ID1 == ID1_fileff) & (ID2 == ID2_fileff) & (ID3 == ID3_fileff) & (ID4 == ID4_fileff):
                bit0 = str()  # Co initially
                bit1 = str()  # C1 initially
                bit2 = str()  # C2 initially
                bit3 = str()  # C3 initially
                bit4 = str()  # C4 initially
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 6) & (flag == 1):
                        bit0 = bit0+dump[it_4]
                    if (count2 == 7) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 9) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                    if (count2 == 10) & (flag == 1):
                        bit4 = bit4+dump[it_4]
                # ---------------------------------------
                # Here I am going to transform the Ci originally informed
                # in the GROMACS .itp file into the Fi values, which are
                # arleady the ones to be used by LAMMPS. Please refer to
                # the GROMACS manual to see this.
                bit4 = float(bit4)
                bit4 = -bit4/4
                bit3 = float(bit3)
                bit3 = -bit3/2
                bit2 = float(bit2)
                bit2 = -(bit2-(4*bit4))
                bit1 = float(bit1)
                bit1 = -(2*bit1-(3*bit3))
                # Once I get here, the variables bit1, bit2, bit3 and bit4
                # will have the due values of K1, K2, K3 and K4, respecti-
                # vely.
                bit1 = "{0:.5f}".format(bit1/4.184)
                bit2 = "{0:.5f}".format(bit2/4.184)
                bit3 = "{0:.5f}".format(bit3/4.184)
                bit4 = "{0:.5f}".format(bit4/4.184)
                bit1 = str(bit1)
                bit2 = str(bit2)
                bit3 = str(bit3)
                bit4 = str(bit4)
                # ---------------------------------------
                out = str(int(dihedrals[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3+str("\t")+bit4
                print(out)
            if (ID4 == ID1_fileff) & (ID3 == ID2_fileff) & (ID2 == ID3_fileff) & (ID1 == ID4_fileff):
                bit1 = str()  # C1 initially
                bit2 = str()  # C2 initially
                bit3 = str()  # C3 initially
                bit4 = str()  # C4 initially
                count2 = 0
                for it_4 in range(0, len(dump)):
                    if (dump[it_4] == ' ') or (dump[it_4] == '\t'):
                        flag = 0
                        continue
                    # --------------------------------------------------
                    if flag == 0:
                        count2 = count2 + 1
                        flag = 1
                    # ---------------------------------------------------
                    if (count2 == 7) & (flag == 1):
                        bit1 = bit1+dump[it_4]
                    if (count2 == 8) & (flag == 1):
                        bit2 = bit2+dump[it_4]
                    if (count2 == 9) & (flag == 1):
                        bit3 = bit3+dump[it_4]
                    if (count2 == 10) & (flag == 1):
                        bit4 = bit4+dump[it_4]
                # ---------------------------------------
                bit4 = float(bit4)
                bit4 = -bit4/4
                bit3 = float(bit3)
                bit3 = -bit3/2
                bit2 = float(bit2)
                bit2 = -(bit2-(4*bit4))
                bit1 = float(bit1)
                bit1 = -(2*bit1-(3*bit3))

                bit1 = "{0:.5f}".format(bit1/4.184)
                bit2 = "{0:.5f}".format(bit2/4.184)
                bit3 = "{0:.5f}".format(bit3/4.184)
                bit4 = "{0:.5f}".format(bit4/4.184)
                bit1 = str(bit1)
                bit2 = str(bit2)
                bit3 = str(bit3)
                bit4 = str(bit4)
                # ---------------------------------------
                out = str(int(dihedrals[it_1,1]))+str("\t")+bit1+str("\t")+bit2+str("\t")+bit3+str("\t")+bit4
                print(out)

# Same comment that I wrote in the end of the section where I print the bond
# potential parameters hold here for the dihedrals.

# This part is only necessary if parameters for a dihedral type is not found: not
# finding it would mean that either it doesnt exist listed on the file (kind of un-
# likely provided that your "OPLs-AA atom type" classification is accurate) OR if it
# is formed by a wildcard, case in which you should check the parameters manually in
# ffdihedral file for the given sequence of "atom names".
print("\n")
print("Dihedral atom names")
for it_1 in range(0, len(dihedrals)):
    ID1 = atom_name[int(dihedrals[it_1,2]-1)]
    ID2 = atom_name[int(dihedrals[it_1,3]-1)]
    ID3 = atom_name[int(dihedrals[it_1,4]-1)]
    ID4 = atom_name[int(dihedrals[it_1,5]-1)]
    print(it_1 +1, ID1, ID2, ID3, ID4)
    
# ----------------------------------------------------------------------
impropers = []
ofi = open("impropers", 'r')
for it_1 in range(0, number_of_impropers):
        dump = ofi.readline()
        #dump = dump[0:32]
        for e,it_2 in zip(dump.split('\t'), range(6)):
            impropers.append(float(e))
impropers = np.array(impropers,float)
impropers = impropers.reshape(number_of_impropers,6)

# Here we simply print the atom names of the atoms forming a given improper in the
# 1-molecule LAMMPS data file that is output by LigParGen.
# This approach is not so bad since the section of dihedrals of OPLS-AA is very small
# (few compounds have impropers within OPLS-AA anyways).
print("\n")
print("Improper atom names")
for it_1 in range(0, len(impropers)):
    ID1 = atom_name[int(impropers[it_1,2]-1)]
    ID2 = atom_name[int(impropers[it_1,3]-1)]
    ID3 = atom_name[int(impropers[it_1,4]-1)]
    ID4 = atom_name[int(impropers[it_1,5]-1)]
    print(it_1 + 1, ID1, ID2, ID3, ID4)

# -----------------------------------------------------------------------------------
# Within the scope of the methodology that I use to find the parameters for each bond,
# angle and dihedral type, I am assuming each bond, angle and dihedral declared in the
# respective sections of the LAMMPS data file as output by LigParGen have a unique type 
# (which in this case would match the ID). Indeed this seems to be the approach and if
# this happens, for sure I wont have any problem in the methodology implied in the code
# to find the parameters. If not, there may be problems. I am going to make therefore a
# sanity check on this here and print a warning in case more than one ID shares a type.
for it_1 in range(0, len(bonds)):
    if bonds[it_1,0] != bonds[it_1,1]:
        print("BOND TYPE DOES NOT MATCH BOND ID")
        
for it_1 in range(0, len(angles)):
    if angles[it_1,0] != angles[it_1,1]:
        print("ANGLE TYPE DOES NOT MATCH ANGLE ID")
        
for it_1 in range(0, len(dihedrals)):
    if dihedrals[it_1,0] != dihedrals[it_1,1]:
        print("DIHEDRAL TYPE DOES NOT MATCH DIHEDRAL ID")