import h5py

with h5py.File('template_bank.hdf', 'r') as bankf:
    temp_mass1 = bankf['mass1'][:]
    temp_mass2 = bankf['mass2'][:]
    temp_s1z = bankf['spin1z'][:]
    temp_s2z = bankf['spin2z'][:]
    
print(len(temp_mass1))
print(temp_mass1)
