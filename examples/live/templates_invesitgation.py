import h5py

with h5py.File('template_bank.hdf', 'r') as bankf:
    temp_mass1 = bankf['mass1'][:]
    temp_mass2 = bankf['mass2'][:]
    temp_s1z = bankf['spin1z'][:]
    temp_s2z = bankf['spin2z'][:]
    
l = len(temp_mass1)
for i in range(l):
    dat = [temp_mass1[i],temp_mass2[i],temp_s1z[i],temp_s2z[i]]
    print(dat)
