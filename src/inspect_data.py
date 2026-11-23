import numpy as np

p = r'C:\Users\OSINOR TRILLIONS\cosmo\data\sn_data'
with open(p + r'\PantheonPlus\Pantheon+SH0ES.dat') as f:
    for i in range(3):
        print('P+ line', i, repr(f.readline()[:300]))
with open(p + r'\DES-Dovekie\DES-Dovekie_HD.csv') as f:
    for i in range(3):
        print('Dovekie line', i, repr(f.readline()[:300]))
z = np.loadtxt(p + r'\DES-Dovekie\DES-Dovekie_HD.csv', delimiter=',', skiprows=1)
print('dovekie shape', z.shape)
npz = np.load(p + r'\DES-Dovekie\covtot_inv_000.npz')
print('dovekie invcov keys', list(npz.keys()), npz[list(npz.keys())[0]].shape)
