# Try to calculate the energy as a function of the external magnetic flux

import numpy as np
import matplotlib.pyplot as plt

j=1j
# Let's consider a 2D lattice model, H in real space, Energy as a function of Magnetic flux

# Lattice constant
a,b=1,1

# The \Phi=q/p*\Phi_0 for each unit cell
# Landau Gauge, B=(0,0,B), Ay=(0,Bx,0), Let's take x=1
# Determine the magnetic supercell
q,p=101,1 # \Phi=p/q*\Phi_0
Num_wann=1  # Number of wannier band
ta, tb, E = 0.2, 0.2, 0 # Hopping term and on-site energy
Ndim = q # Dimension of the H matrix
# Some constant, C=e/hbar, \Phi_0 magnetic flux
Num_Phi = q     # I think this case would allow one period, if it is 2q, then it is two period

Phi=np.zeros(Num_Phi)
for i in range(Num_Phi):
    Phi[i]=i/q

    

H = np.zeros([Ndim,Ndim,Num_Phi],dtype=complex)
k1,k2 = 0,np.pi/2
# Let's define a function so that I can solve the eigenvalue with different flux
    # Diagonal elements
    # Diagonal elements
    # We only consider the change of k2 in this script
for m in range(Num_Phi): # Different Phi 
    for i in range(q):   # Matrix element
        H[i,i,m] = E - 2*ta*np.cos(k1+2*np.pi*Phi[m]*i)

        
    # Upper matrix element
    # Upper matrix element

for m in range(Num_Phi):
    for i in range(q-1):
        H[i,i+1,m] = -tb*np.exp(j*k2)
            
for m in range(Num_Phi):
    H[0,q-1,m] = -tb*np.exp(j*k1*q)
        
    ## Lower Matrix elements

for m in range(Num_Phi):
    for i in range(q-1):
        H[i+1,i,m] = -tb*np.exp(-j*k2)
            
for m in range(Num_Phi):
    H[q-1,0,m] = -tb*np.exp(-j*k1*q)
        
        
    ## Check Hermitian
H_check = np.zeros((Ndim,Ndim))
for m in range(Num_Phi-1):
    H_check=H[:,:,m]-np.transpose(np.conj(H[:,:,m]))
    #print(H_check[:,:])
    #print(H_boor)
    for n in range(Ndim):
        for l in range(Ndim):
            if H_check[n,l]!=0:
                print("Not Hermitian")
                break

print("Hermitian")
        
        #    break
                
    ## Eigenvalue
eigenvalue=np.zeros((Ndim,Num_Phi))
for i in range(Num_Phi):
    eigenvalue[:,i] = np.linalg.eigvals(H[:,:,i])
    eigenvalue[:,i] = np.sort(eigenvalue[:,i])

#print(eigenvalue[:,0])

plt.figure()

for q in range(Ndim):
    plt.plot(Phi,eigenvalue[q,:],'k')

#plt.plot(k2,eigenvalue[2,:])
#print(np.linalg.eigvals(H[:,:,0]))
#print(eigenvalue[:,9])
plt.show()    
