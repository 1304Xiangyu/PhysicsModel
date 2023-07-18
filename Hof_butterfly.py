import numpy as np
import matplotlib.pyplot as plt

j=1j
# Let's consider a 2D lattice model, H in real space, Energy as a function of Magnetic flux

# Lattice constant
a,b=1,1

# The \Phi=q/p*\Phi_0 for each unit cell
# Landau Gauge, B=(0,0,B), Ay=(0,Bx,0), Let's take x=1
# Determine the magnetic supercell
q,p=9,1 # \Phi=p/q*\Phi_0
Num_wann=1  # Number of wannier band
ta, tb, E = 0.2, 0.05, 0 # Hopping term and on-site energy
Ndim = q # Dimension of the H matrix

# Kpath
Num_k = 201
k1=np.zeros(Num_k)
k2=np.zeros(Num_k)

for i in range(Num_k):
    k1[i]=-np.pi+2*np.pi*i/(Num_k+1)
    k2[i]=-np.pi+2*np.pi*i/(Num_k*q+1)

# Some constant, C=e/hbar, \Phi_0 magnetic flux
Phi = p/q

H = np.zeros((Ndim,Ndim,Num_k),dtype=complex)

# Diagonal elements
# Diagonal elements
# We only consider the change of k2 in this script
for m in range(Num_k):
    for i in range(q):
        H[i,i,m] = E - 2*ta*np.cos(k1[m]+2*np.pi*Phi*i)

    
# Upper matrix element
# Upper matrix element

for m in range(Num_k):
    for i in range(q-1):
        H[i,i+1,m] = -tb*np.exp(j*k2[21])
        
for m in range(Num_k):
    H[0,q-1,m] = -tb*np.exp(j*k1[m]*q)
    
## Lower Matrix elements

for m in range(Num_k):
    for i in range(q-1):
        H[i+1,i,m] = -tb*np.exp(-j*k2[21])
        
for m in range(Num_k):
    H[q-1,0,m] = -tb*np.exp(-j*k1[m]*q)
    
    
## Check Hermitian
H_check = np.zeros((Ndim,Ndim))
for m in range(Num_k-1):
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
eigenvalue=np.zeros((Ndim,Num_k))
for i in range(Num_k):
    eigenvalue[:,i] = np.linalg.eigvals(H[:,:,i])
    eigenvalue[:,i] = np.sort(eigenvalue[:,i])

plt.figure()

for q in range(q):
    plt.plot(k1,eigenvalue[q,:],'k')

#plt.plot(k2,eigenvalue[2,:])
#print(np.linalg.eigvals(H[:,:,0]))
#print(eigenvalue[:,9])
plt.show()    


#print(eigenvalue[:,2])
#print(H[:,:,2])


