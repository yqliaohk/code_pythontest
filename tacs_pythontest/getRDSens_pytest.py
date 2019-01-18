# python test for implemented getRDSens
import numpy as np 

def getRDSens(Uaerob, Xaero, Xaerob, pt, vars, varsb, Xpts, Xptsb):
    N = np.empty(([1,number_nodes]))
    Na = np.empty_like(N)
    Nb = np.empty_like(N)
    Nc = np.empty_like(N)

    getshapefunction(pt, N, Na, Nb, Nc)

    X = np.empty(([1,3]))
    Xa = np.empty_like(X)
    Xb = np.empty_like(X)
    Xc = np.empty_like(X)

    solidjacobian(X, Xa, N, Na, Nb, Nc, Xpts)

    J = np.empty([1,9])
    jacobian3d(Xa, J)

    Ud = np.empty([1,9])
    getDeformGradient(Ud, J, Na, Nb, Nc, vars)

    U = np.zeros([1,3])

    for i in range(num_nodes):
        U[0] = U[0] + N[i]*u[3*i+0]
        U[1] = U[1] + N[i]*u[3*i+1]
        U[2] = U[2] + N[i]*u[3*i+2]

    R = np.empty([1,3])
    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]
    
    rot = np.empty([1,3])
    rot[0] = 0.5*(Ud[5]-Ud[7])
    rot[1] = 0.5*(Ud[6]-Ud[2])
    rot[2] = 0.5*(Ud[1]-Ud[3])

    Rb = np.zeros([1,3])
    Ub = np.zeros([1,3])
    rotb = np.zeros([1,3])
    
    Ub[2] = Ub[2] + Uaerob[2]
    rotb[0] = rotb[0] + R[1]*Uaerob[2]
    Rb[1] = Rb[1] + rot[0]*Uaerob[2]
    rotb[1] = rotb[1] - R[0]*Uaerob[2]
    Rb[0] = Rb[0] - rot[1]*Uaerob[2]
    
    Ub[1] = Ub[1] + Uaerob[1]
    rotb[2] = rotb[2] + R[0]*Uaerob[1]
    Rb[0] = Rb[0] + rot[2]*Uaerob[1]
    rotb[0] = rotb[0] - R[2]*Uaerob[1]
    Rb[2] = Rb[2] - rot[0]*Uaerob[1]
    
    Ub[0] = Ub[0] + Uaerob[0]
    rotb[1] = rotb[1] + R[2]*Uaerob[0]
    Rb[2] = Rb[2] + rot[1]*Uaerob[0]
    rotb[2] = rotb[2] - R[1]*Uaerob[0]
    Rb[1] = Rb[1] - rot[2]*Uaerob[0]

    Xaerob[0] = Rb[0] 
    Xb[0] = - Rb[0]
    Xaerob[1] = Rb[1]
    Xb[1] = - Rb[1]
    Xaerob[2] = Rb[2]
    Xb[2] = - Rb[2]

    for i in range(num_nodes):
        ub[3*i] = n[i]*Ub[0]
        ub[3*i+1] = n[i]*Ub[1]
        ub[3*i+2] = n[i]*Ub[2]

    rot[0] = 0.5*(Ud[5] - Ud[7])
    rot[1] = 0.5*(Ud[6] - Ud[2]) 
    rot[2] = 0.5*(Ud[1] - Ud[3])

    Udb = np.zeros([1,9])

    Udb[0] = Udb[0];
    Udb[1] = Udb[1]+0.5*rotb[2];
    Udb[2] = Udb[2]-0.5*rotb[1];
    Udb[3] = Udb[3]-0.5*rotb[2];
    Udb[5] = Udb[4]+0.5*rotb[0];
    Udb[6] = Udb[5]+0.5*rotb[1];
    Udb[7] = Udb[6]-0.5*rotb[0];

    Jb = np.empty([1,9])
    getDeformGradientReverse(Udb, Jb, Na, Nb, Nc, varsb)

    jacobian3dReverse(Xab, Xa, Jb)

    solidjacobianReverse(Xb, Xab, N, Na, Nb, Nc, Xptsb)



def getshapefunction(pt, N, Na, Nb, Nc):
    lagrangeSF(na,dna, pt[0])
    lagrangeSF(nb,dnb, pt[1])
    lagrangeSF(nc,dnc, pt[2])
    
    for k in range(order):
        for j in range(order):
            for i in range(order):
                N = na[i]*nb[j]*nc[k];
                Na = dna[i]*nb[j]*nc[k]
                Nb = na[i]*dnb[j]*nc[k]
                Nc = na[i]*nb[j]*dnc[k]
            

def lagrangeSF(sf, dsf, a):
    sf[0] = 1.0
    dsf[0] = 0.0

def solidjacobian(X, Xa, N, Na, Nb, Nc, Xpts):
    Xa = np.zeros([1,9])
    X = np.zeros([1,3])

    for i in range(num_nodes):
        X[0] = X[0] + Xpts[3*i]*N[i]
        X[1] = X[1] + Xpts[3*i+1]*N[i]
        X[2] = X[2] + Xpts[3*i+2]*N[i]

        Xa[0] = Xa[0] + Xpts[3*i]*Na[i]
        Xa[1] = Xa[1] + Xpts[3*i]*Nb[i]
        Xa[2] = Xa[2] + Xpts[3*i]*Nc[i]

        Xa[3] = Xa[3] + Xpts[3*i+1]*Na[i]
        Xa[4] = Xa[4] + Xpts[3*i+1]*Nb[i]
        Xa[5] = Xa[5] + Xpts[3*i+1]*Nc[i]
        
        Xa[6] = Xa[6] + Xpts[3*i+2]*Na[i]
        Xa[7] = Xa[7] + Xpts[3*i+2]*Nb[i]
        Xa[8] = Xa[8] + Xpts[3*i+2]*Nc[i]

def jacobian3d(Xd, Jinv):
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1])- Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2])+ Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
    hinv = 1.0/h
    
    Jinv[0] =   (Xd[4]*Xd[8] - Xd[5]*Xd[7])*hinv
    Jinv[1] = - (Xd[1]*Xd[8] - Xd[2]*Xd[7])*hinv
    Jinv[2] =   (Xd[1]*Xd[5] - Xd[2]*Xd[4])*hinv
        
    Jinv[3] = - (Xd[3]*Xd[8] - Xd[5]*Xd[6])*hinv
    Jinv[4] =   (Xd[0]*Xd[8] - Xd[2]*Xd[6])*hinv
    Jinv[5] = - (Xd[0]*Xd[5] - Xd[2]*Xd[3])*hinv
        
    Jinv[6] =   (Xd[3]*Xd[7] - Xd[4]*Xd[6])*hinv
    Jinv[7] = - (Xd[0]*Xd[7] - Xd[1]*Xd[6])*hinv
    Jinv[8] =   (Xd[0]*Xd[4] - Xd[1]*Xd[3])*hinv

    return h
        

def getDeformGradient(Ud, J, Na, Nb, Nc, vars):
    Ua = np.zeros([1,9])
    
    for i in range(num_nodes):
        Ua[0] = Ua[0] +vars[3*i]*Na[i]
        Ua[1] = Ua[1] +vars[3*i]*Nb[i]
        Ua[2] = Ua[2] +vars[3*i]*Nc[i]

        Ua[3] = Ua[3] +vars[3*i+1]*Na[i]
        Ua[4] = Ua[4] +vars[3*i+1]*Nb[i]
        Ua[5] = Ua[5] +vars[3*i+1]*Nc[i]
        
        Ua[6] = Ua[6] +vars[3*i+2]*Na[i]
        Ua[7] = Ua[7] +vars[3*i+2]*Nb[i]
        Ua[8] = Ua[8] +vars[3*i+2]*Nc[i]

    Ud[0] = Ua[0]*J[0] + Ua[1]*J[3] + Ua[2]*J[6]
    Ud[3] = Ua[3]*J[0] + Ua[4]*J[3] + Ua[5]*J[6]
    Ud[6] = Ua[6]*J[0] + Ua[7]*J[3] + Ua[8]*J[6]

    Ud[1] = Ua[0]*J[1] + Ua[1]*J[4] + Ua[2]*J[7]
    Ud[4] = Ua[3]*J[1] + Ua[4]*J[4] + Ua[5]*J[7]
    Ud[7] = Ua[6]*J[1] + Ua[7]*J[4] + Ua[8]*J[7]

    Ud[2] = Ua[0]*J[2] + Ua[1]*J[5] + Ua[2]*J[8]
    Ud[5] = Ua[3]*J[2] + Ua[4]*J[5] + Ua[5]*J[8]
    Ud[8] = Ua[6]*J[2] + Ua[7]*J[5] + Ua[8]*J[8]


def getDeformGradientReverse()
        