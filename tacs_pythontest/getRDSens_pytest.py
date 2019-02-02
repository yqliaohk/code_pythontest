# python test for implemented getRDSens
import numpy as np 

def getRDSens(Uaerob, Xaero, pt, vars,Xpts):
    global num_nodes
    

    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)

    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros((3))
    Xa = np.zeros((9))
    Xb = np.zeros_like(X)
    Xab = np.zeros_like(Xa)

    Xa = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    J = np.zeros(9)
    h, J = jacobian3d(Xa)

    Ud = np.zeros(9)
    Ud = getDeformGradient(Ud, J, Na, Nb, Nc, vars)

    U = np.zeros(3)

    u = vars
    for i in range(num_nodes):
        U[0] = U[0] + N[i]*u[3*i+0]
        U[1] = U[1] + N[i]*u[3*i+1]
        U[2] = U[2] + N[i]*u[3*i+2]

    R = np.zeros(3)
    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]
    
    rot = np.zeros(3)
    rot[0] = 0.5*(Ud[5]-Ud[7])
    rot[1] = 0.5*(Ud[6]-Ud[2])
    rot[2] = 0.5*(Ud[1]-Ud[3])

    Rb = np.zeros(3)
    Ub = np.zeros(3)
    rotb = np.zeros(3)
    
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

    Xaerob = np.zeros(3)
    Xaerob[0] = Rb[0] 
    Xb[0] = - Rb[0]
    Xaerob[1] = Rb[1]
    Xb[1] = - Rb[1]
    Xaerob[2] = Rb[2]
    Xb[2] = - Rb[2]

    ub = np.zeros_like(u)

    n = N
    for i in range(num_nodes):
        ub[3*i] = n[i]*Ub[0]
        ub[3*i+1] = n[i]*Ub[1]
        ub[3*i+2] = n[i]*Ub[2]

    varsb = ub

    rot[0] = 0.5*(Ud[5] - Ud[7])
    rot[1] = 0.5*(Ud[6] - Ud[2]) 
    rot[2] = 0.5*(Ud[1] - Ud[3])

    Udb = np.zeros(9)

    Udb[0] = Udb[0]
    Udb[1] = Udb[1]+0.5*rotb[2]
    Udb[2] = Udb[2]-0.5*rotb[1]
    Udb[3] = Udb[3]-0.5*rotb[2]
    Udb[5] = Udb[4]+0.5*rotb[0]
    Udb[6] = Udb[5]+0.5*rotb[1]
    Udb[7] = Udb[6]-0.5*rotb[0]

    Jb = np.zeros(9)
    Jb, varsb = getDeformGradientReverse(Udb, J, Jb, Na, Nb, Nc, vars, varsb)

    Xab = jacobian3dReverse(Xa, Jb)

    Xptsb = np.zeros_like(Xpts)
    Xptsb = solidjacobianReverse(Xptsb, Xb, Xab, N, Na, Nb, Nc)

    return Xaero, varsb, Xptsb

def getRD(Xaero,pt,vars,Xpts):
    global num_nodes
    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)

    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros((3))
    Xa = np.zeros((9))
    Xb = np.zeros_like(X)
    Xab = np.zeros_like(Xa)

    Xa = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    J = np.zeros(9)
    h, J = jacobian3d(Xa)

    Ud = np.zeros(9)
    Ud = getDeformGradient(Ud, J, Na, Nb, Nc, vars)

    U = np.zeros(3)

    u = vars
    for i in range(num_nodes):
        U[0] = U[0] + N[i]*u[3*i+0]
        U[1] = U[1] + N[i]*u[3*i+1]
        U[2] = U[2] + N[i]*u[3*i+2]

    R = np.zeros(3)
    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]
    
    rot = np.zeros(3)
    rot[0] = 0.5*(Ud[5]-Ud[7])
    rot[1] = 0.5*(Ud[6]-Ud[2])
    rot[2] = 0.5*(Ud[1]-Ud[3])


    Uaero = np.zeros(3)
    Uaero[0] = U[0] + rot[1]*R[2] - rot[2]*R[1]
    Uaero[1] = U[1] + rot[2]*R[0] - rot[0]*R[2]
    Uaero[2] = U[2] + rot[0]*R[1] - rot[1]*R[0]

    return Uaero

def getRDXptSens(Xaero,pt,vars,Xpts):
    global num_nodes
    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)

    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros((3))
    Xa = np.zeros((9))

    Xa = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    XSens = np.zeros((3))
    XaSens = np.zeros((9))

    XaSens = solidjacobian(XSens,N, Na, Nb, Nc, Xpts)

    sh = np.zeros(1)
    J = np.zeros(9)
    JSens = np.zeros(9)

    J, JSens, sh, h = jacobian3dSens(Xa, XaSens)

    Ud = np.zeros(9)
    UdSens = np.zeros(9)
    Ud, UdSens = getDeformGradientSens(J, JSens, Na, Nb, Nc, vars)

    R = np.zeros(3)
    RSens = np.zeros(3)

    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]
    RSens[0] = Xaero[0] - XSens[0]
    RSens[1] = Xaero[1] - XSens[1]
    RSens[2] = Xaero[2] - XSens[2]
    
    rot = np.zeros(3)
    rotSens = np.zeros(3)
    rot[0] = 0.5*(Ud[5]-Ud[7])
    rot[1] = 0.5*(Ud[6]-Ud[2])
    rot[2] = 0.5*(Ud[1]-Ud[3])
    rotSens[0] = 0.5*(UdSens[5] - UdSens[7])
    rotSens[1] = 0.5*(UdSens[6] - UdSens[2])
    rotSens[2] = 0.5*(UdSens[1] - UdSens[3])

    UaeroSens = np.zeros(3)
    Uaero[0] = U[0] + rot[1]*R[2] - rot[2]*R[1]
    Uaero[1] = U[1] + rot[2]*R[0] - rot[0]*R[2]
    Uaero[2] = U[2] + rot[0]*R[1] - rot[1]*R[0]

    UaeroSens[0] = (rotSens[1]*R[2] - rotSens[2]*R[1] + 
		  rot[1]*RSens[2] - rot[2]*RSens[1])
    UaeroSens[1] = (rotSens[2]*R[0] - rotSens[0]*R[2] + 
		  rot[2]*RSens[0] - rot[0]*RSens[2])
    UaeroSens[2] = (rotSens[0]*R[1] - rotSens[1]*R[0] + 
		  rot[0]*RSens[1] - rot[1]*RSens[0])

    return UaeroSens

def jacobian3dSens(Xd,sXd):
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1])- Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2])+ Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
    hinv = 1.0/h

    sh= sXd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) + \
    Xd[8]*(Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[3]*sXd[1] - sXd[3]*Xd[1]) \
    - sXd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) - \
    Xd[7]*(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[3]*sXd[2] - sXd[3]*Xd[2]) \
    + sXd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]) + \
    Xd[6]*(Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4])\

    Jinv[0] =   (Xd[4]*Xd[8] - Xd[5]*Xd[7])*hinv
    Jinv[1] = - (Xd[1]*Xd[8] - Xd[2]*Xd[7])*hinv
    Jinv[2] =   (Xd[1]*Xd[5] - Xd[2]*Xd[4])*hinv
        
    Jinv[3] = - (Xd[3]*Xd[8] - Xd[5]*Xd[6])*hinv
    Jinv[4] =   (Xd[0]*Xd[8] - Xd[2]*Xd[6])*hinv
    Jinv[5] = - (Xd[0]*Xd[5] - Xd[2]*Xd[3])*hinv
        
    Jinv[6] =   (Xd[3]*Xd[7] - Xd[4]*Xd[6])*hinv
    Jinv[7] = - (Xd[0]*Xd[7] - Xd[1]*Xd[6])*hinv
    Jinv[8] =   (Xd[0]*Xd[4] - Xd[1]*Xd[3])*hinv

    for i in range(9):
        sJinv[i] = - Jinv[i]*hinv*sh
    
    sJinv[0] +=  (Xd[4]*sXd[8] + sXd[4]*Xd[8] - Xd[5]*sXd[7] - sXd[5]*Xd[7])*hinv
    sJinv[1] += -(Xd[1]*sXd[8] + sXd[1]*Xd[8] - Xd[2]*sXd[7] - sXd[2]*Xd[7])*hinv
    sJinv[2] +=  (Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4])*hinv
        
    sJinv[3] += -(Xd[3]*sXd[8] + sXd[3]*Xd[8] - Xd[5]*sXd[6] - sXd[5]*Xd[6])*hinv
    sJinv[4] +=  (Xd[0]*sXd[8] + sXd[0]*Xd[8] - Xd[2]*sXd[6] - sXd[2]*Xd[6])*hinv
    sJinv[5] += -(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[2]*sXd[3] - sXd[2]*Xd[3])*hinv
        
    sJinv[6] +=  (Xd[3]*sXd[7] + sXd[3]*Xd[7] - Xd[4]*sXd[6] - sXd[4]*Xd[6])*hinv
    sJinv[7] += -(Xd[0]*sXd[7] + sXd[0]*Xd[7] - Xd[1]*sXd[6] - sXd[1]*Xd[6])*hinv
    sJinv[8] +=  (Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[1]*sXd[3] - sXd[1]*Xd[3])*hinv

    return Jinv, sJinv, sh, h

def getDeformGradientSens(J,JSens, Na, Nb, Nc, vars):
    Ua = np.zeros(9)
    
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

    UdSens[0] = Ua[0]*JSens[0] + Ua[1]*JSens[3] + Ua[2]*JSens[6]
    UdSens[3] = Ua[3]*JSens[0] + Ua[4]*JSens[3] + Ua[5]*JSens[6]
    UdSens[6] = Ua[6]*JSens[0] + Ua[7]*JSens[3] + Ua[8]*JSens[6]

    UdSens[1] = Ua[0]*JSens[1] + Ua[1]*JSens[4] + Ua[2]*JSens[7]
    UdSens[4] = Ua[3]*JSens[1] + Ua[4]*JSens[4] + Ua[5]*JSens[7]
    UdSens[7] = Ua[6]*JSens[1] + Ua[7]*JSens[4] + Ua[8]*JSens[7]

    UdSens[2] = Ua[0]*JSens[2] + Ua[1]*JSens[5] + Ua[2]*JSens[8]
    UdSens[5] = Ua[3]*JSens[2] + Ua[4]*JSens[5] + Ua[5]*JSens[8]
    UdSens[8] = Ua[6]*JSens[2] + Ua[7]*JSens[5] + Ua[8]*JSens[8]

    return Ud, UdSens

def getshapefunction(pt, N, Na, Nb, Nc):
    order = 2
    na = np.zeros(order)
    dna = np.zeros_like(na)
    nb = np.zeros_like(na)
    dnb = np.zeros_like(nb)
    nc = np.zeros_like(na)
    dnc = np.zeros_like(nc)
    na, dna = lagrangeSF(na,dna, pt[0])
    nb, dnb = lagrangeSF(nb,dnb, pt[1])
    nc, dnc = lagrangeSF(nc,dnc, pt[2])
    
    
    index = 0
    for k in range(order):
        for j in range(order):
            for i in range(order):
                N[index] = na[i]*nb[j]*nc[k]
                print('dna[i]=',dna[i])
                print('nb[j]=',nb[j])
                print('nc[k]=',nc[k])
                Na[index] = dna[i]*nb[j]*nc[k]
                print('na[i]=',na[i])
                print('dnb[j]=',dnb[j])
                print('nc[k]=',nc[k])
                Nb[index] = na[i]*dnb[j]*nc[k]
                Nc[index] = na[i]*nb[j]*dnc[k]
                index = index +1
    print('in getshape=',N,Na, Nb, Nc)  

    return N, Na, Nb, Nc     

def lagrangeSF(sf, dsf, a):
    sf = np.array([1,1.6,1.4])+a*np.ones(3)
    dsf =a*np.array([2.3,1.1,0.5])
    return sf,dsf

def solidjacobian(X, N, Na, Nb, Nc, Xpts):
    Xa = np.zeros(9)
    X = np.zeros(3)
    print X
    print('Xpts=',Xpts)
    print('N=',N)
    print('Na=',Na)
    print('Nb=',Nb)
    
    for i in range(num_nodes):
        print('i=',i)
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

    print('in solidjacobian, Xd=',Xa)
    return Xa

def jacobian3d(Xd):
    print('Xd=',Xd)
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1])- Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2])+ Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
    print('h=',h)
    hinv = 1.0/h
    Jinv = np.zeros(9)
    Jinv[0] =   (Xd[4]*Xd[8] - Xd[5]*Xd[7])*hinv
    Jinv[1] = - (Xd[1]*Xd[8] - Xd[2]*Xd[7])*hinv
    Jinv[2] =   (Xd[1]*Xd[5] - Xd[2]*Xd[4])*hinv
        
    Jinv[3] = - (Xd[3]*Xd[8] - Xd[5]*Xd[6])*hinv
    Jinv[4] =   (Xd[0]*Xd[8] - Xd[2]*Xd[6])*hinv
    Jinv[5] = - (Xd[0]*Xd[5] - Xd[2]*Xd[3])*hinv
        
    Jinv[6] =   (Xd[3]*Xd[7] - Xd[4]*Xd[6])*hinv
    Jinv[7] = - (Xd[0]*Xd[7] - Xd[1]*Xd[6])*hinv
    Jinv[8] =   (Xd[0]*Xd[4] - Xd[1]*Xd[3])*hinv

    return h, Jinv
        

def getDeformGradient(Ud, J, Na, Nb, Nc, vars):
    Ua = np.zeros(9)
    
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

    return Ud

def getDeformGradientReverse(Udb, J, Jb, Na, Nb, Nc, vars, varsb):

    Ua = np.zeros(9)
    
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

    Uab = np.zeros(9)
    Uab[0] = Uab[0]+Udb[0]*J[0]
    Uab[1] = Uab[1]+Udb[0]*J[3]
    Uab[2] = Uab[2]+Udb[0]*J[6] 

    Jb[0] = Jb[0]+Udb[0]*Ua[0]
    Jb[3] = Jb[3]+Udb[0]*Ua[1]
    Jb[6] = Jb[6]+Udb[0]*Ua[2]
    
    Uab[3] = Uab[3]+Udb[3]*J[0]
    Uab[4] = Uab[4]+Udb[3]*J[3]
    Uab[2] = Uab[2]+Udb[3]*J[6]

    Jb[0] = Jb[0]+Udb[3]*Ua[3]
    Jb[3] = Jb[3]+Udb[3]*Ua[4]
    Jb[6] = Jb[6]+Udb[3]*Ua[5]

    Uab[6] = Uab[6]+Udb[6]*J[0]
    Uab[7] = Uab[7]+Udb[6]*J[3]
    Uab[8] = Uab[8]+Udb[6]*J[6]

    Jb[0] = Jb[0]+Udb[6]*Ua[6]
    Jb[3] = Jb[3]+Udb[6]*Ua[7]
    Jb[6] = Jb[6]+Udb[6]*Ua[8]

    Uab[0] = Uab[0]+Udb[1]*J[1]
    Uab[1] = Uab[1]+Udb[1]*J[4]
    Uab[2] = Uab[2]+Udb[1]*J[7]

    Jb[1] = Jb[1]+Udb[1]*Ua[0]
    Jb[4] = Jb[4]+Udb[1]*Ua[1]
    Jb[7] = Jb[7]+Udb[1]*Ua[2]

    Uab[3] = Uab[3]+Udb[4]*J[1]
    Uab[4] = Uab[4]+Udb[4]*J[4]
    Uab[5] = Uab[5]+Udb[4]*J[7] 

    Jb[1] = Jb[1]+Udb[4]*Ua[3]
    Jb[4] = Jb[4]+Udb[4]*Ua[4]
    Jb[7] = Jb[7]+Udb[4]*Ua[5]

    Uab[6] = Uab[6]+Udb[7]*J[1]
    Uab[7] = Uab[7]+Udb[7]*J[4]
    Uab[8] = Uab[8]+Udb[7]*J[7] 

    Jb[1] = Jb[1]+Udb[7]*Ua[6]
    Jb[4] = Jb[4]+Udb[7]*Ua[7]
    Jb[7] = Jb[7]+Udb[7]*Ua[8]

    Uab[0] = Uab[0]+Udb[2]*J[2]
    Uab[1] = Uab[1]+Udb[2]*J[5]
    Uab[2] = Uab[2]+Udb[2]*J[8]

    Jb[2] = Jb[2]+Udb[2]*Ua[0]
    Jb[5] = Jb[5]+Udb[2]*Ua[1]
    Jb[8] = Jb[8]+Udb[2]*Ua[2]
    
    Uab[3] = Uab[3]+Udb[5]*J[2]
    Uab[4] = Uab[4]+Udb[5]*J[5]
    Uab[5] = Uab[5]+Udb[5]*J[8]

    Jb[2] = Jb[2]+Udb[5]*Ua[3]
    Jb[5] = Jb[5]+Udb[5]*Ua[4]
    Jb[8] = Jb[8]+Udb[5]*Ua[5]

    Uab[6] = Uab[6]+Udb[8]*J[2]
    Uab[7] = Uab[7]+Udb[8]*J[5]
    Uab[8] = Uab[8]+Udb[8]*J[8]

    Jb[2] = Jb[2]+Udb[8]*Ua[6]
    Jb[5] = Jb[5]+Udb[8]*Ua[7]
    Jb[8] = Jb[8]+Udb[8]*Ua[8]

    for i in range(num_nodes):
        varsb[3*i] = varsb[3*i]+Uab[0]*Na[i]
        varsb[3*i] = varsb[3*i]+Uab[1]*Nb[i]
        varsb[3*i] = varsb[3*i]+Uab[2]*Nc[i]

        varsb[3*i+1] = varsb[3*i+1]+Uab[3]*Na[i]
        varsb[3*i+1] = varsb[3*i+1]+Uab[4]*Nb[i]
        varsb[3*i+1] = varsb[3*i+1]+Uab[5]*Nc[i]
        
        varsb[3*i+2] = varsb[3*i+2]+Uab[6]*Na[i]
        varsb[3*i+2] = varsb[3*i+2]+Uab[7]*Nb[i]
        varsb[3*i+2] = varsb[3*i+2]+Uab[8]*Nc[i]

    return Jb, varsb
    
        
def jacobian3dReverse(Xd, Jinvb):
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) 
                  - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) 
                  + Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
    hinv = 1.0/h
    Xdb = np.zeros(9)
    Xdb[0] = Xdb[0] + Xd[4]*hinv*Jinvb[8]
    Xdb[4] = Xdb[4] + Xd[0]*hinv*Jinvb[8]
    Xdb[1] = Xdb[1] - Xd[3]*hinv*Jinvb[8]
    Xdb[3] = Xdb[3] - Xd[1]*hinv*Jinvb[8]
    hinvb = (Xd[0]*Xd[4]-Xd[1]*Xd[3])*Jinvb[8]

    Xdb[0] = Xdb[0] + Xd[7]*(-(hinv*Jinvb[7]))
    Xdb[7] = Xdb[7] + Xd[0]*(-(hinv*Jinvb[7]))
    Xdb[1] = Xdb[1] - Xd[6]*(-(hinv*Jinvb[7]))
    Xdb[6] = Xdb[6] - Xd[1]*(-(hinv*Jinvb[7]))
    hinvb = hinvb - (Xd[0]*Xd[7]-Xd[1]*Xd[6])*Jinvb[7]

    Xdb[3] = Xdb[3] + Xd[7]*hinv*Jinvb[6]
    Xdb[7] = Xdb[7] + Xd[3]*hinv*Jinvb[6]
    Xdb[4] = Xdb[4] - Xd[6]*hinv*Jinvb[6]
    Xdb[6] = Xdb[6] - Xd[4]*hinv*Jinvb[6]
    hinvb = hinvb + (Xd[3]*Xd[7]-Xd[4]*Xd[6])*Jinvb[6]

    Xdb[0] = Xdb[0] + Xd[5]*(-(hinv*Jinvb[5]))
    Xdb[5] = Xdb[5] + Xd[0]*(-(hinv*Jinvb[5]))
    Xdb[2] = Xdb[2] - Xd[3]*(-(hinv*Jinvb[5]))
    Xdb[3] = Xdb[3] - Xd[2]*(-(hinv*Jinvb[5]))
    hinvb = hinvb - (Xd[0]*Xd[5]-Xd[2]*Xd[3])*Jinvb[5]

    Xdb[0] = Xdb[0] + Xd[8]*hinv*Jinvb[4]
    Xdb[8] = Xdb[8] + Xd[0]*hinv*Jinvb[4]
    Xdb[2] = Xdb[2] - Xd[6]*hinv*Jinvb[4]
    Xdb[6] = Xdb[6] - Xd[2]*hinv*Jinvb[4]
    hinvb = hinvb + (Xd[0]*Xd[8]-Xd[2]*Xd[6])*Jinvb[4]

    Xdb[3] = Xdb[3] + Xd[8]*(-(hinv*Jinvb[3]))
    Xdb[8] = Xdb[8] + Xd[3]*(-(hinv*Jinvb[3]))
    Xdb[5] = Xdb[5] - Xd[6]*(-(hinv*Jinvb[3]))
    Xdb[6] = Xdb[6] - Xd[5]*(-(hinv*Jinvb[3]))
    hinvb = hinvb - (Xd[3]*Xd[8]-Xd[5]*Xd[6])*Jinvb[3]

    Xdb[1] = Xdb[1] + Xd[5]*hinv*Jinvb[2]
    Xdb[5] = Xdb[5] + Xd[1]*hinv*Jinvb[2]
    Xdb[2] = Xdb[2] - Xd[4]*hinv*Jinvb[2]
    Xdb[4] = Xdb[4] - Xd[2]*hinv*Jinvb[2]
    hinvb = hinvb + (Xd[1]*Xd[5]-Xd[2]*Xd[4])*Jinvb[2]

    Xdb[1] = Xdb[1] + Xd[8]*(-(hinv*Jinvb[1]))
    Xdb[8] = Xdb[8] + Xd[1]*(-(hinv*Jinvb[1]))
    Xdb[2] = Xdb[2] - Xd[7]*(-(hinv*Jinvb[1]))
    Xdb[7] = Xdb[7] - Xd[2]*(-(hinv*Jinvb[1]))
    hinvb = hinvb - (Xd[1]*Xd[8]-Xd[2]*Xd[7])*Jinvb[1]

    Xdb[4] = Xdb[4] + Xd[8]*hinv*Jinvb[0]
    Xdb[8] = Xdb[8] + Xd[4]*hinv*Jinvb[0]
    Xdb[5] = Xdb[5] - Xd[7]*hinv*Jinvb[0]
    Xdb[7] = Xdb[7] - Xd[5]*hinv*Jinvb[0]
    hinvb = hinvb + (Xd[4]*Xd[8]-Xd[5]*Xd[7])*Jinvb[0]

    hb = - hinvb/(h*h)
    Xdb[8] = Xdb[8] + (Xd[0]*Xd[4]-Xd[3]*Xd[1])*hb
    Xdb[0] = Xdb[0] + Xd[5]*(-(Xd[7]*hb)) + Xd[4]*Xd[8]*hb
    Xdb[4] = Xdb[4] + Xd[0]*Xd[8]*hb - Xd[2]*Xd[6]*hb
    Xdb[3] = Xdb[3] - Xd[2]*(-(Xd[7]*hb)) - Xd[1]*Xd[8]*hb
    Xdb[1] = Xdb[1] + Xd[5]*Xd[6]*hb - Xd[3]*Xd[8]*hb
    Xdb[7] = Xdb[7] - (Xd[0]*Xd[5]-Xd[3]*Xd[2])*hb
    Xdb[5] = Xdb[5] + Xd[1]*Xd[6]*hb + Xd[0]*(-(Xd[7]*hb))
    Xdb[2] = Xdb[2] - Xd[4]*Xd[6]*hb - Xd[3]*(-(Xd[7]*hb))
    Xdb[6] = Xdb[6] + (Xd[1]*Xd[5]-Xd[2]*Xd[4])*hb
        
    return Xdb


def solidjacobianReverse(Xptsb, Xb, Xab, N, Na, Nb, Nc):
    global num_nodes


    for i in range(num_nodes):
        Xptsb[3*i] = Xptsb[3*i]+Xb[0]*N[i]
        Xptsb[3*i+1] = Xptsb[3*i+1]+Xb[1]*N[i]
        Xptsb[3*i+2] = Xptsb[3*i+2]+Xb[2]*N[i]
        
        Xptsb[3*i] = Xptsb[3*i]+Xab[0]*Na[i]
        Xptsb[3*i] = Xptsb[3*i]+Xab[1]*Nb[i]
        Xptsb[3*i] = Xptsb[3*i]+Xab[2]*Nc[i]

        Xptsb[3*i+1] = Xptsb[3*i+1]+Xab[3]*Na[i]
        Xptsb[3*i+1] = Xptsb[3*i+1]+Xab[4]*Nb[i]
        Xptsb[3*i+1] = Xptsb[3*i+1]+Xab[5]*Nc[i]
        
        Xptsb[3*i+2] = Xptsb[3*i+2]+Xab[6]*Na[i]
        Xptsb[3*i+2] = Xptsb[3*i+2]+Xab[7]*Nb[i]
        Xptsb[3*i+2] = Xptsb[3*i+2]+Xab[8]*Nc[i]

    return Xptsb 

num_nodes = 8
Xpts = 0.1*np.array([5,6,3,2,1,4,7,8,6,5,6,3,2,1,4,7,8,6,5,6,3,2,1,4,7,8,6])
vars = 0.1*np.array([2,5,6,6,7,5,4,7,1,2,5,6,6,7,5,4,7,1,2,5,6,6,7,5,4,7,1])
Uaerob = 0.1*np.array([5,4,7])
Xaero = 0.1*np.array([2,9,8])
pt = 0.1*np.array([0.2,0.5,0.25])



# Xaerob, varsb, Xptsb = getRDSens(Uaerob,Xaero, pt, vars, Xpts)
Uaero = getRD(Xaero,pt, vars, Xpts)
print('Uaero = ', Uaero)

print('X_aero = \n',Xaero)
print('varsb = \n',varsb)
print('Xptsb = \n',Xptsb)