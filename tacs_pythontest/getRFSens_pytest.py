import numpy as np
# from getRFSens_pytest import getshapefunction
# from getRFSens_pytest import jacobian3d
# from getRFSens_pytest import jacobian3dSens
# from getRFSens_pytest import jacobian3dReverse
# from getRFSens_pytest import solidjacobian
# from getRFSens_pytest import solidjacobianSens
# from getRFSens_pytest import solidjacobianReverse

def getRF(Faero, Xaero, pt, Xpts):
    global num_nodes
      

    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)

    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros(3)
    Xa = np.zeros(9)
 

    Xa, X = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    J = np.zeros(9)
    h, J = jacobian3d(Xa)


    R = np.zeros(3)
    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]

    
    n = N
    na = Na
    nb = Nb
    nc = Nc

    res = np.zeros(3*num_nodes)
    
    for i in range(num_nodes):
        if USE_RIGID_MOMENT:
            Dn = np.zeros(9)
            Dn[0] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[3] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[6] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            
            Dn[1] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[4] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[7] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            
            Dn[2] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[5] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[8] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            
            rot = np.zeros(3)
            rot[0] = 0.0
            rot[1] = -0.5*Dn[2]
            rot[2] =  0.5*Dn[1]
            res[3*i] = - (Faero[0]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) + 
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))
            
            rot[0] =  0.5*Dn[5]
            rot[1] =  0.0
            rot[2] = -0.5*Dn[3]
            res[3*i+1] = - (Faero[1]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) + 
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))
            
            rot[0] = -0.5*Dn[7]
            rot[1] =  0.5*Dn[6]
            rot[2] = 0.0
            res[3*i+2] = - (Faero[2]*n[i] +
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))
        
        else:
            res[3*i] = -Faero[0]*n[i]
            res[3*i+1] = -Faero[1]*n[i+1]
            res[3*i+2] = -Faero[2]*n[i+2]
      
    return res


def getRFSens(resb,Faero,  Xaero, pt, Xpts):
    global num_nodes
    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)


    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros(3)
    Xa = np.zeros(9)
    Xb = np.zeros_like(X)
    Xab = np.zeros_like(Xa)

    Xa, X = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    J = np.zeros(9)
    h, J = jacobian3d(Xa)


    R = np.zeros(3)
    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]

    
    n = N
    na = Na
    nb = Nb
    nc = Nc

    Faerob = np.zeros(3)
    Jb = np.zeros(9)
    Rb = np.zeros(3)
    
    rotb = np.zeros(3)
    Dnb = np.zeros(9)
    Dn = np.zeros(9)
    res = np.zeros(3*num_nodes)
    print('in reverse mode, Faero=', Faero)
    print('in reverse mode, na=', na)
    print('in reverse mode, R=', R)
    for i in range(num_nodes):
        if USE_RIGID_MOMENT:
            
            Dn[0] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[3] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[6] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            
            Dn[1] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[4] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[7] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            
            Dn[2] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[5] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[8] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            
            rot = np.zeros(3)

            # ===================================================
            rot[0] = -0.5*Dn[7]
            rot[1] =  0.5*Dn[6]
            rot[2] = 0.0
            res[3*i+2] = - (Faero[2]*n[i] +
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            Faerob[0] += - ((rot[1]*R[2] - rot[2]*R[1])) * resb[3*i+2]
            Faerob[1] += - ((rot[2]*R[0] - rot[0]*R[2])) * resb[3*i+2]
            Faerob[2] += - (n[i]+(rot[0]*R[1] - rot[1]*R[0])) * resb[3*i+2]


            rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[3*i+2]
            rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[3*i+2]
            rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[3*i+2]

            Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[3*i+2]
            Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[3*i+2]
            Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[3*i+2]

            Dnb[6] = 0.5*rotb[1]
            Dnb[7] = -0.5*rotb[0]
            # ===================================================


            rot[0] =  0.5*Dn[5]
            rot[1] =  0.0
            rot[2] = -0.5*Dn[3]
            res[3*i+1] = - (Faero[1]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            Faerob[0] += - ((rot[1]*R[2] - rot[2]*R[1])) * resb[3*i+1]
            Faerob[1] += - (n[i]+(rot[2]*R[0] - rot[0]*R[2])) * resb[3*i+1]
            Faerob[2] += - ((rot[0]*R[1] - rot[1]*R[0])) * resb[3*i+1]


            rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[3*i+1]
            rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[3*i+1]
            rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[3*i+1]

            Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[3*i+1]
            Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[3*i+1]
            Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[3*i+1]
            

            Dnb[3] = -0.5*rotb[2]
            Dnb[5] = 0.5*rotb[0]
            # ===================================================


            rot[0] = 0.0
            rot[1] = -0.5*Dn[2]
            rot[2] =  0.5*Dn[1]
            res[3*i] = - (Faero[0]*n[i] +
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            Faerob[0] += - (n[i] + (rot[1]*R[2] - rot[2]*R[1])) * resb[3*i]
            Faerob[1] += - ((rot[2]*R[0] - rot[0]*R[2])) * resb[3*i]
            Faerob[2] += - ((rot[0]*R[1] - rot[1]*R[0])) * resb[3*i]

            
            rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[3*i]
            rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[3*i]
            rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[3*i]

            Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[3*i]
            Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[3*i]
            Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[3*i]
          
            
            Dnb[2] = -0.5*rotb[1]
            Dnb[1] = 0.5*rotb[2]
            # ===================================================

            Jb[0] +=  na[i] *(Dnb[0]+Dnb[3]+Dnb[6])
            Jb[3] +=  nb[i] *(Dnb[0]+Dnb[3]+Dnb[6])
            Jb[6] +=  nc[i] *(Dnb[0]+Dnb[3]+Dnb[6])

            Jb[1] +=  na[i] *(Dnb[1]+Dnb[4]+Dnb[7])
            Jb[4] +=  nb[i] *(Dnb[1]+Dnb[4]+Dnb[7])
            Jb[7] +=  nc[i] *(Dnb[1]+Dnb[4]+Dnb[7])

            Jb[2] +=  na[i] *(Dnb[2]+Dnb[5]+Dnb[8])
            Jb[5] +=  nb[i] *(Dnb[2]+Dnb[5]+Dnb[8])
            Jb[8] +=  nc[i] *(Dnb[2]+Dnb[5]+Dnb[8])
        else:

            Faerob[0] += -n[i]*resb[3*i]
            Faerob[1] += -n[i+1]*resb[3*i+1]
            Faerob[2] += -n[i+2]*resb[3*i+2]

        # // Increment the res pointers
    print('in reverse, res=',res)
    # // Compute the rigid link Reverse
    Xaerob = np.zeros(3)
    Xb = np.zeros(3)
    Xptsb = np.zeros(3*num_nodes)
    Xaerob[0] += Rb[0]
    Xb[0] += -Rb[0]
    Xaerob[1] += Rb[1]
    Xb[1] += -Rb[1]
    Xaerob[2] += Rb[2]
    Xb[2] += -Rb[2]

    Xab = jacobian3dReverse(Xa, Jb)

    Xptsb = solidjacobianReverse(Xb, Xab, N, Na, Nb, Nc)

    return Faerob, Xaerob, Xptsb

def getRFforward(Faero,FaeroSens, Xaero,XaeroSens, pt, Xpts,XptsSens):
    global num_nodes

    N = np.zeros(num_nodes)
    Na = np.zeros(num_nodes)
    Nb = np.zeros(num_nodes)
    Nc = np.zeros(num_nodes)

    N, Na, Nb, Nc = getshapefunction(pt, N, Na, Nb, Nc)

    X = np.zeros(3)
    Xa = np.zeros(9)

    Xa, X = solidjacobian(X,N, Na, Nb, Nc, Xpts)

    XSens = np.zeros(3)
    XaSens = np.zeros(9)

    XaSens, XSens = solidjacobianSens(XSens,N, Na, Nb, Nc, XptsSens)

    sh = np.zeros(1)
    J = np.zeros(9)
    JSens = np.zeros(9)

    J, JSens, sh, h = jacobian3dSens(Xa, XaSens, sh)

    R = np.zeros(3)
    RSens = np.zeros(3)
    res = np.zeros(3*num_nodes)
    resSens = np.zeros(3*num_nodes)

    R[0] = Xaero[0] - X[0]
    R[1] = Xaero[1] - X[1]
    R[2] = Xaero[2] - X[2]
    RSens[0] = XaeroSens[0] - XSens[0]
    RSens[1] = XaeroSens[1] - XSens[1]
    RSens[2] = XaeroSens[2] - XSens[2]


    n = N
    na = Na
    nb = Nb
    nc = Nc

    DnSens = np.zeros(9)
    rot = np.zeros(3)
    rotSens = np.zeros(3)
    resSens = np.zeros(3*num_nodes)

    print('in forward mode, Faero=', Faero)
    print('in forward mode, na=', na)
    print('in forward mode, R=', R)
    for i in range(num_nodes):
        if USE_RIGID_MOMENT:
            Dn = np.zeros(9)
            Dn[0] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[3] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            Dn[6] = na[i]*J[0] + nb[i]*J[3] + nc[i]*J[6]
            
            Dn[1] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[4] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            Dn[7] = na[i]*J[1] + nb[i]*J[4] + nc[i]*J[7]
            
            Dn[2] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[5] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]
            Dn[8] = na[i]*J[2] + nb[i]*J[5] + nc[i]*J[8]

            
            DnSens[0] = na[i]*JSens[0] + nb[i]*JSens[3] + nc[i]*JSens[6]
            DnSens[3] = na[i]*JSens[0] + nb[i]*JSens[3] + nc[i]*JSens[6]
            DnSens[6] = na[i]*JSens[0] + nb[i]*JSens[3] + nc[i]*JSens[6]
            
            DnSens[1] = na[i]*JSens[1] + nb[i]*JSens[4] + nc[i]*JSens[7]
            DnSens[4] = na[i]*JSens[1] + nb[i]*JSens[4] + nc[i]*JSens[7]
            DnSens[7] = na[i]*JSens[1] + nb[i]*JSens[4] + nc[i]*JSens[7]
            
            DnSens[2] = na[i]*JSens[2] + nb[i]*JSens[5] + nc[i]*JSens[8]
            DnSens[5] = na[i]*JSens[2] + nb[i]*JSens[5] + nc[i]*JSens[8]
            DnSens[8] = na[i]*JSens[2] + nb[i]*JSens[5] + nc[i]*JSens[8]
            
            
            rot[0] = 0.0
            rot[1] = -0.5*Dn[2]
            rot[2] =  0.5*Dn[1]

            rotSens[0] = 0.0
            rotSens[1] = -0.5*DnSens[2]
            rotSens[2] =  0.5*DnSens[1]

            res[3*i] = - (Faero[0]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) + 
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            
            resSens[3*i] = - (FaeroSens[0]*n[i]
             + FaeroSens[0]*(rot[1]*R[2] - rot[2]*R[1])
             + Faero[0]*(rotSens[1]*R[2] + rot[1]*RSens[2] - rotSens[2]*R[1] - rot[2]*RSens[1])
             + FaeroSens[1]*(rot[2]*R[0] - rot[0]*R[2])
             + Faero[1]*(rotSens[2]*R[0] + rot[2]*RSens[0] - rotSens[0]*R[2] - rot[0]*RSens[2])
             + FaeroSens[2]*(rot[0]*R[1] - rot[1]*R[0])
             + Faero[2]*(rotSens[0]*R[1] + rot[0]*RSens[1] - rotSens[1]*R[0] - rot[1]*RSens[0]))

            
            rot[0] =  0.5*Dn[5]
            rot[1] =  0.0
            rot[2] = -0.5*Dn[3]

            rotSens[0] =  0.5*DnSens[5]
            rotSens[1] =  0.0
            rotSens[2] = -0.5*DnSens[3]



            res[3*i+1] = - (Faero[1]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) + 
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            resSens[3*i+1] = - (FaeroSens[1]*n[i] + 
            FaeroSens[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[0]*(rotSens[1]*R[2] + rot[1]*RSens[2] - rotSens[2]*R[1] - rot[2]*RSens[1]) +
            FaeroSens[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[1]*(rotSens[2]*R[0] + rot[2]*RSens[0] - rotSens[0]*R[2] - rot[0]*RSens[2]) + 
            FaeroSens[2]*(rot[0]*R[1] - rot[1]*R[0]) +
            Faero[2]*(rotSens[0]*R[1] + rot[0]*RSens[1] - rotSens[1]*R[0] - rot[1]*RSens[0]))

            
            rot[0] = -0.5*Dn[7]
            rot[1] =  0.5*Dn[6]
            rot[2] = 0.0

            rotSens[0] = -0.5*DnSens[7]
            rotSens[1] =  0.5*DnSens[6]
            rotSens[2] = 0.0

            res[3*i+2] = - (Faero[2]*n[i] + 
            Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) + 
            Faero[2]*(rot[0]*R[1] - rot[1]*R[0]))

            resSens[3*i+2] = - (FaeroSens[2]*n[i] + 
            FaeroSens[0]*(rot[1]*R[2] - rot[2]*R[1]) + 
            Faero[0]*(rotSens[1]*R[2] + rot[1]*RSens[2] - rotSens[2]*R[1] - rot[2]*RSens[1]) +
            FaeroSens[1]*(rot[2]*R[0] - rot[0]*R[2]) +
            Faero[1]*(rotSens[2]*R[0] + rot[2]*RSens[0] - rotSens[0]*R[2] - rot[0]*RSens[2]) + 
            FaeroSens[2]*(rot[0]*R[1] - rot[1]*R[0]) +
            Faero[2]*(rotSens[0]*R[1] + rot[0]*RSens[1] - rotSens[1]*R[0] - rot[1]*RSens[0]))


        
        else:
            resSens[3*i] = -FaeroSens[0]*n[i]
            resSens[3*i+1] = -FaeroSens[1]*n[i+1]
            resSens[3*i+2] = -FaeroSens[2]*n[i+2]
    
    print('in forward, res=',res)
    return resSens

def jacobian3dSens(Xd,sXd,sh):
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1])- Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2])+ Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
    hinv = 1.0/h

    sh= sXd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) + \
    Xd[8]*(Xd[0]*sXd[4] + sXd[0]*Xd[4] - Xd[3]*sXd[1] - sXd[3]*Xd[1]) \
    - sXd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) - \
    Xd[7]*(Xd[0]*sXd[5] + sXd[0]*Xd[5] - Xd[3]*sXd[2] - sXd[3]*Xd[2]) \
    + sXd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]) + \
    Xd[6]*(Xd[1]*sXd[5] + sXd[1]*Xd[5] - Xd[2]*sXd[4] - sXd[2]*Xd[4])\

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

    sJinv = np.zeros(9)
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
                Na[index] = dna[i]*nb[j]*nc[k]
                Nb[index] = na[i]*dnb[j]*nc[k]
                Nc[index] = na[i]*nb[j]*dnc[k]
                index = index +1

    return N, Na, Nb, Nc 

def lagrangeSF(sf, dsf, a):
    sf = np.array([1,1.6,1.4])+a*np.ones(3)
    dsf =a*np.array([2.3,1.1,0.5])
    return sf,dsf



def solidjacobian(X, N, Na, Nb, Nc, Xpts):
    Xa = np.zeros(9)
    X = np.zeros(3)
    
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

    return Xa, X

def solidjacobianSens(XSens, N, Na, Nb, Nc, XptsSens):
    XaSens = np.zeros(9)
    XSens = np.zeros(3)
    
    for i in range(num_nodes):
        XSens[0] += XptsSens[3*i]*N[i]
        XSens[1] += XptsSens[3*i+1]*N[i]
        XSens[2] += XptsSens[3*i+2]*N[i]

        XaSens[0] += XptsSens[3*i]*Na[i]
        XaSens[1] += XptsSens[3*i]*Nb[i]
        XaSens[2] += XptsSens[3*i]*Nc[i]

        XaSens[3] += XptsSens[3*i+1]*Na[i]
        XaSens[4] += XptsSens[3*i+1]*Nb[i]
        XaSens[5] += XptsSens[3*i+1]*Nc[i]
        
        XaSens[6] += XptsSens[3*i+2]*Na[i]
        XaSens[7] += XptsSens[3*i+2]*Nb[i]
        XaSens[8] += XptsSens[3*i+2]*Nc[i]

    return XaSens, XSens

def jacobian3d(Xd):

    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1])- Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2])+ Xd[6]*(Xd[1]*Xd[5] - Xd[2]*Xd[4]))
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


def jacobian3dReverse(Xd, Jinvb):
    h = (Xd[8]*(Xd[0]*Xd[4] - Xd[3]*Xd[1]) \
                  - Xd[7]*(Xd[0]*Xd[5] - Xd[3]*Xd[2]) \
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


def solidjacobianReverse(Xb, Xab, N, Na, Nb, Nc):
    Xptsb = np.zeros(3*num_nodes)

    for i in range(num_nodes):
        Xptsb[3*i] = Xptsb[3*i]+Xb[0]*N[i]
        Xptsb[3*i+1] = Xptsb[3*i+1]+Xb[1]*N[i]
        Xptsb[3*i+2] = Xptsb[3*i+2]+Xb[2]*N[i]
        
        Xptsb[3*i] = Xptsb[3*i]+Xab[0]*Na[i]+Xab[1]*Nb[i]+Xab[2]*Nc[i]

        Xptsb[3*i+1] = Xptsb[3*i+1]+Xab[3]*Na[i]+Xab[4]*Nb[i]+Xab[5]*Nc[i]

        Xptsb[3*i+2] = Xptsb[3*i+2]+Xab[6]*Na[i]+Xab[7]*Nb[i]+Xab[8]*Nc[i]


    return Xptsb 


USE_RIGID_MOMENT = 1
num_nodes = 8
Xpts = np.array([5,6,3,2,1,4,7,8,6,5,6,3,2,1,4,7,8,6,5,6,3,2,1,4])
dXptsf = 0.1*np.array([5,6,1,2,4,8,7,8,6,6,8,3,2,1,7,7,8,6,5,6,1,2,1,4])

# XptsSens = 0.1*np.array([4,7,0,3,4,7,5,3,4,7,1,9,3,4,6,8,1,4,6,5,7,3,7,2,3,4,1])
Faero = 1.0*np.array([2,5,6])
dFaerof = 0.1*np.array([2,6,1])
resb = 0.1*np.array([5,7,1,2,4,4,7,5,3,1,8,3,2,1,7,2,8,6,5,3,1,7,3,4])

Xaero = 0.1*np.array([2,9,8])
XaeroSens = 0.1*np.array([6,4,1])
# XaeroSens = np.zeros(3)
pt = 0.1*np.array([0.2,0.8,0.25])


# 1. Start check the overall derivative ==================

res = getRF(Faero, Xaero, pt, Xpts)
print('res = ', res)
step = 1e-7
dXpts = step*dXptsf
Xpts_2 = Xpts+dXpts
dFaero = step*dFaerof
Faero_2 = Faero+dFaero
dXaero = step*XaeroSens
Xaero_2 = Xaero+dXaero
res_2 = getRF(Faero_2, Xaero_2, pt, Xpts_2)
FD = (res_2-res)/step
print('FD_res = ', FD)
resSens = getRFforward(Faero,dFaerof, Xaero,XaeroSens, pt, Xpts,dXptsf)
print('resSens = ', resSens)

Faerob, Xaerob, Xptsb = getRFSens(resb,Faero,  Xaero, pt, Xpts)

print('Xptsb=',Xptsb)
print('Xaerob=',Xaerob)
print('Faerob=',Faerob)

A_dotproc = np.inner(resSens,resb)
y_dotproc = np.inner(dXptsf,Xptsb) + np.inner(dFaerof,Faerob) +np.inner(XaeroSens,Xaerob)
print('A_dotproc = ',A_dotproc)
print('y_dotproc = ',y_dotproc)