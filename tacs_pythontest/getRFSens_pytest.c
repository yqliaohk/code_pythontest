template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRF( TacsScalar res[],
				      const TacsScalar Faero[],
				      const TacsScalar Xaero[],
				      const double pt[], 
				      const TacsScalar Xpts[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the Jacobian
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the rigid link
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  // Compute the contribution from each node
  double *n = N, *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Compute the derivative of the displacement gradient
    // with respect to the nodal variables
    if (USE_RIGID_MOMENT){
      TacsScalar Dn[9];
      Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      
      Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      
      Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      
      TacsScalar rot[3];
      rot[0] = 0.0;
      rot[1] = -0.5*Dn[2];
      rot[2] =  0.5*Dn[1];
      res[0] = - (Faero[0]*n[0] +
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
      
      rot[0] =  0.5*Dn[5];
      rot[1] =  0.0;
      rot[2] = -0.5*Dn[3];
      res[1] = - (Faero[1]*n[0] + 
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
      
      rot[0] = -0.5*Dn[7];
      rot[1] =  0.5*Dn[6];
      rot[2] = 0.0;
      res[2] = - (Faero[2]*n[0] +
		  Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));
    }
    else {
      res[0] = -Faero[0]*n[0];
      res[1] = -Faero[1]*n[1];
      res[2] = -Faero[2]*n[2];
    }

    // Increment the res pointers
    res += 3;    
    n++; na++; nb++; nc++;
  }
}




template <int NUM_NODES>
void TACS3DElement<NUM_NODES>::getRFSens(const TacsScalar resb[],
				      const TacsScalar Faero[], TacsScalar Faerob[],
				      const TacsScalar Xaero[], TacsScalar Xaerob[],
				      const double pt[], 
				      const TacsScalar Xpts[], TacsScalar Xptsb[] ){
  // Compute the element shape functions
  double N[NUM_NODES];
  double Na[NUM_NODES], Nb[NUM_NODES], Nc[NUM_NODES];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  TacsScalar X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the Jacobian
  TacsScalar J[9];
  FElibrary::jacobian3d(Xa, J);

  // Compute the rigid link
  TacsScalar R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  TacsScalar Faerob[3]
  TacsScalar Rb[3]
  TacsScalar Jb[9];

  // Compute the contribution from each node
  double *n = N, *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < NUM_NODES; i++ ){
    // Compute the derivative of the displacement gradient
    // with respect to the nodal variables
    if (USE_RIGID_MOMENT){
      TacsScalar Dn[9];
      Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      
      Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      
      Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      
      TacsScalar rot[3];
      rot[0] = 0.0;
      rot[1] = -0.5*Dn[2];
      rot[2] =  0.5*Dn[1];
      // res[0] = - (Faero[0]*n[0] +
		  // Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  // Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  // Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));

      Faerob[0] += - (n[0] + (rot[1]*R[2] - rot[2]*R[1])) * resb[0];
      Faerob[1] += - ((rot[2]*R[0] - rot[0]*R[2])) * resb[0];
      Faerob[2] += - ((rot[0]*R[1] - rot[1]*R[0])) * resb[0];

      TacsScalar rotb[3];
      rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[0];
      rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[0];
      rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[0];

      Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[0];
      Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[0];
      Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[0];
      
      TacsScalar Dnb[9];
      Dnb[2] = -0.5*rotb[1];
      Dnb[1] = 0.5*rotb[2];
      

      
      rot[0] =  0.5*Dn[5];
      rot[1] =  0.0;
      rot[2] = -0.5*Dn[3]

      // res[1] = - (Faero[1]*n[0] + 
		  // Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  // Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  // Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));

      Faerob[0] += - ((rot[1]*R[2] - rot[2]*R[1])) * resb[1];
      Faerob[1] += - (n[0]+(rot[2]*R[0] - rot[0]*R[2])) * resb[1];
      Faerob[2] += - ((rot[0]*R[1] - rot[1]*R[0])) * resb[1];


      rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[1];
      rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[1];
      rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[1];

      Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[1];
      Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[1];
      Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[1];
      

      Dnb[3] = -0.5*rotb[2];
      Dnb[5] = 0.5*rotb[0];

      
      rot[0] = -0.5*Dn[7];
      rot[1] =  0.5*Dn[6];
      rot[2] = 0.0;
      // res[2] = - (Faero[2]*n[0] +
		  // Faero[0]*(rot[1]*R[2] - rot[2]*R[1]) +
		  // Faero[1]*(rot[2]*R[0] - rot[0]*R[2]) +
		  // Faero[2]*(rot[0]*R[1] - rot[1]*R[0]));

      Faerob[0] += - ((rot[1]*R[2] - rot[2]*R[1])) * resb[2];
      Faerob[1] += - ((rot[2]*R[0] - rot[0]*R[2])) * resb[2];
      Faerob[2] += - (n[0]+(rot[0]*R[1] - rot[1]*R[0])) * resb[2];


      rotb[0] = - (Faero[2]*R[1]-Faero[1]*R[2]) * resb[2];
      rotb[1] = - (Faero[0]*R[2]-Faero[2]*R[0]) * resb[2];
      rotb[2] = - (Faero[1]*R[0]-Faero[0]*R[1]) * resb[2];

      Rb[0] += - (Faero[1]*rot[2]-Faero[2]*rot[1]) * resb[2];
      Rb[1] += - (Faero[2]*rot[0]-Faero[0]*rot[2]) * resb[2];
      Rb[2] += - (Faero[0]*rot[1]-Faero[1]*rot[0]) * resb[2];
      

      Dnb[6] = 0.5*rotb[1];
      Dnb[7] = -0.5*rotb[0];

      
      Jb[0] +=  na[0] *(Dnb[0]+Dnb[3]+Dnb[6]);
      Jb[3] +=  nb[0] *(Dnb[0]+Dnb[3]+Dnb[6]);
      Jb[6] +=  nc[0] *(Dnb[0]+Dnb[3]+Dnb[6]);

      Jb[1] +=  na[0] *(Dnb[1]+Dnb[4]+Dnb[7]);
      Jb[4] +=  nb[0] *(Dnb[1]+Dnb[4]+Dnb[7]);
      Jb[7] +=  nc[0] *(Dnb[1]+Dnb[4]+Dnb[7]);

      Jb[2] +=  na[0] *(Dnb[2]+Dnb[5]+Dnb[8]);
      Jb[5] +=  nb[0] *(Dnb[2]+Dnb[5]+Dnb[8]);
      Jb[8] +=  nc[0] *(Dnb[2]+Dnb[5]+Dnb[8]);

    }
    else {
      // res[0] = -Faero[0]*n[0];
      // res[1] = -Faero[1]*n[1];
      // res[2] = -Faero[2]*n[2];

      Faerob[0] += -n[0]*resb[0];
      Faerob[1] += -n[1]*resb[1];
      Faerob[2] += -n[2]*resb[2];
    }

    // Increment the res pointers
    resb += 3;    
    n++; na++; nb++; nc++;
  }

  // Compute the rigid link Reverse
  TacsScalar Xb[3];
  Xaerob[0] = Rb[0];
  Xb[0] = -Rb[0];
  Xaerob[1] = Rb[1];
  Xb[1] = -Rb[1];
  Xaerob[2] = Rb[2];
  Xb[2] = -Rb[2];

  FElibrary::jacobian3dReverse(Xab, Xa, Jb);

  solidJacobianReverse(Xb, Xab, N, Na, Nb, Nc, Xptsb);

              }

  