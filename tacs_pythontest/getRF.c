void getRF();

int main(){
    float Xaero[9], Faero[3], res[24], Xpts[24], pt[3];

    float res = getRF(res, Faero, Xaero, pt, Xpts);
    return 0;
}

float getRF( float res[],
				      const float Faero[],
				      const float Xaero[],
				      const double pt[], 
				      const float Xpts[] ){
  // Compute the element shape functions
  double N[8];
  double Na[8], Nb[8], Nc[8];
  getShapeFunctions(pt, N, Na, Nb, Nc);

  // Compute the derivative of the position w.r.t. the parametric
  // coordinates 
  float X[3], Xa[3];
  solidJacobian(X, Xa, N, Na, Nb, Nc, Xpts);
  
  // Compute the inverse of the Jacobian
  float J[9];
  jacobian3d(Xa, J);

  // Compute the rigid link
  float R[3];
  R[0] = Xaero[0] - X[0];
  R[1] = Xaero[1] - X[1];
  R[2] = Xaero[2] - X[2];

  // Compute the contribution from each node
  double *n = N, *na = Na, *nb = Nb, *nc = Nc;
  for ( int i = 0; i < 8; i++ ){
    // Compute the derivative of the displacement gradient
    // with respect to the nodal variables
    if (1){
      float Dn[9];
      Dn[0] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[3] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      Dn[6] = na[0]*J[0] + nb[0]*J[3] + nc[0]*J[6];
      
      Dn[1] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[4] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      Dn[7] = na[0]*J[1] + nb[0]*J[4] + nc[0]*J[7];
      
      Dn[2] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[5] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      Dn[8] = na[0]*J[2] + nb[0]*J[5] + nc[0]*J[8];
      
      float rot[3];
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