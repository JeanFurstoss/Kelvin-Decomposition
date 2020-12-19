#include <iostream>
#include <math.h>
#include <vector>
using namespace std;

// Construc the different Kelvin bases associated to different crystal
// symmetries mainly based on the work :
// EigenTensors of linear anisotropic elastic materials. M.M.
// Mehrabadi and S.C. Cowin. 1990.

class KelvinBase{
protected:
  vector<double*> *EigenTensors;
  double* EigenValues;
  double* ElasticConstants;
  int **myClassicToVoigt;
  double * SOTensor_tmp;
  double * FOTensor_tmp;
public:
  KelvinBase(double* ElasticConstants):
    ElasticConstants(ElasticConstants) {
    EigenValues = new double[6];
    EigenTensors = new vector<double*>();
    for (int i = 0; i < 6; ++i) {
      EigenValues[i] = 0;
      EigenTensors->push_back(new double[9]);
      for (int j = 0; j < 9; ++j) {
	(*EigenTensors)[i][j] = 0.;
      }
    }

    myClassicToVoigt = new int*[3];
    for (int i = 0; i < 3; ++i)
      myClassicToVoigt[i] = new int[3];

    myClassicToVoigt[0][0] = 0;
    myClassicToVoigt[1][1] = 1;
    myClassicToVoigt[2][2] = 2;
    myClassicToVoigt[0][1] = 3;
    myClassicToVoigt[1][0] = 3;
    myClassicToVoigt[0][2] = 4;
    myClassicToVoigt[2][0] = 4;
    myClassicToVoigt[1][2] = 5;
    myClassicToVoigt[2][1] = 5;

    SOTensor_tmp = new double[9];
    FOTensor_tmp = new double[81];
    
  }
  ~KelvinBase() {
    for (unsigned int i = 0; i < (*EigenTensors).size(); i++)
      delete[] EigenTensors->at(i);

    delete EigenTensors;
    delete[] EigenValues;

    for (int i = 0; i < 3; ++i)
      delete[] myClassicToVoigt[i];
    delete[] myClassicToVoigt;

    delete[] SOTensor_tmp;
    delete[] FOTensor_tmp;
  }
  double* getEigenValues() { return EigenValues; }
  vector<double*>* getEigenTensors() { return EigenTensors; }
  void printFullBase();
  void printBase();
  virtual void constructBase() = 0;
  void EigenValues3x3(const double* A, double * EigenValues);
  void FOVoigt2Tensor(const double* FOVoigt, double* tensor);
  void FOTensorDotSOTensor(const double* FOTensor,
			   const double* SOTensor, double* product);
  void CheckBase();
};

void KelvinBase::EigenValues3x3(const double* A, double * EigenValues) {
	
  double TraceA = 0;
  for (int i = 0; i < 3; i++) {
    TraceA += A[i*3 + i];
  }

  double *Adev = SOTensor_tmp;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
	Adev[i*3 + j] = A[i*3 + j];
    }
    Adev[i*3 + i] -= (TraceA/3);
  }

  double J2 = sqrt ( sqrt ( pow ( ( ( ( Adev[0] * Adev[4] ) + ( Adev[4] * Adev[8] ) + ( Adev[0] * Adev[8] )
    - ( Adev[1] * Adev[3] ) - ( Adev[5] * Adev[7] ) - ( Adev[2] * Adev[6] ) ) / 3 ) ,2 ) ) );

  for (int i = 0; i<9; i++)
    Adev[i] /= J2;
	

  double J3;
  J3 = ( Adev[0] * ( (Adev[4]*Adev[8]) - (Adev[5]*Adev[7]) ) )
    -  ( Adev[1] * ( (Adev[3]*Adev[8]) - (Adev[5]*Adev[6]) ) )
    +  ( Adev[2] * ( (Adev[3]*Adev[7]) - (Adev[4]*Adev[6]) ) );
  
  double alpha;
  alpha = acos(J3/2)/3;
  for (int i = 0; i<3; i++) {
    EigenValues[i] = ( J2 * 2 * cos(alpha) ) + (TraceA/3);
    alpha += (2 * M_PI / 3);
  }
	
}

void KelvinBase::FOVoigt2Tensor(const double* FOVoigt,
				double* tensor) {
 
  int dim3 = 27;
  int dim2 = 9;
  int dim = 3;

  for (int i = 0; i < dim; ++i) {
    for (int j = 0; j < dim; ++j) {
      for (int k = 0; k < dim; ++k) {
	for (int l = 0; l < dim; ++l) {
	  tensor[i * dim3 + j * dim2 + k * dim + l] =
	    FOVoigt[myClassicToVoigt[i][j] * 6 +
		    myClassicToVoigt[k][l]];
	}
      }
    }
  }

}

void KelvinBase::FOTensorDotSOTensor(const double* FOTensor,
				     const double* SOTensor,
				     double* product) {
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      SOTensor_tmp[i*3+j] = SOTensor[i*3 + j];
      product[i*3 + j]=0;
    }
  }
	
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
	for(int l = 0; l < 3; l++) {
	  product[i*3 + j] +=
	    FOTensor[i*27+j*9+k*3+l]*SOTensor_tmp[k*3 + l];
	}
      }
    }
  }
	
}

void KelvinBase::CheckBase() {
  double *C = FOTensor_tmp;
  double *tmp = SOTensor_tmp;
  FOVoigt2Tensor(ElasticConstants, C);
  std::cout << "Fourth order elastic contants" << std::endl;
  for (int i=0;i<3;i++) {
    for (int j=0;j<3;j++) {
      for (int k=0;k<3;k++) {
	for (int l=0;l<3;l++) {
	  std::cout << C[i*27 + j*9 + k*3 + l] << " " ;
	}
	std::cout << std::endl;
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
  std::cout << "We check that the base tensors are actual "
	    << "eigentensors of the fourth order elastic tensor"
	    << std::endl;

  for (int ti=0; ti<6; ti++) {
    FOTensorDotSOTensor(C, (*EigenTensors)[ti], tmp);
    std::cout << "T"<< ti << " Tensor before multipliaction: " << std::endl;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
	cout << (*EigenTensors)[ti][i*3+j] << " ";
      }
      cout << endl;
    }
    std::cout << "T"<< ti << " Tensor after multiplication: " << std::endl;
    for(int i = 0; i < 3; i++) {
      for(int j = 0; j < 3; j++) {
	cout << tmp[i*3+j] << " ";
      }
      cout << endl;
    }
    std::cout << " Expected eigenvalue: " << EigenValues[ti]
	      << std::endl;
    for (int i=0;i<3;i++) {
      for (int j=0;j<3;j++) {
	std::cout << tmp[i*3+j]/((*EigenTensors)[ti][i*3+j]) << " ";
      }
    }
    std::cout << std::endl;
  }
}

//This method prints the tensors of the Kelvin base, their norms and all deltas
void KelvinBase::printFullBase() {	
  std::cout << "$$$$   Complete Kelvin base    $$$$" << std::endl;

  for (int ti=0; ti<6; ti++) {
    double NormB_ti(0);
    for (int i = 0; i<9; i++)
      NormB_ti += pow((*EigenTensors)[ti][i],2);
    
    double delta_ti1(0), delta_ti2(0), delta_ti3(0), delta_ti4(0),
      delta_ti5(0), delta_ti6(0);
    for (int i = 0; i<3; i++) {
      for (int j = 0; j<3; j++) {
	delta_ti1 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[0][j*3+i];
	delta_ti2 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[1][j*3+i];
	delta_ti3 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[2][j*3+i];
	delta_ti4 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[3][j*3+i];
	delta_ti5 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[4][j*3+i];
	delta_ti6 +=
	  (*EigenTensors)[ti][i*3+j]*(*EigenTensors)[5][j*3+i];
      }
    }
    std::cout << "Delta" << ti+1 << "1: " << delta_ti1
	      << "Delta" << ti+1 << "2: " << delta_ti2
	      << "Delta" << ti+1 << "3: " << delta_ti3
	      << "Delta" << ti+1 << "4: " << delta_ti4
	      << "Delta" << ti+1 << "5: " << delta_ti5
	      << "Delta" << ti+1 << "6: " << delta_ti6
	      << std::endl;
	
    std::cout << "K" << ti+1 << "1, eigenvalue = " << EigenValues[ti]
	      << ", eigen tensor :" << std::endl;
	
    for (int i = 0; i<3; i++) {
      for (int j = 0; j<3; j++) {
	std::cout << (*EigenTensors)[ti][i * 3 + j] << " " ;
      }
      std::cout << std::endl;
    }
    std::cout << "Norm tensor "<< ti+1 << " : "
	      << NormB_ti << std::endl;
  }

}

//This method print the base only
void KelvinBase::printBase() {
  for (int ti=0; ti<6; ti++) {
    std::cout << "K" << ti+1 << ", eigenvalue = " << EigenValues[ti]
	      << ", eigen tensor :" << std::endl;
	
    for (int i = 0; i<3; i++) {
      for (int j = 0; j<3; j++) {
	std::cout << (*EigenTensors)[ti][i * 3 + j] << " " ;
      }
      std::cout << std::endl;
    }
  }
}

class IsotropicKelvinBase : public KelvinBase{
public:
  void constructBase();
  IsotropicKelvinBase(double* ElasticConstants):KelvinBase(ElasticConstants) {
    constructBase();
  };
  ~IsotropicKelvinBase() {
  };
};

void IsotropicKelvinBase::constructBase() {
// First five eigenvalues (multiplicity 5) = 2 * C44
  EigenValues[0] = 2 * ElasticConstants[21];
  EigenValues[1] = 2 * ElasticConstants[21];
  EigenValues[2] = 2 * ElasticConstants[21];
  EigenValues[3] = 2 * ElasticConstants[21];
  EigenValues[4] = 2 * ElasticConstants[21];
// Last eigenvalue associated with dilatational strain C11 + 2 * C12
  EigenValues[5] = ElasticConstants[0] + 2 * ElasticConstants[1];

  // b1
  (*(EigenTensors))[0][0] = -1.0/sqrt(6.0);
  (*(EigenTensors))[0][1] = 0;
  (*(EigenTensors))[0][2] = 0;
  (*(EigenTensors))[0][3] = 0;
  (*(EigenTensors))[0][4] = -1.0/sqrt(6.0);
  (*(EigenTensors))[0][5] = 0;
  (*(EigenTensors))[0][6] = 0;
  (*(EigenTensors))[0][7] = 0;
  (*(EigenTensors))[0][8] = 2.0/sqrt(6.0);
	
  // b2
  (*(EigenTensors))[1][0] = -1.0/sqrt(2.0);
  (*(EigenTensors))[1][1] = 0;
  (*(EigenTensors))[1][2] = 0;
  (*(EigenTensors))[1][3] = 0;
  (*(EigenTensors))[1][4] = 1.0/sqrt(2.0);
  (*(EigenTensors))[1][5] = 0;
  (*(EigenTensors))[1][6] = 0;
  (*(EigenTensors))[1][7] = 0;
  (*(EigenTensors))[1][8] = 0;
	
  // b3
  (*(EigenTensors))[2][0] = 0;
  (*(EigenTensors))[2][1] = 0;
  (*(EigenTensors))[2][2] = 0;
  (*(EigenTensors))[2][3] = 0;
  (*(EigenTensors))[2][4] = 0;
  (*(EigenTensors))[2][5] = 1.0/sqrt(2.0);
  (*(EigenTensors))[2][6] = 0;
  (*(EigenTensors))[2][7] = 1.0/sqrt(2.0);
  (*(EigenTensors))[2][8] = 0;
	
  // b4
  (*(EigenTensors))[3][0] = 0;
  (*(EigenTensors))[3][1] = 0;
  (*(EigenTensors))[3][2] = 1.0/sqrt(2.0);
  (*(EigenTensors))[3][3] = 0;
  (*(EigenTensors))[3][4] = 0;
  (*(EigenTensors))[3][5] = 0;
  (*(EigenTensors))[3][6] = 1.0/sqrt(2.0);
  (*(EigenTensors))[3][7] = 0;
  (*(EigenTensors))[3][8] = 0;
	
  // b5
  (*(EigenTensors))[4][0] = 0;
  (*(EigenTensors))[4][1] = 1.0/sqrt(2.0);
  (*(EigenTensors))[4][2] = 0;
  (*(EigenTensors))[4][3] = 1.0/sqrt(2.0);
  (*(EigenTensors))[4][4] = 0;
  (*(EigenTensors))[4][5] = 0;
  (*(EigenTensors))[4][6] = 0;
  (*(EigenTensors))[4][7] = 0;
  (*(EigenTensors))[4][8] = 0;
	
  // b6
  (*(EigenTensors))[5][0] = 1.0/sqrt(3.0);
  (*(EigenTensors))[5][1] = 0;
  (*(EigenTensors))[5][2] = 0;
  (*(EigenTensors))[5][3] = 0;
  (*(EigenTensors))[5][4] = 1.0/sqrt(3.0);
  (*(EigenTensors))[5][5] = 0;
  (*(EigenTensors))[5][6] = 0;
  (*(EigenTensors))[5][7] = 0;
  (*(EigenTensors))[5][8] = 1.0/sqrt(3.0);

}

class OrthorhombicKelvinBase : public KelvinBase{
public:
  void constructBase();
  OrthorhombicKelvinBase(double* ElasticConstants):KelvinBase(ElasticConstants) {
    constructBase();
  };
  ~OrthorhombicKelvinBase() {
  };
};

void OrthorhombicKelvinBase::constructBase() {
  double* ReducedElasticConstants = SOTensor_tmp;
  for (int i = 0; i<3; i++) {
    for (int j = 0; j<3; j++) {
      ReducedElasticConstants[i*3+j] = ElasticConstants[i*6+j];
    }
  }
// First three eigenvalues are eigenvalues of (C11, C12, C13; C21, C22, C23; C31, C32, C33)
  EigenValues3x3(ReducedElasticConstants, EigenValues);
// Last three eigenvalues 2 * C44, 2 * C55, 2 * C66
  EigenValues[3] = 2 * ElasticConstants[21];
  EigenValues[4] = 2 * ElasticConstants[28];
  EigenValues[5] = 2 * ElasticConstants[35];

  double lambda1 = EigenValues[0],
    lambda2 = EigenValues[1],
    lambda3 = EigenValues[2];
  double Ktilde1, Ktilde2, Ktilde3;
  
  Ktilde1 = 1 / ( (lambda1 - lambda3) * (lambda1 - lambda2) * 
		  ( ( ElasticConstants[1] * ( pow(ElasticConstants[2], 2) - pow(ElasticConstants[8], 2) ) ) - 
		    ( ElasticConstants[2] * ElasticConstants[8] * ( ElasticConstants[0] - ElasticConstants[7] ) ) ) );

  Ktilde2 = 1 / ( (lambda2 - lambda1) * (lambda2 - lambda3) * 
		  ( ( ElasticConstants[1] * ( pow(ElasticConstants[2], 2) - pow(ElasticConstants[8], 2) ) ) - 
		    ( ElasticConstants[2] * ElasticConstants[8] * ( ElasticConstants[0] - ElasticConstants[7] ) ) ) );
		
  Ktilde3 = 1 / ( (lambda3 - lambda2) * (lambda3 - lambda1) * 
		  ( ( ElasticConstants[1] * ( pow(ElasticConstants[2], 2) - pow(ElasticConstants[8], 2) ) ) - 
		    ( ElasticConstants[2] * ElasticConstants[8] * ( ElasticConstants[0] - ElasticConstants[7] ) ) ) );	

  double E11, E22, E33;		

  //Tensor 1
	
  E11 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ( ElasticConstants[0] - lambda3 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ( ElasticConstants[7] - lambda3 ) ) );

  E33 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[8] ) ); 

  (*(EigenTensors))[0][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow(E11,2) + pow(E22,2) + pow(E33,2) );
  (*(EigenTensors))[0][1] = 0;
  (*(EigenTensors))[0][2] = 0;	
	
  E11 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ( ElasticConstants[0] - lambda3 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[1] ) );
		
  E22 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ( ElasticConstants[7] - lambda3 ) ) );	
		
  E33 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda1 ) ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[0][3] = 0;
  (*(EigenTensors))[0][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow(E11,2) + pow(E22,2) +	pow(E33,2) ); 
  (*(EigenTensors))[0][5] = 0;

  E11 = ( ( ( ElasticConstants[0] - lambda1 ) * ( ElasticConstants[7] - lambda1 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ( ElasticConstants[0] - lambda3 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ( ElasticConstants[0] - lambda1 ) * ( ElasticConstants[7] - lambda1 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ( ElasticConstants[7] - lambda3 ) ) );

  E33 = ( ( ( ElasticConstants[0] - lambda1 ) * ( ElasticConstants[7] - lambda1 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde1 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda2 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[0][6] = 0;
  (*(EigenTensors))[0][7] = 0;
  (*(EigenTensors))[0][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
	
  //Tensor 2
	
  E11 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ( ElasticConstants[0] - lambda1 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[1] ) );
	
  E22 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ( ElasticConstants[7] - lambda1 ) ) );
	
  E33 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[1][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow(E11,2) + pow(E22,2) + pow(E33,2) );
  (*(EigenTensors))[1][1] = 0;
  (*(EigenTensors))[1][2] = 0;

  E11 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ( ElasticConstants[0] - lambda1 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ( ElasticConstants[7] - lambda1 ) ) );
	
  E33 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda2 ) ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[1][3] =  0;
  (*(EigenTensors))[1][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow(E11,2) + pow(E22,2) +	pow(E33,2) );
  (*(EigenTensors))[1][5] = 0;

  E11 = ( ( ( ElasticConstants[0] - lambda2 ) * ( ElasticConstants[7] - lambda2 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ( ElasticConstants[0] - lambda1 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ( ElasticConstants[0] - lambda2 ) * ( ElasticConstants[7] - lambda2 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ( ElasticConstants[7] - lambda1 ) ) );

  E33 = ( ( ( ElasticConstants[0] - lambda2 ) * ( ElasticConstants[7] - lambda2 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde2 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda3 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[1][6] = 0;
  (*(EigenTensors))[1][7] = 0;
  (*(EigenTensors))[1][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );

  //Tensor 3
	
  E11 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ( ElasticConstants[0] - lambda2 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[1] ) );
	
  E22 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ( ElasticConstants[7] - lambda2 ) ) );
	
  E33 = ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[8] ) );
	
  (*(EigenTensors))[2][0] =  ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[2][1] = 0;
  (*(EigenTensors))[2][2] = 0;
						
  E11 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ( ElasticConstants[0] - lambda2 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ( ElasticConstants[7] - lambda2 ) ) );

  E33 = ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] * ( ElasticConstants[0] - lambda3 ) ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[2][3] = 0;
  (*(EigenTensors))[2][4] =  ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[2][5] = 0;
					
  E11 = ( ( ( ElasticConstants[0] - lambda3 ) * ( ElasticConstants[7] - lambda3 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ( ElasticConstants[0] - lambda2 ) )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[1] ) );

  E22 = ( ( ( ElasticConstants[0] - lambda3 ) * ( ElasticConstants[7] - lambda3 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[1] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ( ElasticConstants[7] - lambda2 ) ) );

  E33 = ( ( ( ElasticConstants[0] - lambda3 ) * ( ElasticConstants[7] - lambda3 ) ) - pow( ElasticConstants[1], 2 ) )
    * Ktilde3 * ( ( ( ( ElasticConstants[1] * ElasticConstants[2] ) - ( ElasticConstants[8] *  ( ElasticConstants[0] - lambda1 ) ) ) * ElasticConstants[2] )
		  - ( ( ( ElasticConstants[1] * ElasticConstants[8] ) - ( ElasticConstants[2] * ( ElasticConstants[7] - lambda1 ) ) ) * ElasticConstants[8] ) );

  (*(EigenTensors))[2][6] = 0;
  (*(EigenTensors))[2][7] = 0;
  (*(EigenTensors))[2][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
	
  //Tensor 4
  (*(EigenTensors))[3][0] = 0;
  (*(EigenTensors))[3][1] = 1/sqrt(2);
  (*(EigenTensors))[3][2] = 0;
  (*(EigenTensors))[3][3] = 1/sqrt(2);
  (*(EigenTensors))[3][4] = 0;
  (*(EigenTensors))[3][5] = 0;
  (*(EigenTensors))[3][6] = 0;
  (*(EigenTensors))[3][7] = 0;
  (*(EigenTensors))[3][8] = 0;
	
  //Tensor 5
  (*(EigenTensors))[4][0] = 0;
  (*(EigenTensors))[4][1] = 0;
  (*(EigenTensors))[4][2] = 1/sqrt(2);
  (*(EigenTensors))[4][3] = 0;
  (*(EigenTensors))[4][4] = 0;
  (*(EigenTensors))[4][5] = 0;
  (*(EigenTensors))[4][6] = 1/sqrt(2);
  (*(EigenTensors))[4][7] = 0;
  (*(EigenTensors))[4][8] = 0; 
	
  //Tensor 6
  (*(EigenTensors))[5][0] = 0;
  (*(EigenTensors))[5][1] = 0;
  (*(EigenTensors))[5][2] = 0;
  (*(EigenTensors))[5][3] = 0;
  (*(EigenTensors))[5][4] = 0;
  (*(EigenTensors))[5][5] = 1/sqrt(2);
  (*(EigenTensors))[5][6] = 0;
  (*(EigenTensors))[5][7] = 1/sqrt(2);
  (*(EigenTensors))[5][8] = 0;

}

class HexagonalKelvinBase : public KelvinBase{
public:
  void constructBase();
  HexagonalKelvinBase(double* ElasticConstants):KelvinBase(ElasticConstants) {
    constructBase();
  };
  ~HexagonalKelvinBase() {
  };
};

void HexagonalKelvinBase::constructBase() {
  double alpha = atan( sqrt(2.) * ( ElasticConstants[0] + ElasticConstants[1] - ElasticConstants[14] ) / ( 4. * ElasticConstants[2] ) );
// First two eigenvalues associated with dilatational strain C33+-sqrt(2)C13(+-tan(alpha) + sec(alpha))
  EigenValues[0] = ElasticConstants[14] + ( sqrt(2.) * ElasticConstants[2] * ( tan(alpha) + ( 1. / cos(alpha) ) ) );
  EigenValues[1] = ElasticConstants[14] - ( sqrt(2.) * ElasticConstants[2] * ( - tan(alpha) + ( 1. / cos(alpha) ) ) );
// Last eigenvalues of multiplicity 2 : C11 - C12 and 2 * C44
  EigenValues[2] = (ElasticConstants[0] - ElasticConstants[1]);
  EigenValues[3] = 2 * ElasticConstants[21];
  EigenValues[4] = (ElasticConstants[0] - ElasticConstants[1]);
  EigenValues[5] = 2 * ElasticConstants[21];

  double E11, E22, E33;
  //Tensor 1
  E11 = 0.25 * ( 1 + sin(alpha) );
  E22 = 0.25 * ( 1 + sin(alpha) );
  E33 = sqrt(2.) * cos(alpha) / 4.;
  (*(EigenTensors))[0][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[0][1] = 0;
  (*(EigenTensors))[0][2] = 0;
  (*(EigenTensors))[0][3] = 0;
  (*(EigenTensors))[0][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[0][5] = 0;
  (*(EigenTensors))[0][6] = 0;
  (*(EigenTensors))[0][7] = 0;
  E11 = sqrt(2.) * cos(alpha) / 4.;
  E22 = sqrt(2.) * cos(alpha) / 4.;
  E33 = 0.5 * (1 - sin(alpha));
  (*(EigenTensors))[0][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  //Tensor 2
  E11 = 0.25 * ( 1 - sin(alpha) );
  E22 = 0.25 * ( 1 - sin(alpha) );
  E33 = - sqrt(2.) * cos(alpha) / 4.;
  (*(EigenTensors))[1][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[1][1] = 0;
  (*(EigenTensors))[1][2] = 0;
  (*(EigenTensors))[1][3] = 0;
  (*(EigenTensors))[1][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[1][5] = 0;
  (*(EigenTensors))[1][6] = 0;
  (*(EigenTensors))[1][7] = 0;
  E11 = - sqrt(2.) * cos(alpha) / 4.;
  E22 = - sqrt(2.) * cos(alpha) / 4.;
  E33 = 0.5 * (1 + sin(alpha));
  (*(EigenTensors))[1][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );

  //Tensor 3
  (*(EigenTensors))[2][0] = - 1. / sqrt(2);
  (*(EigenTensors))[2][1] = 0.;
  (*(EigenTensors))[2][2] = 0.;
  (*(EigenTensors))[2][3] = 0;
  (*(EigenTensors))[2][4] = 1. / sqrt(2);
  (*(EigenTensors))[2][5] = 0;
  (*(EigenTensors))[2][6] = 0;
  (*(EigenTensors))[2][7] = 0;
  (*(EigenTensors))[2][8] = 0;
	
  //Tensor 4
  (*(EigenTensors))[3][0] = 0;
  (*(EigenTensors))[3][1] = 1. / sqrt(2);
  (*(EigenTensors))[3][2] = 0;
  (*(EigenTensors))[3][3] = 1. / sqrt(2);
  (*(EigenTensors))[3][4] = 0;
  (*(EigenTensors))[3][5] = 0;
  (*(EigenTensors))[3][6] = 0;
  (*(EigenTensors))[3][7] = 0;
  (*(EigenTensors))[3][8] = 0;

  //Tensor 5
  (*(EigenTensors))[4][0] = 0;
  (*(EigenTensors))[4][1] = 0;
  (*(EigenTensors))[4][2] = 0;
  (*(EigenTensors))[4][3] = 0;
  (*(EigenTensors))[4][4] = 0;
  (*(EigenTensors))[4][5] = 1. / sqrt(2.);
  (*(EigenTensors))[4][6] = 0.;
  (*(EigenTensors))[4][7] = 1. / sqrt(2.);
  (*(EigenTensors))[4][8] = 0; 
	
  //Tensor 6
  (*(EigenTensors))[5][0] = 0;
  (*(EigenTensors))[5][1] = 0;
  (*(EigenTensors))[5][2] = 1. / sqrt(2.);
  (*(EigenTensors))[5][3] = 0;
  (*(EigenTensors))[5][4] = 0;
  (*(EigenTensors))[5][5] = 0;
  (*(EigenTensors))[5][6] = 1. / sqrt(2.);
  (*(EigenTensors))[5][7] = 0.;
  (*(EigenTensors))[5][8] = 0;

}

class TetragonalKelvinBase : public KelvinBase{
public:
  void constructBase();
  TetragonalKelvinBase(double* ElasticConstants):KelvinBase(ElasticConstants) {
    constructBase();
  };
  ~TetragonalKelvinBase() {
  };
};

void TetragonalKelvinBase::constructBase() {
  double alpha = atan( sqrt(2.) * ( ElasticConstants[0] + ElasticConstants[1] - ElasticConstants[14] ) / ( 4. * ElasticConstants[2] ) );
// First two eigenvalues associated with dilatational strain C33+-sqrt(2)C13(+-tan(alpha) + sec(alpha))
  EigenValues[0] = ElasticConstants[14] + ( sqrt(2.) * ElasticConstants[2] * ( tan(alpha) + ( 1. / cos(alpha) ) ) );
  EigenValues[1] = ElasticConstants[14] - ( sqrt(2.) * ElasticConstants[2] * ( - tan(alpha) + ( 1. / cos(alpha) ) ) );
// Third eigenvalue multiplicity 1 = C11 - C12
  EigenValues[2] = (ElasticConstants[0] - ElasticConstants[1]);
// Last three eigenvalues multiplicity = 2 * C44
  EigenValues[3] = 2 * ElasticConstants[35];
  EigenValues[4] = 2 * ElasticConstants[21];
  EigenValues[5] = 2 * ElasticConstants[21];

  double E11, E22, E33;
  //Tensor 1
  E11 = 0.25 * ( 1 + sin(alpha) );
  E22 = 0.25 * ( 1 + sin(alpha) );
  E33 = sqrt(2.) * cos(alpha) / 4.;
  (*(EigenTensors))[0][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[0][1] = 0;
  (*(EigenTensors))[0][2] = 0;
  (*(EigenTensors))[0][3] = 0;
  (*(EigenTensors))[0][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[0][5] = 0;
  (*(EigenTensors))[0][6] = 0;
  (*(EigenTensors))[0][7] = 0;
  E11 = sqrt(2.) * cos(alpha) / 4.;
  E22 = sqrt(2.) * cos(alpha) / 4.;
  E33 = 0.5 * (1 - sin(alpha));
  (*(EigenTensors))[0][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  //Tensor 2
  E11 = 0.25 * ( 1 - sin(alpha) );
  E22 = 0.25 * ( 1 - sin(alpha) );
  E33 = - sqrt(2.) * cos(alpha) / 4.;
  (*(EigenTensors))[1][0] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[1][1] = 0;
  (*(EigenTensors))[1][2] = 0;
  (*(EigenTensors))[1][3] = 0;
  (*(EigenTensors))[1][4] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );
  (*(EigenTensors))[1][5] = 0;
  (*(EigenTensors))[1][6] = 0;
  (*(EigenTensors))[1][7] = 0;
  E11 = - sqrt(2.) * cos(alpha) / 4.;
  E22 = - sqrt(2.) * cos(alpha) / 4.;
  E33 = 0.5 * (1 + sin(alpha));
  (*(EigenTensors))[1][8] = ( -1 * (((E11<0)&(E22<0))||((E22<0)&(E33<0))||((E11<0)&(E33<0))) + 1 * (((E11>0)&(E22>0))||((E22>0)&(E33>0))||((E11>0)&(E33>0))) ) 
    * sqrt( pow (E11,2) + pow (E22,2) + pow (E33,2) );

  //Tensor 3
  (*(EigenTensors))[2][0] = - 1. / sqrt(2);
  (*(EigenTensors))[2][1] = 0.;
  (*(EigenTensors))[2][2] = 0.;
  (*(EigenTensors))[2][3] = 0;
  (*(EigenTensors))[2][4] = 1. / sqrt(2);
  (*(EigenTensors))[2][5] = 0;
  (*(EigenTensors))[2][6] = 0;
  (*(EigenTensors))[2][7] = 0;
  (*(EigenTensors))[2][8] = 0;
	
  //Tensor 4
  (*(EigenTensors))[3][0] = 0;
  (*(EigenTensors))[3][1] = 1. / sqrt(2);
  (*(EigenTensors))[3][2] = 0;
  (*(EigenTensors))[3][3] = 1. / sqrt(2);
  (*(EigenTensors))[3][4] = 0;
  (*(EigenTensors))[3][5] = 0;
  (*(EigenTensors))[3][6] = 0;
  (*(EigenTensors))[3][7] = 0;
  (*(EigenTensors))[3][8] = 0;

  //Tensor 5
  (*(EigenTensors))[4][0] = 0;
  (*(EigenTensors))[4][1] = 0;
  (*(EigenTensors))[4][2] = 0;
  (*(EigenTensors))[4][3] = 0;
  (*(EigenTensors))[4][4] = 0;
  (*(EigenTensors))[4][5] = 1. / sqrt(2.);
  (*(EigenTensors))[4][6] = 0.;
  (*(EigenTensors))[4][7] = 1. / sqrt(2.);
  (*(EigenTensors))[4][8] = 0; 
	
  //Tensor 6
  (*(EigenTensors))[5][0] = 0;
  (*(EigenTensors))[5][1] = 0;
  (*(EigenTensors))[5][2] = 1. / sqrt(2.);
  (*(EigenTensors))[5][3] = 0;
  (*(EigenTensors))[5][4] = 0;
  (*(EigenTensors))[5][5] = 0;
  (*(EigenTensors))[5][6] = 1. / sqrt(2.);
  (*(EigenTensors))[5][7] = 0.;
  (*(EigenTensors))[5][8] = 0;


}

// Compute the decompostion of an elastic constant tensor along the different symetries
// Monoclinic, orthorhombic, tetragonal, hexagonal, isotropic
// Formalism based on J.T. Browaeys and S. Chevrot Geophys. J. Int (2004)
// all matrix for projection are defined in the article
// Jean Furstoss 26/09/19

class Projector{
protected:
  double* matrix;
public:
  Projector() { matrix = new double[21 * 21]; }
  ~Projector() { delete[] matrix; }
  virtual void ConstructMatrixProjection() =0;
  void project(const double* ToProject, double* Projected);
};

void Projector::project(const double* ToProject, double* Projected) {
  for (int i=0; i<21; i++) {
    Projected[i] = 0;
    for (int j=0; j<21; j++) {
      Projected[i] += matrix[i*21 + j] * ToProject[j];
    }
  }
}

class MonoclinicProjection : public Projector{
public:
  MonoclinicProjection():Projector() {	ConstructMatrixProjection(); }
  void ConstructMatrixProjection();
};

void MonoclinicProjection::ConstructMatrixProjection() {
  for (int i=0; i<21; i++) {
    for (int j=0; j<21; j++) {
      if( i==10 || i == 11 || i==13 || i==14 || i==16
	  || i==17 || i==19 || i==20 || i!=j ) 
	matrix[i*21 +j] = 0;
      else 
	matrix[i*21 +j] = 1;
    }
  }
}

class OrthorhombicProjection : public Projector{
public:
  OrthorhombicProjection():Projector() {	ConstructMatrixProjection(); }
  void ConstructMatrixProjection();
};

void OrthorhombicProjection::ConstructMatrixProjection() {
  for (int i=0; i<21; i++) {
    for (int j=0; j<21; j++) {
      if( i>=10 || i!=j ) 
	matrix[i*21 +j] = 0;
      else 
	matrix[i*21 +j] = 1;
    }
  }
}

class TetragonalProjection : public Projector{
public:
  TetragonalProjection():Projector() {	ConstructMatrixProjection(); }
  void ConstructMatrixProjection();
};

void TetragonalProjection::ConstructMatrixProjection() {
  for (int i=0; i<21; i++) {
    for (int j=0; j<21; j++) {
      matrix[i*21 +j] = 0;
    }
  }
  matrix[0*21 + 0] = 0.5;
  matrix[0*21 + 1] = 0.5;
  matrix[1*21 + 0] = 0.5;
  matrix[1*21 + 1] = 0.5;
  matrix[2*21 + 2] = 1;
  matrix[3*21 + 3] = 0.5;
  matrix[3*21 + 4] = 0.5;
  matrix[4*21 + 3] = 0.5;
  matrix[4*21 + 4] = 0.5;
  matrix[5*21 + 5] = 1;
  matrix[6*21 + 6] = 0.5;
  matrix[6*21 + 7] = 0.5;
  matrix[7*21 + 7] = 0.5;
  matrix[7*21 + 6] = 0.5;
  matrix[8*21 + 8] = 1;
}

class HexagonalProjection : public Projector{
public:
  HexagonalProjection():Projector() {	ConstructMatrixProjection(); }
  void ConstructMatrixProjection();
};

void HexagonalProjection::ConstructMatrixProjection() {
  for (int i=0; i<21; i++) {
    for (int j=0; j<21; j++) {
      matrix[i*21 +j] = 0;
    }
  }
  matrix[0*21 + 0] = 0.375;
  matrix[0*21 + 1] = 0.375;
  matrix[1*21 + 1] = 0.375;
  matrix[1*21 + 0] = 0.375;
  matrix[0*21 + 5] = 1./(4* sqrt(2));
  matrix[1*21 + 5] = 1./(4* sqrt(2));
  matrix[0*21 + 8] = 1./4;
  matrix[1*21 + 8] = 1./4;
  matrix[2*21 + 2] = 1;
  matrix[3*21 + 3] = 0.5;
  matrix[3*21 + 4] = 0.5;
  matrix[4*21 + 3] = 0.5;
  matrix[4*21 + 4] = 0.5;
  matrix[5*21 + 0] = 1./(4* sqrt(2));
  matrix[5*21 + 1] = 1./(4* sqrt(2));
  matrix[5*21 + 5] = 3./4;
  matrix[5*21 + 8] = - 1./ (2* sqrt(2));
  matrix[6*21 + 6] = 0.5;
  matrix[6*21 + 7] = 0.5;
  matrix[7*21 + 6] = 0.5;
  matrix[7*21 + 7] = 0.5;
  matrix[8*21 + 0] = 0.25;
  matrix[8*21 + 1] = 0.25;
  matrix[8*21 + 5] = - 1./(2* sqrt(2));
  matrix[8*21 + 8] = 0.5;
}

class IsotropicProjection : public Projector{
public:
  IsotropicProjection():Projector() {	ConstructMatrixProjection(); }
  void ConstructMatrixProjection();
};

void IsotropicProjection::ConstructMatrixProjection() {
  for (int i=0; i<21; i++) {
    for (int j=0; j<21; j++) {
      matrix[i*21 +j] = 0;
    }
  }
  for (int i=0; i<3; i++) {
    for (int j=0; j<3; j++) {
      matrix[i*21 +j] = 3./15;
    }
  }
  for (int i=0; i<3; i++) {
    for (int j=3; j<6; j++) {
      matrix[i*21 +j] = sqrt(2) /15;
      matrix[j*21 +i] = sqrt(2) /15;
    }
  }
  for (int i=0; i<3; i++) {
    for (int j=6; j<9; j++) {
      matrix[i*21 +j] = 2. /15;
      matrix[j*21 +i] = 2. /15;
    }
  }
  for (int i=3; i<6; i++) {
    for (int j=3; j<6; j++) {
      matrix[i*21 +j] = 4. /15;
    }
  }
  for (int i=3; i<6; i++) {
    for (int j=6; j<9; j++) {
      matrix[i*21 +j] = -sqrt(2) /15;
      matrix[j*21 +i] = -sqrt(2) /15;
    }
  }
  for (int i=6; i<9; i++) {
    for (int j=6; j<9; j++) {
      matrix[i*21 +j] = 1. /5;
    }
  }
}

class ElasticConstantsProjection{
private:
  IsotropicProjection iso;
  HexagonalProjection hexa;
  TetragonalProjection tetra;
  OrthorhombicProjection ortho;
  MonoclinicProjection mono;
  double NormRelIso, NormRelHexa, NormRelTetra, NormRelOrtho, NormRelMono;
protected:
  double* Vector21;
  double* Vector21Projected;
  double* CVoigtProjected;
  double Norm;
public:
  double getNormIso() { return NormRelIso; }
  double getNormHexa() { return NormRelHexa; }
  double getNormTetra() { return NormRelTetra; }
  double getNormOrtho() { return NormRelOrtho; }
  double getNormMono() { return NormRelMono; }
  ElasticConstantsProjection(const double* VoigtConstants) {
    Vector21 = new double[21];
    Vector21Projected = new double[21];
    CVoigtProjected = new double[36];
    Voigt2Vector21(VoigtConstants);
    decomposeTensor();
  }
  ~ElasticConstantsProjection() { 
    delete[] Vector21;
    delete[] Vector21Projected;
    delete[] CVoigtProjected;
  }
  void Voigt2Vector21(const double* CVoigt);
  void Vector212Voigt(const double* _vector21);
  double computeNorm(const double* vector21);
  double* getElasticConstants() { return Vector21; }
  void ProjectOnIso() { 
    iso.project(Vector21, Vector21Projected);
    Vector212Voigt(Vector21Projected);
  }
  void ProjectOnHexa() { 
    hexa.project(Vector21, Vector21Projected);
    Vector212Voigt(Vector21Projected);
  }
  void ProjectOnTetra() { 
    tetra.project(Vector21, Vector21Projected);
    Vector212Voigt(Vector21Projected);
  }
  void ProjectOnOrtho() { 
    ortho.project(Vector21, Vector21Projected);
    Vector212Voigt(Vector21Projected);
  }
  void ProjectOnMono() { 
    mono.project(Vector21, Vector21Projected);
    Vector212Voigt(Vector21Projected);
  }
  void decomposeTensor();
  void printProjected() {
    for(int i=0;i<6;i++) {
      for(int j=0;j<6;j++) {
	cout << CVoigtProjected[i*6 + j] << " ";
      }
      cout << endl;
    }
  }
		
};

double ElasticConstantsProjection::computeNorm(const double* vector21) {
  double norm=0;
  for(int i=0;i<21;i++) norm+= vector21[i]*vector21[i];
  return sqrt(norm);
}
	

void ElasticConstantsProjection::decomposeTensor() {
  double normIni = computeNorm(Vector21);
  ProjectOnIso();
  double normIso = computeNorm(Vector21Projected) / normIni;
  cout << "Isotropic Projection : " << endl;
  printProjected();
  for(int i=0;i<21;i++) Vector21[i] -= Vector21Projected[i];
  ProjectOnHexa();
  double normHexa = computeNorm(Vector21Projected) / normIni;
  cout << "Hexagonal Projection : " << endl;
  printProjected();
  for(int i=0;i<21;i++) Vector21[i] -= Vector21Projected[i];
  ProjectOnTetra();
  double normTetra = computeNorm(Vector21Projected) / normIni;
  cout << "Tetragonal Projection : " << endl;
  printProjected();
  for(int i=0;i<21;i++) Vector21[i] -= Vector21Projected[i];
  ProjectOnOrtho();
  double normOrtho = computeNorm(Vector21Projected) / normIni;
  cout << "Orthorhombic Projection : " << endl;
  printProjected();
  for(int i=0;i<21;i++) Vector21[i] -= Vector21Projected[i];
  ProjectOnMono();
  double normMono = computeNorm(Vector21Projected) / normIni;
  cout << "Monoclinic Projection : " << endl;
  printProjected();
  for(int i=0;i<21;i++) Vector21[i] -= Vector21Projected[i];

  double normIni1 = (normIso + normHexa + normTetra + normOrtho + normMono);
  cout << "Isotropic tensor magnitude / Total tensor magnitude : " << normIso /( normIni1 )  << endl;
  NormRelIso = normIso / normIni1;
  cout << "Hexagonal tensor magnitude / Total tensor magnitude : " << normHexa /( normIni1 )  << endl;
  NormRelHexa = normHexa / normIni1;
  cout << "Tetragonal tensor magnitude / Total tensor magnitude : " << normTetra /(  normIni1 )  << endl;
  NormRelTetra = normTetra / normIni1;
  cout << "Orthorhombic tensor magnitude / Total tensor magnitude : " << normOrtho /(  normIni1 )  << endl;
  NormRelOrtho = normOrtho / normIni1;
  cout << "Monoclinic tensor magnitude / Total tensor magnitude : " << normMono /(  normIni1 ) << endl;
  NormRelMono = normMono / normIni1;
}

void ElasticConstantsProjection::Voigt2Vector21(const double* CVoigt) {
  Vector21[0] = CVoigt[0];
  Vector21[1] = CVoigt[7];
  Vector21[2] = CVoigt[14];
  Vector21[3] = sqrt(2) * CVoigt[8];
  Vector21[4] = sqrt(2) * CVoigt[2];
  Vector21[5] = sqrt(2) * CVoigt[1];
  Vector21[6] = 2 * CVoigt[21];
  Vector21[7] = 2 * CVoigt[28];
  Vector21[8] = 2 * CVoigt[35];
  Vector21[9] = 2 * CVoigt[3];
  Vector21[10] = 2 * CVoigt[10];
  Vector21[11] = 2 * CVoigt[17];
  Vector21[12] = 2 * CVoigt[15];
  Vector21[13] = 2 * CVoigt[4];
  Vector21[14] = 2 * CVoigt[11];
  Vector21[15] = 2 * CVoigt[9];
  Vector21[16] = 2 * CVoigt[16];
  Vector21[17] = 2 * CVoigt[5];
  Vector21[18] = ( 2. / sqrt(2) ) * CVoigt[29];
  Vector21[19] = ( 2. / sqrt(2) ) * CVoigt[23];
  Vector21[20] = ( 2. / sqrt(2) ) * CVoigt[22];
}

void ElasticConstantsProjection::Vector212Voigt(const double* _vector21) {
  CVoigtProjected[0] = _vector21[0];
  CVoigtProjected[7] = _vector21[1];
  CVoigtProjected[14] = _vector21[2];
  CVoigtProjected[8] = _vector21[3] / sqrt(2);
  CVoigtProjected[13] = _vector21[3] / sqrt(2);
  CVoigtProjected[2] = _vector21[4] / sqrt(2);
  CVoigtProjected[12] = _vector21[4] / sqrt(2);
  CVoigtProjected[1] = _vector21[5] / sqrt(2);
  CVoigtProjected[6] = _vector21[5] / sqrt(2);
  CVoigtProjected[21] = _vector21[6] / 2;
  CVoigtProjected[28] = _vector21[7] / 2;
  CVoigtProjected[35] = _vector21[8] / 2;
  CVoigtProjected[3] = _vector21[9] / 2;
  CVoigtProjected[18] = _vector21[9] / 2;
  CVoigtProjected[25] = _vector21[10] / 2;
  CVoigtProjected[10] = _vector21[10] / 2;
  CVoigtProjected[17] = _vector21[11] / 2;
  CVoigtProjected[32] = _vector21[11] / 2;
  CVoigtProjected[15] = _vector21[12] / 2;
  CVoigtProjected[20] = _vector21[12] / 2;
  CVoigtProjected[4] = _vector21[13] / 2;
  CVoigtProjected[24] = _vector21[13] / 2;
  CVoigtProjected[11] = _vector21[14] / 2;
  CVoigtProjected[31] = _vector21[14] / 2;
  CVoigtProjected[9] = _vector21[15] / 2;
  CVoigtProjected[19] = _vector21[15] / 2;
  CVoigtProjected[16] = _vector21[16] / 2;
  CVoigtProjected[26] = _vector21[16] / 2;
  CVoigtProjected[5] = _vector21[17] / 2;
  CVoigtProjected[30] = _vector21[17] / 2;
  CVoigtProjected[29] = _vector21[18] / ( 2. / sqrt(2) );
  CVoigtProjected[34] = _vector21[18] / ( 2. / sqrt(2) );
  CVoigtProjected[23] = _vector21[19] / ( 2. / sqrt(2) );
  CVoigtProjected[33] = _vector21[19] / ( 2. / sqrt(2) );
  CVoigtProjected[22] = _vector21[20] / ( 2. / sqrt(2) );
  CVoigtProjected[27] = _vector21[20] / ( 2. / sqrt(2) );
}

int main() {

  // TESTS FOR DIFFERENT SYMMETRIES
  //	double* ElasticConstantsIso = new double[36];
  //	for(int i=0; i<36; i++) ElasticConstantsIso[i] = 0;
  //	
  //	ElasticConstantsIso[0] = 150; //C11
  //	ElasticConstantsIso[1] = 100; //C12
  //	ElasticConstantsIso[2] = 100; //C13
  //	ElasticConstantsIso[6] = 100; //C21
  //	ElasticConstantsIso[7] = 150; //C22
  //	ElasticConstantsIso[8] = 100; //C23
  //	ElasticConstantsIso[12] = 100; //C31
  //	ElasticConstantsIso[13] = 100; //C32
  //	ElasticConstantsIso[14] = 150; //C33
  //	ElasticConstantsIso[21] = 130; //C44
  //	ElasticConstantsIso[28] = 130; //C55
  //	ElasticConstantsIso[35] = 130; //C66
  //	cout << "Isotrop elastic constants in Voigt notation :" << endl;
  //	for(int i = 0; i < 6; i++) {
  //		for(int j = 0; j < 6; j++) {
  //			cout << ElasticConstantsIso[i*6+j] << " ";
  //		}
  //		cout << endl;
  //	}
  //
  //	IsotropicKelvinBase IsoBase(ElasticConstantsIso);
  //	IsoBase.printBase();
  //
  //	double* ElasticConstantsOrtho = new double[36];
  //	for(int i=0; i<36; i++) ElasticConstantsOrtho[i] = 0;
  //	
  //	ElasticConstantsOrtho[0] = 290; //C11
  //	ElasticConstantsOrtho[1] = 55; //C12
  //	ElasticConstantsOrtho[2] = 60; //C13
  //	ElasticConstantsOrtho[6] = 55; //C21
  //	ElasticConstantsOrtho[7] = 170; //C22
  //	ElasticConstantsOrtho[8] = 65; //C23
  //	ElasticConstantsOrtho[12] = 60; //C31
  //	ElasticConstantsOrtho[13] = 65; //C32
  //	ElasticConstantsOrtho[14] = 210; //C33
  //	ElasticConstantsOrtho[21] = 55; //C44
  //	ElasticConstantsOrtho[28] = 70; //C55
  //	ElasticConstantsOrtho[35] = 68; //C66
  //	cout << "Orthorhombic elastic constants in Voigt notation :" << endl;
  //	for(int i = 0; i < 6; i++) {
  //		for(int j = 0; j < 6; j++) {
  //			cout << ElasticConstantsOrtho[i*6+j] << " ";
  //		}
  //		cout << endl;
  //	}
  //
  //	OrthorhombicKelvinBase OrthoBase(ElasticConstantsOrtho);
  //	OrthoBase.printBase();
  //
  //	double* ElasticConstantsHexa = new double[36];
  //	for(int i=0; i<36; i++) ElasticConstantsHexa[i] = 0;
  //	
  //	ElasticConstantsHexa[0] = 26.91; //C11
  //	ElasticConstantsHexa[1] = 9.6; //C12
  //	ElasticConstantsHexa[2] = 6.61; //C13
  //	ElasticConstantsHexa[6] = 9.6; //C21
  //	ElasticConstantsHexa[7] = 26.91; //C22
  //	ElasticConstantsHexa[8] = 6.61; //C23
  //	ElasticConstantsHexa[12] = 6.61; //C31
  //	ElasticConstantsHexa[13] = 6.61; //C32
  //	ElasticConstantsHexa[14] = 23.61; //C33
  //	ElasticConstantsHexa[21] = 6.53; //C44
  //	ElasticConstantsHexa[28] = 6.53; //C55
  //	ElasticConstantsHexa[35] = 8.65; //C66
  //	cout << "Hexagonal elastic constants in Voigt notation :" << endl;
  //	for(int i = 0; i < 6; i++) {
  //		for(int j = 0; j < 6; j++) {
  //			cout << ElasticConstantsHexa[i*6+j] << " ";
  //		}
  //		cout << endl;
  //	}
  //
  //	HexagonalKelvinBase HexaBase(ElasticConstantsHexa);
  //	HexaBase.printBase();
  //	//HexaBase.CheckBase();
  //
  //	double* ElasticConstantsTetra = new double[36];
  //	for(int i=0; i<36; i++) ElasticConstantsTetra[i] = 0;
  //	
  //	ElasticConstantsTetra[0] = 8.391; //C11
  //	ElasticConstantsTetra[1] = 4.87; //C12
  //	ElasticConstantsTetra[2] = 2.81; //C13
  //	ElasticConstantsTetra[6] = 4.87; //C21
  //	ElasticConstantsTetra[7] = 8.391; //C22
  //	ElasticConstantsTetra[8] = 2.81; //C23
  //	ElasticConstantsTetra[12] = 2.81; //C31
  //	ElasticConstantsTetra[13] = 2.81; //C32
  //	ElasticConstantsTetra[14] = 9.665; //C33
  //	ElasticConstantsTetra[21] = 1.754; //C44
  //	ElasticConstantsTetra[28] = 1.754; //C55
  //	ElasticConstantsTetra[35] = 0.7407; //C66
  //	cout << "Tetragonal elastic constants in Voigt notation :" << endl;
  //	for(int i = 0; i < 6; i++) {
  //		for(int j = 0; j < 6; j++) {
  //			cout << ElasticConstantsTetra[i*6+j] << " ";
  //		}
  //		cout << endl;
  //	}
  //
  //	TetragonalKelvinBase TetraBase(ElasticConstantsTetra);
  //	TetraBase.printBase();
  //
  //	delete[] ElasticConstantsIso;
  //	delete[] ElasticConstantsOrtho;
  //	delete[] ElasticConstantsHexa;
  //	delete[] ElasticConstantsTetra;


  double* ElasticConstants = new double[36];
  for(int i=0; i<36; i++) ElasticConstants[i] = 0;

  cout << "Enter the elastic constants of your material in "
       << "Voigt notation" << endl;
  
  cout << "C11 :" << endl;
  cin >> ElasticConstants[0];
  cout << "C22 :" << endl;
  cin >> ElasticConstants[7];
  cout << "C33 :" << endl;
  cin >> ElasticConstants[14];
  cout << "C44 :" << endl;
  cin >> ElasticConstants[21];
  cout << "C55 :" << endl;
  cin >> ElasticConstants[28];
  cout << "C66 :" << endl;
  cin >> ElasticConstants[35];
  cout << "C12 :" << endl;
  cin >> ElasticConstants[1];
  cout << "C13 :" << endl;
  cin >> ElasticConstants[2];
  cout << "C23 :" << endl;
  cin >> ElasticConstants[8];

  ElasticConstants[6] = ElasticConstants[1];
  ElasticConstants[12] = ElasticConstants[2];
  ElasticConstants[13] = ElasticConstants[8];

  cout << "Elastic constants in Voigt notation :" << endl;
  for(int i = 0; i < 6; i++) {
    for(int j = 0; j < 6; j++) {
      cout << ElasticConstants[i*6+j] << " ";
    }
    cout << endl;
  }

  ElasticConstantsProjection Elas(ElasticConstants);

  if(Elas.getNormHexa() < 1e-10) {
    cout << "$$$ The lowest symmetry is isotrop $$$" << endl; 
    cout << "The associated Kelvin base is : " << endl;
    IsotropicKelvinBase IsoBase(ElasticConstants);
    IsoBase.printBase();
  }
  if(Elas.getNormTetra() < 1e-10) {
    cout << "$$$ The lowest symmetry is hexagonal $$$" << endl; 
    cout << "The associated Kelvin base is : " << endl;
    HexagonalKelvinBase HexaBase(ElasticConstants);
    HexaBase.CheckBase();
    HexaBase.printBase();
  }
  if(Elas.getNormOrtho() < 1e-10) {
    cout << "$$$ The lowest symmetry is tetragonal $$$" << endl; 
    cout << "The associated Kelvin base is : " << endl;
    TetragonalKelvinBase TetraBase(ElasticConstants);
    TetraBase.printBase();
  }
  if(Elas.getNormMono() < 1e-10) {
    cout << "$$$ The lowest symmetry is orthorhombic $$$" << endl; 
    cout << "The associated Kelvin base is : " << endl;
    OrthorhombicKelvinBase OrthoBase(ElasticConstants);
    OrthoBase.printBase();
  }
  else cout << "The lowest symmetry is mono(or tri)-clinic for "
	    << "which the Kelvin base is not implemented, or "
	    << "those elastic constants may not represent a "
	    << "crystal symmetry" << endl;

  delete[] ElasticConstants;
  return 0;
}
