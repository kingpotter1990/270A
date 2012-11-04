//
//  FIXED COROTATED.h
//  270AHW2
//
//  Created by Franklin Fang on 10/22/12.
//  Copyright (c) 2012 Franklin Fang. All rights reserved.
//

#ifndef _70AHW2_FIXED_COROTATED_h
#define _70AHW2_FIXED_COROTATED_h

#include "DEFORMABLE_OBJECTS.h"

template <class T>
class FIXED_COROTATED_3D:public HYPERELASTICITY_CONSTITUTIVE_MODEL_3D<T>{
    T lambda;
	T mu;
	
public:
	FIXED_COROTATED_3D(const T youngs_modulus,const T poisson_ratio){
		lambda=youngs_modulus*poisson_ratio/(((T)1+poisson_ratio)*((T)1-(T)2*poisson_ratio));
		mu=youngs_modulus/((T)2*((T)1+poisson_ratio));
	}
	
	MATRIX_3X3<T> P(const MATRIX_3X3<T>& F){
        MATRIX_3X3<T> U,V,F_IT;
        VECTOR_2D<T> sig;
        F.SVD(U,sig,V);
        T J = F.Determinant();
        F_IT = F;
        F_IT.Invert().Transpose();
        return 2*mu*F - 2*mu*U*V.Transposed() + lambda*(J - 1)*J*F_IT;
	}
	
	MATRIX_3X3<T> dP(const MATRIX_3X3<T>& F, const MATRIX_3X3<T>& dF){
		MATRIX_3X3<T> delta_JF_IT, F_IT, F_I,delta_R;
        delta_JF_IT(0,0) = dF(1,1)*F(2,2) + F(1,1)*dF(2,2) - (dF(1,2)*F(2,1) + F(1,2)*dF(2,1));
        delta_JF_IT(0,1) = dF(1,2)*F(2,0) + F(1,2)*dF(2,0) - (dF(1,0)*F(2,1) + F(1,0)*dF(2,1));
        delta_JF_IT(0,2) = dF(1,0)*F(2,1) + F(1,0)*dF(2,1) - (dF(1,1)*F(2,0) + F(1,1)*dF(2,0));
        
        delta_JF_IT(1,0) = dF(0,2)*F(2,1) + F(0,2)*dF(2,1) - (dF(0,1)*F(2,2) + F(0,1)*dF(2,2));
        delta_JF_IT(1,1) = dF(0,0)*F(2,2) + F(0,0)*dF(2,2) - (dF(0,2)*F(2,0) + F(0,2)*dF(2,0));
        delta_JF_IT(1,2) = dF(0,1)*F(2,0) + F(0,1)*dF(2,0) - (dF(0,0)*F(2,1) + F(0,0)*dF(2,1));
        
        delta_JF_IT(2,0) = dF(1,2)*F(0,1) + F(1,2)*dF(0,1) - (dF(0,2)*F(1,1) + F(0,2)*dF(1,1));
        delta_JF_IT(2,1) = dF(0,2)*F(1,0) + F(0,2)*dF(1,0) - (dF(0,0)*F(1,2) + F(0,0)*dF(1,2));
        delta_JF_IT(2,2) = dF(0,0)*F(1,1) + F(0,0)*dF(1,1) - (dF(0,1)*F(1,0) + F(0,1)*dF(1,0));
        
        F.Delta_R(dF,delta_R);
        F_IT = F;
        F_I = F_IT.Invert();
        F_IT = F_I.Transpose();
        T J = F.Determinant();
        
        return 2*mu*dF - 2*mu*delta_R + lambda*(F_IT*dF*J*F_I + J*delta_JF_IT) - lambda*delta_JF_IT;
        
	}
	
	void Element_Stifness_Matrix(const int element,const MATRIX_3X3<T>& Dm_Inverse,MATRIX_MXN<T>& element_stiffness){
        
	}

}
;

#endif
