// ------------------------------------------------------------------------ -
// Projekt JA
// Dominik R¹czka(JA_2016_gr1)
// Projekt: JA_D.Raczka_Rozwiazywanie_URL
// Wersja: 1.1
// ------------------------------------------------------------------------ -

#include "stdafx.h"

extern "C" 
{
	/////////////////////////////////////////////////////////////////////
	//
	//
	//
	//	Function that calculates set of equations using Cholesky method with Gauss elimination
	// 
	//	*************  Parameters:  ***********************
	//	/param size - size of matrixes
	//	/param matL - pointer to first element of matrix L
	//	/param matL - pointer to first element of matrix U
	//	/param matL - pointer to first element of vector Y
	//  /param matL - pointer to first element of vector B
	//  /param matL - pointer to first element of vector X
	//  
	//	
	//	*************  Function Returns:  ******************
	//	> Vecrtor X - vector of result of set of equations. It is not returned by parameter, because
	//  > function operates on original data.
	// 
	//
	//
	//
	///////////////////////////////////////////////////////////////////////
	__declspec(dllexport) void calculateC(int size, double* matL, double* matU, double* vecY, double* vecB, double* vecX)
	{
		//E = sigma; n - size; indexes in code are all reduced by 1 to not waste memory
		vecY[0] = vecB[0];																//vecY[1] = vecB[1]
																						//						i-1		
		for (int i = 1; i < size; i++)													//vecY[i] = vecB[i] - E(matL[i,k]*vecY[k]	i = 2,1,3,...,n
		{																				//						k=1
			double sum = 0.0;

			for (int k = 0; k < i; k++)
				sum += matL[i*size + k] * vecY[k];

			vecY[i] = vecB[i] - sum;
		}

		vecX[size - 1] = vecY[size - 1] / matU[(size - 1) * size + (size - 1)];		//vecX[n] = (vecY[n]/matU[n][n])
																					//						n
		for (int i = size - 2; i >= 0; i--)											//vecX[i] = (vecY[i] - E(matU[i][k]*vecX[k]))/(matU[i][i])	i=n-1,n-2,...,1
		{																			//						k=i+1
			double sum = 0.0;

			for (int k = i + 1; k < size; k++)
				sum += matU[i*size + k] * vecX[k];

			vecX[i] = (vecY[i] - sum) / matU[i*size + i];
		}
	}
}