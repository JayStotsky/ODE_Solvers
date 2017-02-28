#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>


template<class T> Matrix<T>::Matrix()
{
	nrows = 0;
	ncols = 0;
	v = 0;
};

template<class T> Matrix<T>::Matrix(int n)
{
	nrows = 1;
	ncols = n;
	v = new T*[nrows];
	for (int i =0; i<nrows; i++)
		v[i] = new T[ncols];
		
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			v[i][j] = static_cast<T>(0.);		
}

template<class T> Matrix<T>::Matrix(const int n1, const int n2)
{
	nrows = n1;
	ncols = n2;
	v = new T*[nrows];
	for (int i =0; i<nrows; i++)
		v[i] = new T[ncols];
		
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			v[i][j] = static_cast<T>(0.);		
};

template<class T> Matrix<T>::Matrix(const Matrix<T>& A)
{
	if (this != &A)
	{
		A.GetSize(nrows,ncols);
		v = new T*[nrows];
		for (int i =0; i<nrows; i++)
			v[i] = new T[ncols];
		
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncols; j++)
				v[i][j] = A(i,j);	
	}
};

template<class T> Matrix<T>::~Matrix()
{
	for (int i =0; i<nrows; i++)
		delete [] v[i];
	delete [] v;
	nrows = 0;
	ncols = 0;
};

template<class T> void Matrix<T>::SetSize(const int n1, const int n2)
{
	nrows = n1; 
	ncols = n2;
	
	v = new T*[nrows];
	for (int i =0; i<nrows; i++)
		v[i] = new T[ncols];
		
	for (int i=0; i<nrows; i++)
		for (int j=0; j<ncols; j++)
			v[i][j] = static_cast<T>(0.);		
};

template<class T> void Matrix<T>::GetSize(int& n1, int& n2) const
{
	n1 = nrows;
	n2 = ncols;
};

template<class T> Matrix<T> Matrix<T>::operator+(const Matrix<T>& A) const
{
	int nrowsA, ncolsA;
	A.GetSize(nrowsA,ncolsA);
	
	if (nrowsA==nrows && ncolsA==ncols)
	{
		Matrix<T> B(nrows,ncols);
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncols; j++)
			B(i,j) = A(i,j)+v[i][j];
		return B;
	}
	else 
	{
		printf("Warning, inconsistent dimensions in matrix addition, returning A\n");
		return A;
	}
};

template<class T> Matrix<T> Matrix<T>::operator-(const Matrix<T>& A) const
{
	int nrowsA, ncolsA;
	A.GetSize(nrowsA,ncolsA);
	
	if (nrowsA==nrows && ncolsA==ncols)
	{
		Matrix<T> B(nrows,ncols);
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncols; j++)
			B(i,j) = A(i,j)-v[i][j];
		return B;
	}
	else 
	{
		printf("Warning, inconsistent dimensions in matrix subtraction, returning A\n");
		return A;
	}
};

template<class T> Matrix<T> Matrix<T>::operator*(const Matrix<T>& A) const
{
	int nrowsA, ncolsA;
	A.GetSize(nrowsA,ncolsA);
	if (nrowsA==ncols)
	{
		Matrix<T> B(nrows,ncolsA);
		for (int i=0; i<nrows; i++)
			for (int j=0; j<ncolsA; j++)
				for (int k=0; k<nrowsA; k++)
					B(i,j) += v[i][k]*A(k,j);
		return B;
	}
	else
	{	
		printf("Warning, dimenions of A and B inconsistent in matrix muliplication, returning A\n");	
		return A;
	}
};

template<class T> Matrix<T> Matrix<T>::operator*(const T r) const
{
	Matrix<T> B(nrows,ncols);
	for (int i = 0; i<nrows; i++)
		for (int j = 0; j<ncols; j++)
		B(i,j) = v[i][j]*r;
	return B;
};

template<class T> Vector<T> Matrix<T>::operator*(const Vector<T>& X) const
{
	Vector<T> B(nrows);
	if (X.GetDim() == ncols)
	{
		for (int i = 0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				B(i)+=v[i][j]*X(j);
	}
	else
		printf("Warning, inconsistent dimensions on matrix-vector multiplication. Returning vector of 0's\n");
	return B;
}


template<class T> Matrix<T>& Matrix<T>::operator=(const Matrix<T>& A)
{
	if (this != &A)
	{
		int nrowsA, ncolsA;
		A.GetSize(nrowsA,ncolsA);
		
		if (nrows!=nrowsA || ncols!=ncolsA)
		{
			for (int i = 0; i<nrows; i++)
				delete [] v[i];
			delete [] v;
			
			nrows = nrowsA;
			ncols = ncolsA;
			v = new T*[nrowsA];
			for (int i = 0; i<nrows; i++)
				v[i] = new T[ncolsA];
		}
		for (int i=0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				v[i][j] = A(i,j);
	}
	return *this;
};

template<class T> T Matrix<T>::operator()(const int row, const int col) const
{return v[row][col];};

template<class T> T& Matrix<T>::operator()(const int row, const int col)
{return v[row][col];};

template<class T> bool Matrix<T>::operator==(const Matrix<T> &A) const
{
	bool isequal=true;
	for (int i = 0; i<nrows; i++)
		for (int j = 0; j<ncols; j++)
			if (v[i][j]!=A(i,j))
				isequal=false;
	return isequal;
};

template<class T> Matrix<T> Matrix<T>::operator-()
{
	for (int i =0; i<nrows; i++)
		for (int j =0; j<ncols; j++)
			v[i][j] = -v[i][j];
	return *this;
};

template<class T> Matrix<T> operator*(double r, Matrix<T>& A)
{
	Matrix<T> B = A*r;
	return B;
}

template<class T> Vector<T> operator*(Vector<T>& v, Matrix<T>& A)
{
	Vector<T> b = A*v;
	return b;
}

template<class T> T Dot(const Matrix<T> &A, const Matrix<T> &B)
{
	T d = static_cast<T>(0.);
	int nrowsA, ncolsA, nrowsB, ncolsB;
	A.GetSize(nrowsA,ncolsA);
	B.GetSize(nrowsB,ncolsB);
	
	if (ncolsA==ncolsB && nrowsA==nrowsB)
	{
		for (int i=0; i<nrowsA; i++)
			for (int j = 0; j<ncolsA; j++)
				d+=A(i,j)*B(j,i);
	}
	else
		printf("Warning A and B are not the same size in matrix double dot product, returning 0 as output\n");
	return d;
};

template<class T> T Max(const Matrix<T> &A)
{
	int nrows, ncols;
	A.GetSize(nrows,ncols);
	
	T r = static_cast<T>(0.);
	if (nrows>0)
	{
		r = A(0,0);
		for (int i = 0; i<nrows; i++)
			for(int j = 0;j<ncols; j++)
				if (A(i,j)>r)
					r = A(i,j);
	}
	return r;
};

template<class T> T Min(const Matrix<T> &A)
{
	int nrows, ncols;
	A.GetSize(nrows,ncols);
	
	T r = static_cast<T>(0.);
	if (nrows>0)
	{
		r = A(0,0);
		for (int i = 0; i<nrows; i++)
			for(int j = 0;j<ncols; j++)
				if (A(i,j)<r)
					r = A(i,j);
	}
	return r;
};

template<class T> void ArgMax(const Matrix<T> &A, int& row, int& col)
{
	int nrows, ncols; 
	A.GetSize(nrows,ncols);
	row = 0;
	col = 0;
	if (row>0)
	{
		for (int i=0 ; i<nrows; i++)
			for (int j=0; j<ncols; j++)
				if (A(i,j)>A(row,col))
				{
					row = i;
					col = j;
				}
	}	
};

template<class T> void ArgMin(const Matrix<T> &A, int& row, int& col)
{
	int nrows, ncols; 
	A.GetSize(nrows,ncols);
	row = 0;
	col = 0;
	if (row>0)
	{
		for (int i=0 ; i<nrows; i++)
			for (int j=0; j<ncols; j++)
				if (A(i,j)<A(row,col))
				{
					row = i;
					col = j;
				}
	}	
};


template<class T> Matrix<T> Abs(const Matrix<T> &A)
{
	int nrows, ncols; 
	A.GetSize(nrows,ncols);
	Matrix<T> B(nrows, ncols);
	if (nrows>0)
	{
		for (int i = 0; i<nrows; i++)
			for (int j = 0; j<ncols; j++)
				B(i,j) = abs(A(i,j)); 
	}
	return B;
};

template<class T> T FrobeniusNorm(const Matrix<T> &A)
{
	int nrows, ncols; 
	A.GetSize(nrows,ncols);
	double fn = 0;
	for (int i = 0; i<nrows; i++)
		for (int j = 0; j<ncols; j++)
			fn += A(i,j)*A(i,j);
		return sqrt(fn);
};
