#ifndef _MATRIX_H_
#define _MATRIX_H_
#include "Vector.H"


template <class T> class Matrix {
	public:
		Matrix<T>();
		Matrix<T>(const int);
		Matrix<T>(const int,const int);
		Matrix<T>(const Matrix&);
		~Matrix<T>();
		
		
		 //basic functionality and operator overloading
		 
		void SetSize(const int,const int);
		void GetSize(int&,int&) const;
		
		Matrix<T> operator+(const Matrix<T>&) const;
		Matrix<T> operator-(const Matrix<T>&) const;
		Matrix<T> operator*(const Matrix<T>&) const;
		Matrix<T> operator*(const T) const;
		Vector<T> operator*(const Vector<T>&) const;
		Matrix<T>& operator=(const Matrix<T>&);
		T operator()(const int,const int) const;
		T& operator()(const int,const int);
		
		bool operator==(const Matrix<T> &) const;
		Matrix<T> operator-();

	private:
		T **v;
		int nrows=0;
		int ncols=0;


};

#include "MatrixImplem.H"

#endif
