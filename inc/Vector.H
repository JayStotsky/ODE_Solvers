#ifndef _VECTOR_H_
#define _VECTOR_H_
#include "UtilityFunctions.H"

template <class T> class Vector {
	public:
		Vector<T>();
		Vector<T>(int);
		Vector<T>(const Vector&);
		~Vector<T>();
		
		
		 //basic functionality and operator overloading
		 
		void SetDim(const int);
		int GetDim() const;
		void Set(const int, const T);
		T Get(const int) const;
		void Get(T*) const;
		void Set(const T*);
		
		Vector<T> operator+(const Vector<T>&);
		Vector<T> operator-(const Vector<T>&);
		Vector<T> operator*(const Vector<T>&);
		Vector<T> operator*(const T);
		Vector<T> operator/(const Vector<T>&);
		Vector<T> operator/(const double);
		Vector<T>& operator=(const Vector<T>&);
		T operator()(const int) const;
		T& operator()(const int);
		T operator[](const int) const;
		T& operator[](const int);
		
		bool operator==(const Vector<T> &) const;
		bool operator!=(const Vector<T> &) const;
		Vector<T> operator-=(const Vector<T>&);
		Vector<T> operator+=(const Vector<T>&);
		Vector<T> operator*=(const Vector<T>&);
		Vector<T> operator/=(const Vector<T>&);
		Vector<T> operator-();
		
		void Append(const int);

	protected:
		T *v;
		int n=0;


};

#include "VectorImplem.H"

#endif
