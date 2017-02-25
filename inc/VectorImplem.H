#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>

//member functions:
/*
template<class T> Vector<T>::Vector();
template<class T> Vector<T>::Vector(const int n1);
template<class T> Vector<T>::Vector(const Vector<T>& A);
template<class T> Vector<T>::~Vector();
template<class T> void Vector<T>::SetDim(const int k);
template<class T> int Vector<T>::GetDim() const;
template<class T> void Vector<T>::Set(const int k, const T vk);
template<class T> void Vector<T>::Set(const T* vk);
template<class T> T Vector<T>::Get(const int k) const;
template<class T> void Vector<T>::Get(T* vk) const;
template<class T> Vector<T> Vector<T>::operator+(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator-(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator*(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator/(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator=(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator*(const T r);
template<class T> Vector<T> Vector<T>::operator/(const T r);
template<class T> T Vector<T>::operator()(const int k) const;
template<class T> T& Vector<T>::operator()(const int k);
template<class T> T Vector<T>::operator[](const int k) const;
template<class T> T& Vector<T>::operator[](const int k);
template<class T> bool Vector<T>::operator==(const Vector<T> &A) const;
template<class T> bool Vector<T>::operator!=(const Vector<T> &A) const;
template<class T> Vector<T> Vector<T>::operator-=(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator+=(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator*=(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator/=(const Vector<T>& A);
template<class T> Vector<T> Vector<T>::operator-();
template<class T> void Vector<T>::Append(const int m);

//non-member functions
template<class T> Vector<T> operator*(double r, Vector<T>& A);
template<class T> T Dot(const Vector<T> &A, const Vector<T> &B);
template<class T> T Max(const Vector<T> &A);
template<class T> T Min(const Vector<T> &A);
template<class T> int ArgMax(const Vector<T> &A);
template<class T> int ArgMin(const Vector<T> &A);
template<class T> Vector<T> Abs(const Vector<T> &A);
template<class T> T Norm(const Vector<T> &A);
template<class T> T ANorm(const Vector<T> &A);
template<class T> T Mean(const Vector<T> &A);
template<class T> T Sum(const Vector<T> &A);
template<class T> Vector<T> Sort(const Vector<T> &A);
void print(const Vector<double> &A);
void print(const Vector<int> &A);
*/




template<class T> Vector<T>::Vector()
{
	n = 0;
	v = 0;
};

template<class T> Vector<T>::Vector(const int n1)
{
	n = n1;
	v = new T[n];
	for (int i=0; i<n; i++)
		v[i] = static_cast<T>(0.);		
};

template<class T> Vector<T>::Vector(const Vector<T>& A)
{
	if (this != &A)
	{
		n = A.GetDim();
		v = new T[n];
		for (int i=0; i<n; i++)
			v[i] = A(i);
	}
};

template<class T> Vector<T>::~Vector()
{
	delete [] v;
	n = 0;
};

template<class T> void Vector<T>::SetDim(const int k)
{
	n=k;
	v = new T[n];
	for (int i =0; i<n; i++)
		v[i] = static_cast<T>(0.);
};

template<class T> int Vector<T>::GetDim() const
{
	return n;
};

template<class T> void Vector<T>::Set(const int k, const T vk)
{
	if (k<n)
		v[k] = vk;
	else
		printf("Error in SET, k>=n\n");
};

template<class T> void Vector<T>::Set(const T* vk)
{
	if (n>0)
	{
		for (int i=0; i<n; i++)
			v[i] = vk[i];
	}
};

template<class T> T Vector<T>::Get(const int k) const
{
	double vk = 0.;
	if (k<n)
		vk = v[k];
	else
		printf("Error in GET, k>=n\n");
	return vk;
};

template<class T> void Vector<T>::Get(T* vk) const
{
	if (n>0)
	{
		for (int i=0; i<n; i++)
			vk[i] = v[i];
	}
};

template<class T> Vector<T> Vector<T>::operator+(const Vector<T>& A)
{
	Vector<T> B(n);
	for (int i=0; i<n; i++)
		B(i) = A(i)+v[i];
	return B;
};

template<class T> Vector<T> Vector<T>::operator-(const Vector<T>& A)
{
	Vector<T> B(n);
	for (int i = 0; i<n; i++)
		B(i) = v[i] - A(i);
	return B;
};

template<class T> Vector<T> Vector<T>::operator*(const Vector<T>& A)
{
	Vector<T> B(n);
	for (int i = 0; i<n; i++)
		B(i) = v[i]*A(i);
	return B;
};

template<class T> Vector<T> Vector<T>::operator*(const T r)
{
	Vector<T> B(n);
	for (int i = 0; i<n; i++)
		B(i) = v[i]*r;
	return B;
};

template<class T> Vector<T> Vector<T>::operator/(const Vector<T>& A)
{
	Vector<T> B(n);
	for (int i = 0; i<n; i++)
		B(i) = v[i]/A(i);
	return B;
};

template<class T> Vector<T> Vector<T>::operator/(const double r)
{
	Vector<T> B(n);
	if (r!=0)
	{
		for (int i = 0; i<n; i++)
			B(i) = v[i]/r;
	}
	return B;
};

template<class T> Vector<T>& Vector<T>::operator=(const Vector<T>& A)
{
	if (this != &A)
	{
		if (n!=A.GetDim())
		{
			delete [] v;
			n = A.GetDim();
			v = new T[n];
		}
		for (int i=0; i<n; i++)
			v[i] = A(i);
	}
	return *this;
};

template<class T> T Vector<T>::operator()(const int k) const
{return v[k];};

template<class T> T& Vector<T>::operator()(const int k)
{return v[k];};

template<class T> T Vector<T>::operator[](const int k) const
{return v[k];};

template<class T> T& Vector<T>::operator[](const int k)
{return v[k];};


template<class T> bool Vector<T>::operator==(const Vector<T> &A) const
{
	bool isequal=true;
	for (int i = 0; i<n; i++)
		if (v[i]!=A(i))
			isequal=false;
	return isequal;
};

template<class T> bool Vector<T>::operator!=(const Vector<T> &A) const
{
	bool isequal=true;
	for (int i = 0; i<n; i++)
		if (v[i]==A(i))
			isequal=false;
	return isequal;
};

template<class T> Vector<T> Vector<T>::operator-=(const Vector<T>& A)
{
	for (int i =0; i<n; i++)
		v[i] = v[i]-A(i);
	return *this;
};

template<class T> Vector<T> Vector<T>::operator+=(const Vector<T>& A)
{
	for (int i =0; i<n; i++)
		v[i] = v[i]+A(i);
	return *this;
};

template<class T> Vector<T> Vector<T>::operator*=(const Vector<T>& A) 
{
	for (int i =0; i<n; i++)
		v[i] = v[i]*A(i);
	return *this;
};

template<class T> Vector<T> Vector<T>::operator/=(const Vector<T>& A) 
{
	for (int i =0; i<n; i++)
		v[i] = v[i]/A(i);
	return *this;
};

template<class T> Vector<T> Vector<T>::operator-()
{
	for (int i =0; i<n; i++)
		v[i] = -v[i];
	return *this;
};

template<class T> void Vector<T>::Append(const int m)
{
	T* v_temp = new T[n+m];
	for (int i=0; i<n+m; i++)
	{	
		if (i<n)
			v_temp[i] = v[i];
		else	
			v_temp[i] = static_cast<T>(0.);		
	}
	n += m;
	delete [] v;
	for (int i=0; i<n; i++)
		v[i] = v_temp[i];
	delete [] v_temp;
};

template<class T> Vector<T> operator*(double r, Vector<T>& A)
{
	Vector<T> B = A*r;
	return B;
}

template<class T> T Dot(const Vector<T> &A, const Vector<T> &B)
{
	T d = static_cast<T>(0.);
	for (int i=0; i<A.GetDim(); i++)
		d+=A(i)*B(i);
	return d;
};

template<class T> T Max(const Vector<T> &A)
{
	int n = A.GetDim();
	T r = static_cast<T>(0.);
	if (n>0)
	{
		r = A(0);
		for (int i = 0; i<n; i++)
			if (A(i)>r)
				r = A(i);
	}
	return r;
};

template<class T> T Min(const Vector<T> &A)
{
	int n = A.GetDim();
	T r = static_cast<T>(0.);
	if (n>0)
	{
		r = A(0);
		for (int i = 0; i<n; i++)
			if (A(i)<r)
				r = A(i);
	}
	return r;
};

template<class T> int ArgMax(const Vector<T> &A)
{
	int n = A.GetDim();
	int r = 0;
	if (n>0)
	{
		for (int i = 1; i<n; i++)
			if (A(i)>A(r))
				r = i;
	}
	return r;
};

template<class T> int ArgMin(const Vector<T> &A)
{
	int n = A.GetDim();
	int r = 0;
	if (n>0)
	{
		for (int i = 1; i<n; i++)
			if (A(i)<A(r))
				r = i;
	}
	return r;
};



template<class T> Vector<T> Abs(const Vector<T> &A)
{
	int n = A.GetDim();
	Vector<T> B(n);
	if (n>0)
	{
		for (int i = 0; i<n; i++)
			B(i) = abs(A(i)); 
	}
	return B;
};

template<class T> T Norm(const Vector<T> &A)
{
	int n = A.GetDim();
	T r = static_cast<T>(0.);
	if (n>0)
	{
		for (int i = 0; i<n; i++)
			r = r+A(i)*A(i);
	}
	return sqrt((double)r);
};

template<class T> T ANorm(const Vector<T> &A)
{
	int n = A.GetDim();	
	T r = static_cast<T>(0.);
	if (n>0)
	{
		for (int i = 0; i<n; i++)
			r = r+abs(A(i));
	}
	return r;
};

template<class T> T Mean(const Vector<T> &A)
{
	int n = A.GetDim();
	T r = static_cast<T>(0.);
	if (n>0)
	{
		for (int i = 0; i<n; i++)
			r = r+A(i)/n;
	}
	return r;
};

template<class T> T Sum(const Vector<T> &A)
{
	int n = A.GetDim();
	T r = static_cast<T>(0.);
	if (n>0)
	{
		for (int i = 0; i<n; i++)
			r = r+A(i);
	}
	return r;
};

template<class T> Vector<T> Sort(const Vector<T> &A)
{
	//insert-sorting algorithm - effective for small vectors, not for large
	int n = A.GetDim();
	Vector<T> B(n);
	B = A;
	
	int j;
	for (int i=1; i<n; i++)
	{
		j=i;
		while (j-1>=0 && B(j)<B(j-1)) 
		{
			std::swap(B(j),B(j-1));
			j--;
		}
	}	
	return B;
};

void print(const Vector<double> &A)
{
	int n = A.GetDim();
	for (int i = 0; i < n; i++)
		printf("%f\n", A(i));
}

void print(const Vector<int> &A)
{
	int n = A.GetDim();
	for (int i = 0; i < n; i++)
		printf("%d\n", A(i));
}

/*
template<class T, class S> bool IsVector(const S& A)
{
	if (typeid(A)==typeid(Vector<T>))
		return 1;
	else 
		return 0;
}


*/
