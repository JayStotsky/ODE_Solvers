#ifndef _UTILITYFUNCTIONS_H_
#define _UTILITYFUNCTIONS_H_

#include <cmath>

template <typename T> T min(T A, T B)
{
	if (A<B)
		return A;
	else
		return B;
};

template <typename T> T max(T A, T B)
{
	if (A>B)
		return A;
	else
		return B;
};

double Dot(double A, double B)
{
	return A*B;
};

double Norm(double A)
{
	return fabs(A);
};

#endif
