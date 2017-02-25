#ifndef _UTILITYFUNCTIONS_H_
#define _UTILITYFUNCTIONS_H_

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

#endif
