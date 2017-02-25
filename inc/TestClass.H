#ifndef _TESTCLASS_H_
#define _TESTCLASS_H_

#include <cstdlib>
#include <iostream>
#include <cstdio>
#include <string>
#include <cmath>
#include <vector>
#include <memory>
#include "Vector.H"



template<class S, int T> 
class Vector2:public Vector<S>{
	public:
		Vector2<S,T>()
		{
			this->v = new S[T];
			this->n = T;
		};
		
};




#endif
