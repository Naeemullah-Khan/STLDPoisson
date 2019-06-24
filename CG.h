#include <math.h>
#include <stdio.h>
#include <iostream>

using namespace std;

/*
 * Solves the system Ax = b using the conjugate gradient method.  
 * Assumes A is linear, symmetric, positive definite.
 *
 * The matrix A is not stored, instead the class uses a function pointer
 * such that when the function pointed to is called, Ax is computed.
 */

template <class T>
inline int ComputeAx(T *_x, T *_Ax,char *region,int YSize,int XSize,int alpha)
{

		// assuming two components
T Gpx, Gpy, Gmx, Gmy;

for (int x2 = 0, p = 0; x2<YSize; x2++) 
{
	for (int x1 = 0; x1<XSize; x1++, p++) 
    {
		if (region[p]) 
        {
			Gpx = x1 < XSize - 1 && region[p + 1] == region[p] ? _x[p + 1] - _x[p] : 0;
			Gmx = x1 > 0 && region[p - 1] == region[p] ? _x[p - 1] - _x[p] : 0;
			Gpy = x2 < YSize - 1 && region[p + XSize] == region[p] ? _x[p + XSize] - _x[p] : 0;
			Gmy = x2 > 0 && region[p - XSize] == region[p] ? _x[p - XSize] - _x[p] : 0;

			_Ax[p] = _x[p] - alpha*(Gpx + Gmx + Gpy + Gmy);
		}
	}
		
}
  return 1;
}



template <class T>
inline void add(T *a, float scalea, T *b, float scaleb, T*c,int N,char* active)
{
  for (int i=0; i<N; i++)
    if (active[i]) c[i]=a[i]*scalea + b[i]*scaleb;

  return;
}


template <class T>
inline void copy(T *a, const T* b,int N,char* active)
{
  for (int i=0; i<N; i++)
    if (active[i]) a[i]=b[i];

  return;
}


template <class T>
inline float inner(const T* a, const T* b,int N,char* active)
{
  float ret=0;
 
  for (int i=0; i<N; i++) {
	   if (active[i]) ret+=a[i]*b[i];
  }
  
  return ret;
}

template <class T>
inline float error(int N,char* active)
{
  float maxerr=0;
  float err;
  
  for (int i=0; i<N; i++) {
    if (active[i]) {
      err=0;
      maxerr=err > maxerr ? err : maxerr;
    }
  }

  return maxerr;
}


template <class T>
int ComputeCGSolution(float errorTol,T *x,T *b,char *active,int alpha,int width, int height)
{

  int N=width*height;
  int XSize=width;
  int YSize=height;

  T *r, *p, *Ap;
  if ( !(r  = new T[N]) ||
       !(p  = new T[N]) ||
       !(Ap = new T[N])  ) {
    return -1;
  }

  float _alpha, beta;
  T *Ax=r;
  float tol=1;
  
  float bb=inner(b,b,N,active);
  tol=bb>0 ? errorTol*errorTol*bb : tol;
  
  for (int i=0; i<N; i++) if (!active[i]) x[i]=0;
  int numbActive=0;
  for (int i=0; i<N; i++) if (active[i]&&numbActive<100) numbActive++;
  
  ComputeAx(x,Ax,active,YSize,XSize,alpha);                      // Ax = A(x)
  add(Ax, -1, b, 1, r,N,active);          // r = b - Ax
  copy(p, r,N,active);                    // p = r
  float rr=inner(r,r,N,active), rrnew;
  int count=0;

  //printf("rr start=%f\n", rr);
  
  while (rr > tol && count < 2*numbActive) {

    ComputeAx(p, Ap,active,YSize,XSize,alpha);                   // Ap = A(p)
    _alpha=rr/inner(p,Ap,N,active);       // _alpha = r.*r / p .* Ap
    //printf("iter.=%d/%d, err=%f, errorTol=%f, _alpha=%f\n", count, numbActive, rr, tol, _alpha);
    add( x, 1,  p,  _alpha, x,N,active);   // x <= x + _alpha *  p
    add( r, 1, Ap, -_alpha, r,N,active);   // r <= r - _alpha * Ap

    //err=error();

    rrnew=inner(r,r,N,active);
    beta=rrnew/rr;
    add( r, 1, p, beta, p,N,active);      // p <= r + beta*p
    rr=rrnew;

    count++;
  }

  //printf("iter.=%d/%d, err=%f, errorTol=%f\n", count, numbActive, rr, tol);
 
  delete[]  r;  r=0;
  delete[]  p;  p=0;
  delete[] Ap; Ap=0;

  return 1;

}


