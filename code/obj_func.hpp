#include <cmath>
#include <cstdlib>

#define NUM_IND 200

#define K_NNZ 3

#define LIVE_SIZE 10

#define K_SAC 10

//#define PRINT_RESULTS

#define DEF_TOL (0.0000001)


void get_initials(double *x, int n){
  int i;
  for(i = 0; i < n; i++){
    x[i] = 1.0 - 1.0 / (i + 10.0);
  }
}


template <typename T>
T eval_random_identity(T& rhs) {
  int M = 6;
  int l = rand() % M;
  switch (l) {
    case 0:
      return sqrt(rhs*rhs);
    case 1:
      return 2.0 + rhs - 2.0;
    case 2:
      return rhs * 2.0 * 0.5;
    case 3:
      return log(exp(rhs));
    case 4:
      return 1.0 / (1.0 / rhs);
    case 5:
      // we guarantee that 0 <= rhs < 1, so this is identity
      return sin(asin(rhs));
    default:
      return rhs;
  }
  return rhs;
}

template <typename T>
T eval_func(T * x, int n){
  srand(12345);
  int  i, j, k;
  int r;
  T fad=0;
  T* live = new T[LIVE_SIZE - K_NNZ];
  for(int i = 0; i < n; ++i) {
    T xc = x[i];
    for (int j = 0; j < K_SAC; j++) {
      xc = eval_random_identity(xc);
    }
    T t = 0.0;
    for(int j = 0; j < K_NNZ; j++) {
      int r = rand() % NUM_IND;
      t = t + x[r];
    }
    for (int j = 0; j < LIVE_SIZE - K_NNZ; j++) {
      live[j] = x[i];
      t = t + live[j];    
    }
    fad = fad + xc * t;
  }
  delete[] live;
  return(fad);
}
