#ifndef hyperParams_h
#define hyperParams_h 1

// hyper parameters are simply a pointer to an array of double arrays,
// the number of hyper params and the length of each.

class HyperParam {
 public:
  double **val;
  int nArgs;
  int *nEachArg;
  HyperParam(double **val_);
}

#endif
