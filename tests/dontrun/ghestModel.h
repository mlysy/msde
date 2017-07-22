/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

/////////////////////////////////////////

// Generalized Heston model:
// dXt = (alpha - .5 * Vt)dt + Vt^{1/2} dB_Xt
// dVt = -gamma * (Vt - mu)dt + sigma*Vt^lambda * dB_Vt
// cor(B_Xt, B_Vt) = rho

// sde model object
class sdeModel {
 public:
  static const int nParams = 6;
  static const int nDims = 2;
  static const bool sdDiff = true;
  static const bool diagDiff = false;
  void sdeDr(double *dr, double *x, double *theta);
  void sdeDf(double *df, double *x, double *theta);
  bool isValidData(double *x, double *theta);
  bool isValidParams(double *theta);
};


// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0] - .5 * x[1]; // x
  dr[1] = -theta[1] * (x[1] - theta[2]); // v
  return;
}

// diffusion function
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = sqrt(x[1]);
  df[2] = theta[3] * pow(x[1], theta[4]);
  df[3] = sqrt(1.0-theta[5]*theta[5]) * df[2];
  df[2] *= theta[5];
  return;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return(x[1] > 0.0);
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool isValid;
  isValid = (theta[1] > 0.0) && (theta[2] > 0.0) && (theta[3] > 0.0) && (theta[4] > 0.0); // gamma,mu,sigma,lambda > 0
  isValid = isValid && (2*theta[1]*theta[2] > theta[3]*theta[3]); // V_t > 0
  isValid = isValid && (-1.0 < theta[5]) && (1.0 > theta[5]); // -1 < rho < 1
  return(isValid);
}

#endif
