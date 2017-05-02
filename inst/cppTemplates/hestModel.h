/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

// sdeModel class: constains functions unique to a given sde

class sdeModel {
  // put private storage variables here
  // end private storage
 public:
  static const int nParams = 5;
  static const int nDims = 2;
  // TODO: indicate whether the information is on the var or sd scale
  // also whether to use a diagonal var/sd.
  void sdeDr(double *dr, double *x, double *theta);
  void sdeDf(double *df, double *x, double *theta);
  bool isValidData(double *x, double *theta);
  bool isValidParams(double *theta);
  sdeModel();
  ~sdeModel();
};

// constructor
inline sdeModel::sdeModel() {
}

// destructor
inline sdeModel::~sdeModel() {} // do nothing

// TODO: drop t and sqrtT in these functions
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = (theta[0] - .125 * x[1]*x[1]); // x
  dr[1] = (theta[2]/x[1] - .5 * theta[1]*x[1]); // z
  return;
}

inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = .5 * x[1];
  df[2] = theta[3];
  df[3] = sqrt(1.0-theta[4]*theta[4]) * df[2];
  df[2] *= theta[4];
  return;
}

inline bool sdeModel::isValidData(double *x, double *theta) {
  return(x[1] > 0.0);
}

inline bool sdeModel::isValidParams(double *theta) {
  int isValid;
  isValid = (theta[1] > 0) && (theta[3] > 0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  isValid = isValid && (theta[2] > 0.5 * theta[3] * theta[3]);
  return(isValid);
}

#endif
