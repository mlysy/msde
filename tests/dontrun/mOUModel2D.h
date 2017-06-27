/////////////////////////////////////////

#ifndef sdeModel_h
#define sdeModel_h 1

/////////////////////////////////////////

// mOU model:
//	dYt = (Gamma * Yt + Lambda) * dt + Psi dBt

// sde model object
class sdeModel {
  // put private storage variables here
  // end private storage
 public:
  static const int nParams = 1;
  static const int nDims = 2;
	static const bool sdDiff = true;
	static const bool diagDiff = false;
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
inline sdeModel::sdeModel() {} // do nothing

// destructor
inline sdeModel::~sdeModel() {} // do nothing

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = (-5 * x[0] * theta[0] + -1.378 * x[1] * theta[0] + 1 * theta[0]);
  dr[1] = (-0.122 * x[0] * theta[0] + -2.2 * x[1] * theta[0] + 0.467 * theta[0]);
  return;
}

// diffusion function
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = theta[0] * 1.853;
  df[1] = 0.0;
  df[2] = theta[0] * 0.483;
  df[3] = theta[0] * 1;
  return;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return true;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  return 1;
}

#endif
