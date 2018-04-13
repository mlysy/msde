#ifndef sdeModel_h
#define sdeModel_h 1

//exponential OU model
// d xt = (alpha - exp(2*vt)/2) dt + exp(vt) d Bxt
// d vt = gamma*(mu - vt) dt + sigma d Bzt
// d Bxt d Bzt = rho dt

class sdeModel {
public:
	static const int nParams = 5;
	static const int nDims = 2;
	static const bool sdDiff = true; 
	static const bool diagDiff = false; 
	void sdeDr(double *dr, double *x, double *theta); 
        void sdeDf(double *df, double *x, double *theta); 
        bool isValidParams(double *theta);
        bool isValidData(double *x, double *theta);
};

inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0] - 0.5*exp(x[1]*2);  
  dr[1] = theta[1]*(theta[2] - x[1]);
  return;
}

inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  df[0] = exp(x[1]); 
  df[2] = theta[3];  
  df[3] = sqrt(1.0 - theta[4]*theta[4])*df[2];  
  df[2] *= theta[4];
  return;
}

inline bool sdeModel::isValidParams(double *theta) {
  bool isValid;
  isValid = (theta[1] > 0.2) && (theta[3] > 0.0);
  isValid = isValid && (-1.0 < theta[4]) && (1.0 > theta[4]);
  return(isValid);
}

inline bool sdeModel::isValidData(double *x, double *theta) {
  return(1);
}

#endif
