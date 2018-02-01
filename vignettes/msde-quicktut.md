---
title: "Inference for Multivariate Stochastic Differential Equations with **`msde`**"
author: "Martin Lysy, JunYong Tong, Nigel Delaney"
date: "2018-01-31"
output:
  html_vignette:
    toc: true
bibliography: references.bib
csl: taylor-and-francis-harvard-x.csl
link-citations: true
vignette: >
  %\VignetteIndexEntry{Getting started with msde}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

\newcommand{\bm}[1]{\boldsymbol{#1}}
\newcommand{\tth}{\bm{\theta}}
\newcommand{\pphi}{\bm{\phi}}
\newcommand{\eeta}{\bm{\eta}}
\newcommand{\Y}{\bm{Y}}
\newcommand{\y}{\bm{y}}
\newcommand{\B}{\bm{B}}
\renewcommand{\a}{\alpha}
\renewcommand{\b}{\beta}
\newcommand{\g}{\gamma}
\newcommand{\rv}[3][1]{#2_{#1},\ldots,#2_{#3}}
\newcommand{\dr}{\bm{\Lambda}}
\newcommand{\df}{\bm{\Sigma}}
\newcommand{\dt}{\Delta t}
\newcommand{\ud}{\mathrm{d}}
\newcommand{\var}{\mathrm{var}}
\newcommand{\cov}{\mathrm{cov}}
\newcommand{\N}{\mathcal{N}}
\newcommand{\Ym}{\Y_{\mathrm{miss}}}
\newcommand{\iid}{\stackrel{\mathrm{iid}}{\sim}}
\newcommand{\ind}{\stackrel{\mathrm{ind}}{\sim}}

## Introduction

A $d$-dimensional stochastic differential equation (SDE) $\Y_t = (\rv [1t] Y {dt})$ is written as
$$
\ud \Y_t = \dr_{\tth}(\Y_t)\,\ud t + \df_{\tth}(\Y_t)^{1/2}\,\ud \B_t,
$$
where $\dr_{\tth}(\y)$ and $\df_{\tth}(\y)$ are the drift and diffusion functions,
<!-- $$ -->
<!-- \dr_{\tth}(\y) = \lim_{\dt \to 0} \frac{E[\Y_{t+\dt} - \Y_t \mid \Y_t = \y]}{\dt}, \qquad \df_{\tth}(\y) = \lim_{\dt \to 0} \frac{\var(\Y_{t+\dt} - \Y_t \mid \Y_t = \y)}{\dt}, -->
<!-- $$ -->
and $\B_t = (\rv [1t] B {dt})$ is $d$-dimensional Brownian motion.  The **`msde`** package implements a Markov Chain Monte Carlo (MCMC) algorithm to sample from the posterior distribution $p(\tth \mid \Y)$ of the parameters given discrete observations $\Y = (\rv [0] {\Y} N)$ recorded at times $\rv [0] t N$, with some of the $d$ components of $\Y_t$ possibly latent.  To do this efficiently, **`msde`** requires on-the-fly **C++** compiling of user-specified models.  Instructions for setting up **R** to compile **C++** code are provided in the [Installation](#install) section.

## Creating an `sde.model` object

The SDE model used throughout this vignette is the so-called Lotka-Volterra predator-prey model.  Let $H_t$ and $L_t$ denote the number of Hare and Lynx at time $t$ coexisting in a given habitat.  The Lotka-Volterra SDE describing the interactions between these two animal populations is given by [@golightly-wilkinson10]:
$$
\begin{bmatrix} \mathrm{d} H_t \\ \mathrm{d} L_t \end{bmatrix} = \begin{bmatrix} \a H_t - \b H_tL_t \\ \b H_tL_t - \g L_t \end{bmatrix}\, \mathrm{d} t + \begin{bmatrix} \a H_t + \b H_tL_t & -\b H_tL_t \\ -\b H_tL_t & \b H_tL_t + \g L_t\end{bmatrix}^{1/2} \begin{bmatrix} \mathrm{d} B_{1t} \\ \mathrm{d} B_{2t} \end{bmatrix}.
$$
Thus we have $d = 2$, $\Y_t = (H_t, L_t)$, and $\tth = (\a, \b, \g)$.

### The `sdeModel` class definition

In order to build this model in **C++**, we create a header file `lotvolModel.h` containing the class definition for an `sdeModel` object.  The basic structure of this class is given below:

```cpp
// sde model object
class sdeModel {
 public:
  static const int nParams = 3; // number of model parameters
  static const int nDims = 2; // number of sde dimensions
  static const bool sdDiff = false; // whether diffusion function is on sd or var scale
  static const bool diagDiff = false; // whether diffusion function is diagonal
  void sdeDr(double *dr, double *x, double *theta); // drift function
  void sdeDf(double *df, double *x, double *theta); // diffusion function
  bool isValidParams(double *theta); // parameter validator
  bool isValidData(double *x, double *theta); // data validator
};
```
The meaning of each class member is as follows:

* `nParams`: The number of model parameters.  For the Lotka-Volterra model we have $\tth = (\a, \b, \g)$, such that `nParams = 3`.
* `nDims`: The number of dimensions in the multivariate SDE.  In this case we have $\Y_t = (H_t, L_t)$, such that `nDims = 2`.
* `sdeDr`: The SDE drift function.  In **R**, this function would be implemented as
    
    ```r
    sde.drift <- function(x, theta) {
      dr <- c(theta[1]*x[1] - theta[2]*x[1]*x[2], # alpha * H - beta * H*L
          theta[2]*x[1]*x[2] - theta[3]*x[2]) # beta * H*L - gamma * L
      dr
    }
    ```
    In **C++** the same thing is accomplished with
    
    ```cpp
    void sdeDr(double *dr, double *x, double *theta) {
      dr[0] = theta[0]*x[0] - theta[1]*x[0]*x[1]; // alpha * H - beta * H * L
      dr[1] = theta[1]*x[0]*x[1] - theta[2]*x[1]; // beta * H * L - gamma * L
      return;
    }
    ```
* `sdeDf`: The SDE diffusion function.  This can be specified on the standard deviation scale, or on the variance scale as above.  In this case, an **R** implementation would be
    
    ```r
    sde.diff <- function(x, theta) {
      df <- matrix(NA, 2, 2)
      df[1,1] <- theta[1]*x[1] + theta[2]*x[1]*x[2] # alpha * H + beta * H*L
      df[1,2] <- -theta[2]*x[1]*x[2] # -beta * H*L
      df[2,1] <- df[1,2] # -beta * H*L
      df[2,2] <- theta[2]*x[1]*x[2] + theta[3]*x[2] # beta * H*L + gamma * L
      df
    }
    ```
    In **C++** the specification is slightly different.  First we set `sdDiff = false` in order to tell **`msde`** to use the variance scale.  Next, the diffusion function is coded as
    
    ```cpp
    void sdeDf(double *df, double *x, double *theta) {
      df[0] = theta[0]*x[0] + theta[1]*x[0]*x[1]; // matrix element (1,1)
      df[2] = -theta[1]*x[0]*x[1]; // element (1,2)
      df[3] = theta[1]*x[0]*x[1] + theta[2]*x[1]; // element (2,2)
      return;
    }
    ```
    Thus there are two major differences with the **R** version.  The first is that the `df` matrix is stored as a vector (or "array" in **C++**).  Its elements are stored in column-major order, i.e., by stacking the columns one after the other into one long vector.  The second difference is that **`msde`** only uses the *upper triangular* portion of the (symmetric) matrix $\df_{\tth}(\y)$.  The elements below the diagonal can be set to any value or not set at all without affecting the computations.

	**`msde`** internally computes the Cholesky decomposition of the diffusion function when it is specified on the variance scale.  For small problems such as this one, it is more efficient to specify the Cholesky decomposition directly, i.e., specify the diffusion on the standard deviation scale.  In this case, the Cholesky decomposition of the diffusion function is
$$
\begin{bmatrix} \a H_t + \b H_tL_t & -\b H_tL_t \\ -\b H_tL_t & \b H_tL_t + \g L_t\end{bmatrix}^{1/2} = \begin{bmatrix} \sqrt{\a H_t + \b H_tL_t} & -\frac{\b H_tL_t}{\sqrt{\a H_t + \b H_tL_t}} \\ 0 & \sqrt{\b H_tL_t + \g L_t - \frac{(\b H_tL_t)^2}{\a H_t + \b H_tL_t}} \end{bmatrix},
$$
which is passed to **`msde`** by setting `sdDiff = true` and
    
    ```cpp
    void sdeDf(double *df, double *x, double *theta) {
      double bHL = theta[1]*x[0]*x[1]; // beta * H*L
      df[0] = sqrt(theta[0]*x[0] + bHL); // sqrt(alpha * H + bHL)
      df[2] = -bHL/df[0];
      df[3] = sqrt(bHL + theta[2]*x[1] - df[2]*df[2]);
      return;
    }
    ```
* `diagDiff`: A logical specifying whether or not the diffusion matrix $\df_\tth(\y)$ is diagonal.  For the Lotka-Volterra SDE model above, we set `diagDiff = false`.  However, if the diffusion function were of the form
$$
\df_\tth(\Y_t) = \begin{bmatrix} \a H_t + \b H_tL_t & 0 \\ 0 & \b H_tL_t + \g L_t\end{bmatrix},
$$
then we would set `diagDiff = true`, **and importantly**, treat the `df` argument of `sdeDf` directly as the diagonal elements of the diffusion matrix.  That is, the diffusion (on the variance scale) is encoded as
    
    ```cpp
    void sdeDf(double *df, double *x, double *theta) {
      double bHL = theta[1]*x[0]*x[1]; // beta * H*L
      df[0] = theta[0]*x[0] + bHL; // alpha * H + bHL
      // note the assignment to df[1] and not df[3]
      df[1] = bHL + theta[2]*x[1]; // bHL + gamma * L
      return;
    }
    ```
* `isValidParams`: A logical used to specify the parameter support.  In this case we have $\a, \b, \g > 0$, such that
    
    ```cpp
    bool isValidParams(double *theta) {
      bool val = theta[0] > 0.0;
      val = val && theta[1] > 0.0;
      val = val && theta[2] > 0.0;
      return val;
    }
    ```
* `isValidData`: A logical used to specify the SDE support, which can be parameter-dependent.  In this case we simply have $H_t, L_t > 0$, such that
    
    ```cpp
    bool isValidData(double *x, double *theta) {
      return (x[0] > 0.0) && (x[1] > 0.0);
    }
    ```

Thus the whole file `lotvolModel.h` is given below:

```cpp
#ifndef sdeModel_h
#define sdeModel_h 1

// Lotka-Volterra Predator-Prey model

// class definition
class sdeModel {
 public:
  static const int nParams = 3; // number of model parameters
  static const int nDims = 2; // number of sde dimensions
  static const bool diagDiff = false; // whether diffusion function is diagonal
  static const bool sdDiff = true; // whether diffusion is on sd or var scale
  void sdeDr(double *dr, double *x, double *theta); // drift function
  void sdeDf(double *df, double *x, double *theta); // diffusion function
  bool isValidParams(double *theta); // parameter validator
  bool isValidData(double *x, double *theta); // data validator
};

// drift function
inline void sdeModel::sdeDr(double *dr, double *x, double *theta) {
  dr[0] = theta[0]*x[0] - theta[1]*x[0]*x[1]; // alpha * H - beta * H*L
  dr[1] = theta[1]*x[0]*x[1] - theta[2]*x[1]; // beta * H*L - gamma * L
  return;
}

// diffusion function (sd scale)
inline void sdeModel::sdeDf(double *df, double *x, double *theta) {
  double bHL = theta[1]*x[0]*x[1]; // beta * H*L
  df[0] = sqrt(theta[0]*x[0] + bHL); // sqrt(alpha * H + bHL)
  df[2] = -bHL/df[0];
  df[3] = sqrt(bHL + theta[2]*x[1] - df[2]*df[2]);
  return;
}

// parameter validator
inline bool sdeModel::isValidParams(double *theta) {
  bool val = theta[0] > 0.0;
  val = val && theta[1] > 0.0;
  val = val && theta[2] > 0.0;
  return val;
}

// data validator
inline bool sdeModel::isValidData(double *x, double *theta) {
  return (x[0] > 0.0) && (x[1] > 0.0);
}

#endif
```
The additions to the previous code sections are:

1. The header include guards (`#ifndef`/`#define`/`#endif`).
2. The `sdeModel::` identifier is prepended to the class member definitions when these are written outside of the class declaration itself.
3. The `inline` keyword before the class member definitions, both of which ensure that only one instance of these functions is passed to the **C++** compiler.

### Compiling and checking the `sde.model` object

One the `sdeModel` class is created as the **C++** level, it is compiled in **R** using the following commands:
































