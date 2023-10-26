# HarmLAB
This toolbox offers an implementation of the Generalised Harmonic Balance method in MATLAB, supporting up to 2 base frequencies. This library can find periodic and quasi-periodic solutions of non-linear ODEs expressed in the following general form:

$$M_x \ddot{x} + C_x \dot{x} + K_x x + f_0 + f_{nl}\left({x},\dot{{x}},\ddot{{x}},{u},\dot{{u}},\ddot{{u}},\omega\right)={f}_e (t) = {M}_u\ddot{{u}} + {C}_u\dot{{u}}+{K}_u{u}$$

with an output given by:

$$y = g(x, \dot{x}, \ddot{x}, u, \dot{u}, \ddot{u} )$$

This toolbox has been developed to solve problems in structural dynamics, so much of the notation and terminology used stems from this field. However, it could equally be applied to problems in a variety of other disciplines.

The input is assumed to be of the following multi-harmonic form:

$${u} = A \sum_{k} \Re\left({U}_{k} e^{k\omega t}\right)$$

where the variable $A$ controls the excitation level.

The output is assumed to contain content at the same frequencies, so can be expressed as:

$${x} = \sum_{k} \Re\left({X}_{k} e^{k\omega t}\right)$$

This toolbox supports both Standard HBM with a single base frequency $\omega$ or Generalised HBM with 2 base frequencies which are interpendent:

$$ \begin{bmatrix}\omega_1 &  \omega_2\end{bmatrix} = \begin{bmatrix}\lambda_1 & \lambda_2\end{bmatrix} \Omega $$

Analytical jacbobians have been implemented throughout to speed up the execution. 

# Problem definition
## `problem` structure
The `problem` structure defines the system being simulated. This must have the following fields for the linear parts of the system:
- The mass, stiffness and damping matrices `K`, `C`, `M`
- The constant term `F0` (which defaults to `zeros`)
- The matrices for the excitation `Ku`, `Cu`, `Mu`

The number of DOF and inputs is inferred from the dimensions of the matrices.

The non-linear parts must be specified via the `model` callback which must have the form:

``` MATLAB
function varargout = test_model(part,States,hbm,problem)
switch part
    case 'nl'
        %f_nl
    case 'output'
        %y
    ...
end
end
```

where `part` can be one of:
  - `nl` which means this function should return $f_{nl}$
  - `nl_x`, `nl_xdot`, `nl_xddot` which means this function should return $\frac{\partial f_{nl}}{\partial x}$ etc
  - `nl_u`, `nl_udot`, `nl_uddot` which means this function should return $\frac{\partial f_{nl}}{\partial u}$ etc
  - `output` which means this function should return the output $y$

The excitation must be specified via the `excite` callback which must has the form:

``` MATLAB
function U = test_excite(hbm, problem, w0)
...
end
```
and returns the $U_k$.

You can also store other useful information in the `problem` structure which is needed by your model (eg model parameters).

### `problem.res` structure
For resonance problems (using `hbm_res` or `hbm_bb`), it is necessary to configure the term to be maximised. This can be of the following form:

$$ \frac{\|X_k|}{|F_{e,k}|} $$

The `res` structure should have the following fields to set the numerator and denominator:
- `iHarm` which sets harmonic to use the terms from
- `output` which can be either `x` (or its derivatives),`fnl`, or `none` (to have no output).
- `input` which can be either `u` (or its derivatives), `fe` or `unity` (to have no denominator).
- `NInput` which selects which DOF from `x` to use
- `NOutput` which selects which index from `fe` to use

## HBM structure
The other key structure needed to solve problems is the `hbm` structure, which stores information about the harmonics and other options. This has the following fields:
- `harm` contains information about harmonics
- `dependence` contains information about the form of $f_{nl}$
- `cont` contains settings for the continuation algorithm
- `options` contains settings for the continuation algorithm

### `harm` structure
This must have the following fields:

- `rFreqRatio` is the ratio of the base harmonics to the base frequency. ie $\lambda$. For HBM standard problems `rFreqRatio` should be set to 1.0. For GHBM problems this should be vector. 
- `NHarm` is the number of harmonics to include for each base frequency. This should have the same length as `rFreqRatio`.
- `Nfft` is the number of FFT points to use for the AFT transformation. This should be set to a power of 2 and have the same length as `rFreqRatio`.
- `iHarmPlot` sets which harmonic to plot on the FRF when using `hbm_frf` or `hbm_bb`. Typically you want to see the first harmonic so this should be set to 1.

### `dependence` structure
This should have the following fields:
- `x`, `xdot`, `xddot`  if the non-linearity has a dependence on the state and its derivatives
- `u`, `udot`, `uddot`  if the non-linearity has an explicit dependence on the state and its derivatives
- `w` if the non-linearity has an explicit frequency dependence

The value should be set to `1.0` or `True` where there is a dependence, and `0.0` or `False` otherwise
### `cont` structure
This is used to configure settings of the continuation algorithm. This should have the following fields:
- `method`: this should be one of the following values:
  - `predcorr` for a predictor-corrector algorithm (default)
  - `none` to simply step in frequency for `hbm_frf` or amplitude for `hbm_amp` and `hbm_bb`
If method is `predcorr`, you can set the following additional settings in the `cont.predcorr` structure:
- `solver` to choose which non-linear solver to use (`fipopt` to use IPOPT or `fsolve` to use the inbuilt solver)
- `step0`  which sets the initial step size
- `min_step`/`max_step`  which sets the minimum and maximum step size
- `num_iter_increase`/ `num_iter_reduce` which sets the threshold for increasing/decreasing the step size depending on the number of iterations
- `C` and `c` which sets the ratio to increase/decrease stepsize after a successful/unsuccessful step

### `options` structure
This sets other more general options and can have the following fields:
- `bUseStandardHBM`: force the solver to use the standard HBM when there is a single base frequency.
- `bVerbose`: toggle whether to suppress output to console
- `bPlot`: toggle whether to suppress plots

## Usage
Once `problem` and `hbm` have been defined, you must call the setup code:

``` MATLAB
[hbm,problem] = setuphbm(hbm,problem);
```
You can then call one of the following functions to solve a problem using HBM.

### One-off problems
The most simple case is to find a periodic solution for a given excitation. Two different functions are provided to solve such problems, depending on the asumptions made about the base frequency:

* `hbm_solve` : this assumes the base frequency $\Omega$ is fixed to a known value. This is the simplest and most common case. This can be called as follows:
  ``` MATLAB
  sol = hbm_solve(hbm,problem,w0,A);
  ```
  where `sol` contains information about the solution including the components of $x$, $f_{nl}$ and $u$ at each harmonic
* `hbm_res` : this assumes that the base frequency $\Omega$ can vary in order to find the maximum response at resonance. This is configured by the `problem.res` field:
  ``` MATLAB
  sol = hbm_res(hbm,problem,w0,A,X0);
  ```
If you do not have initial guess for `X0`, then this can be omitted or set to `[]`.

### Continuation problems
Each of the different types of one-off problem can then be solved over a range of frequencies. If the base frequency is assumed fixed, there are two potential continuation parameters:

* ```hbm_frf```: The base frequency $\Omega$ is used as the the continuation parameter, keeping the amplitude $A$ fixed. This yields a non-linear frequency response.
  ``` MATLAB
  sol = hbm_frf(hbm,problem,A,w0,X0,wEnd,XEnd);
  ```
* ```hbm_amp```: The base frequency $\Omega$ is kept fixed, and the amplitude $A$ is used as the the continuation parameter.
  ``` MATLAB
  sol = hbm_frf(hbm,problem,A,w0,X0,wEnd,XEnd);
  ```
If the base frequency is allowed to vary, and chosen to satisfy an objective, then there is only one potential continuation parameter:

* ```hbm_bb```: The excitation amplitude $A$ is the continuation parameter, and the base frequency $\Omega$ is chosen to satisfy the objective. This yields the "backbone" curve.
    ``` MATLAB        
    sol = hbm_bb(bb,problem,A0,w0,X0,Aend,wEnd,XEnd);
    ```
If you not have an intial guess for `X0` or `XEnd`, then this can be set to `[]` for any of these functions.
