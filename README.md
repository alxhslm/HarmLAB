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

## Types of problem

### One-off problems
The most simple case is to find a periodic solution for a given excitation. Two different functions are provided to solve such problems, depending on the asumptions made about the base frequency:

* ```hbm_solve``` : this assumes the base frequency is fixed to a known value. This is the simplest and most common case.
* ```hbm_res``` : this assumes that the base frequency can vary. To make the problem solvable, it is necessary to choose an objective to be maximised or minimised, therby finding the 'optimal' solution. For example, this function could be used to find the maximum response at resonance.

### Continuation problems
Each of the different types of one-off problem can then be solved over a range of frequencies. If the base frequency is assumed fixed, there are two potential continuation parameters:

* ```hbm_frf```: The base frequency is used as the the continuation parameter, keeping the amplitude fixed. This yields a non-linear frequency response.
* ```hbm_amp```: The base frequency is kept fixed, and the amplitude is used as the the continuation parameter.

If the base frequency is allowed to vary, and chosen to satisfy an objective, then there is only one potential continuation parameter:

* ```hbm_bb```: The excitation amplitude is the continuation parameter, and the base frequency is chosen to satisfy the objective. This yields the 'optimal' response over a range of amplitudes. For example, this can used to compute the backbone curve.

## Solvers
This toolbox can interface with two different optimisation toolboxes:

* IPOPT through a custom wrapper ```fipopt```
* The inbuilt MATLAB optimisers ```fsolve``` and ```fmincon```

For the continuation problems, two different continuation toolboxes can be used:

* Predictor-corrector algorithm from this toolbox
* The ```ep``` toolbox from the third-party continuation toolbox ```coco```

For the predictor-corrector algorithm, two different types of corrector can be used:

* Pseudo-arc length with Moore-penrose Inverse corrections
* Arc-length using either ```fipopt``` or ```fsolve``` to find the solution

Analytical jacbobians have been implemented throughout to speed up the execution.

# Usage
## Problem structure
The `problem` structure defines the system being simulated. This must have the following fields for the linear parts of the system:
- The mass, stiffness and damping matrices `K`, `C`, `M`
- The constant term `F0` (which defaults to `zeros`)
- The matrices for the excitation `Ku`, `Cu`, `Mu`

The non-linear parts must be specified via 2 additional fields:

- `model` is a callback specified by the user which must have the form:
    ``` MATLAB
    function varargout = test_model(part,States,hbm,problem)
    ...
    end
    ```
    where `part` can be one of:
  - `nl` which means this function should return $f_{nl}$
  - `nl_x`, `nl_xdot`, `nl_xddot` which means this function should return $\frac{\partial f_{nl}}{\partial x}$ etc
  - `nl_u`, `nl_udot`, `nl_uddot` which means this function should return $\frac{\partial f_{nl}}{\partial u}$ etc
  - `output` which means this function should return the output $y$
- `excite` which is another callback which must has the form:
    ``` MATLAB
    function U = test_excite(hbm, problem, w0)
    ...
    end
    ```
    and returns the `U_k`.

You can also store other useful information in the `problem` structure which is needed by your model (eg model parameters).