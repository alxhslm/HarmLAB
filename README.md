# matlab-hbm
This toolbox offers and implementation of the Generalised Harmonic Balance method in MATLAB, supporting up to 2 base frequencies. This library can find periodic and quasi-periodic solutions of non-linear ODEs expressed in the following general form:

<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbf{M}_x\ddot{\mathbf{x}}&space;&plus;&space;\mathbf{C}_x\dot{\mathbf{x}}&plus;\mathbf{K}_x\mathbf{x}&space;&plus;&space;\mathbf{f}_{nl}\left(\right&space;\mathbf{x}&space;,\dot{\mathbf{x}},\ddot{\mathbf{x}},\mathbf{u},\dot{\mathbf{u}},\ddot{\mathbf{u}},\omega)=\mathbf{f}_e&space;(t)&space;=&space;\mathbf{M}_u\ddot{\mathbf{u}}&space;&plus;&space;\mathbf{C}_u\dot{\mathbf{u}}&plus;\mathbf{K}_u\mathbf{u}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbf{M}_x\ddot{\mathbf{x}}&space;&plus;&space;\mathbf{C}_x\dot{\mathbf{x}}&plus;\mathbf{K}_x\mathbf{x}&space;&plus;&space;\mathbf{f}_{nl}\left(\right&space;\mathbf{x}&space;,\dot{\mathbf{x}},\ddot{\mathbf{x}},\mathbf{u},\dot{\mathbf{u}},\ddot{\mathbf{u}},\omega)=\mathbf{f}_e&space;(t)&space;=&space;\mathbf{M}_u\ddot{\mathbf{u}}&space;&plus;&space;\mathbf{C}_u\dot{\mathbf{u}}&plus;\mathbf{K}_u\mathbf{u}" title="\mathbf{M}_x\ddot{\mathbf{x}} + \mathbf{C}_x\dot{\mathbf{x}}+\mathbf{K}_x\mathbf{x} + \mathbf{f}_{nl}\left(\right \mathbf{x} ,\dot{\mathbf{x}},\ddot{\mathbf{x}},\mathbf{u},\dot{\mathbf{u}},\ddot{\mathbf{u}},\omega)=\mathbf{f}_e (t) = \mathbf{M}_u\ddot{\mathbf{u}} + \mathbf{C}_u\dot{\mathbf{u}}+\mathbf{K}_u\mathbf{u}" /></a>

This toolbox has been developed to solve problems in structural dynamics, so much notation and terminology used stems from this field. However, it could equally be applied to problems in a variety of other disciplines.

## Types of problem

### One-off problems
The most simple case is to find a periodic solution for a given excitation. Two different functions are provided to solve such problems, depending on the asumptions made about the base frequency:

* ```hbm_solve``` : this assumes the base frequency is fixed to a known value. This is the simplest and most common case.
* ```hbm_res``` : this assumes that the base frequency can vary. To make the problem solvable, it is necessary to choose an objective to be maximised or minimised, therby finding the 'optimal' solution. For example, this function could be used to find the maximum response at resonance.

### Continuation problems
Each of the different types of one-off problem can then be solved over a range of frequencies. Two different functions are provided:

* ```hbm_frf```: If the base frequency is used as the the continuation parameter, this yields a non-linear frequency response. This can be used to compute the non-linear frequeny response function (FRF).
* ```hbm_bb```: If the base frequency is allowed to vary, the excitation amplitude is the continuation parameter. This yields the 'optimal' response over a range of amplitudes. For example, this can used to compute the backbone curve.

## Solvers
This toolbox can interface with two different optimisation toolboxed:

* IPOPT through a custom wrapper ```fipopt```
* The inbuilt MATLAB optimisers ```fsolve``` and ```fmincon```

For the continuation problems, three different continuation algorithm can be used:

* Pseudo-arc length continuation
* Arc-length continutation
* Using the 'ep' toolbox from ```coco```

Analytical jacbobians have been provided throughout to speed up the execution.