# Effects of a rigid spherical scatterer in spatial audio reproduction fields

#### A library for analysing the effects of a rigid spherical scatterer in virtual loudspeaker setups

This Matlab library was developed during master's thesis research in the Department of Signal Processing and Acoustics, Aalto University, Finland. The thesis [1] is available through Aalto University's websites.

This library also complements the related convention paper [2].

## Description

This Matlab library can be used to calculate and analyse pressure and velocity fields around a rigid spherical scatterer in different virtual loudspeaker setups.

### Quick guide

* Calculate and save pressure/velocity values into matrices with functions in [save_fields/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/save_fields/)
* Analyse the pressure fields and recreate the plots presented in [2] with [pressure_errors.m](https://github.com/mapu8/Spherical-scatterer-effects/blob/master/src/analysis/pressure_errors.m)

### Prerequisites

Matlab libraries from [Politis](https://github.com/polarch) should be installed and added to the Matlab path before using this library. The necessary libraries are:

* [Spherical-Harmonic-Transform](https://github.com/polarch/Spherical-Harmonic-Transform)
* [Spherical-Array-Processing](https://github.com/polarch/Spherical-Array-Processing)
* [Higher-Order-Ambisonics](https://github.com/polarch/Higher-Order-Ambisonics)
* [Vector-Base-Amplitude-Panning](https://github.com/polarch/Vector-Base-Amplitude-Panning)
* [Array-Response-Simulator](https://github.com/polarch/Array-Response-Simulator)


### Calculating and saving pressure/velocity fields

To calculate and save pressure/velocity fields into files, functions in the [save_fields/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/save_fields/) folder should be used. During the master's thesis process, the calculations were executed on Aalto's high-performance computing cluster, [Triton](https://scicomp.aalto.fi/triton/). As a remnant, [savePressureFieldsTriton.m](https://github.com/mapu8/Spherical-scatterer-effects/blob/master/src/save_fields/savePressureFieldsTriton.m) function contains a collection of different variations for calculations, including pressure/velocity fields in different loudspeaker setups in diffuse and direct sound fields. The function can be used as it is, or as a reference on how to call the functions in [save_fields/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/save_fields/). 

As an example, the function could be called as:

```
savePressureFieldsTriton("/some/root/folder/","data/simulated_responses/test_output.mat",0,"dodeca",1234)
```
Inputs:
- "/some/root/folder/" is the path to this library
- "data/simulated_responses/test_output.mat" is the path to the output file
- 0 is the number of loops used in diffuse field calculations (use 0 for direct sound field calculations)
- "dodeca" determines the virtual loudspeaker setup
- 1234 is the rng seed for the random number generator inside the function (NOTE: this should be different for each diffuse field calculation with the same setup. For example, if multiple function calls are used to calculate diffuse field responses, each function call should have it's own rng seed. Otherwise the randomisation of the phases fails, and each function outputs the same results, because each function uses the same random phases.)

The above call would calculate the pressure field for a dodecahedron loudspeaker setup with a scatterer in the field. The frequency band, DOAs, points of calculations, and the scatterer size and position are determined inside the function [savePressureFieldsTriton.m](https://github.com/mapu8/Spherical-scatterer-effects/blob/master/src/save_fields/savePressureFieldsTriton.m). These should be changed according to what is the desired output. The calculated pressure matrix is then stored to test_output.mat.

Alternatively, the actual calculation function could be called:

```
savePressureFieldsReproduction(file,head_diameter,head_pos,X,Y,D,d,f_vec,dop_vec,decoding_method,speaker_dirs,source_distance,indexes,max_dist_p)
```
In the above function call, the arguments should be filled with meaningful values. See savePressureFieldsTriton.m and the function descriptions in [save_fields/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/save_fields/) for more detail.

NOTE: the calculations can be heavy depending on the resolution used in the grid, the amount of measurement points calculated, the frequency band used, and the amount of different source directions taken into account. For lighter calculations, only single points in the field can be calculated by using indexed values in the function call.


### Analysing the results

When the calculations are done and the matrices are obtained, the results can be analysed. This can be done with scripts in the [analysis/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/analysis/) folder.

[pressure_errors.m](https://github.com/mapu8/Spherical-scatterer-effects/blob/master/src/analysis/pressure_errors.m) contains only the pressure error calculations. It can be used to reproduce the plots found in [2].

## Gifs & Images

The [gifs/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/gifs/) folder contains various gifs, that illustrate the effects of a rigid sphere to a single point source. The [images/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/images/) folder contains the pressure error plots used in [2].

![Couldn't find image](https://github.com/mapu8/Spherical-scatterer-effects/blob/master/gifs/pressure_f1000_r0.50.gif "Glorienttes gif")



## Authors

* **[Lauros Pajunen](https://github.com/mapu8)** - *Initial work*
* **[Archontis Politis](https://github.com/polarch)** - *Implementations of pressure and velocity field calculations*


## References

[1] Pajunen, L., 2019,
    Effects of a rigid spherical scatterer in spatial audio reproduction fields.
    Master's thesis, Aalto University, Espoo, Finland.

[2] Pajunen, L., Politis, A., Pulkki, V., Vaalgamaa, M., Strömmer, S., 2020,
    Effects of rigid spherical scatterer on spatial audio reproduction quality.
    Submitted for Audio Engineering Society Convention 148.





