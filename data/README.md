## Simulated responses

The [simulated_responses/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/data/simulated_responses/) folder contains the pressure and velocity responses that were simulated around a rigid spherical scatterer. The simulated data is stored into multi-dimensional matrices, that are constructed most of the time as:

* 1D - point in the field, 8 points with indexes:  
      1: center (right ear position when head in the field)  
      2: 2cm front  
      3: 5cm front  
      4: 10cm front  
      5: 2cm right side  
      6: 5cm right side  
      7: 10cm right side  
      8: 10cm left side (left ear position when head in the field)
* 2D - just 1
* 3D - used DOP (361 angles accross the horizontal plane, for example)
* 4D - frequencies (200 frequencies, if 50:50:10000)
* 5D - loudspeaker setup/reproduction (reference, dodeca, t-design,...,etc.)

In total the matrices are sized something like M(8,1,361,200,6).
Example: M(1,1,91,50,2) --> outputs the calculated pressure values for the center point (next to right ear) from 91-1=90 DOP angle (DOA=270) in unit circle notation, with frequency of 2500 Hz, using dodeca reproduction setup. Check savePressureFieldsXXX.m in [src/save_fields/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/src/save_fields/) for more detail on how the matrices are constructed.

NOTE: in the pre-calculated simulations without a scatterer, the 8th measurement point (i.e. 10cm to the left) was erroneously set to 20cm left. If you use this point to some calculations, it should be re-simulated.

### Naming of the simulated matrices in [simulated_responses/](https://github.com/mapu8/Spherical-scatterer-effects/tree/master/data/simulated_responses/)

* Some files have a pre-fix "diffuse" or "non-diffuse". This means if the file contains diffuse or non-diffuse pressure calculations. NOTE: by default, the files are non-diffuse. If a file does not have a pre-fix "diffuse", it contains non-diffuse data.
* Each file name contains the loudspeaker setup used, e.g. "dodeca", "reference", etc.
* Text "center" indicates that the file contains data, where the calculations were conducted in a setting where the scattering sphere is located in the center of the reproduction field. NOTE: if there is no "center" text within a file, the file contains data where the measurement points are in the center of the field (next to the right ear position), and the sphere is slightly off-center towards the left.
* Text "noscat" indicates that the file contains data made in a field without a scatterer. In this case, the measurement positions are positioned in the center of field.
* Text "delay" indicates that the file contains data, where the center measurement point has been delayed by 0.058 ms (corresponding to a 2cm travel length). These are used in the delay-and-sum beamformer calculations.
* Text "full" indicates that responses were calculated in every point in the horizontal plane. NOTE: by default, only single points are calculated next to the sphere. If a file name does not contain "full", the data was calculated with only single points.
* Some files have texts "d02", d05", d10", "0.02", "0.05", "0.10". These indicate the diameter of the sphere used in the calculations. NOTE: by default a sphere of diameter 0.02 m is used in the calculations. If some of the previously mentioned texts don't appear in the file name, 0.02 m diameter is used.
* Some files have numbers 1-6 in their filename, like ""non_diffuse_dodeca3_full_203.mat". These were used in the development to separate different versions of files, these can be ignored.
* All files have some number at the end. These indicate the array number in Triton server, when the matrices were calculated. These can be ignored.
* "Velocity" in a file name indicates that the file contains velocity data, instead of pressure data. In this case, the saved matrices are 6D, where the last dimension relates to the axis (1=x, 2=y, z=3). NOTE: by default the files contain pressure data.  
IMPORTANT: the coordinate system follows the right-handed coordinate system, where x is to the right, y is to the front, and z is up. Figure 1 in [1] is misleading in this sense, because x is set to the front. This was only made to make the graph more readable, as most publications use 0 degree as the frontal direction.


## References

[1] Pajunen, L., Politis, A., Pulkki, V., Vaalgamaa, M., Strömmer, S., 2020,
    Effects of rigid spherical scatterer on spatial audio reproduction quality.
    Submitted for Audio Engineering Society Convention 148.




