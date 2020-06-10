# vec_spherical_harmonic
Generate plots of vector spherical harmonics

[![View vec_spherical_harmonic on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://au.mathworks.com/matlabcentral/fileexchange/72096-vec_spherical_harmonic)

Use VSphHarmonic() and SumVSphHarmonic() only.

Use:

> vharm = VSphHarmonic(1,-1,Parity.Complex);

> vharm.plot_MRadiation;

> svh = SumVSphHarmonic([1 2], [0 2], [Parity.Even Parity.Odd], [true false], 16);

> svh.plot_MNangVec;
