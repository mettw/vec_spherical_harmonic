# vec_spherical_harmonic
Generate plots of vector spherical harmonics

Use VSphHarmonic() and SumVSphHarmonic() only.
Use:

> vharm = VSphHarmonic(1,-1,Parity.Complex);
> vharm.plot_MRadiation;

> svh = SumVSphHarmonic([1 2], [0 2], [Parity.Even Parity.Odd], [true false], 16);
> svh.plot_MNangVec;
