# Land-cover maps based on pollen observations for Europe

The excel files in Land-cover Maps
1. land-cover-maps_1900CE,
2. land-cover-maps_11725CE,
3. land-cover-maps_4000BCE,

consist of six models described in details in (https://arxiv.org/abs/1703.06719)

	a) constant
	b) elevation
	c) K-L_ESM
	d) K-L-RCA3
	e) H-L_ESM
	f) H-L-RCA3

`K`: Kapplan human-land use (KK10),
`L`: LPJ-GUESS dynamic vegetation model forced by climate data from:
`ESM` :  Earth System model
`RCA3`: Rossby Centre Regional Climate Model

and

4. land-cover-maps_1425CE,
5. land-cover-maps_1000BCE,

only consist of `K-L_ESM` based on the results in
Pirzamanbein et al. (2018) (https://doi.org/10.1016/j.spasta.2018.03.005)

The columns represent

* Lon: longitude
* Lat: latitude
* C: coniferous forest
* B: Broadleaved forest
* U: Unforested land

These results obtained by running the Demo.m function in the src folder in MATLAB.
Demo.m sets the structures, dependency and the names. It calls Main.m function to run the MCMC sampling. The results of the Main.m functions are based on the papers mentioned above. However, some of the functions can be used in general cases, dealing with Gaussian Markov Random fields, Compositional Data and Dirichlet distribution.

[Behnaz's website](www.behnaz.pirzamanbin.name)
