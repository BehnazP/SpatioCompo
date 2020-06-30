# Land-cover maps based on pollen observations for Europe
### The usage of the land-cover maps is free, however attribution to the publication is required.

The `.csv` files in Land-cover Maps folder
1. land-cover-maps_1900CE,
2. land-cover-maps_1725CE,
3. land-cover-maps_4000BCE,

consist of six models described in details in (https://arxiv.org/abs/1703.06719)

a) constant  
b) elevation  
c) K-L_ESM  
d) K-L_RCA3  
e) H-L_ESM  
f) H-L_RCA3

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

### Reference:
* Pirzamanbein, Behnaz, Johan Lindström, Anneli Poska, and Marie-José Gaillard. "Modelling Spatial Compositional Data: Reconstructions of past land cover and uncertainties." Spatial statistics 24 (2018): 14-31.
* Pirzamanbein, Behnaz, Anneli Poska, and Johan Lindström. "Bayesian reconstruction of past land-cover from pollen data: model robustness and sensitivity to auxiliary variables." arXiv preprint arXiv:1703.06719 (2017).

### Licenses
#### Code
MATLAB code for modelling spatial compositional data, specifically for reconstructions of past land cover and uncertainties.
Copyright (C) 2018  Pirzamanbein Pirzamanbein

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
#### Data
<a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-sa/4.0/88x31.png" /></a><br /><span xmlns:dct="http://purl.org/dc/terms/" href="http://purl.org/dc/dcmitype/Dataset" property="dct:title" rel="dct:type">Land-cover maps based on pollen observations for Europe</span> by <span xmlns:cc="http://creativecommons.org/ns#" property="cc:attributionName">Behnaz Pirzamanbein</span> is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-sa/4.0/">Creative Commons Attribution-ShareAlike 4.0 International License</a>.<br />Based on a work at <a xmlns:dct="http://purl.org/dc/terms/" href="https://github.com/BehnazP/SpatioCompo" rel="dct:source">https://github.com/BehnazP/SpatioCompo</a>.
