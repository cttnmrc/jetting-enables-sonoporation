# Cyclic jetting enables microbubble-mediated drug delivery

MATLAB codes used in the paper *Cyclic jetting enables microbubble-mediated drug delivery*, submitted to Nature Physics.

## Scripts description
* `BubbleRadiusDisplacementImpulseStress.m`: Script to numerically simulate the time evolution of the bubble radius, bubble displacement, dimensionless impulse and stress generated. Corresponding to Fig. 2a-e and Fig. 4e.  
* `ModifiedRayleighPlesset_MarmottantZhou.m`: Function to solve the modified Rayleigh-Plesset equation with the Marmottant model as shell model and the Zhou model as the gas core model.  
- `ExperimentalPressureEnvelope.m`: Function to extract the envelope of the experimentally recorded ultrasound driving signal.  
- `UltrasoundPulse.dat`: Experimentally recorded ultrasound driving signal.  
- `BubbleFeaturesExtraction.m`: Script to extract bubble radius and displacement from video recordings.

## Environment to run to scripts 

All the codes are written in MATLAB and tested on the release R2023a. 

## Authors

**Marco Cattaneo<sup>1\*</sup> , Giulia Guerriero<sup>1</sup>, Gazendra Shakya<sup>1</sup>, Lisa A. Krattiger<sup>2</sup>, Lorenza G. Paganella<sup>3,4</sup>, Maria L. Narciso<sup>4,5</sup> and Outi Supponen<sup>1\*</sup>**

<sup>1</sup> Institute of Fluid Dynamics, ETH Zürich, Sonneggstrasse 3, Zürich,8092, Switzerland.  
<sup>2</sup> Department of Obstetrics, University Hospital Zürich, University of Zürich, Schmelzbergstrasse 12, Zürich, 8091, Switzerland.  
<sup>3</sup> Institute of Energy and Process Engineering, ETH Zürich, Leonhardstrasse 21, Zürich, 8092, Switzerland.  
<sup>4</sup> Institute for Mechanical Systems, ETH Zürich, Leonhardstrasse 21, Zürich, 8092, Switzerland.  
<sup>5</sup> Swiss Federal Laboratories for Materials Science and Technology (EMPA), Überlandstrasse 129, Dübendorf, 8600, Switzerland.  

<sup>\*</sup> Corresponding authors: mcattaneo@ethz.ch, outis@ethz.ch

## Funding 

We acknowledge ETH Zürich for the financial support (Research Grant 1-010206-000).

## License

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

This library is a free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details http://www.gnu.org/licenses.


