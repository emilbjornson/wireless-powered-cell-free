Joint Power Control and LSFD for Wireless-Powered Cell-Free Massive MIMO
==================

This is a code package is related to the following scientific article:

Özlem Tuğfe Demir and Emil Björnson, “[Joint Power Control and LSFD for Wireless-Powered Cell-Free Massive MIMO](https://ieeexplore.ieee.org/document/9258425
),” IEEE Transactions on Wireless Communications, vol. 20, no. 3, pp. 1756-1769, March 2021, doi: 10.1109/TWC.2020.3036281.

The package contains a simulation environment, based on Matlab, that reproduces some of the numerical results and figures in the article. We encourage you to also perform reproducible research!


## Abstract of Article

Wireless communication systems have almost exclusively operated in the far-field of antennas and antenna arrays,
which is conventionally characterized by having propagation
distances beyond the Fraunhofer distance. This is natural since
the Fraunhofer distance is normally only a few wavelengths.
With the advent of active arrays and passive reconfigurable
intelligent surfaces (RIS) that are physically large, it is plausible
that the transmitter or receiver is located in between the
Fraunhofer distance of the individual array/surface elements and
the Fraunhofer distance of the entire array. An RIS then can be
configured to reflect the incident waveform towards a point in
the radiative near-field of the surface, resulting in a beam with
finite depth, or as a conventional angular beam with infinity
focus, which only results in amplification in the far-field. To
understand when these different options are viable, an accurate
characterization of the near-field behaviors is necessary. In this
paper, we revisit the motivation and approximations behind the
Fraunhofer distance and show that it is not the right metric for
determining when near-field focusing is possible. We obtain the
distance range where finite-depth beamforming is possible and
the distance where the beamforming gain tapers off

## Content of Code Package

The article contains 12 simulation figures, numbered 2-13. Figures 2 and 3 are generated by the Matlab script Fig2_3.m. Figures 4, 5, 12, and 13 are generated by the Matlab script Fig4_5_12_13.m. Figures 6 and 7 are generated by the Matlab script Fig6_7.m. Figures 8, 9, 10, and 11 are generated respectively by the Matlab scripts Fig8.m, Fig9.m, Fig10.m, and Fig11.m. The package also contains the Matlab scripts functionExampleSetup.m, functionExpectations_lmmse.m, functionExpectations_lmmse_random_pilot_sequence.m, functionExpectations_ls.m, optimize_power_coh_lin.m, optimize_power_coh_nonlin.m, optimize_power_noncoh_lin.m, and optimize_power_noncoh_nonlin.m that are MATLAB functions used by some of the scripts.

For the optimization problems, the convex programming solver CVX http://cvxr.com/cvx/ is used. Please adjust the default solver of CVX to SDPT3 not to encounter any issues. The version we have tested was SDPT3 4.0. 

See each file for further documentation.

## Acknowledgements

The work of Ö. T. Demir and E. Björnson was partially supported by ELLIIT and the Wallenberg AI, Autonomous Systems and Software Program (WASP) funded by the Knut and Alice Wallenberg Foundation. 

## License and Referencing

This code package is licensed under the GPLv2 license. If you in any way use this code for research that results in publications, please cite our original article listed above.
