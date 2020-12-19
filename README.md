# Kelvin-Decomposition
This code is provided with the article "Handling tensors using tensorial Kelvin bases : application to olivine polycrystal deformation modeling using elastically anistropic CPFEM" Jean Furstoss · David Alejandro Ruiz Sarrazola · Marc Bernacki · Daniel Pino Muñoz

It permits, after giving the elastic constants of a material, to decompose the elastic constant matrix into different symmetries : isotropic, hexagonal, tetragonal, orthorhombic and finally monoclinic. The norms of the different projections are then computed and the last non-null projection gives the lowest symmetry of the material. This first part of the code is based on the work of J.T. Browaeys and S. Chevrot, Geophys. J. Int., 2004.

After identifying the lowest symmetry, the Kelvin base associated with the elastic constants provided is computed and then printed. The computation of the Kelvin bases is mainly based on the work of M.M. Mehrabadi and S.C. Cowin, Mech. appl. Math., 1990. The Kelvin base for monoclinic and triclinic symmetries are not considered here.

The Voigt notation used here assumes : (Sxx, Syy, Szz, Sxy, Syz, Szx) = CVoigt (Exx, Eyy, Ezz, 2Exy, 2Eyz, 2Ezx).

The code compilation can be done using gcc.
