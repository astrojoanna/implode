# implode
IMPLODE stays for "planetesImal forMation via Pebble cLOuD collapsE"

This is a Fortran code to calculate collapse of a bound pebble cloud driven by inelastic collisions of pebbles and gas drag written by Joanna Drążkowska and Rico Visser.

The representative particle approach is used throught the code. The motion of particles is integrated with a Runge-Kutta Fehlberg variable step scheme. A Monte Carlo algorithm is used to solve the collisional evolution.

## Code usage

### Setup

Compile with `make`

Run with `./implode parameters`
where `parameters` is a parameter file. A couple of example files are included in the repository, in particular `virialtest.par` is a setup with no gas and no collisions considered where the particles are initialized with random velocities corresponding to virial equilibrium, `settlingtest.par` is a test in which there is gas but no collisions and the particles are not given any initial velocity dispersion, and `setup.par` corresponds to one of the models presented in Visser et al. paper.

WARNING: the code is OPENMP parallel and thus will by default use all availble resources. To set the number of threads use "set OMP_NUM_THREADS=n" or "export OMP_NUM_THREADS=n" depending on your operating system.

Before re-compiling, clean the directory with "make clean" - WARNING: it will remove the output if you do not re-name it. When you just change the parameter file and not the source code, there is no need to re-compile. The simulation will not run if the output is overwritten.

Parameters that can be set in the parameter file and changed without re-compiling the code are:
1. `gas`                      is gas drag taken into account? (logical)
2. `collisions`               are collisions taken into account? (logical)
3. `nr_parts`                 number of representative particles (integer)
4. `nr_parts_per_zone`        minimum number of representative particles per zone (integer, note that doe to Monte Carlo method convergence it should not be below 200)
5. `radial_distance_AU`       radial distance of the clump (float, [AU])
6. `gas_surfde_cgs`           surface density of the gas (at the clump location, float [g cm-2])
7. `final_radius_km`          final radius of planetesimal (float, [km])
8. `max_time_yr`              maximum time of the simulation (float, [years])
9. `delta_t_omegaK`           basic timestep for advection solver (float, [1/OmegaK the local Keplerian frequency])
10. `softening_length_Rsolid` softening length for the advection solver (float, [in the units of the final planetesimal radius])
11. `init_min_size`           minimum size of pebbles in the initial size distribution (float, [cm])
12. `init_max_size`           maximum size of pebbles in the initial size distribution (float, [cm])
13. `initial_size_distr`      is the initial size distribution considered? (logical, if false, all grains start at `init_max_size`)
14. `size_distr_kappa`        slope of the initial size distribution (float, default is 0.1666666667 corresponding to the MRN distribution)
15. `monomer_size_cm`         monomer size (sets the minimum grain size in the case of fragmentation, float, [cm])
16. `material_density_cgs`    material denisty of pebbles and the final planetesimal (float, [g cm-3])
17. `bouncing`                is bouncing included as a collision outcome? (logical)
18. `COR`                     coefficient of restitution for bouncing collisions (float)
19. `fragmentation`           is fragmentation included as a collision outcome? (logical)
20. `initial_vel_disp_factor` multiplier of the initial velocities (float, 1 corresponds to the virial equilibrium)
21. `init_free_fall`          initialize the radial velocity component with free-fall/terminal speed velocity? (logical)
22. `output_freq`             number of timesteps between outputs (integer)
23. `output_dt_year`          maximum time between outputs (float [years])
24. `screen_out`              write extended screen output? (logical)
25. `full_output`             should full output be written? (see below)

### Output

The full output is stored in `output...dat` file and it is written every `output_freq` timesteps, the frequency can be set in the parameter file. Additionally, there is a timestep criterion `output_dt_year` to avoid sparse output, the maximum timestep between outputs can also be set in the parameter file. The output format is:
1. The first row has 3 columns: these are time in years, physical mass represented by one superparticle `mswarm` (in gram, constant), initial radius of the cloud `R0` (in cm)
The other rows have NN(=total number of superparticles used) columns each and are:
2. x position [cm]
3. y position [cm]
4. z position [cm]
5. mass of the physical particle represented by each superparticle, normalized by the initial minimum pebble mass specified in the parameter file
6. x velocity component [cm/s]
7. y velocity component [cm/s]
8. z velocity component [cm/s]
9. number of collision the superparticle underwent (integer)
10. is the particle part of the core (not evolving anymore) (logical)
11. number of bouncing collisions the superparticle underwent (integer)
12. number of fragmenting collisions the superparticle underwent (integer)

Additionally, there is the `energy...dat` file, which, for every timestep, saves the 1. time [yrs], 2. total kinetic energy T [erg], 3. potential energy abs(U) [erg], 4. core mass fraction, 5. total number of collisions, 6. total number of bouncing collisions, 7. total number of fragmentation collisions, 8. optical depth of active particles containing the inner 1% of mass, 9. settling threshold radius location [Hill radii]. 

The `sizedistr-init...dat` contains the initial size distribution of particles and `sizedistr...dat` contains the final size distribution. The columns are 1. mass in the units of monomer mass (specified in the parameter file), 2. size [cm], 3. the distribution function f(m)m^2.

## License

This code is licensed under the GNU General Public License v3. Among other things, this means that if you modified this code for a publication, it needs to be open sourced under the same license.

## Credits

Please cite Visser, Drążkowska, & Dominik (2020).
