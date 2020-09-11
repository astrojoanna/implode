# implode
IMPLODE stays for "planetesImal forMation via Pebble cLOuD collapsE"

This is a Fortran code to calculate collapse of a bound pebble cloud driven by inelastic collisions of pebbles and gas drag written by Joanna Drążkowska and Rico Visser.

## Code usage

Compile with `make`

Run with `./implode parameters`
where `parameters` is a parameter file. A couple of example files are included in the repository, in particular `virialtest.par` is a setup with no gas and no collisions considered where the particles are initialized with random velocities corresponding to virial equilibrium, `settlingtest.par` is a test in which there is gas but no collisions and the particles are not given any initial velocity dispersion, and `setup.par` corresponds to one of the models presented in Visser et al. paper.

WARNING: the code is OPENMP parallel and thus will by default use all availble resources. To set the number of threads use "set OMP_NUM_THREADS=n" or "export OMP_NUM_THREADS=n" depending on your operating system.

Before re-compiling, clean the directory with "make clean" - WARNING: it will remove the output if you don't re-name it. When you just change the parameter file and not the source code, there is no need to re-compile. The simulation will not run if the output would be overwritten.

Output is stored in `output...dat` file and it is written every n timesteps, the frequency can be set in the parameter file. Additionally, there is a time step criterion to avoid sparse output, the maximum time step between outputs can also be set in the parameter file. The output format is:
1. The first row has 3 columns: these are time in years, physical mass represented by one superparticle `mswarm` (in gram, constant), initial radius of the cloud `R0` (in cm)
The other rows have NN(=total number of superparticles used) columns each and are:
2. x position [cm]
3. y position [cm]
4. z position [cm]
5. mass of the physical particle represented by each superparticle, normalized by its initial mass (=1.0 if the mass didn’t change during collisions)
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
