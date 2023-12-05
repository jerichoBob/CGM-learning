# CGM-learning

A repo for my research on CGM (Circumgalactic Medium)

I am beginning work with Rongmon Bordoloi's research group @ NC State.  This will be where I keep my notes. Will open it up if that seems appropriate.

## First things first - Create your environment

We use conda, so install miniconda by opening up the terminal of your choice do the following:

```sh
# Download the things
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-MacOSX-arm64.sh
bash Miniforge3-MacOSX-arm64.sh

# make sure that "our" conda is active when a new shell starts up
~/miniforge3/bin/conda config --set auto_activate_base true

# create the environment
conda create -n astroresearch python=3.9.5

# activate the environment
conda activate astroresearch
```

<!-- I have switch allegences from conda over to simple venv (built into python distro) to minimize the number of competing environments I have on my machine.

```sh
python3 -m venv ~/venvs/astroresearch
source ~/venvs/astroresearch/bin/activate
``` -->

Next, download the things (if you haven't done that already)

```sh
cd ~/School/github_repos # or wherever you have your local github repos

# clone the needful
git clone https://github.com/rongmon/rbcodes.git
git clone https://github.com/pypeit/kcwitools.git
git clone https://github.com/linetools/linetools.git (NOTE: there is a conda installer for linetools now!!!)
```

And finally, install all the things ...

```sh
cd  ~/School/github_repos/rbcodes; pip install -r requirements_simple.txt

# and add these to fill in the gaps - don't know what I could do if online tools depend on these libraries... they need to be installable via pip/conda but aren't yet
cd ~/School/github_repos/kcwitools/; python setup.py install
cd ~/School/github_repos/linetools; python setup.py develop

# now you can run the extraction tool, rbspecgui and Analyze-2.18133.py (which depends on astropy, kcwitools and linetools)
```

## Useful Links

* [rbcodes](https://github.com/rongmon/rbcodes) - various python tools for Rongmon's group
* [linetools](https://github.com/linetools/linetools) - also [https://linetools.readthedocs.io/](https://linetools.readthedocs.io/)
* [kcwitools](https://github.com/pypeit/kcwitools) - Tools for KCWI Reduction and Analysis
* [astropy-learning](https://github.com/jerichoBob/astropy-learning) - my documentation on coming up to speed on astropy 

## Papers

* [2017-The Circumgalactic Medium.pdf](https://www.annualreviews.org/doi/abs/10.1146/annurev-astro-091916-055240)
* [2022-Nature-Resolving the H i in damped Lyman α systems that power star formation-.pdf](https://www.nature.com/articles/s41586-022-04616-1)
--
* [2007-Lilly_2007_ApJS_172_70.pdf](https://iopscience.iop.org/article/10.1086/516589)
* [2008-What is L*?: Anatomy of the Galaxy Luminosity Function
](https://arxiv.org/abs/astro-ph/0504580)
* [2011-Bordoloi-THE RADIAL AND AZIMUTHAL PROFILES OF Mg ii ABSORPTION AROUND 0.5 < z < 0.9 zCOSMOS GALAXIES OF DIFFERENT COLORS, MASSES, AND ENVIRONMENTS](http://iopscience.iop.org/article/10.1088/0004-637X/743/1/10/pdf)
* [2011-Into the central 10 pc of the most distant known radio quasar](https://www.aanda.org/articles/aa/abs/2011/07/aa17341-11/aa17341-11.html)
* [2011-CO (2–1) LINE EMISSION IN REDSHIFT 6 QUASAR HOST GALAXIES](https://iopscience.iop.org/article/10.1088/2041-8205/739/1/L34)
* [2018-A Budget and Accounting of Metals at z=0 -- Results from the COS-HALOS Survey.pdf](https://ui.adsabs.harvard.edu/abs/2014ApJ...786...54P/abstract)

## Some New Terms (for me)

* **Column Density** - [here](https://astronomy.swin.edu.au/cosmos/C/Column+Density) - a measure of the amount of intervening matter between an observer and the object being observed.
* **Lyman alpha forest** - [here](https://en.wikipedia.org/wiki/Lyman-alpha_forest) - The Lyman-alpha absorption lines in the quasar spectra result from intergalactic gas through which the galaxy or quasar's light has traveled. Since neutral hydrogen clouds in the intergalactic medium are at different degrees of redshift (due to their varying distance from Earth), their absorption lines are observed at a range of wavelengths. Each individual cloud leaves its fingerprint as an absorption line at a different position in the observed spectrum.
* **Virial Theorem** - [here](https://en.wikipedia.org/wiki/Virial_theorem) -   relates the average over time of the total kinetic energy of a stable system of discrete particles, bound by potential forces (forces characterized exclusively by potential difference), with that of the total potential energy of the system.
* **Virial Radius** - here - The virial radius of a gravitationally bound astrophysical system is the radius within which the virial theorem applies. It is defined as the radius at which the density is equal to the critical density .
  * how is this different from [Hill Sphere](https://en.wikipedia.org/wiki/Hill_sphere) / [Roche Lobe](https://en.wikipedia.org/wiki/Roche_lobe)?
    * nice explanation of Hill Sphere vs Roche Lobe [here](https://astronomy.stackexchange.com/questions/47907/whats-the-difference-between-the-roche-lobe-and-roche-sphere).
* **Critical Density** - [here](https://en.wikipedia.org/wiki/Friedmann_equations#Density_parameter) - The density parameter Ω is defined as the ratio of the actual (or observed) density ρ to the critical density ρc of the Friedmann universe. The relation between the actual density and the critical density determines the overall geometry of the universe
* **Galactic Mass/Virial Mass/Virial Radius** - nice commentary [here](https://physics.stackexchange.com/questions/406867/intuitive-understanding-of-the-virial-radius-mass) - 
  * excerpt: The basic problem of defining the mass of a galaxy is that the density decreases from the center and out, but never reaches zero. Hence, there is no well-defined radius at which to stop counting particles, so when you say "the actual mass", this concept doesn't really exist.
  * **Virial mass**
  * **Overdensity mass**
* **Half-light mass**
* odd turn of phrase: **3 decades in stellar mass** - 
* **L\*, sub-\*, - What is L⋆?**: Anatomy of the Galaxy Luminosity Function (2005) - Asantha Cooray AND Milos Milosavljevic 
* **DLA** - Damped Lyman Alpha
  * A ground-based imaging study of galaxies causing damped Lyman α (DLA), sub-DLA and Lyman limit system absorption in quasar spectra*
  * Column Density: N ~ 20
* **LLS** - Lyman Limit System
* **z-index / redshift**
  * https://en.wikipedia.org/wiki/Redshift#/media/File:Distance_compared_to_z.png
* **MOC** - Multi-Order Coverage map - [here](https://ivoa.net/documents/MOC/)
* [QSO - Quasi-Stellar Object](http://www.stargazing.net/david/qso/index.html#:~:text=Quasi%2DStellar%20Objects%20(QSO)%20is%20class%20of%20objects%20beyond,radio%20sources%20that%20are%20starlike.)
  * Quasi-Stellar Objects (QSO) is class of objects beyond our Milky Way Galaxy that have a starlike visual appearance except that the optical spectrum has a large redshift. Quasars (quasi-stellar radio sources) are strong radio sources that are starlike. A QSO may have a strong radio source or may not have a strong radio source.
* [Voigt Profile](https://en.wikipedia.org/wiki/Voigt_profile)
  * The Voigt profile (named after Woldemar Voigt) is a probability distribution given by a convolution of a Cauchy-Lorentz distribution and a Gaussian distribution. It is often used in analyzing data from spectroscopy or diffraction.
  * **More info**
    * [VoIgt profile Parameter Estimation Routine (VIPER): H I photoionization rate at z < 0.5](https://academic.oup.com/mnras/article/467/3/3172/3062531)
