# readBT
A set of helper functions for reading in and working with BT-Settl synthetic spectra.
PHOENIX BT-Settl models can be downloaded from [France Allard's website](http://perso.ens-lyon.fr/france.allard/) or my grids are available on [my website](http://people.bu.edu/mveyette/phoenix/index.html).

The most useful functions are:

## readBT
Function to read in synthetic spectra. Returns two arrays: wavelenght in Angstroms and flux in erg s^-1 cm^-2 A^-1

To read in the full model as is:
```python
from readBT import readBT
wave, flam = readBT('lte3000-5.28-M-1.00-A+0.40-C+0.35-O+0.89.txt.xz')
print(wave)
>>> [1.000e+00 1.050e+00 1.100e+00 ... 9.985e+06 9.990e+06 9.995e+06]
print(flam)
>>> [4.12952405e-95 4.12952405e-95 4.12952405e-95 ... 6.09817664e-06 6.08555231e-06 6.07295412e-06]
```

This function can also convolved the model down to a desired resolution. However, make sure to set a wavelength range in microns withe the `waverange` keyword or else it will take a VERY long time to process the whole spectum. You can also use `npix` toset the desired sampling per resolution element of the output.
```python
from readBT import readBT
wave, flam = readBT('lte3000-5.28-M-1.00-A+0.40-C+0.35-O+0.89.txt.xz', R=2500, waverange=[2.0,2.5], npix=10)
```

## parseBT
This function parses the filename of a model (or list of filenames) and returns a record array with the parameters of the model(s).

```python
from readBT import parseBT
files = ['lte3000-5.28-M-1.00-A+0.40-C+0.35-O+0.89.txt.xz','lte3000-5.28-M-1.00-A+0.00-C+0.35-O+0.89.txt.xz']
params = parseBT(files)
params
>>> rec.array([(3000., 5.28, -1., 0.4, 0.35, 0.89, 0.),
               (3000., 5.28, -1., 0. , 0.35, 0.89, 0.)],
              dtype=[('teff', '<f8'), ('logg', '<f8'), ('mh', '<f8'), ('am', '<f8'), ('c', '<f8'), ('o', '<f8'), ('ano', '<f8')])
params['teff']
>>> array([3000., 3000.])
params['am']
>>> array([0.4, 0. ])
```
