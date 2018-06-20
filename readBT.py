import math
import numpy as np
import scipy.interpolate
import scipy.signal
import time
from scipy.ndimage.filters import maximum_filter, uniform_filter
from findFile import findFile
import os
import astropy.io.fits as fits
import re
   
def parseBT(files):
    """
    Parses BT-settl file name.
    Input - filename (string), or array of files
    Output - recarray with parameters: teff, logg, mh, am, c, o, ano (if avail)
    """
    
    if np.size(files) == 1:
        files = [files]
        
    files = [os.path.basename(file) for file in files]        
    params = np.recarray(len(files), dtype=[('teff',float), ('logg',float), ('mh',float), ('am',float), ('c',float), ('o',float), ('ano',float)])

    for i, file in enumerate(files):
        match = re.match(r'lte(\d*)-(\d*.\d*)-M([+?,-?]\d*.\d*)-A([+?,-?]\d*.\d*)-C([+?,-?]\d*.\d*)-O([+?,-?]\d*.\d*)-ANO([+?,-?]\d*.\d*)\..*', file)
        if match:
            params[i] = np.array(match.expand(r'\1 \2 \3 \4 \5 \6 \7').split(), dtype=float)
        else:
            match = re.match(r'lte(\d*)-(\d*.\d*)-M([+?,-?]\d*.\d*)-A([+?,-?]\d*.\d*)-C([+?,-?]\d*.\d*)-O([+?,-?]\d*.\d*)\..*', file)
            params[i] = np.array(match.expand(r'\1 \2 \3 \4 \5 \6 0').split(), dtype=float)
        
    return params

def parseBT_Allard(file):
    """
    Parses France Allard's BT-settl file name.
    Input - filename (string)
    Output - float array with values for [teff, logg, MH, A]
    """
    
    name  = os.path.split(file)[-1]
    
    fteff = name.find('-')
    flogg = np.min([name.find('-', fteff+1), name.find('+', fteff+1)])
    fMH   = name.find('a', flogg+1)
    fA    = name.find('.BT', fMH+1)
    
    iteff = name.find('lte') + 3
    ilogg = fteff + 1
    iMH   = flogg
    iA    = fMH   + 1
    
    return np.array([name[iteff:fteff]+'00', name[ilogg:flogg], name[iMH:fMH], name[iA:fA]], dtype='float')
    
def parseBT_old(files):
    """
    Parses BT-settl file name. DEPRECATED - specific to older versions of PHOENIX output 
    Input - filename (string), or array of files
    Output - recarray with parameters: teff, logg, mh, am, c, o
    """
    
    if np.size(files) == 1:
        files = [files]
    
    params = np.recarray(len(files),
               dtype=[('teff',float), ('logg',float), ('mh',float), ('am',float), ('c',float), ('o',float)])
    for i, file in enumerate(files):
    
        name  = os.path.split(file)[-1]
        
        fteff = name.find('-')
        flogg = name.find('-M')
        fMH   = name.find('-A')
        fA    = name.find('-C')
        fC    = name.find('-O')
        
        Odot  = name.find('.', fC+2)
        nexti = np.array([name.find('.', Odot+1),name.find('_', Odot+1)], dtype='float')
        nexti[nexti < 0] = np.nan
        fO    = int(np.nanmin(nexti))
        
        iteff = name.find('lte') + 3
        ilogg = fteff + 1
        iMH   = flogg + 2
        iA    = fMH   + 2
        iC    = fA    + 2
        iO    = fC    + 2
        
        params[i] = np.array([name[iteff:fteff], name[ilogg:flogg], name[iMH:fMH], name[iA:fA], name[iC:fC], name[iO:fO]], dtype='float')
    
    return params
    
def readBT(file, R=None, npix=3.0, samp=2.5E7, waverange=None, air=False, verbose=False,
           wave=None, flam=None):
    """
    Returns wavelength in Angstroms and flux in erg s^-1 cm^-2 A^-1
    
    Arguments
    file - name of file to read from
    R    - resolution to convolve spectrum to, default to not change resolution
    npix - number of pixels per resolution element of output spectrum, most spectrographs use 2-4
    samp - sampling to interpolate wavelength grid to for convolution, must be higher than original sampling
    waverange - two-element array with start and end wavelengths of output ***IN MICRONS***
    air  - set to True to convert output wavelengths from vacuum to air
    verbose - set to True to print timing info
    wave - input wavelength array to process instead of reading in from file
    flam - input flux array to process instead of reading from file
    """
    
    #Default to no broadening
    if R is None:
        asIs = True
    else:
        asIs = False
    
    if waverange is None and not asIs:
        raise ValueError("You really don't want to run this without setting a waverange or reading in asIs. "
              "Your computer will probably freeze.")
    
    if verbose:
        stime = time.time()
        print('Starting File: ' + str(file))
    
    #NOT IMPLEMENTED YET
    # if vsini == 0:
        # vsini = 3.0 # km/s
        # rotBroad = 0
    # else:
        # rotbroad = 1
    
    #Read in data unless passed
    if verbose: stime1 = time.time()
    if wave is None or flam is None:
    
        readin = []
        if os.path.splitext(file)[1] == '.xz':
            import lzma
            openit = lzma.open
        else:
            openit = open
        with openit(file, 'rt') as file1:
            for line in file1:
                subbed = re.sub(r'(\d)-(\d)', r'\1 -\2', line).replace('D','E')
                readin.append(subbed.split())
            
        wave = []
        flam = []    
        for line in readin:
            wave.append(float(line[0]))
            flam.append(float(line[1]))
            
        DF = -8.0
            
        wave = np.array(wave)
        flam = 10.0**(np.array(flam) + DF)

        wave, uwave = np.unique(wave, return_index=True)
        flam = flam[uwave]
        
    if verbose: print('Reading file took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Return without processing if asIs is True
    if asIs:
        if waverange is not None:
            keep = (waverange[0] <= wave/1e4) & (wave/1e4 <= waverange[1])
            wave = wave[keep]
            flam = flam[keep]
    
        if verbose:
            etime = time.time()
            print('Finished in ' + str(int(round(etime-stime))) + ' seconds.')
        return wave, flam
    
    ss = 10.0 * npix
    
    #Default to full data range
    if waverange is None: waverange = np.array([wave[ss+1], wave[-1*ss - 2]]) * 1e-4
    
    #Interpolate to higher sampling
    if verbose: stime1 = time.time()
    sampBump = math.ceil(samp / R / npix)
    iWaveStart = waverange[0] * 1e4
    iWaveEnd = waverange[1] * 1e4
    #if np.min(wave) > iWaveStart or np.max(wave) < iWaveEnd:
    #    print('Oh no! Not enough wavelength range. File: ' + file)
    iWaveFull = iWaveStart * np.exp( np.arange( samp * np.log(iWaveEnd/iWaveStart) ) / samp )
    w = (iWaveFull > np.min(wave)) & (iWaveFull < np.max(wave))
    iWaveOffset = np.where(w)[0][0]
    iWave = iWaveFull[w]
    ww = (wave >= 0.9*iWaveStart) & (wave <= 1.1*iWaveEnd)
    tck = scipy.interpolate.splrep(wave[ww], flam[ww], s=0)
    iflam = np.array(scipy.interpolate.splev(iWave, tck, der=0))
    if verbose: print('Interpolation took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Convolve with spectrograph resolution
    if verbose: stime1 = time.time()
    fwhm = samp / R
    stdv = fwhm / ( 2.0 * math.sqrt( 2.0 * math.log(2.0) ) )
    #Create Kernel and normalize
    lsf  = np.array(scipy.signal.gaussian(int(4.0 * fwhm), stdv))
    lsf /= lsf.sum()
    #Convolve
    specflam = scipy.signal.fftconvolve(iflam, lsf, mode='same')
    if verbose: print('Convolution took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    #Integrate over sampBump to get npix per resolution element
    if verbose: stime1 = time.time()
    i1Full  = ( float(sampBump) * np.arange( np.floor(
              (np.size(iWaveFull)-1.0) / sampBump ) ) )        # beginning index of each integral
    i1sampFull = i1Full + sampBump
    w = (iWaveFull[i1Full.astype('int').tolist()] > np.min(wave)) & (iWaveFull[i1sampFull.astype('int').tolist()] < np.max(wave))
    i1     = (i1Full.astype('int')[w] - iWaveOffset).tolist()
    i1samp = (i1sampFull.astype('int')[w] - iWaveOffset).tolist()
    trapA   = ( ( iWave[1:] - iWave[0:-1] ) *
              0.5 * ( specflam[0:-1] + specflam[1:] ) )         # trapezoidal area under each samp point
    intFlux = np.sum((trapA[int(i1[0]):int(i1[-1]+sampBump)]
                      ).reshape(len(i1), sampBump), axis=1)     # sum over sampBump points for each pixel
    #import pdb ; pdb.set_trace()
    delLam  = iWave[i1samp] - iWave[i1]   # delta lambda
    intWave = ( (iWave[i1samp]**2  - iWave[i1]**2)
              / ( 2.0 * delLam ) )                              # mean wavelength of each pixel
    intflam = intFlux / delLam
    nSpec   = len(intWave)
    if verbose: print('Integration took ' + str(int(round(time.time() - stime1))) + ' seconds')
    
    if air: intWave /= (1.0 + 2.735182e-4 + 131.4182 / intWave**2 + 2.76249e8 / intWave**4)
    
    if verbose:
        etime = time.time()
        print('Finished in ' + str(int(round(etime-stime))) + ' seconds.')
    
    return intWave, intflam
    
def readHusser(waveFile, flamFile, **kwargs):
    """
    Runs readBT on the Husser files.
    """
    return readBT('', wave=fits.getdata(waveFile, 0), flam=1e-8*fits.getdata(flamFile, 0), **kwargs)
    
def getCont(flam, ss=65.0):
    """
    Returns continuum of spectrum.
    
    Arguments
    flam - spectrum
    ss - scan size in pixels. Usually 32 delLam works. Best to make it odd
    """
    
    #Scanning max filter
    cont = maximum_filter(flam, ss)
    
    #Scanning mean filter of twice the width
    cont = uniform_filter(cont, 2.0*ss)
    
    return cont

def getMag(lam, flam):
    """
    Returns 2MASS magnitudes
    
    Arguments
    lam - wavelength in angstroms
    flam - spectrum in erg s^-1 cm^-2 A^-1
    """
    
    #2MASS Zero Points (W cm^-2 um^-1 to erg s^-1 cm^2 A^-1)
    twomass0 = np.array([3.129e-13, 1.133e-13, 4.283e-14]) * 1e-4 * 1e7
    
    filterFiles = np.array(['2MASS_J.txt', '2MASS_H.txt', '2MASS_K.txt'])
    
    mags = [-99.0, -99.0, -99.0]
    
    for i, file in enumerate(filterFiles):
        
        #Read in filter data
        filtData = np.genfromtxt(findFile(file))
        filtLam = filtData[:,0]
        filtR   = filtData[:,1]
        
        #Change wavelength to A
        filtLam *= 1e4
        
        #Only use wavelengths within band
        good = (lam > np.min(filtLam)) & (lam < np.max(filtLam))
        
        #Interpolate filter response to match input
        intFiltR = np.interp(lam[good], filtLam, filtR)
      
        #Flux in the band
        filtFlux = np.trapz(flam[good]*intFiltR, lam[good]) / np.trapz(intFiltR, lam[good]) 
      
        #Return mag of input spectrum
        mags[i] = -2.5 * np.log10(filtFlux/twomass0[i])
        
    return mags
    