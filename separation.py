import astropy.io.fits as fits
import numpy as np
import matplotlib.pyplot as plt
import figsizemargins
import matplotlib.patches as patches
import astropy.io.fits as io
import glob

def plotSpec(wave, spec, fit, err):
    ratio = fit/spec
    notnan = ~np.isnan(ratio)
    fig, ax = figsizemargins.makeFigureInstance(x = 1, y=2, figureSize=(32,2))
    ax[0].plot(wave, spec, lw=0.5, label='data', zorder=10)
    ax[0].plot(wave, fit, lw=0.5, label='model')
    ax[1].plot(wave, err, lw=0.5, label='err')
    ax[1].plot(wave[notnan], ratio[notnan], lw=0.5, label='model/data')

    for axis in ax:
        axis.set_xlabel('wavelength')
        axis.legend(loc='best')
        axis.set_xlim(np.min(wave), np.max(wave))

    ax[0].set_ylim(0, 1.1)
    ax[1].set_ylim(np.min(ratio[notnan]), np.max(ratio[notnan]), 1.5)
    #plt.yscale('log')
    fig.savefig('test.pdf', dpi=400)
    plt.close(fig)


def getSpectra(objID, wave=True, nPixel=8575):
    dataDir = '/mnt/ceph/users/landerson/l30e.2/'
    filename = 'aspcapStar-r6-l30e.2-' + objID + '.fits'
    out = glob.glob(dataDir + '*/' + filename)
    try:
        filen =  out[0]
        hdu = fits.open(filen)
        spec = hdu[1].data
        err = hdu[2].data
        fit = hdu[3].data
        nPixel = len(spec)
        if wave:
            wave = 10.**(hdu[1].header['CRVAL1'] + np.arange(nPixel)*hdu[1].header['CDELT1'])
            return spec, err, fit, wave
        else:
            return spec, err, fit
    except IndexError:
        print filename, out
        foo = np.zeros(nPixel)
        foo.fill(np.nan) 
        return foo, foo, foo

plot=False
data = fits.getdata('allStar-l30e.2.fits')

good = data['LOGG'] != -9999.0

#define inital box
Lmin = 2
Lmax = 3.2
Tmin = 4200
Tmax = 5500
rc = (data['TEFF'] > Tmin) & (data['TEFF'] < Tmax) & (data['LOGG'] > Lmin) & (data['LOGG'] < Lmax)

#define zoom of inital box
LminZ = 2.38
LmaxZ = 2.42
TminZ = 4650
TmaxZ = 4700
rcZoom = (data['TEFF'] > TminZ) & (data['TEFF'] < TmaxZ) & (data['LOGG'] > LminZ) & (data['LOGG'] < LmaxZ)

if plot:
    fig, ax = figsizemargins.makeFigureInstance(x=3, y=1, wspace=0.75)
    ax[0].plot(data['TEFF'][good], data['LOGG'][good], 'ko', ms=1, rasterized=True, alpha=0.05, zorder=0)
    ax[0].add_patch(patches.Rectangle((Tmin, Lmin), Tmax-Tmin, Lmax-Lmin, fill=False, linewidth=1, zorder=1, color='#1f77b4'))
    ax[1].scatter(data['TEFF'][rc], data['LOGG'][rc], c='k', s=1, rasterized=True, alpha=0.05, zorder=0)
    ax[1].add_patch(patches.Rectangle((TminZ, LminZ), TmaxZ-TminZ, LmaxZ-LminZ, fill=False, linewidth=1, zorder=1, color='#1f77b4'))
    ax[2].scatter(data['TEFF'][rcZoom], data['LOGG'][rcZoom], c = 'k', s=1, alpha=0.1, zorder=0)

    for axis in ax:
        axis.invert_xaxis()
        axis.invert_yaxis()
        axis.set_xlabel('Teff [K]')
        axis.set_ylabel('log g')
        #axis.grid(zorder=2)
    fig.savefig('apogeeTeffLogg.pdf', dpi=400)
    plt.close(fig)

try: 
    specMatrix = np.load('specMatrix.npy')
    print 'read in npy matrix'

except IOError:

    nstars = np.sum(rcZoom)
    nPixels = 8575
    specMatrix = np.zeros((nstars, nPixels))

    for i, index in enumerate(np.where(rcZoom)[0]):
        objID = data['APSTAR_ID'][index].split('.')[-1]
        spec, err, fit = getSpectra(objID, wave=False)
        specMatrix[i] = spec
        if plot:
            plotSpec(wave, spec, fit, err)

    np.save('specMatrix', specMatrix)

print 'calculating SVD'
U ,s, V = np.linalg.svd(specMatrix)

print 'length of rcZoom: ', np.sum(rcZoom)
