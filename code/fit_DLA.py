"""fit_DLA.py
    Rongmon's UI/script to fit a DLA profile to a spectrum.
"""

import matplotlib
matplotlib.use('TkAgg')
from linetools.spectra.xspectrum1d import XSpectrum1D
import numpy as np
import matplotlib.pyplot as plt 
from astropy.io import ascii
from scipy.interpolate import splrep,splev
import sys
import os
from vfit import rb_voigt as v
import pandas
# Fitting with a lsqcurve fitting routine
from scipy.optimize import curve_fit

def model_profile(theta,wave,FWHM):
    theta=[theta[0],theta[1],theta[2],1.]
    model_const=pandas.DataFrame(columns=('wave0','gamma','f','FWHM'))
    # Enter for the First clump
    #         lam0,     gamma,      fval,    FWHM[km/s]
    model_const.loc[0]=[1215.6701,6.265E8,0.416400 ,FWHM]
    flx=v.full_voigt_profile(theta,wave,model_const.loc[:].values)
    return flx

def model_test(wave,*params):
    theta=np.zeros((len(params)))
    for i in range(0, len(params)):
        theta[i] = params[i]
    FWHM=75.
    flx=model_profile(theta,wave,FWHM)
    return flx


def onclick(event):
    toolbar = plt.get_current_fig_manager().toolbar
    if event.button==1 and toolbar.mode=='':
        window = ((event.xdata-5)<=wave) & (wave<=(event.xdata+5))
        y = np.median(flux[window])
        plt.plot(event.xdata,y,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
    plt.draw()

def onpick(event):
    # when the user clicks right on a continuum point, remove it
    if event.mouseevent.button==3:
        if hasattr(event.artist,'get_label') and event.artist.get_label()=='cont_pnt':
            event.artist.remove()

def ontype(event):

    if event.key=='enter':
        cont_pnt_coord = []
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='cont_pnt':
                cont_pnt_coord.append(artist.get_data())
            elif hasattr(artist,'get_label') and artist.get_label()=='continuum':
                artist.remove()
        cont_pnt_coord = np.array(cont_pnt_coord)[...,0]
        sort_array = np.argsort(cont_pnt_coord[:,0])
        x,y = cont_pnt_coord[sort_array].T
        spline = splrep(x,y,k=3)
        continuum = splev(wave,spline)
        plt.plot(wave,continuum,'r-',lw=2,label='continuum')
        theta=[logN[0],dopper_b[0],vel[0]]
        #flx_mdl=model_profile(theta,wave,FWHM)
        #plt.step(wave,flx_mdl*continuum,'g--')
        #lb=[20.,10.,-10.]
        #ub=[22.,50.,10.]
        #bounds=[lb,ub]     
        #popt, pcov = curve_fit(model_test, wave, flux/continuum, p0=theta,sigma=error/continuum,bounds=bounds)
        #fit = model_profile(popt,wave,FWHM)
        fit = model_profile(theta,wave,FWHM)
        plt.plot(wave,fit,'k--',linewidth=1)
        plt.plot(wave,fit*continuum,'r-')
        #plt.title('Best Fit logN ='+np.str(popt[0]))
        plt.title('Best Fit logN ='+str(logN[0]))

    elif event.key=='N':
        logN[0] +=0.05
        print('Best Fit logN ='+str(logN[0]))
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                # Extract fitted continuum
                continuum = artist.get_data()[1]
                plt.cla()
                plt.plot(wave,continuum,'r-',lw=2,label='continuum')
                plt.step(wave,flux,'b-',linewidth=1)
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.ylim([yr[0],yr[1]])
                plt.xlim([xr[0], xr[1]])
                theta=[logN[0],dopper_b[0],vel[0]]
                fit = model_profile(theta,wave,FWHM)
                plt.plot(wave,fit,'k--',linewidth=1)
                plt.plot(wave,fit*continuum,'r-')
                #plt.title('Best Fit logN ='+str(popt[0]))
                plt.title('Best Fit logN ='+str(logN[0]))

  
    elif event.key=='n':
        logN[0] -=0.05
        print('Best Fit logN ='+str(logN[0]))
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                # Extract fitted continuum
                continuum = artist.get_data()[1]
                plt.cla()
                plt.plot(wave,continuum,'r-',lw=2,label='continuum')
                plt.step(wave,flux,'b-',linewidth=1)
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.ylim([yr[0],yr[1]])
                plt.xlim([xr[0], xr[1]])
                theta=[logN[0],dopper_b[0],vel[0]]
                fit = model_profile(theta,wave,FWHM)
                plt.plot(wave,fit,'k--',linewidth=1)
                plt.plot(wave,fit*continuum,'r-')
                #plt.title('Best Fit logN ='+np.str(popt[0]))
                plt.title('Best Fit logN ='+str(logN[0]))



    elif event.key=='V':
        dopper_b[0] +=5.
        print('Doppler b ='+str(dopper_b[0]))
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                # Extract fitted continuum
                continuum = artist.get_data()[1]
                plt.cla()
                plt.plot(wave,continuum,'r-',lw=2,label='continuum')
                plt.step(wave,flux,'b-',linewidth=1)
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.ylim([yr[0],yr[1]])
                plt.xlim([xr[0], xr[1]])
                theta=[logN[0],dopper_b[0],vel[0]]
                fit = model_profile(theta,wave,FWHM)
                plt.plot(wave,fit,'k:',linewidth=1)
                plt.plot(wave,fit*continuum,'r-')
                #plt.title('Best Fit logN ='+np.str(popt[0]))
                plt.title('Doppler b ='+str(dopper_b[0]))

  
    elif event.key=='v':
        dopper_b[0] -=5.0
        print('Doppler b ='+str(dopper_b[0]))
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':
                # Extract fitted continuum
                continuum = artist.get_data()[1]
                plt.cla()
                plt.plot(wave,continuum,'r-',lw=2,label='continuum')
                plt.step(wave,flux,'b-',linewidth=1)
                plt.xlabel('Wavelength')
                plt.ylabel('Flux')
                plt.ylim([yr[0],yr[1]])
                plt.xlim([xr[0], xr[1]])
                theta=[logN[0],dopper_b[0],vel[0]]
                fit = model_profile(theta,wave,FWHM)
                plt.plot(wave,fit,'k--',linewidth=1)
                plt.plot(wave,fit*continuum,'r-')
                #plt.title('Best Fit logN ='+np.str(popt[0]))
                plt.title('Doppler b ='+str(dopper_b[0]))




    # when the user hits 'r': clear the axes and plot the original spectrum
    elif event.key=='r':
        plt.cla()
        plt.step(wave,flux,'b-',linewidth=1)
        plt.xlabel('Wavelength')
        plt.ylabel('Flux')
        plt.ylim([yr[0],yr[1]])
        plt.xlim([xr[0], xr[1]])



    # when the user hits 'b': selects a handpicked x,y value
    elif event.key=='b':
        plt.plot(event.xdata,event.ydata,'ro',ms=5,picker=5,label='cont_pnt',markeredgecolor='k')
        plt.draw()

        # when the user hits 'w': if the normalised spectrum exists, write it to a
    # file.
    elif event.key=='w':
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='continuum':#'normalised':
                continuum = artist.get_data()[1]
                outfilename=filepath+os.path.splitext(filename)[0]+'_norm.fits'
                from astropy.table import Table, Column, MaskedColumn
                table=Table([wave,flux,continuum],names=['wave','flux','cont'])
                table.write(outfilename,format='fits')
                break

  


    # At any time pressing q means graceful exit
    elif event.key=='q':
        quit_index=0;
        for artist in plt.gca().get_children():
            if hasattr(artist,'get_label') and artist.get_label()=='normalised':
                quit_index=1
            if quit_index==1:
                plt.close()
                print('Interactive Contunuum Normalization Done.')
                print('Hope you remembered to save the fit by pressing w!')
                print('Good Bye!')
                break
            else:
                plt.close()
                print('Quitting without normalizing. Moving along.....')
                break




    plt.draw()




if __name__ == "__main__":
    # Get the filename of the spectrum from the command line, and plot it
   
    from astropy.io import fits
    #dat=fits.open("/Users/bordoloi/Dropbox/Mage_COSMOS/Data/QSO-819402_full.fits")
    #data=dat[1].data
    #wave=np.array(data['wave'][0])
    #flux=np.array(data['flux'][0])
    #error=np.array(data['error'][0])
    #hdu = fits.open('/Users/bordoloi/Dropbox/KCWI/LensedArc/Extracted_Spectra/arclens_16_3_1_1.fits')

    #filename='SDSSJ140505.77+470441.1_coadd_G140L_final_lpALL.fits'#'J1148_HIRES.fits' #'test3.fits'
    #filepath='/Users/bordoloi/Dropbox/Proposals/2019/Keck/2019B/LensedArc/j2222_data/'
    #filepath='/Users/bordoloi/Dropbox/COS-Blue/Targets/SDSSJ140505.77+470441.1_target/'#'/Users/bordoloi/Dropbox/Proposals/2019/Keck/2019A/J1148_Observing/J1148/'
    

    filepath='/Users/bordoloi/Dropbox/Research/COS_GC_outflow/rb_spec_test/1H1613-097/BW_reduction/'
    filename='1H1613-097-G130M'
    sp = XSpectrum1D.from_file(filepath+filename)  
    #hdu = fits.open(filepath+filename)

    '''
    hdu_hdr=hdu[0].header
    crval3 = hdu_hdr['CRVAL1']
    crpix3 = hdu_hdr['CRPIX1']
    cd3_3 = hdu_hdr['CDELT1']
    wavedim = hdu_hdr['NAXIS1']
    # Do it
    wave = 10.**(crval3 + (crpix3 + np.arange(0, wavedim, 1.0)) * cd3_3)
    



    flux = hdu[1].data[0]
    error=hdu[1].data[0]
    wave=hdu[1].data[0]
    hdu.close()
    '''
    #data=hdu[1].data
    #wave=data['wave']
    #flux=data['flux']
    #error=data['error']

    wave=sp.wavelength.value
    flux=sp.flux.value 
    error=sp.sig.value
    zabs=0.#0.357#6.419#2.3000527#2.54366#2.41273#2.208483
    wave=wave/(1.+zabs)

    xr=[1190.,1240.]
    yr=[0,3]
    q=np.where((wave >xr[0]) & (wave <xr[1]) & np.isfinite(wave) )
    cont1=np.median(np.median(flux[q]))
    wave=wave[q]
    flux=flux[q]/cont1#/np.median(flux[q])
    #error=error[q]/cont1#/np.median(flux[q])

    FWHM=np.round(120.)
    logN=np.array([20.])
    dopper_b=np.array([30.])
    vel=np.array([0.])



    print("FWHM of the spectrum ="+str(FWHM)+" km/s")





    spectrum, = plt.step(wave,flux,'b-',label='spectrum',linewidth=1)
    plt.xlabel('Wavelength')
    plt.ylabel('Flux')
    plt.ylim([yr[0],yr[1]])
    plt.xlim([xr[0], xr[1]])


    # Connect the different functions to the different events
    plt.gcf().canvas.mpl_connect('key_press_event',ontype)
    plt.gcf().canvas.mpl_connect('button_press_event',onclick)
    plt.gcf().canvas.mpl_connect('pick_event',onpick)
    plt.show() # show the window
