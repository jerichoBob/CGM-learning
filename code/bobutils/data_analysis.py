def create_gaussian_fitter(wv, fx, sig, xlims):
    ''' fit the absorption line, bounded by the evaluation limits xlims, with a gaussian '''
    import numpy as np
    from astropy.modeling import models, fitting

    delwv  = np.double(xlims[1])-np.double(xlims[0])
    amp_guess = fx.min() # this should find where the flux bottoms out
    amp_disp = fx.max() # this should approximate the continuum
    mean_guess = np.mean(xlims) # this should be the center of the absorber
    stddev_guess = 0.5*(delwv) # this should be half of the width of the absorber
    
    # since we have an absorber, the gaussian fit drops down from the continuum
    g_init = (models.Const1D(amp_disp) +
              models.Gaussian1D(amplitude=(amp_guess - amp_disp), mean=mean_guess, stddev=stddev_guess))
    
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, wv, fx)

    fit_info = fit_g.fit_info # fit_info contains the results from the last fit performed
    print(f"fit_info: {fit_info}")
    
    if fit_info['param_cov'] is None:
        raise ValueError('gaussian_ew: The fit is not sensible! Check initial_guesses')

    return g, fit_g

def caculate_ew_original(g, wv, fit_g):
    """ comes from linetools/analysis/utils.py - gaussian_ew() -- assuming it is correct :D """
    import numpy as np
    # Area under curve of Gaussian is [amplitude*stddev*sqrt(2*pi)]
    EW = g.amplitude.value * g.stddev.value * np.sqrt(2 * np.pi) #unitless
    EW = EW * wv.unit #add the same unit as wv
    
    #error estimation
    cov = fit_g.fit_info['param_cov'] #covariance matrix
    x = g.parameters[0] # amplitude
    y = g.parameters[2] # stddev
    sigEW = EW * np.sqrt(cov[0,0] / x**2 + cov[2,2] / y**2 + 2 * cov[0,2] / (x*y))
    return EW, sigEW

def calculate_ew_from_gaussian(g, xlims):
    """ a simple integration of the gaussian fit, but doesn't yet give us the error for EW """

    import numpy as np
    from scipy.integrate import quad
    
    ew = quad(g, xlims[0], xlims[1])[0]
    return ew