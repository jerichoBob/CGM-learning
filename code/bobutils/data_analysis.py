def create_gaussian_fitter(wv, fx, xlims):
    ''' fit the absorption line, bounded by the evaluation limits xlims, with a gaussian '''
    import numpy as np
    from astropy.modeling import models, fitting

    delwv  = np.double(xlims[1])-np.double(xlims[0])
    amp_min = fx.min() # this should find where the flux bottoms out
    amp_cont = fx.max() # this should approximate the continuum
    mean_guess = np.mean(xlims) # this should be the center of the absorber
    stddev_guess = 0.5*(delwv) # this should be half of the width of the absorber
    
    # since we have an absorber, the gaussian fit drops down from the continuum
    g_init = (models.Const1D(amp_cont) +
              models.Gaussian1D(amplitude=(amp_min - amp_cont), mean=mean_guess, stddev=stddev_guess))
    
    fit_g = fitting.LevMarLSQFitter()
    g = fit_g(g_init, wv, fx)

    fit_info = fit_g.fit_info # fit_info contains the results from the last fit performed
    # print(f"fit_info: {fit_info}")
    
    if fit_info['param_cov'] is None:
        raise ValueError('gaussian_ew: The fit is not sensible! Check initial_guesses')

    return g, fit_g

def calculate_ew_original(g, wv, fit_g):
    """ 
    this implementations comes from linetools/analysis/utils.py - gaussian_ew() -- 
    I am assuming (without verification) that it produces correct results :D 
    """

    import numpy as np
    # Components: 
    #     [0]: <Const1D(amplitude=1.00002032)>
    #     [1]: <Gaussian1D(amplitude=-0.89579374, mean=3838.41347354, stddev=1.42832398)>
    # Parameters:
    #        amplitude_0         amplitude_1           mean_1            stddev_1     
    #     ------------------ ------------------- ------------------ ------------------
    #     1.0000203235364862 -0.8957937376523996 3838.4134735367384 1.4283239796305052
    # Area under curve of Gaussian is [amplitude*stddev*sqrt(2*pi)]

    amplitude = np.abs(g.parameters[1]) # amplitude of the gaussian
    stddev = g.parameters[3] # stddev of the gaussian
    print("Inside calculate_ew_original  ----- ")
    print(f"      amplitude: {amplitude}")
    print(f"         stddev: {stddev}")
    EW = amplitude * stddev * np.sqrt(2 * np.pi) #unitless
    # EW = EW * wv.unit #add the same unit as wv

    #error estimation
    x = amplitude
    y = stddev
    cov = fit_g.fit_info['param_cov'] #covariance matrix
    # still don't know if these are the correct indices into the covariance matrix
    EWsig = EW * np.sqrt(cov[0,0] / x**2 + cov[2,2] / y**2 + 2 * cov[0,2] / (x*y))
    return EW, EWsig

def calculate_ew_from_gaussian(g, xlims):
    """ a simple integration of the gaussian fit, but doesn't yet give us the error for EW """

    import numpy as np
    from scipy.integrate import quad
    
    ew = quad(g, xlims[0], xlims[1])[0]
    return ew