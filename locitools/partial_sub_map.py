import numpy as np
from scipy import ndimage
from scipy import interpolate
import pyfits as pyf

def partial_sub_map(r1, r2, theta1, theta2, sub_coefs, pa):

    """
    Function partial_sub_map turns the fractional flux loss coefficients
    computed by loci_sub into a 2D map.  The function averages values
    azimuthally and linearly interpolates radially.  It takes six
    arguments:
    1.  r1:  1D array, the minimum radius of each subtraction region
    2.  r2:  1D array, the maximum radius of each subtraction region
    3.  theta1:  1D array, the minimum angle of each subtraction
        region.  In the range [0, 2\pi).
    4.  theta2:  1D array, the maximum angle of each subtraction
        region.  In the range [0, 2\pi).
    5.  sub_coefs:  2D array, the fractional flux loss coefficients.
        The second dimension runs over frame ID.  
    6.  pa:  1D array, the position angles of the frames.

    partial_sub_map returns a 2D array, a map, of fractional flux loss.    
    
    """

    nframes = pa.shape[0]
    rindex = np.ones(r1.shape, int)
    for i in range(1, r1.shape[0]):
        if r1[i] == r1[i - 1]:
            rindex[i] = rindex[i - 1]
        else:
            rindex[i] = rindex[i - 1] + 1
    rmax = r2[-1]
    
    ###################################################################
    # We'll define our new array to have rmax*pi theta elements.
    # This equates to every-other-pixel spacing at the outer edge.
    ###################################################################

    ntheta = int(rmax * np.pi)
    vals = np.zeros((rindex[-1] + 1, ntheta), np.float32)
    newtheta = np.linspace(0, 2 * np.pi, ntheta)
    new_r = np.zeros(rindex[-1])

    i1 = 0
    i2 = 0

    ###################################################################
    # Go over all of the regions at a given radius, then step forward.
    ###################################################################

    for i in range(rindex[-1]):
        i1 = i2
        i2 = max(np.extract(rindex == rindex[i1], np.arange(r1.shape[0]))) + 1

        phi = np.zeros(((i2 - i1) * 3))
        y = np.zeros(phi.shape)
        dphi = theta2[i1] - theta1[i1]
        new_r[i] = 0.5 * (r2[i1] + r1[i1])
        
        ###################################################################
        # Add the closest flux loss point to each position in theta.
        # Add the position angle and shift by multiples of 2pi to make the
        # result fall in the range [0, 2pi).  Finally, normalize to get
        # the average fractional flux loss at each position.
        ###################################################################

        for j in range(nframes):
            di = i2 - i1
            phi[di:di * 2] = 0.5 * (theta1[i1:i2] + theta2[i1:i2]) + pa[j]
            
            offset = phi[di] - (phi[di] % (2 * np.pi))
            phi[di:di * 2] -= offset
            #print phi[di], phi[-1]
            phi[:di] = phi[di:di * 2] - 2 * np.pi
            phi[di * 2:] = phi[di:di * 2] + 2 * np.pi

            y[:di] = sub_coefs[i1:i2, j]
            y[di:2 * di] = y[:di]
            y[2 * di:] = y[:di]

            tck = interpolate.interp1d(phi, y, kind='nearest')
            vals[i] += tck(newtheta)

        vals[i] /= nframes

    ###################################################################
    # Now interpolate onto a 2D array to compute a map.  The radial
    # component will be linearly interpolated between region midpoints.
    ###################################################################

    x = np.arange(2 * int(rmax) + 1) - int(rmax)
    #x = np.arange(1501) - 750
    x, y = np.meshgrid(x, x)
    r = np.sqrt(x**2 + y**2)
    theta = np.arctan2(y, x) + np.pi

    indx_theta = theta * (ntheta - 1) / (2 * np.pi)
    tck = interpolate.interp1d(new_r, np.arange(new_r.shape[0]), 
                               bounds_error=False, fill_value=rindex[-1])
    indx_r = tck(r)

    out_im = ndimage.map_coordinates(vals, [indx_r, indx_theta], order=1)
    np.putmask(out_im, out_im == 0, np.nan)

    return out_im


