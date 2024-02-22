import numpy as np
import scipy.integrate as integrate
import scipy.interpolate as interpolate
import numpy.random as random
import scipy.stats as stats

class PDF(object):
    def __init__(self, pdf, low, high):
        self.pdf = pdf

        x = np.linspace(low, high, 150)

        y_pdf = interpolate.UnivariateSpline(x, pdf(x))
        
        y = np.zeros(len(x))

        for i in range(1, len(y)):
            #y[i] = integrate.quad(pdf, low, x[i])[0]
            y[i] = y_pdf.integral(low, x[i])
            y[i] = max(y[i], y[i-1])
        y /= y[-1]
        
        self.low, self.high = low, high
        self.cdf = interpolate.interp1d(x, y, "linear")
        ok = np.ones(len(x))
        ok[1:] = y[1:] != y[:-1]
        self.inverse_cdf = interpolate.interp1d(y, x, "linear")

    def mean(self):
        return (integrate.quad(lambda x: self.pdf(x)*x, self.low, self.high)[0]/
                integrate.quad(lambda x: self.pdf(x), self.low, self.high)[0])

    def sigma(self):
        x = self.mean()
        x2 = (integrate.quad(lambda x: self.pdf(x)*x*x, self.low, self.high)[0]/
              integrate.quad(lambda x:self.pdf(x), self.low, self.high)[0])
        return np.sqrt(x2 - x*x)

    def sample(self, n):
        y = random.random(n)
        return self.inverse_cdf(y)

def gaussian_coupala_sample(x, y_cdf, rho):
    n = len(x)

    r = random.multivariate_normal([0, 0], [[1, rho], [rho, 1]], n)
    xr, yr = r[:,0], r[:,1]
    i_xr = stats.rankdata(xr, "ordinal")-1
    i_yr = stats.rankdata(yr, "ordinal")-1

    x_order = np.argsort(x)
    x_idx = np.arange(n, dtype=int)
    x_orig_idx = x_idx[x_order]

    y_low, y_high = y_cdf[1,0], y_cdf[1,-1]
    y_range = y_high - y_low
    inverse_y_cdf = interpolate.interp1d(y_cdf[1,:], y_cdf[0,:])
    y_rand = inverse_y_cdf(random.random(len(x))*y_range + y_low)
    y_sort = np.sort(y_rand)
    y = y_sort[i_yr]

    xr_order = np.argsort(xr)

    out = np.zeros(n)
    out[x_orig_idx] = y[xr_order]

    return out
