import numpy as np
import scipy
from scipy import stats

def estimate_pi0(
    pv,
    lambdas=None,
    pi0=None,
    df=3,
    method="smoother",
    smooth_log_pi0=False,
    verbose=True,
):
    """Estimate pi0 based on the pvalues"""
    try:
        pv = np.array(pv)
    except:
        pv = pv.copy()
    assert pv.min() >= 0 and pv.max() <= 1, "p-values should be between 0 and 1"
    if lambdas is None:
        epsilon = 1e-8
        lambdas = scipy.arange(0, 0.9 + 1e-8, 0.05)
    if len(lambdas) > 1 and len(lambdas) < 4:
        raise ValueError(
            """if length of lambda greater than 1, you need at least 4 values"""
        )
    if len(lambdas) >= 1 and (min(lambdas) < 0 or max(lambdas) >= 1):
        raise ValueError("lambdas must be in the range[0, 1[")
    m = float(len(pv))

    pv = pv.ravel()  # flatten array
    if pi0 is not None:
        pass
    elif len(lambdas) == 1:
        pi0 = np.mean(pv >= lambdas[0]) / (1 - lambdas[0])
        pi0 = min(pi0, 1)
    else:
        # evaluate pi0 for different lambdas
        pi0 = [np.mean(pv >= this) / (1 - this) for this in lambdas]
        # in R
        # lambda = seq(0,0.09, 0.1)
        # pi0 = c(1.0000000, 0.9759067, 0.9674164, 0.9622673, 0.9573241,
        #         0.9573241 0.9558824, 0.9573241, 0.9544406, 0.9457901)
        # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F, spar=0)
        # predict(spi0, x=max(lambda))$y  --> 0.9457946
        # spi0 = smooth.spline(lambda, pi0, df=3, all.knots=F)
        # predict(spi0, x=max(lambda))$y  --> 0.9485383
        # In this function, using pi0 and lambdas, we get 0.9457946
        # this is not too bad, the difference on the v17 data set
        # is about 0.3 %
        if method == "smoother":
            if smooth_log_pi0:
                pi0 = np.log(pi0)
            # In R, the interpolation is done with smooth.spline
            # within qvalue. However this is done with default
            # parameters, and this is different from the Python
            # code. Note, however, that smooth.spline has a parameter
            # called spar. If set to 0, then we would get the same
            # as in scipy. It looks like scipy has no equivalent of
            # the smooth.spline function in R if spar is not 0
            tck = scipy.interpolate.splrep(lambdas, pi0, k=df)
            pi0 = scipy.interpolate.splev(lambdas[-1], tck)
            if smooth_log_pi0:
                pi0 = np.exp(pi0)
            pi0 = min(pi0, 1.0)
        elif method == "bootstrap":
            raise NotImplementedError
            """minpi0 = min(pi0)
            mse = rep(0, len(lambdas))
            pi0.boot = rep(0, len(lambdas))
            for i in range(1,100):
                p.boot = sample(p, size = m, replace = TRUE)
                for i in range(0,len(lambdas)):
                    pi0.boot[i] <- mean(p.boot > lambdas[i])/(1 - lambdas[i])
                mse = mse + (pi0.boot - minpi0)^2
            pi0 = min(pi0[mse == min(mse)])
            pi0 = min(pi0, 1)"""
        if pi0 > 1:
            if verbose:
                print("got pi0 > 1 (%.3f), setting it to 1" % pi0)
            pi0 = 1.0
    assert pi0 >= 0 and pi0 <= 1, "pi0 is not between 0 and 1: %f" % pi0
    return pi0


s = list(np.random.uniform(0, 1, 800))
t = []
for i in range(0, 200):
    a = np.random.randn(3) + 2
    b = np.random.randn(3)
    t.append(stats.ttest_ind(a, b).pvalue)

print(estimate_pi0(s + t))
