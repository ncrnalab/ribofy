import numpy as np
import pandas as pd

from scipy.stats import wilcoxon, binomtest, f
import statsmodels.formula.api as smf 
import statsmodels.api as sm
import statsmodels.tools.sm_exceptions as sme

from scipy.special import digamma,polygamma
from scipy.stats import nbinom

libmtspec = True
try:
    from mtspec import mtspec
except ModuleNotFoundError:
    libmtspec = False



def get_2D_matrix (psites):
    """converts from 1D array to 2D matrix"""
    
    mat = np.reshape (psites, (int (len(psites)/3), 3))
    
    return (mat)



def get_counts (psites):
    """counts number of p-sites in each frame"""

    mat = get_2D_matrix (psites)

    return ({
        'total' : np.sum(mat),
        'frame0' : np.sum(mat[:,0]),
        'frame1' : np.sum(mat[:,1]),
        'frame2' : np.sum(mat[:,2])
    })


def get_taper (psites, time_bandwidth = 3, ntapers = "default", nfft = "default"):
    """Performs multitaper analysis (as in ribotaper) with Ftest statistics for 1/3 frequency
    psites: 1D array (ORF length) with P-site counts ncodons
    returns: p-value
    """

    if sum (psites) == 0:
        return (np.nan)

    if nfft == "default":
        nfft = int(2 * 2**np.ceil(np.log2(len(psites))))

    if ntapers == "default":
        ntapers = int(2*time_bandwidth) - 1

    # Calculate the spectral estimation.
    spec, freq, jackknife, fstatistics, _ = mtspec(data=np.array(psites), delta = 1, time_bandwidth = time_bandwidth, number_of_tapers=ntapers, nfft=nfft, statistics=True, rshape=0)

    m = int(np.round (nfft/3))
    sf = f.sf (fstatistics[m],dfn=2,dfd=(2*ntapers)-2)
    return (sf)


def get_wilcox (mat):
    """
    Paired wilcoxon-test for frame0 > mean (frame1, frame2)
    mat: 2D matrix with shape (3, ncodons)
    returns: p-value
    """

    frame0 = mat[:,0]
    frame12 = np.mean (mat[:,1:3], axis=1)

    #wilcox_stat, wilcox_p = wilcoxon(frame0, frame12, alternative="greater") if not np.all (frame0-frame12==0) else (np.nan, np.nan)
    wilcox_stat, wilcox_p = wilcoxon(frame0 - frame12, alternative="greater") if not np.all (frame0-frame12==0) else (np.nan, np.nan)
    return (wilcox_p)



def get_binom (mat):
    """
    Perform binomial-test for n(frame0 > frame1 and frame0 > frame2). Adding random noise to reduce draw-bias, otherwise on-frame is max on draw
    mat: 2D matrix with shape (3, ncodons)
    returns: p-value
    """


    mat = mat + np.random.uniform(low=0.0, high=0.99, size=mat.shape)

    index_max = np.argmax (mat, axis=1)
    binom_p = binomtest (k=np.sum (index_max == 0), n=len(index_max), p=1/3, alternative="greater").pvalue if len (index_max) > 0 else np.nan
    return (binom_p)


def get_theta_md (y, limit=20, eps = np.finfo(float).eps**.25):
    """estimates theta for nb GLM - adapted from theta.md (MASS package, R)"""

    y = np.array (y)
    mu = np.mean (y)
    dfr = y.shape[0] - 2

    weights = np.ones (len(y))
    n    = np.sum(weights)    
    t0   = n/np.sum(weights * (y/mu - 1)**2)
    nmax = [np.max ([1,p]) for p in y]
    a    = 2 * np.sum(weights * y * np.log(nmax/mu)) - dfr

    it = 0
    idel = 1
    while (it + 1 < limit and np.abs(idel) > eps and not np.isnan (t0)):
        it = it+1
        t0 = np.abs(t0)
        tmp = np.log((y + t0)/(mu + t0))
        top = a - 2 * np.sum(weights * (y + t0) * tmp)
        bot = 2 * np.sum(weights * ((y - mu)/(mu + t0) - tmp))
        idel = top/bot
        t0 = t0 - idel
    
    if t0 <= 0 or np.isnan (t0) or np.isinf (t0):
        t0 = 1 # default alpha in statsmodels nb glm
            
    return (t0)


def get_theta_ml (y, limit = 10, eps = np.finfo(float).eps**.25, trace = False): 
    """estimates theta for nb GLM - adapted from theta.ml (MASS package, R)"""
    
    def score (n, th, mu, y, w):
        return (sum(w * (digamma(th + y) - digamma(th) + np.log(th) + 1 - np.log(th + mu) - (y + th)/(mu + th))))
    
    def info (n, th, mu, y, w):
        return (sum(w * (-polygamma(1, th + y) + polygamma(1, th) - 1/th + 2/(mu + th) - (y + th)/(mu + th)**2)))

    try:
        mu = np.mean (y)

        weights = np.ones ((len (y)))
        n = np.sum(weights)
            
        t0 = n/sum(weights * (y/mu - 1)**2)
        it = 0
        idel = 1
        
        if (trace): 
            print ("theta.ml: iter", it, "theta", t0)
        
        while (it < limit and abs(idel) > eps):
            
            t0 = abs (t0)
            i = info (n, t0, mu, y, weights)
            idel = score(n, t0, mu, y, weights) / i
            t0 = t0 + idel
            it = it+1
        
        if t0 <= 0 or np.isnan (t0) or np.isinf (t0):
            t0 = 1
            
        if it == limit and trace: 
            print ("iteration limit reached")
    
        return (t0)

    except ZeroDivisionError:
        return (1)

    


def convert_params(mu, theta):
    """
    Convert mean/dispersion parameterization of a negative binomial to the ones scipy supports
    See https://en.wikipedia.org/wiki/Negative_binomial_distribution#Alternative_formulations
    """
    r = theta
    var = mu + 1 / r * mu ** 2
    p = (var - mu) / var
    return r, 1 - p


def pmf(counts, mu, theta):
    """
    """
    return nbinom.pmf(counts, *convert_params(mu, theta))



def get_glm (mat, remove_outliers = False):
    """
    Fits a negative binomial GLM to the p-sites with a two-class frame feature (on or off-frame) and extracts the parameter for the frame coefficient. 
    mat: 2D matrix with shape (ncodons, 3)
    returns: p-value
    """

    df_glm = pd.DataFrame ({
        'counts' : mat.reshape (-1),
        'frame' : ['onframe', 'offframe', 'offframe'] * mat.shape[0]
    })

    try:

        if remove_outliers:

            theta_g = df_glm.groupby ("frame").agg ([np.mean, get_theta_ml])
            
            df_glm['pmf'] = pmf (df_glm.counts.values, theta_g.loc[df_glm.frame, ('counts','mean')], theta_g.loc[df_glm.frame, ('counts','get_theta_ml')])

            df_glm['adj_pmf'] = p_adjust_bh (df_glm.pmf)

            df_glm = df_glm[df_glm.adj_pmf >= 0.01]


        theta = get_theta_ml (df_glm.counts.values)        
        model = smf.glm(formula = "counts ~ frame", data=df_glm, family=sm.families.NegativeBinomial(alpha=1/theta)).fit()        
        glm_p = model.pvalues[1] # glm_ttest.pvalue

        # converting to one-tailed
        if model.params[1] > 0: #== max (model.params):
            glm_p_onetailed = glm_p/2
        else:
            glm_p_onetailed = 1-glm_p/2

        return (glm_p_onetailed)


    except sme.PerfectSeparationError:
        return (np.nan)
    except ValueError:        
        return (np.nan)
    except IndexError:
        return (np.nan)





def p_adjust_bh (p):
    """
    Benjamini-Hochberg p-value correction for multiple hypothesis testing.
    adapted from here: https://stackoverflow.com/questions/7450957/how-to-implement-rs-p-adjust-in-python to allow NaNs
    """
    p = np.asfarray(p)
    
    nna = ~np.isnan (p)
    q = np.empty ((len(p)))
    q[:] = np.nan
    pnna = p[nna]

    by_descend = pnna.argsort()[::-1]
    by_orig = by_descend.argsort()

    n = len(pnna) #[~np.isnan (p)])
    i = np.arange(len(pnna), 0, -1)
    q[nna] = np.minimum(1, np.fmin.accumulate((float (n)/i) * pnna[by_descend]))[by_orig]
    return q
        



def get_filtered_padj (s, pcol="p_glm", name="filtered_padj"):
    """
    Adapted from DESeq2; filtering by expression, the BH padjustment if performed solely on ORFs exceeding the expression threshold. 
    Then, the threshold that maximized number of rejections (i.e. significant ORFs) are used. In contrast to DESeq2, the maximization is
    not based on lowess regression, but simply the cutoff with max rejections (lowess implementation TODO).
    """
    
    filter=np.array (s['n'])
    p=np.array(s[pcol])
    nrows = s.shape[0]

    if nrows < 50:        
        s[name] = p_adjust_bh(p)   
        return (s)

    lq = np.mean(filter == 0)
    uq = .95 if lq < .95 else 1

    r = np.array (np.linspace(start=lq, stop=uq, num=50))

    cutoffs = np.quantile (filter, r)

    result = np.empty((nrows,len(cutoffs)))
    result[:] = np.nan

    for i in range (len (cutoffs)):
        
        use = filter >= cutoffs[i]    
        
        if (np.any(use)):
            
            use_p = p[use]        
            result[use, i] = p_adjust_bh(use_p)        


    numRej = np.sum (result < 0.05, axis=0)
    j = np.argmax(numRej)

    s[name] = result[:,j]
    return (s)
