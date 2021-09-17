from base_cwas import *
import code
from pytest_bdd import scenario, given, when, then

@given('a connectir-based MDMR folder "{mdir}" looking at the factor "{factor}"')
def step(context, mdir, factor):
    # Distances
    sfile         = op.join(op.dirname(mdir), "subdist.desc")
    bigmemory     = importr('bigmemory')
    dmats         = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % sfile))
    context.nobs  = np.sqrt(dmats.shape[0])
    context.nvoxs = dmats.shape[1]
    context.dmats = dmats.T.reshape(context.nvoxs, context.nobs, context.nobs)
    # Fperms
    ffile         = op.join(mdir, "fperms_%s.desc" % factor)
    context.r_F_perms = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % ffile))
    # Fstats
    context.r_Fs  = context.r_F_perms[0,:]
    # Pvalues
    pfile         = op.join(mdir, "pvals.desc")
    pvals         = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % pfile))
    if len(pvals.shape) > 1:
        mfile     = op.join(mdir, "modelinfo.rda")
        robjects.r("load(%s)" % mfile)
        fnames    = np.array(robjects.r("names(modelinfo$H2s)"))
        fi        = (fnames == factor).nonzero()[0]
        pvals     = pvals[:,fi]
    context.r_ps  = pvals
    # Permutations
    pfile         = op.join(mdir, "perms_%s.desc" % factor)
    perms         = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % pfile))
    context.perms = (perms.T - 1)[1:,:]

@given('{nperms} permutations')
def step(context, nperms):
    context.perms = int(nperms)

@given('regressors from "{rfile}" with columns of interest {cols}')
def step(context, rfile, cols):
    context.regressors = np.loadtxt(rfile)
    context.cols = eval(cols)

@when('we compute many mdmrs in python')
def step(context):
    context.voxs   = range(12)    
    context.Fs, context.ps = calc_mdmrs(context.dmats[context.voxs], context.regressors, 
                                        context.cols, context.perms)

@when('we compute a single MDMR in python')
def step(context):
    context.voxs = [0]
    ps, Fs, F_perms, _ = mdmr(context.dmats[0].reshape(context.nobs**2,1), context.regressors,  
                              context.cols, context.perms)
    context.ps = ps
    context.Fs = Fs
    context.F_perms = F_perms

@then('the many pseudo-F values should be like R')
def step(context):
    comp = np.allclose(context.Fs, context.r_Fs[context.voxs])
    assert_that(comp, "pseudo-F stats")

@then('the many F permutations should be like R')
def step(context):
    comp = np.allclose(context.F_perms, context.r_F_perms[:,context.voxs])
    assert_that(comp, "F permutations")

@then('the many p-values should be like R')
def step(context):
    comp = np.allclose(context.ps, context.r_ps[context.voxs])
    assert_that(comp, "p-values")

@then('the many p-values should be similar to R')
def step(context):
    comp = np.corrcoef(context.ps, context.r_ps[context.voxs])[0,1]
    assert_that(comp, greater_than(0.98), "p-values")
