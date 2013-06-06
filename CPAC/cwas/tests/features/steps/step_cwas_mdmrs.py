from base_cwas import *
import code
 
@given('a connectir-based MDMR folder "{mdir}" looking at the factor "{factor}"')
def step(context, mdir, factor):
    # Distances
    sfile         = op.join(op.dirname(mdir), "subdist.desc")
    dmats         = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % sfile))
    context.nobs  = np.sqrt(dmats.shape[0])
    context.dmats = dmats
    # Fstats
    ffile         = op.join(mdir, "fperms_%s.desc" % factor)
    context.r_Fs  = np.array(robjects.r("attach.big.matrix('%s')[1,]" % ffile))    
    # Pvalues
    mfile         = op.join(op.dirname(mdir), "mask.nii.gz")
    mask          = nib.load(mfile).get_data().nonzero()
    pfile         = op.join(mdir, "one_minus_significance_%s.nii.gz" % factor)
    context.r_ps  = 1 - nib.load(pfile).get_data()[mask]

@given('regressors from "{rfile}" with columns of interest {cols}')
def step(context, rfile, cols):
    context.regressors = np.loadtxt(rfile)
    context.cols = eval(cols)

@when('we compute many mdmrs in python')
def step(context):
    context.nperms = 1000
    context.Fs, context.ps = calc_mdmrs(context.dmats, context.regressors, context.cols, context.nperms)

@then('the many pseudo-F values should be like R')
def step(context):
    comp = np.allclose(context.Fs, context.r_Fs)
    assert_that(comp, "pseudo-F stats")

@then('the many p-values should be similar to R')
def step(context):
    comp = np.corrcoef(context.ps, context.r_ps)[0,1]
    code.interact(local=locals())
    assert_that(comp, greater_than(0.98))
