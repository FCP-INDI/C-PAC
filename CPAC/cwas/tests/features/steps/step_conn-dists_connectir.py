from base_cwas import *
from pytest_bdd import scenario, given, when, then

@given('a connectir-based distance folder "{sdir}"')
def step(context, sdir):
    bigmemory = importr('bigmemory')
    # Distances
    sfile           = op.join(sdir, "subdist.desc")
    dmats           = np.array(robjects.r("as.matrix(attach.big.matrix('%s'))" % sfile))
    context.nobs    = np.sqrt(dmats.shape[0])
    context.nvoxs   = dmats.shape[1]
    context.r_dmats = dmats.T.reshape(context.nvoxs, context.nobs, context.nobs)
    # Mask
    mfile           = op.join(sdir, "mask.nii.gz")
    mask            = nib.load(mfile).get_data().astype('bool')
    mask_indices    = np.where(mask)
    context.mask    = mask
    context.indices = mask_indices
    # Convert from column to row major format (aka take R data to python)
    mfile           = op.join(sdir, "mask_r2py.nii.gz")
    inds_r2py       = nib.load(mfile).get_data().astype('int')
    inds_r2py       = inds_r2py[inds_r2py.nonzero()] - 1
    context.r_dmats = context.r_dmats[inds_r2py,:,:]
    
## Done in step_connectivity_multiple_subjects.py
#
#@given('the subject list "{sfile}"')
#def step(context, sfile):
#    context.slist = [ l.strip().strip('"') for l in open(sfile).readlines() ]
#    context.nsubs = len(context.slist)
#
#@given('the subjects data are masked')
#def step(context):
#    context.sdata = [ nib.load(sfile).get_data().astype('float64')[context.indices].T
#                        for sfile in context.slist ]
#    context.sdata = np.array(context.sdata)
#    context.ntpts = context.sdata.shape[1]
#    context.nvoxs = context.sdata.shape[2]

@when('we calculate distances for {nseeds} seeds with cpac')
def step(context, nseeds):
    context.nseeds  = int(nseeds)
    context.seeds   = range(context.nseeds)
    context.dmats   = calc_subdists(context.sdata, (0,context.nseeds))

@when('we calculate distances for {nseeds} seeds with numpy')
def step(context, nseeds):
    context.nseeds      = int(nseeds)
    context.seeds       = range(context.nseeds)
    context.ref_dmats   = np.zeros((context.nseeds, context.nsubs, context.nsubs))
    Smaps = np.zeros((context.nsubs, context.nseeds, context.nvoxs))
    for i in range(context.nsubs):
        x = custom_corrcoef(context.sdata[i][:,context.seeds].T, context.sdata[i])
        Smaps[i,:,:] = x
    for j in range(context.nseeds):
        S   = Smaps[:,context.seeds[j],:]
        S0  = np.delete(S, context.seeds[j], 1)
        S0  = np.arctanh(S0)
        D   = 1 - np.corrcoef(S0)
        context.ref_dmats[j] = D

@then('the cpac distances should be the same as numpy')
def step(context):
    comp = np.allclose(context.dmats, context.ref_dmats)
    assert_that(comp, "cpac vs numpy distances")

@then('the cpac distances should be like connectir')
def step(context):
    import code
    code.interact(local=locals())
    
    dmats   = context.dmats.reshape((context.nseeds, context.nsubs**2)).T
    r_dmats = context.r_dmats.reshape((context.nvoxs, context.nsubs**2))[context.seeds,:].T
    
    corrs = np.diag(custom_corrcoef(dmats, r_dmats))
    comp  = (corrs>0.99).all()
    assert_that(comp, "cpac vs R distances")
    
