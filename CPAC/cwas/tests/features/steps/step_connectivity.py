from base_cwas import *

@given('simulated time-series data')
def step(context):
    context.ntpts = 100
    context.nvoxs = 200
    context.sdata = np.random.random((context.ntpts, context.nvoxs))

@given('subject data from "{sfile}"')
def step(context, sfile):
    context.sdata = nib.load(sfile).get_data().astype('float64')

@given('mask data from "{mfile}"')
def step(context, mfile):
    context.mask = nib.load(mfile).get_data().astype('bool')
    context.indices = np.where(context.mask)

@given('the subject data is masked')
def step(context):
    context.sdata = context.sdata[context.indices].T

@when('we norm the data')
def step(context):
    context.sdata_normed = norm_cols(context.sdata)

@when('we compute the connectivity on normed data')
def step(context):
    context.vox_inds = range(12)
    context.cmat = ncor(context.sdata_normed, context.vox_inds)

@then('the correlation values should be like the standard correlation function')
def step(context):
    import code
    #code.interact(local=locals())
    context.ref_cmat = np.corrcoef(context.sdata.T)[context.vox_inds,:]
    comp = np.allclose(context.cmat, context.ref_cmat)
    assert_that(comp, "subject correlation")