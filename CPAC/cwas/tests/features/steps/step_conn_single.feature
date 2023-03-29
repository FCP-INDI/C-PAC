from base_cwas import *
from pytest_bdd import scenario, given, when, then

@given('simulated time-series data')
def step(context):
    context.ntpts = 100
    context.nvoxs = 200
    context.sdata = np.random.random((context.ntpts, context.nvoxs))

@given('subject data from "{sfile}"')
def step(context, sfile):
    context.sdata = nib.load(sfile).get_fdata().astype('float64')

@given('mask data from "{mfile}"')
def step(context, mfile):
    context.mask = nib.load(mfile).get_fdata().astype('bool')
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
    context.cmats = ncor(context.sdata_normed, context.vox_inds)

@when('we compute the connectivity with numpy')
def step(context):
    context.ref_cmats = custom_corrcoef(context.sdata[:,context.vox_inds].T, context.sdata)

@then('the correlation values for cpac should be like numpy')
def step(context):
    comp = np.allclose(context.cmats, context.ref_cmats)
    assert_that(comp, "subject correlations")