from base_cwas import *

@given('simulated subjects time-series data')
def step(context):
    context.nsubs = 10
    context.ntpts = 100
    context.nvoxs = 200
    context.sdata = np.random.random((context.nsubs, context.ntpts, context.nvoxs))

@given('the subject list "{sfile}"')
def step(context, sfile):
    context.slist = [ l.strip().strip('"') for l in open(sfile).readlines() ]
    context.nsubs = len(context.slist)

## Defined in step_connectivity.py
#@given('mask data from "{mfile}"')
#def step(context, mfile):
#    context.mask = nib.load(mfile).get_data().astype('bool')
#    context.indices = np.where(context.mask)

@given('the subjects data are masked')
def step(context):
    context.sdata = [ nib.load(sfile).get_data().astype('float64')[context.indices].T
                        for sfile in context.slist ]
    context.ntpts = context.sdata[0].shape[0]
    context.nvoxs = context.sdata[0].shape[1]

@when('we norm the subjects data')
def step(context):
    context.sdata_normed = norm_subjects(context.sdata)

@when('we compute the connectivity for multiple subjects on normed data')
def step(context):
    context.vox_inds = range(12)
    context.cmats    = ncor_subjects(context.sdata_normed, context.vox_inds)

@when('we compute the connectivity for multiple subjects with numpy')
def step(context):
    context.ref_cmats = np.zeros_like(context.cmats)
    for i in range(context.nsubs):
        context.ref_cmats[i] = custom_corrcoef(context.sdata[i][:,context.vox_inds].T, context.sdata[i])

## Defined in conn_single.feature
#@then('the correlation values for cpac should be like numpy')
#def step(context):
#    comp = np.allclose(context.cmats, context.ref_cmats)
#    assert_that(comp, "subject correlation")