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
    context.sdata = np.array(context.sdata)
    context.ntpts = context.sdata.shape[1]
    context.nvoxs = context.sdata.shape[2]

@when('we norm the subjects data')
def step(context):
    context.sdata_normed = norm_subjects(context.sdata)

@when('we compute the connectivity on normed subjects data')
def step(context):
    context.vox_inds = range(12)
    context.carr     = ncor_subjects(context.sdata_normed, context.vox_inds)

@then('the correlation values for multiple subjects should be like the standard correlation function')
def step(context):
    context.ref_carr = np.zeros_like(context.carr)
    for i in range(context.nsubs):
        context.ref_carr[i] = np.corrcoef(context.sdata[i].T)[context.vox_inds,:]
    comp = np.allclose(context.carr, context.ref_carr)
    assert_that(comp, "subject correlation")