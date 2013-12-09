from base_cwas import *

# FROM step_conn_multiple.py
#@given('simulated subjects time-series data')
#def step(context):
#    context.nsubs = 12
#    context.ntpts = 100
#    context.nvoxs = 200
#    context.sdata = np.random.random((context.nsubs, context.ntpts, context.nvoxs))

@when('we compute the distances from functional data with cpac for {nvoxels} voxel block')
def step(context, nvoxels):
    context.dmats = calc_subdists(context.sdata, [0,context.nvoxs], int(nvoxels))

@when('we compute the distances from functional data with cpac using a voxel range')
def step(context, nvoxels):
    dmats           = np.zeros((context.nvoxs, context.nsubs, context.nsubs))
    vranges         = [[0,np.floor(context.nvoxs)],[np.floor(context.nvoxs),context.nvoxs]]
    dmats[vranges[0]] = calc_subdists(context.sdata, vranges[0])
    dmats[vranges[0]] = calc_subdists(context.sdata, vranges[1])
    context.dmats   = dmats

@when('we compute the distances from functional data with numpy')
def step(context):    
    # Calculate connectivity
    Smaps = np.zeros((context.nsubs, context.nvoxs, context.nvoxs))
    for i in range(context.nsubs):
        Smaps[i,:,:] = np.corrcoef(context.sdata[i].T)
    
    # Calculate distances
    context.ref_dmats   = np.zeros((context.nvoxs, context.nsubs, context.nsubs))
    for j in range(context.nvoxs):
        S   = Smaps[:,j,:]
        S0  = S.copy()
        S0[:,j] = 0.99
        S0  = np.arctanh(S0)
        D   = 1 - np.corrcoef(S0)
        context.ref_dmats[j] = D

# FROM step_dists.py
#@then('the cpac distances should be like numpy')
#def step(context):
#    comp = np.allclose(context.dmats, context.ref_dmats)
#    assert_that(comp, "distances")