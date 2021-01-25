from base_cwas import *
from pytest_bdd import scenario, given, when, then

@given('simulated subjects connectivity data for {nseeds} seed voxels')
def step(context, nseeds):
    context.nsubs   = 10
    context.nseeds  = int(nseeds)
    context.nvoxs   = 200
    context.cmats   = np.random.random((context.nsubs, context.nseeds, context.nvoxs))
    context.cmats   = 2*context.cmats - 1

@when('we compute the distances with cpac')
def step(context):
    context.dmats   = np.zeros((context.nseeds, context.nsubs, context.nsubs))
    for i in range(context.nseeds):
        context.dmats[i] = compute_distances(context.cmats[:,i,:])

@when('we compute the distances with numpy')
def step(context):
    context.ref_dmats   = np.zeros((context.nseeds, context.nsubs, context.nsubs))
    for i in range(context.nseeds):
        context.ref_dmats[i] = 1 - np.corrcoef(context.cmats[:,i,:])

@then('the cpac distances should be like numpy')
def step(context):
    comp = np.allclose(context.dmats, context.ref_dmats)
    assert_that(comp, "distances")