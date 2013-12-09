from base_cwas import *
import code

@given("image file '{filepath}'")
def step(context, filepath):
    context.image_file = filepath

@given("mask file '{filepath}'")
def step(context, filepath):
    context.mask_file = filepath

@given('voxel threshold of {thresh}')
def step(context, pval):
    context.thresh = float(thresh)

@when('we calculate the cluster sizes using CPAC')
def step(context):
    labels = cluster_file(context.image_file, context.thresh, context.mask_file)
    context.cpac_sizes = cluster_sizes(labels)

@when('we calculate the cluster sizes using FSL')
def step(context):
    sizes = fsl_cluster_sizes(context.image_file, context.thresh)
    context.fsl_sizes = sizes[sizes>1] # for now our cluster script misses the isolated voxels (cluster of size 1)

@then('the cluster sizes derived from CPAC should be like FSL')
def step(context):
    comp = np.allclose(context.cpac_sizes, context.r_clust)
    assert_that(comp, "cluster sizes")
