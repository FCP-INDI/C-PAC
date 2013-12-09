from base_cwas import *
import code

def fsl_cluster_sizes(img_file, thr):
    
def fsl_cluster_perms(Fperms, fstats, pvals, mask_file):
    import os
    from tempfile import mkstemp
    
    # Load mask (and header data)
    mask_img    = nib.load(mask_file)
    mask        = mask_img.get_data()
    mask_inds   = np.where(mask>0)
    
    # Pseudo-F Threshold
    f_thr       = fstat_given_pval(fstats, pvals, 0.05)
    
    # Setup
    nperms      = Fperms.shape[0]
    max_sizes   = np.zeros(nperms)
    tmpfile     = mkstemp()
    
    # Loop through each permutation to get the maximum cluster size
    for i in range(nperms):
        # Save this particular whole-brain value
        fperms_brain            = np.zeros(mask.shape)
        fperms_brain[mask_inds] = Fperms[i,:]
        img = nib.Nifti1Image(fperms_brain, mask_img.affine(), mask_img.header())
        img.to_filename(tmpfile)
        
        # Cluster with FSL and get max size
        max_sizes[i] = fsl_cluster_sizes(tmpfile, f_thr).max()
        
        # Remove old file
        os.remove(tmpfile)
    
    
    clust_inds  = np.array(range(nclust)) + 1
    table       = np.hstack((clust_inds, clust_sizes[subset], clust_pvals[subset]))
    
    # Save the fstats image
    fstats_brain            = np.zeros(mask.shape)
    fstats_brain[mask_inds] = fstats
    img                     = nib.Nifti1Image(fstats_brain, mask_img.affine(), mask_img.header())
    img.to_filename(tmpfile)
    
    # Get the cluster image of the thresholded fstats image
    from subprocess import Popen, PIPE
    clust_file  = mkstemp()
    cmd = "cluster -i %s -t %f -o %s" % (tmpfile, f_thr, clust_file)
    os.system(cmd)
    clust       = nib.load(clust_file).get_data()[mask>0]
    
    # Get the cluster sizes and p-values
    clust_sizes = fsl_cluster_sizes(tmpfile, f_thr)
    clust_pvals = np.array([ (cz>=max_sizes).mean() for cz in clust_sizes ])
    
    # Threshold a little bit
    subset      = clust_pvals < 0.05
    nclust      = subset.sum()
    
    # Create the thresholded cluster
    clust_thr   = np.zeros_like(clust)
    for i,val in enumerate(subset.nonzero()):
        clust_thr[clust==val] = i+1
    
    # Get the final table
    clust_inds  = np.array(range(nclust)) + 1
    table       = np.hstack((clust_inds, clust_sizes[subset], clust_pvals[subset]))
    
    # Remove temporary folders
    os.remove(tmpfile)
    os.remove(clust_file)
    
    return (table, clust_thr)


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



cluster_correct_perms(Fperms, fstats, pvals, mask, clust_thr=0.05, vox_thr=0.05)