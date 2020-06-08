# -*- coding: utf-8 -*-


def fsl_aff_to_rigid(in_xfm, out_name):
    out_mat = os.path.join(os.getcwd(), out_name)
    cmd = ['aff2rigid', in_xfm, out_mat]
    retcode = subprocess.check_output(cmd)
    return out_mat


def create_3dskullstrip_arg_string(shrink_fac, var_shrink_fac,
                                   shrink_fac_bot_lim, avoid_vent, niter,
                                   pushout, touchup, fill_hole, avoid_eyes,
                                   use_edge, exp_frac, smooth_final,
                                   push_to_edge, use_skull, perc_int,
                                   max_inter_iter, blur_fwhm, fac, monkey, mask_vol):
    """
    Method to return option string for 3dSkullStrip
    
    Parameters
    ----------
    shrink_fac : float
        Parameter controlling the brain vs non-brain intensity threshold (tb)

    var_shrink_fac : boolean
        Vary the shrink factor with the number of iterations

    shrink_fac_bot_lim : float
        Do not allow the varying SF to go below SFBL

    avoid_vent : boolean
        Avoid ventricles

    niter : float
        Number of iterations

    pushout : boolean
        Consider values above each node in addition to values below the node when deciding on expansion

    touchup : boolean
        Perform touchup operations at end to include areas not covered by surface expansion

    fill_hole : float
         Fill small holes that can result from small surface intersections caused by the touchup operation

    avoid_eyes : boolean
        Avoid eyes

    use_edge : boolean
        Use edge detection to reduce leakage into meninges and eyes

    exp_frac : float
        Speed of expansion

    smooth_final : float
        Perform final surface smoothing after all iterations

    push_to_edge : boolean
        Perform aggressive push to edge at the end

    use_skull : boolean
        Use outer skull to limit expansion of surface into the skull due to very strong shading artifacts

    perc_int : float
        Percentage of segments allowed to intersect surface

    max_inter_iter : float
        Number of iteration to remove intersection problems

    blur_fwhm : float
        Blur dset after spatial normalization

    fac : float
         Multiply input dataset by FAC if range of values is too small

    monkey : boolean
        Use monkey option in SkullStripping
    
    mask_vol : boolean
        Output a mask volume instead of a skull-stripped volume.

    Returns
    -------
    opt_str : string
        Command args
    
    """

    expr = ''
    defaults = dict(
        fill_hole=10 if touchup else 0,
        shrink_fac=0.6,
        shrink_fac_bot_lim=0.4 if use_edge else 0.65,
        niter=250,
        exp_frac=0.1,
        smooth_final=20,
        perc_int=0,
        max_inter_iter=4,
        blur_fwhm=0,
        fac=1.0,
        monkey=False,
        mask_vol=False
    )

    if float(shrink_fac) != defaults['shrink_fac']:
        expr += ' -shrink_fac {0}'.format(shrink_fac)

    if not var_shrink_fac:
        expr += ' -no_var_shrink_fac'
    
    if mask_vol:
        expr += ' -mask_vol'

    if monkey:
        expr += ' -monkey'

    if float(shrink_fac_bot_lim) != defaults['shrink_fac_bot_lim']:
        expr += ' -shrink_fac_bot_lim {0}'.format(shrink_fac_bot_lim)

    if not use_edge:
        expr += ' -no_use_edge'

    if not avoid_vent:
        expr += ' -no_avoid_vent'

    if int(niter) != defaults['niter']:
        expr += ' -niter {0}'.format(niter)

    if not pushout:
        expr += ' -no_pushout'

    if not touchup:
        expr += ' -no_touchup'

    if int(fill_hole) != defaults['fill_hole']:
        expr += ' -fill_hole {0}'.format(fill_hole)

    if not avoid_eyes:
        expr += ' -no_avoid_eyes'

    if float(exp_frac) != defaults['exp_frac']:
        expr += ' -exp_frac {0}'.format(exp_frac)

    if int(smooth_final) != defaults['smooth_final']:
        expr += ' -smooth_final {0}'.format(smooth_final)

    if push_to_edge:
        expr += ' -push_to_edge'

    if use_skull:
        expr += ' -use_skull'

    if float(perc_int) != defaults['perc_int']:
        expr += ' -perc_int {0}'.format(perc_int)

    if int(max_inter_iter) != defaults['max_inter_iter']:
        expr += ' -max_inter_iter {0}'.format(max_inter_iter)

    if float(blur_fwhm) != defaults['blur_fwhm']:
        expr += ' -blur_fwhm {0}'.format(blur_fwhm)

    if float(fac) != defaults['fac']:
        expr += ' -fac {0}'.format(fac)

    return expr
