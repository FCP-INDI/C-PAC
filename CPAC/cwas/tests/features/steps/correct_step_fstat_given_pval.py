from base_cwas import *
import code

@given('simulated fstats')
def step(context):
    context.fstats = np.random.f(100, 1, 42)

@given('simulated pvals')
def step(context):
    context.pvals = np.random.random(42)

@given('a pval of {pval}')
def step(context, pval):
    context.given_p = float(pval)

@when('we infer an fstat given a pval using CPAC')
def step(context):
    context.cpac_fstat = fstat_given_pval(context.fstats, context.pvals, context.given_p)

@when('we infer an fstat given a pval using R')
def step(context):
    robjects.globalenv["fstats"]    = context.fstats
    robjects.globalenv["pvals"]     = context.pvals
    fit     = robject.r.lm("fstats ~ pvals")
    new_dat = robjects.DataFrame({'pvals': robjects.IntVector([context.given_p])})
    context.r_fstat = np.array(robjects.r.predict(fit, new_dat))[0]

@then('the fstat derived from CPAC should be like R')
def step(context):
    comp = np.allclose(context.cpac_fstat, context.r_fstat)
    assert_that(comp, "pseudo-F stat")
