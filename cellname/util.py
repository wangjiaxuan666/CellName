import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


def p_adjust_bh(p):
    """The Benjamini-Hochberg p-value correction for multiple hypothesis testing.

    Parameters
    ----------
    p : `list`
        A list of p-values.

    References
    ----------
    .. [1] Benjamini & Hochberg (1995) Controlling the false discovery rate: a practical
        and powerful approach to multiple testing.  J Royal Statistical Society, Series B

    Returns
    -------
    int
        Adjusted p-values.
    """
    p = np.asfarray(p)
    p = np.ma.array(p, mask=np.isnan(p))  # to handle nan
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(float(len(p)), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    q = q[by_orig]
    q[np.isnan(p.data)] = np.nan  # put nan back
    return q.data


def _guess_cell_type(x, median_expr, dd, weights):
    rr = median_expr.loc[:, median_expr.columns == x.name].values.flatten()
    # genes expressed in this cell cluster
    genes_exp = set(x.index[rr > 0])
    # genes _not_ expressed in this cell cluster
    genes_not_exp = set(x.index[rr == 0])
    res = list()
    for ct in dd:
        s = dd[ct]
        x_ss = x[x.index.isin(s)]
        if len(x_ss) == 0:
            continue
        gene_weights = weights[weights.index.isin(x_ss.index)]
        gene_weights = pd.Series(gene_weights, x_ss.index)
        activity_score = sum(x_ss * gene_weights) / len(x_ss) ** 0.3
        # how many expressed genesets are found in the geneset?
        ct_exp = len(genes_exp & s)
        # how many _non_ expressed genes are found in the geneset?
        ct_non_exp = len(genes_not_exp & s)
        # how many expressed genes are NOT found in the geneset?
        ct_exp_not_found = len(genes_exp - s)
        # how many _non_ expressed genes are NOT found in the geneset?
        not_exp_not_found_in_geneset = len(genes_not_exp - s)
        # one sided fisher
        contigency_tbl = [[ct_exp, ct_non_exp],
                          [ct_exp_not_found, not_exp_not_found_in_geneset]]
        odds_ratio, pval = fisher_exact(
            contigency_tbl, alternative='greater')
        markers_found = ','.join(list(genes_exp & s))
        if markers_found == '':
            markers_found = 'NA'
        res.append({'activity_score': activity_score,
                    'ct': ct,
                    'pvalue': pval,
                    'markers': markers_found})
    res = sorted(res, key=lambda k: k['activity_score'], reverse=True)
    return res
