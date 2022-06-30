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

def standard(data, scaling_factor=10000):
    """Performs a standard normalization by scaling with the total
    read depth per cell and then multiplying with a scaling factor.
    Parameters
    ----------
    data : :class:`pandas.DataFrame`
        A pandas data frame object containing raw read counts
        (rows=genes, columns=cells).
    scaling_factor : `int`
        Scaling factor used to multiply the scaled counts
        with. Default: 10000
    References
    ----------
    .. [1] Evans et al. (2018) Briefings in Bioinformatics
           https://academic.oup.com/bib/article/19/5/776/3056951
    .. [2] https://github.com/oscar-franzen/adobo/blob/ee3362eaef1e311a81ceb6e4a272721de420bae6/adobo/normalize.py#L256
    Returns
    -------
    :class:`pandas.DataFrame`
        A normalized data matrix with same dimensions as before.
    """
    col_sums = data.sum(axis=0).values
    data_norm = (data / col_sums) * scaling_factor
    return data_norm
