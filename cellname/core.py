"""Main module."""
import os
import numpy as np
import pandas as pd
import collections
from sklearn.preprocessing import scale as sklearn_scale
from .util import p_adjust_bh, _guess_cell_type,standard


def predict_celltype(obj, cell_type_markers=None, clusters="leiden", q=0.75, qv=0.1):
    """
    change the code within adobo python packages:https://github.com/oscar-franzen/adobo/tree/master/adobo
    But the orgin code can't run sucess,and don't fit the BGI stereo-seq,So I have change some for bettter use!
    :param qv: the Q value thresholds.
    :param q: the percent gene is marker gene.
    :param obj: the anndata which scanpy use
    :param cell_type_markers: the pandas Dataframe which including marker gene and cell tyoe,two columns.
    :param clusters:which cluster you want to use in obj.obs
    :return a cell type annotation pandas Dataframe
    """

    if isinstance(cell_type_markers, pd.DataFrame):
        # custom cell type markers were provided
        ma_ss = cell_type_markers
        ma_ss.columns = ['official gene symbol', 'cell type']
    elif isinstance(cell_type_markers,str):
        ma_ss = pd.read_csv(cell_type_markers,sep="\t")
        ma_ss.columns = ['official gene symbol', 'cell type']
    else:
        ma = pd.read_csv('%s/../data/markers.tsv' %
                         os.path.dirname(__file__), sep='\t')
        # restrict to mouse
        # ma = ma[ma.species.str.match('Hs')]
        markers = ma
        ui = ma.iloc[:, ma.columns == 'ubiquitousness index']
        ma = ma[np.array(ui).flatten() < 0.05]
        ma_ss = ma.iloc[:, ma.columns.isin(['official gene symbol',
                                            'cell type'])]

    marker_freq = ma_ss[ma_ss.columns[0]].value_counts()
    markers = ma_ss

    dd = collections.defaultdict(list)
    for item in markers.groupby('cell type'):
        dd[item[0]] = set(item[1][item[1].columns[0]])
    # down-weighting overlapping genes improves gene set analysis
    # Tarca AL, Draghici S, Bhatti G, Romero R; BMC Bioinformatics 2012 13:136

    weights = 1 + np.sqrt(((max(marker_freq) - marker_freq) /
                           (max(marker_freq) - min(marker_freq))))
    # Get the expression matrixs and the clusters.
    t = obj.raw.X.toarray()
    t_pd = pd.DataFrame(data=t, index=obj.raw.obs_names, columns=obj.raw.var_names)
    df = t_pd.transpose()
    norm = standard(df)
    X = np.log2(norm + 1)
    min_cluster_size = 10
    cl = obj.obs[clusters]
    ret = X.groupby(cl.values, axis=1).quantile(q)
    q = pd.Series(cl).value_counts()
    cl_remove = q[q < min_cluster_size].index
    ret = ret.iloc[:, np.logical_not(ret.columns.isin(cl_remove))]
    median_expr = ret
    median_expr.index = median_expr.index.str.upper()
    # s = np.sum(median_expr.index.str.match('^(.+)_.+'))
    # if median_expr.shape[0] == s:
    #     input_symbols = median_expr.index.str.extract(
    #         '^(.+)_.+')[0]
    #     input_symbols = input_symbols.str.upper()
    #     median_expr.index = input_symbols
    # # (1) centering is done by subtracting the column means
    # # (2) scaling is done by dividing the (centered) by their standard
    # # deviations
    scaled = sklearn_scale(median_expr, with_mean=True, axis=0)
    median_expr_Z = pd.DataFrame(scaled)
    median_expr_Z.index = median_expr.index
    median_expr_Z.columns = median_expr.columns
    ret = median_expr_Z.apply(func=_guess_cell_type, axis=0, args=(median_expr, dd, weights))
    # restructure
    bucket = []
    for i, kk in enumerate(ret):
        lines = [line for line in ret[kk]]
        _df = pd.DataFrame.from_dict(lines, orient='columns')
        _df['cluster'] = [i] * _df.shape[0]
        cols = _df.columns.tolist()
        _df = _df[cols[-1:] + cols[:-1]]
        bucket.append(_df)
    final_tbl = pd.concat(bucket)
    if final_tbl.shape[0] == 0:
        raise Exception('Final table is empty. Check gene symbols of input data.')
    padj = p_adjust_bh(final_tbl['pvalue'])
    final_tbl['padj_BH'] = padj
    final_tbl.columns = ['cluster',
                         'activity score',
                         'cell type',
                         'p-value',
                         'markers',
                         'adjusted p-value BH']
    # save the best scoring for each cluster
    res_pred = final_tbl.groupby('cluster').nth(0)
    _a = res_pred['adjusted p-value BH'] > qv
    res_pred.loc[_a, 'cell type'] = 'Unknown'
    return res_pred