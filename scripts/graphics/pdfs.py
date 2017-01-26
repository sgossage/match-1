from __future__ import print_function, absolute_import
import matplotlib.pyplot as plt
import numpy as np
from .graphics import corner_setup, fix_diagonal_axes, key2label


def pdf_plot(SSP, xattr, yattr=None, ax=None, sub=None, save=False,
             truth=None, cmap=None, plt_kw=None, X=None, prob=None,
             logp=True):
    """Plot -2 ln P vs marginalized attributes

    SSP : SSP class instance

    xattr, yattr : str
        column names to marginalize and plot

    ax : plt.Axes
        add plot to this axis

    sub : str
        if save, add this string to the filename

    save : bool
        save the plot to file with axes lables.

    truth : dict
        truth dictionary with attributes as keys and truth values as values
        overplot as hline and vline on axes.

    cmap : cm.cmap instance
        if yattr isn't None, call plt.pcolor with this cmap.

    plt_kw : dict
        if yattr is None, pass these kwargs to plt.plot

    Returns
    ax : plt.Axes
        new or updated axes instance
    """
    plt_kw = plt_kw or {'lw': 4, 'color': 'k'}
    cmap = cmap or plt.cm.viridis_r
    sub = sub or ''
    truth = truth or {}
    do_cbar = False
    pstr = ''
    if log:
        pstr = 'Log '
    if ax is None:
        _, ax = plt.subplots()
        if yattr is not None:
            do_cbar = True

    if yattr is None:
        # plot type is marginal probability. Attribute vs -2 ln P
        if X is None and prob is None:
            X, prob = SSP.marginalize(xattr, log=logp)
            if not SSP.vdict[xattr]:
                return ax
            SSP.build_posterior(xattr, X, prob)
        l = ax.plot(X, prob, **plt_kw)
        ax.set_xlim(X.min(), X.max())
        if save:
            ptype = 'marginal'
            # ax.set_ylabel(key2label('fit'))
            ax.set_ylabel(key2label(pstr+'Probability'))
    else:
        # plot type is joint probability.
        # Attribute1 vs Attribute2 colored by fit
        [X, Y], prob = SSP.marginalize(xattr, yattr=yattr, log=logp)
        if not SSP.vdict[xattr] or not SSP.vdict[yattr]:
            return ax
        l = ax.pcolor(X, Y, prob, cmap=cmap)
        ax.set_xlim(X.min(), X.max())
        ax.set_ylim(Y.min(), Y.max())

        if do_cbar:
            cbar = plt.colorbar(l)
            # cbar.set_label(key2label('fit'))
            cbar.set_label(key2label(pstr+'Probability'))

        if save:
            ptype = '{}_joint'.format(yattr)
            ax.set_ylabel(key2label(yattr, gyr=SSP.gyr))

        if yattr in truth:
            ax.axhline(truth[yattr], color='k', lw=3)

    if xattr in truth:
        ax.axvline(truth[xattr], color='k', lw=3)

    if save:
        ax.set_xlabel(key2label(xattr, gyr=SSP.gyr))
        # add subdirectory to filename
        if len(sub) > 0:
            sub = '_' + sub
        outfmt = '{}_{}{}_{}{}'
        outname = outfmt.format(SSP.name.replace('.csv', ''),
                                xattr, sub, ptype, EXT)
        plt.savefig(outname, bbox_inches='tight')
        print('wrote {}'.format(outname))
        plt.close()
    return ax


def pdf_plots(SSP, marginals=None, sub=None, twod=False, truth=None,
              text=None, cmap=None, fig=None, axs=None, frompost=False):
    """Call pdf_plot for a list of xattr and yattr"""
    text = text or ''
    sub = sub or ''
    truth = truth or {}
    marginals = marginals or SSP._getmarginals()

    ndim = len(marginals)
    if twod:
        fig, axs = corner_setup(ndim)
        raxs = []
        for c, mx in enumerate(marginals):
            for r, my in enumerate(marginals):
                ax = axs[r, c]
                if r == c:
                    # diagonal
                    # my = 'fit'  # my is reset for ylabel call
                    my = 'Probability'
                    raxs.append(SSP.pdf_plot(mx, ax=ax, truth=truth))
                else:
                    # off-diagonal
                    raxs.append(SSP.pdf_plot(mx, yattr=my, ax=ax, truth=truth,
                                             cmap=cmap))

                if c == 0:
                    # left most column
                    ax.set_ylabel(key2label(my))

                if r == ndim - 1:
                    # bottom row
                    ax.set_xlabel(key2label(mx))
        [ax.locator_params(axis='y', nbins=6) for ax in axs.ravel()]
        fix_diagonal_axes(raxs, ndim)
    else:
        raxs = []
        if fig is None and axs is None:
            fig, axs = plt.subplots(ncols=ndim, figsize=(15, 3))
        [ax.tick_params(left='off', labelleft='off', right='off', top='off')
         for ax in axs]
        X = None
        prob = None
        for i in marginals:
            if frompost:
                pattr = '{:s}prob'.format(i)
                X = SSP.data[i][np.isfinite(SSP.data[i])]
                prob = SSP.data[pattr][np.isfinite(SSP.data[pattr])]
            ax = axs[marginals.index(i)]
            ax = SSP.pdf_plot(i, truth=truth, ax=ax, X=X, prob=prob)
            ax.set_xlabel(key2label(i, gyr=SSP.gyr))
            raxs.append(ax)

        if text:
            axs[-1].text(0.90, 0.10, '${}$'.format(text), ha='right',
                         transform=axs[-1].transAxes)
        fig.subplots_adjust(bottom=0.2, left=0.05)

        raxs[0].set_ylabel(key2label('Probability'))
    [ax.locator_params(axis='x', nbins=6) for ax in axs.ravel()]
    return fig, raxs
