"""Plot that mimics MATCH/pgcmd use with CMD class for enhancements."""
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import numpy as np
from .graphics import square_aspect, zeroed_cmap, add_inner_title
import matplotlib.transforms as mtransforms

import matplotlib as mpl
import os

def setup_imgrid(figsize=[12, 3], nrows=1, ncols=4):
    """Default settings to initialize match_plot"""
    from mpl_toolkits.axes_grid1 import ImageGrid
    igkw = {'nrows_ncols': (nrows, ncols),
            'axes_pad': .7,
            'label_mode': "all",
            'share_all': True,
            'cbar_location': "top",
            'cbar_mode': "each",
            'cbar_size': "7%",
            'cbar_pad': "2%"}

    fig = plt.figure(figsize=figsize)
    grid = ImageGrid(fig, 111, **igkw)
    return grid

def mpl_hack(ax):
    import matplotlib
    ver = float(matplotlib.__version__[:3])
    if ver <= 1.3:
        ax.set_adjustable('box-forced')
        ax.autoscale(False)
    return

def match_plot(hesslist, extent, mag, color, bins, labels=None, twobytwo=True, sig=True,
               xlabel=None, ylabel=None, cmap=None, logcounts=False,
               photf_pts=None, mist_pts=None, best_list=None, savedir=""):
    '''
    Plot four hess diagrams with indivdual color bars using ImageGrid
    hesslist : list
        list of length 3 or 4 of Nx2 arrays to pass to plt.imshow
    extent : list
        extent passed to plt.imshow
    labels=None : list
        labels corresponding to hesslist
    twobytwo : bool [True]
        image grid is 2x2 (or one row)
    sig : bool [True]
        include significance (or 4th hess in the hesslist)
    xlabel, ylabel : str, str
        axes labels
    cmap : colormap or list of colormaps [None]
        colormap to use for each or all hess plots.
        default color map uses zeroed_cmap for diff and sig plots to set 0 to
        be white.
    logcounts: bool [False]
        pass np.log10(hess) to imshow instead of hess
    photf_pts: tuple [None]
        points for overplotting data from a photometry file w/
        format (v-i, v).
    '''
    if twobytwo:
        figsize = [15, 15]
        nrows = 2
        ncols = 2
    else:
        nrows = 1
        ncols = 4
        if not sig:
            ncols = 3
        figsize = [ncols * 3, 3]

    if not sig and len(hesslist) == 4:
        hesslist = hesslist[:-1]

    grid = setup_imgrid(figsize=figsize, nrows=nrows, ncols=ncols)

    for i, (ax, hess) in enumerate(zip(grid, hesslist)):
        if i < 2:
            vmin = np.min(hesslist[0])
            vmax = np.max(hesslist[0])
        elif i == 2:
            vmin = -np.max(hesslist[0])
            vmax = np.max(hesslist[0])
        else:
            vmin = np.min(hesslist[i])
            vmax = np.max(hesslist[i])

        ax = hessimg(ax=ax, hess=hess, mag=mag, color=color, bins=bins, 
                     vmin=vmin, vmax=vmax, extent=extent, labels=labels,
                     photf_pts=photf_pts, mist_pts=mist_pts,
                     best_list=best_list, cmap=cmap, logcounts=logcounts,
                     ax_i=i, mode='series', savedir=savedir)

    if xlabel is not None:
        ind = 0
        if twobytwo:
            ind = 1
        _ = [ax.set_xlabel(xlabel) for ax in grid.axes_row[ind]]
        grid.axes_all[0].xaxis.label.set_visible(True)
    if ylabel is not None:
        _ = [ax.set_ylabel(ylabel) for ax in grid.axes_column[0]]

    return grid

def hessimg(ax, hess, extent, mag, color, bins, vmin, vmax, labels=None, 
            photf_pts=None, mist_pts=None, best_list=None, cmap=None, 
            ax_i=0, logcounts=False, xlabel=None, ylabel=None, 
            mode='single', cbar=True, ymag='V', savedir=""):

    """
       Draws a hess diagrams to a given axis.
    """

    i = ax_i

    if cmap is None:
        if i > 1:
            # bottom row: diff, sig
            colors = zeroed_cmap(hess)
        else:
            # first row: data, model. White will be on the left of color bar
            if i == 0:
                colors = plt.cm.Blues#copper#Blues
            if i == 1:
                colors = plt.cm.Oranges#Greys_r#Greys
            # colors = plt.cm.get_cmap('binary', 11)
    else:
        if isinstance(cmap, list):
            colors = cmap[i]
        else:
            colors = cmap
    if logcounts:
        hess = np.log10(hess)

    
    txtout_name = os.path.join(savedir, 
                               'hess{:d}_{:s}_bins{:d}_{:d}.txt'.format(i, 
                                                                 labels[i].replace("-2 ln P = ", "m2lnP"),
                                                                 bins[0],
                                                                 bins[1]))
    print("Saving ", txtout_name)
    np.savetxt(txtout_name,
                X=np.c_[color, mag, hess], delimiter='\t',fmt="%s")
    print("bins: ", bins)
    print(labels[i])

    #img = ax.imshow(hess, origin='upper', extent=extent,
    #                interpolation="nearest", cmap=colors)

    #nmagbin = len(mag)
    #ncolbin = len(color)

    print("i = ", i)
    print("vmin = ", vmin)
    print("vmax = ", vmax)

    if i > 1:
        h = ax.hist2d(color, mag, bins=bins, weights=hess,
                          cmap=colors, vmin=vmin, vmax=vmax)
    elif i <= 1:
        #if vmin == 0.0:
        #    vmin = 1e-6
        h = ax.hist2d(color, mag, bins=bins, weights=hess,
                          cmap=colors, norm=mpl.colors.LogNorm())


    # trying an affine transform to correct skew when ymag is I, not V.
   
#    if ymag == 'I':
        # skew angle in degrees; need to actually determine this
#        nmagbin = hess.shape[0]
#        ncolbin = hess.shape[1]

#        xpos = np.linspace(start = extent[0], stop = extent[1], num=ncolbin, endpoint=False)
#        for xp in xpos:
#            ax.axvline(x=xp)
#        ypos = np.linspace(start = extent[3], stop = extent[2], num=nmagbin, endpoint=False)
#        for yp in ypos:
#            ax.axhline(y=yp)

        #ax.scatter(max(xpos), max(ypos), s=500, c='r')
        #ax.scatter(min(xpos), max(ypos), s=500, c='r')
        #ax.scatter(max(xpos), min(ypos[hess[-1][:] == 0]), s=200, c='g')
#        ax.scatter(extent[1], extent[2], s=500, c='r')
#        ax.scatter(extent[0], extent[2], s=500, c='r')
#        ax.scatter(extent[1], extent[3], s=500, c='g')       
 
#        aside = np.abs(max(xpos) - min(xpos))
#        oside = np.abs(min(ypos[hess[:][-1] == 0]) - max(ypos))
#        hside = np.sqrt((min(ypos[hess[:][-1] == 0]) - max(ypos))**2 + (max(xpos) - min(xpos))**2)

#        yangle = -skewang#0.0#-np.arccos(aside/hside)*180/np.pi
#        xangle = 0.0

#        transform = mtransforms.Affine2D().skew_deg(xangle, yangle)
#        trans_data = transform + ax.transData
#        img.set_transform(trans_data)
#        x1, x2, y1, y2 = img.get_extent()

#        ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], 
#                "y--", transform=trans_data)       

    mpl_hack(ax)
    #ax.cax.colorbar(img)

    if mode == 'series' and cbar:
        #ax.cax.colorbar(img)
        ax.cax.colorbar(h[3])
    elif mode == 'single':
        if cbar:
            #plt.colorbar(img, ax=ax)
            plt.colorbar(h[3], ax=ax)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if xlabel is not None:
            ax.set_xlabel(xlabel)

    #mpl.rc('text',usetex=True)
 #   if labels is not None:
 #       _ = add_inner_title(ax, labels[i], loc=1)

        # SSG: Also adding best age & logz (to model plot):
 #       if i == 1 and mist_pts is not None:
            #print(best_list)
 #           assert best_list is not None, 'Need best_list for labels'
 #           _ = add_inner_title(ax,
 #                               r"Age = {:.3e}".format(10 ** best_list[0]),
 #                               loc=2)
 #           _ = add_inner_title(ax,
 #                               r"LogZ = {:.2f}".format(best_list[1]), loc=3)

    if labels is not None:
        _ = add_inner_title(ax, labels[i], loc=1)

        # SSG: Also adding best age & logz (to model plot):
        if (i == 1 | (i > 0 and mode == 'single')) & (best_list is not None):
            #print(best_list)
            assert best_list is not None, 'Need best_list for labels'
            _ = add_inner_title(ax,  "Age = {:.3e}".format(10**best_list[0]), loc=2)
            _ = add_inner_title(ax,  "LogZ = {:.2f}".format(best_list[1]), loc=3)

#    mpl.rc('text',usetex=True)
#    if labels is not None:
#        _ = add_inner_title(ax, labels[i], loc=1)

    # SSG: Also adding best age & logz (to model plot):
#    if i == 1 and mist_pts is not None:
        #print(best_list)
#        assert best_list is not None, 'Need best_list for labels'
#        _ = add_inner_title(ax,
#                            "Age = {:.3e}".format(10 ** best_list[0]),
#                            loc=2)
#        _ = add_inner_title(ax,
#                            "LogZ = {:.2f}".format(best_list[1]), loc=3)

#    mpl.rc('text', usetex=False)
    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax = square_aspect(ax)

    # Can add overplotting stuff here for the hess plots.
    if photf_pts is not None:
        ax.scatter(*photf_pts, color='r', lw=0, alpha=0.4, s=10, zorder=9999)
    # print(mist_pts)
    if mist_pts is not None:
        # if there are multiple sets of points to plot for the
        # mist model...
        #print(mist_pts)
        try:
            mist_pts[0][0]
            for n, ptset in enumerate(mist_pts):
                ls = '-' if n%2 else '--'
                ax.plot(*ptset, color='g', alpha=0.8, ls = ls)
        except IndexError:
            ax.plot(*mist_pts, color='g', alpha=0.8)

    return ax
