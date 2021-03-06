"""Class for reading the output.cmd file from calcsfh"""
from __future__ import print_function, absolute_import
import argparse
import os
import sys

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from match.scripts.config import EXT
from match.scripts.fileio import read_match_cmd
from match.scripts.graphics.match_plot import match_plot, hessimg
from match.scripts.utils import parse_pipeline
from match.scripts.wrappers.stats import call_stats
from match.scripts.likelihood import stellar_prob

import matplotlib as mpl

__all__ = ['CMD']


def splitcmds(filename, overwrite=False):
    name, ext = os.path.splitext(filename)

    headlen = 4

    with open(filename, 'r') as inp:
        lines = [l.strip() for l in inp.readlines()]

    filters = lines[2]
    _, x, y = np.array(lines[1].split(), dtype=int)

    nfilelines = len(lines) - headlen
    ncmdlines = x * y

    if ncmdlines != nfilelines:
        cmd0 = lines[:ncmdlines + headlen]
        cmd1 = [lines[0]]
        # more than one cmd.
        cmd1.extend(lines[ncmdlines + headlen:])
        filters1 = lines[ncmdlines + headlen:][1]
        cmd0fn = '{0:s}_{1:s}{2:s}'.format(name, filters, ext)
        cmd1fn = '{0:s}_{1:s}{2:s}'.format(name, filters1, ext)
        if overwrite or not os.path.isfile(cmd0fn):
            with open(cmd0fn, 'w') as outp:
                outp.write('\n'.join(cmd0))
        if overwrite or not os.path.isfile(cmd1fn):
            with open(cmd1fn, 'w') as outp:
                outp.write('\n'.join(cmd1))
    else:
        print('{0:s} appears to be a single cmd, doing nothing.'.format(filename))
        return filename, filename
    return cmd0fn, cmd1fn


class CMD(object):
    """
    An object to read the MATCH CMD file and hold paramters to
    automatically make plots with the same color scale as other MATCH CMD
    files.
    """

    def __init__(self, filename=None, onlyheader=False, params=None, ymag='V'):
        if filename is not None:
            self.base, self.name = os.path.split(os.path.abspath(filename))
            (self.cmd, self.fit, self.colors, self.yfilter,
             self.ncmd, self.nmagbin, self.ncolbin) = \
                read_match_cmd(filename, onlyheader=onlyheader, ymag=ymag)
            if not onlyheader:
                # self.cmd will be an empty array if onlyheader
                self.load_match_cmd()
            if params is not None:
                self.params = params
                dlage = np.array(
                    sorted(np.unique(params['log10 Age'].values))).diff()[0]
                dt = 10**(params['log10 Age'] + dlage) - \
                    10**params['log10 Age']
                self.mass = params['SFR'] * dt

    def set_labels(self):
        """Set up list of labels for pgpro"""
        #mpl.rc('text', usetex=True)
        strfmt = '{:s}'#r'${{\rm {:s}}}$'
        labels = [strfmt.format(i) for i in ['data', 'model', 'd-m']]
        try:
            target, _ = parse_pipeline(self.name)
            labels[0] = strfmt.format(target)
        except:
            pass
        #labels.append(r'$-2\ln P = {:g}$'.format(self.fit))
        labels.append('-2 ln P = {:g}'.format(self.fit))
        #mpl.rc('text', usetex=False)
        return labels

    def set_axis_labels(self, ax=None):
        xlabel = r'$\rm{{{}}}$'.format(
        #xlabel = '{}'.format(
            self.colors.replace('WFC', 'F').replace('UVIS', 'F').replace('Tycho_B', 'B_T').replace('Tycho_V', 'V_T'))
        ylabel = r'$\rm{{{}}}$'.format(
        #ylabel = '{}'.format(
            self.yfilter.replace('WFC', 'F').replace('UVIS', 'F').replace('Tycho_B', 'B_T').replace('Tycho_V', 'V_T'))
        if ax is not None:
            ax.set_xlabel(xlabel)
            ax.set_ylabel(ylabel)
        return xlabel, ylabel

    def load_cmd_fromfake(self, filename, dcol=0.05, dmag=0.1, ymag='V',
                          yfilter=None, colors=None, xlim=None, ylim=None):
        self.data = None
        self.diff = None
        self.sig = None
        self.colors = colors
        self.yfilter = yfilter

        mag1, mag2 = np.loadtxt(filename, unpack=True)
        inds, = np.nonzero((mag1 < 30) & (mag2 < 30))
        mag1 = mag1[inds]
        mag2 = mag2[inds]
        color = mag1 - mag2
        if ymag.upper() == 'V':
            mag = mag1
        elif ymag.upper() == 'I':
            mag = mag2
        else:
            print('Error: ymag needs to be V or I')
            return

        if xlim is None:
            cmin = np.min(color)
            cmax = np.max(color)
        else:
            cmin, cmax = xlim
        if ylim is None:
            mmax = np.max(mag)
            mmin = np.min(mag)
        else:
            mmin, mmax = np.sort(ylim)

        cbins = np.arange(cmin, cmax, dcol)
        mbins = np.arange(mmin, mmax, dmag)
        H, ce, me = np.histogram2d(color, mag, bins=[cbins, mbins])

        self.model = H.T
        self.extent = [cmin, cmax, mmax, mmin]
        self.nmagbin = len(mbins)
        self.ncolbin = len(cbins)

    def load_match_cmd(self):
        """
        pgcmd needs hesses and extent. Optional are max_* which set the vmins
        and vmaxs.
        """
        assert self.ncmd == 1, \
            '{0:s} Mutliple CMDs not supported, use splitcmds'.format(self.name)
        self.data=self.cmd['Nobs'].reshape(self.nmagbin, self.ncolbin)
        self.model=self.cmd['Nsim'].reshape(self.nmagbin, self.ncolbin)
        self.diff=self.cmd['diff'].reshape(self.nmagbin, self.ncolbin)
        self.sig=self.cmd['sig'].reshape(self.nmagbin, self.ncolbin)
        #self.hesses=[self.data, self.model, self.diff, self.sig]
        self.hesses=[self.cmd['Nobs'], self.cmd['Nsim'], self.cmd['diff'], self.cmd['sig']]

        self.extent=[np.min(self.cmd['color']), np.max(self.cmd['color']),
                       np.max(self.cmd['mag']), np.min(self.cmd['mag'])]
        # self.max_counts = np.nanmax(np.concatenate([self.data, self.model]))
        # self.max_diff = np.nanmax(np.abs(self.diff))
        # self.max_sig = np.nanmax(np.abs(self.sig))

    def recalc(self):
        """
            Recalculates the difference and signifigance (e.g. if CMD model values are changed).
        """

        # re-calc the difference (data - model):
        self.cmd['diff'] = self.cmd['Nobs'] - self.cmd['Nsim']
        # the significance (a bit heavy-handed):
        m2lnP, P, self.cmd['sig'] = stellar_prob(self.cmd['Nobs'], self.cmd['Nsim'])
        # replace the hess list with new hesses:
        self.hesses = [self.cmd['Nobs'], self.cmd['Nsim'], self.cmd['diff'], self.cmd['sig']]
        # replace fit (-2lnP) with newly calculated value:
        self.fit = np.sum(m2lnP)

    def resamp(self, hess_name):
        """
            Re-samples the hess diagram using a Poisson distribution. The argument 
        hess_name should be either 'Nobs' or 'Nsim'.
        """
        
        # re-sample:
        new_hess = np.random.poisson(self.cmd[hess_name].flat).reshape(self.cmd[hess_name].shape)
        # replace
        self.cmd[hess_name] = new_hess
        # re-calc hesses:
        self.recalc()

    def pgcmd(self, labels=None, outdir=None, logcounts=False, figname=None,
              twobytwo=True, sig=True, photf_pts=None, mist_pts=None,
              best_list=None, ymag='V'):
        '''produce the image that pgcmd.pro makes
        enhances graphics.match_plot.match_plot:
            automatic titles for each panel
            automatic axes labels
            add exclude/include gates
            logcounts only applies to data, model, not diff and sig.
        '''
        labels = labels or self.set_labels()

        if "/" in figname:
            savedir = "/".join(figname.split("/")[0:-1])
        else:
            print("No leading path in savedir.")
            savedir = ""
        print(savedir)
        
        if figname is None:
            base = outdir or self.base
            assert os.path.isdir(base), \
                '{} directory not found'.format(base)
            figname = os.path.join(base, os.path.split(self.name)[1] + EXT)

        hesses = self.hesses
        if logcounts:
            hesses[0] = np.log10(hesses[0])
            hesses[1] = np.log10(hesses[1])

        xlabel, ylabel = self.set_axis_labels()

        grid = match_plot(hesses, self.extent, self.cmd['mag'], self.cmd['color'], 
                          (self.ncolbin, self.nmagbin), labels=labels, ylabel=ylabel,
                          xlabel=xlabel, twobytwo=twobytwo, sig=sig,
                          photf_pts=photf_pts, mist_pts=mist_pts,
                          best_list=best_list, savedir=savedir)

        gates = self.cmd['gate']
        ugates = np.unique(gates)
        if len(ugates) > 1:
            dinds = np.digitize(gates, bins=np.unique(gates), right=True)
            _ = [grid.axes_all[0].plot(self.cmd['color'][dinds == i],
                                       self.cmd['mag'][dinds == i],
                                       '.', alpha=0.3)
                 for i in range(len(dinds)) if i == 0]

        for ax in grid.axes_all:
            ax.locator_params(axis='x', nbins=6)

        # save each hess diagram individually
        #f = plt.gcf()
        #for i, ax in enumerate(grid.axes_all):
    
            # Save just the portion _inside_ the axis' boundaries
        #    axbounds = ax.get_window_extent().transformed(f.dpi_scale_trans.inverted())

            # preserve any preceding path in the save file name
        #    svpath = os.path.join(*(figname.split('/')[:-1]), 'hess{:d}.png'.format(i))
            # Pad the saved area by 10% in the x-direction and 20% in the y-direction
        #    f.savefig(svpath, bbox_inches=axbounds.expanded(1.1, 1.2))

        plt.savefig(figname)
        plt.close()
        print('wrote {}'.format(figname))
        return grid

    def plthess(self, ax=None,labels=None, outdir=None, logcounts=False, figname=None, cmap=None,
               sig=True, photf_pts=None, mist_pts=None, best_list=None, hess_i=0, ymag='V', save=True,
               vmin=None, vmax=None):
        '''produce the image that pgcmd.pro makes
        enhances graphics.match_plot.match_plot:
            automatic titles for each panel
            automatic axes labels
            add exclude/include gates
            logcounts only applies to data, model, not diff and sig.
        '''
        hess_names = ['data', 'model', 'residual', 'sig']

        labels = labels or self.set_labels()

        if figname is None:
            base = outdir or self.base
            assert os.path.isdir(base), \
                '{} directory not found'.format(base)
            figname = os.path.join(base, os.path.split(self.name)[1]+ '_{:s}'.format(hess_names[hess_i]) + EXT)

        hess = self.hesses[hess_i]
        if logcounts & (hess_i == 0 or hess_i == 1):
            hess = np.log10(hess)

        xlabel, ylabel = self.set_axis_labels()

        if ax is None:
            fig = plt.figure(figsize=(16,9))
            ax = fig.add_subplot(111)

        if vmin == None:
            vmin = np.min(hess)
        if vmax == None:
            vmax = np.max(hess)

        #ax = hessimg(ax = ax, hess = hess, extent=self.extent, labels = labels, photf_pts = photf_pts,
        #                   mist_pts = mist_pts, best_list = best_list, cmap = cmap, logcounts = logcounts, 
        #                   ax_i = hess_i, mode = 'single', ylabel = ylabel, xlabel = xlabel)
        ax = hessimg(ax=ax, hess=hess, mag=self.cmd['mag'], color=self.cmd['color'], 
                     bins=(self.ncolbin, self.nmagbin), 
                     vmin=vmin, vmax=vmax, extent=self.extent, labels=None,
                     photf_pts=photf_pts, mist_pts=mist_pts, xlabel=xlabel, ylabel=ylabel,
                     best_list=best_list, cmap=cmap, logcounts=logcounts,
                     ax_i=hess_i, mode='single', cbar=False)

        gates = self.cmd['gate']
        ugates = np.unique(gates)
        if len(ugates) > 1:
            dinds = np.digitize(gates, bins=np.unique(gates), right=True)
            _ = [ax.plot(self.cmd['color'][dinds == i],
                                       self.cmd['mag'][dinds == i],
                                       '.', alpha=0.3)
                 for i in range(len(dinds)) if i == 0]

        #for ax in grid.axes_all:
        ax.locator_params(axis='x', nbins=6)

        if save:
            plt.savefig(figname)
            plt.close()
            print('wrote {}'.format(figname))

        return

def sortbyfit(cmdfns, onlyheader=False):
    """
    Order CMDs fit parameter given a list of CMDs or their filenames

    Parameters:

    cmdfns: list (of string filesnames or CMDs)
        if list of strings will initialize and return list of CMDs

    onlyheader : bool (False)
        passed to CMD.__init__() if cmdfns in list of strings.
        Only read the header of the CMD file (which contains fit parameter)

    Returns:
        CMDs or CMD filenames sorted by increasing fit parameter
    """
    assert isinstance(cmdfns, list), 'Need a list to sort. {}'.format(cmdfns)

    retv = np.array(cmdfns)
    if isinstance(cmdfns[0], str):
        cmds = np.array([CMD(cmdfn, onlyheader=onlyheader)
                         for cmdfn in cmdfns])
    else:
        cmds = cmdfns

    try:
        icmd = np.argsort(np.concatenate([cmd.fit for cmd in cmds]))
    except ValueError:
        icmd = np.argsort([cmd.fit for cmd in cmds])
    return retv[icmd]


def call_pgcmd_byfit(cmdfns, nmax=5, outdir=None, logcounts=False):
    """
    Call pgcmd with filename order them by increasing best fit value.
    cmdfns : string or list
        .cmd filename or list of filenames.

    nmax : int
        make best nmax plots.
    """
    if not isinstance(cmdfns, list):
        cmd = CMD(cmdfns)
        cmd.pgcmd(outdir=outdir, logcounts=logcounts)
        return

    cmds = sortbyfit(cmdfns)

    for j, cmd in enumerate(cmds):
        if j > nmax:
            break
        jstr = ('{}'.format(j)).zfill(4)
        cmd = CMD(cmd)
        cmd.pgcmd(figname='{}{}{}'.format(cmd.name, jstr, EXT),
                  outdir=outdir, logcounts=logcounts)
        # twobytwo=False, sig=False)
    return


def main(argv):
    """main function for cmd"""
    parser = argparse.ArgumentParser(description="plot cmd file")

    parser.add_argument('-o', '--outdir', type=str,
                        help='directory to place images')

    parser.add_argument('-f', '--byfit', action='store_true',
                        help='number filenames by best fit')

    parser.add_argument('-c', '--stats', action='store_true',
                        help='call match/bin/stats and exit')

    parser.add_argument('-s', '--split', action='store_true',
                        help='split cmd file if it contains more than one CMD')

    parser.add_argument('--logcounts', action='store_true',
                        help='use log binning for data and model')

    parser.add_argument('--nmax', type=int,
                        help='max best cmds to plot with --byfit')

    parser.add_argument('cmdfiles', type=str, nargs='*',
                        help='.cmd files to plot')

    args = parser.parse_args(argv)

    if args.outdir is not None:
        assert os.path.isdir(args.outdir), \
            'directory {} not found'.format(args.outdir)

    if args.stats:
        call_stats(args.cmdfiles, outdir=args.outdir)
        sys.exit()

    if args.split:
        cmds = np.unique(np.concatenate([splitcmds(c) for c in args.cmdfiles]))
        args.cmdfiles = cmds

    if args.byfit:
        if args.nmax is None:
            args.nmax = len(args.cmdfiles)
        call_pgcmd_byfit(args.cmdfiles, nmax=args.nmax,
                         outdir=args.outdir, logcounts=args.logcounts)
    else:
        for cmdfile in args.cmdfiles:
            print(cmdfile)
            cmd = CMD(filename=cmdfile)
            cmd.pgcmd(outdir=args.outdir, logcounts=args.logcounts)


if __name__ == "__main__":
    main(sys.argv[1:])
