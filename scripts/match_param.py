"""
Make a match param file
Note -- will have to edit the param file by hand to insert dmod and Av.
"""
from __future__ import print_function
import argparse
import os
import sys
import time

import numpy as np
import matplotlib.pylab as plt

from .config import PARAMEXT
from .fileio import read_calcsfh_param, calcsfh_input_parameter, match_filters
from .utils import replaceext, parse_pipeline
from .match_phot import make_phot
from .graphics.match_diagnostic import match_diagnostic


def move_on(okay, msg='0 to move on: '):
    """read and return raw input"""
    okay = int(input(msg))
    time.sleep(1)
    return okay


def within_limits(params, fakefile, offset=1.):
    """
    Cull param cmd limits to that of the fake file
    params : dict
        calcsfh_input_parameter dictionary (only need CMD limits)
    fakefile : string
        match AST file
    offset : float
        mag below
    """
    vimin = params['vimin']
    vimax = params['vimax']
    vmin = params['vmin']
    imin = params['imin']
    vmax = params['vmax']
    imax = params['imax']

    mag1in, mag2in, _, _ = np.loadtxt(fakefile, unpack=True)
    colin = mag1in - mag2in
    msg = 'Overwrote'
    if vimin < colin.min():
        vimin = colin.min()
        msg += ' vimin'
    if vimax > colin.max():
        vimax = colin.max()
        msg += ' vimax'
    if vmin < mag1in.min():
        vmin = mag1in.min()
        msg += ' vmin'
    if vmax > mag1in.max() - 1.:
        vmax = mag1in.max() - 1.
        msg += ' vmax'
    if imin < mag2in.min():
        imin = mag2in.min()
        msg += ' imin'
    if imax > mag2in.max() - 1:
        imax = mag2in.max() - 1.
        msg += ' imax'
    msg += ' with values from matchfake'
    print(msg)
    params['vimin'] = vimin
    params['vimax'] = vimax
    params['vmin'] = vmin
    params['imin'] = imin
    params['vmax'] = vmax
    params['imax'] = imax
    return params


def find_match_limits(mag1, mag2, comp1=90., comp2=90., color_only=False,
                      xlim=None, ylim=None):
    """
    click color limits on a cmd and mag1 mag2 limits on a plot of mag1 vs mag2
    """
    col = mag1 - mag2

    _, ax = plt.subplots()
    ax.plot(col, mag2, 'o', color='k', ms=3, alpha=0.3, mec='none')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ax.get_ylim()[::-1])

    if comp1 < 90.:
        ax.hlines(comp2, *ax.get_xlim())

    okay = 1
    while okay == 1:
        print('click color extrema')
        pts = plt.ginput(2, timeout=-1)
        colmin, colmax = [pts[i][0] for i in range(2)]
        if colmin > colmax:
            colmin, colmax = colmax, colmin
        ax.vlines(colmin, *ax.get_ylim())
        ax.vlines(colmax, *ax.get_ylim())
        plt.draw()
        okay = move_on(0)

    plt.close()

    inds, = np.nonzero((col < colmax) & (col > colmin))
    data = (colmin, colmax)
    if not color_only:
        _, ax = plt.subplots()
        ax.plot(mag1, mag2, '.', color='k')
        okay = 1
        while okay == 1:
            print('click mag extrema')
            pts = plt.ginput(2, timeout=-1)
            mag1max, mag2max = pts[0]
            mag1min, mag2min = pts[1]
            if mag1min > mag1max:
                mag1min, mag1max = mag1max, mag1min
            if mag2min > mag2max:
                mag2min, mag2max = mag2max, mag2min

            ax.plot(mag1max, mag2max, 'o', color='r')
            ax.plot(mag1min, mag2min, 'o', color='r')
            plt.draw()
            okay = move_on(okay)

        plt.close()

        inds, = np.nonzero((mag1 < mag1min) & (mag1 > mag1max) &
                           (mag2 < mag2min) & (mag2 > mag2max) &
                           (col < colmax) & (col > colmin))

    _, ax = plt.subplots()
    ax.plot(col, mag2, '.', color='k')
    ax.plot(col[inds], mag2[inds], '.', color='r')
    ax.set_ylim(ax.get_ylim()[::-1])
    if comp2 < 90.:
        ax.hlines(comp2, *ax.get_xlim(), lw=2)
    ax.vlines(colmin, *ax.get_ylim(), lw=2)
    ax.vlines(colmax, *ax.get_ylim(), lw=2)
    if not color_only:
        ax.hlines(mag2max, *ax.get_xlim(), lw=2)
        data = (colmin, colmax, mag1min, mag1max, mag2min, mag2max)

    plt.draw()

    print(data)
    return data


def find_gates(mag1, mag2, param):
    """Click 4 points to make an exclude gate -- does not work in calcsfh!"""
    print('not supported')
    sys.exit()
    col = mag1 - mag2

    lines = open(param, 'r').readlines()
    colmin, colmax = map(float, lines[4].split()[3:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])
    # mag2min, mag2max = map(float, lines[5].split()[:-1])
    # click around
    _, ax = plt.subplots()
    ax.plot(col, mag2, ',', color='k', alpha=0.2)
    ax.set_ylim(mag1max, mag1min)
    ax.set_xlim(colmin, colmax)

    okay = 1
    while okay != 0:
        print('click ')
        pts = np.asarray(plt.ginput(n=4, timeout=-1))
        exclude_gate = '1 {} 0 \n'.format(' '.join(['%.4f' % p
                                                    for p in pts.flatten()]))
        pts = np.append(pts, pts[0]).reshape(5, 2)
        ax.plot(pts[:, 0], pts[:, 1], color='r', lw=3, alpha=0.3)
        plt.draw()
        okay = move_on(0)
    lines[7] = exclude_gate
    # not so simple ... need them to be parallelograms.
    # PASS!

    # write new param file with exclude/include gate
    os.system('mv {0} {0}_bkup'.format(param))
    with open(param, 'w') as outp:
        _ = [outp.write(l) for l in lines]
    print('wrote %s' % param)


def check_filters(filters, hstflag=None):
    """
    Assign filters, check against templates/match_filters.json, and
    add HST instrument if format is from UW pipline filename.

    Paramters
    ---------
    filters : str list
        str list of vfilter, ifilter -- only two filter are currently supported

    hstflag : str (wfc, uvis, or hrc)
        UW pipeline filenames follow format of F???W no matter WFC3, HRC, or
        UVIS.
        However in MATCH, F???W specifically selects WFPC2. Pass this flag
        if the filter is in format F???W but NOT from WFPC2 instrument.
        Otherwise, pass WFC???W, UVIS???W, HRC???W format.

    Returns
    -------
    vfilter, ifilter : (str, str)
        filters selected for calcsfh
    """
    hstflag = hstflag or ''
    vfilter = filters[0]
    ifilter = filters[1]

    if 'wfc' in hstflag.lower():
        cam = 'WFC'
    elif 'uvis' in hstflag.lower():
        cam = 'UVIS'
    elif 'hrc' in hstflag.lower():
        cam = 'HRC'
    else:
        cam = ''

    if len(cam) > 0:
        vfilter = vfilter.replace('F', cam)
        ifilter = ifilter.replace('F', cam)

    if vfilter.startswith('F') and vfilter.endswith('W'):
        cam = 'WFPC2'

    print('filters set to {0:s}, {1:s}'.format(vfilter, ifilter))
    if len(cam) > 0:
        print('assuming HST instrument: {0:s}'.format(cam))

    return vfilter, ifilter


def match_param(mag1, mag2, filters, param, interactive=False, fake=None,
                comp_frac=0.5, param_kw=None, power_law_imf=False, zinc=False,
                bright_lim=20., overwrite=False, hstflag=None, max_tbins=100):
    """
    Make match param file

    Will check filter list against templates/match_filters.json
    Will check the CMD limits against the AST limits (if fake is passed)

    Parameters
    ----------
    mag1, mag2 : array, array
        v mag and i mag (extrema used for CMD limits)

    filters : list of strings
        v and i filter names

    param : string
        template parameter file or if clobber parameter file name

    interactive : bool
        choose cmd limits interactively [not tested in py3+]

    fake : string
        matchfake filename if using comp_frac or want to check CMD limits
        are within fake limits (see FK overflow in MATCH README)

    comp_frac : float
        completeness fraction to set faint mag limit

    param_kw : dict
        parameters of template/calcsfh_input_parameter.json to write

    power_law_imf : bool
        passed to fileio.calcsfh_input_parameter

    zinc : bool
        passed to fileio.calcsfh_input_parameter

    bright_lim : float
        passed to asts.ast.get_completeness_fraction

    overwrite : bool
        overwrite param file if exisiting
    """
    def write_param(fn, pstr):
        with open(fn, 'w') as out:
            out.write(pstr)
        print('wrote {}'.format(fn))
        return

    param_kw = param_kw or {}

    if os.path.isfile(param) and not overwrite:
        template = read_calcsfh_param(param)
        if param_kw['tbin'] is not None:
            print(template['ntbins'])
            del template['ntbins']
        template.update(param_kw)
    else:
        template = param_kw

    template['v'], template['i'] = check_filters(filters, hstflag=hstflag)

    if interactive:
        vimin, vimax, vmin, vmax, imin, imax = find_match_limits(mag1, mag2)
    else:
        color = mag1 - mag2
        vimin = np.min(color)
        vimax = np.max(color)
        vmin = np.min(mag1)
        imin = np.min(mag2)
        vmax = template['vmax']
        imax = template['imax']
        dove = 'user input'
        if vmax is None or imax is None:
            if fake is not None:
                from .asts import ASTs
                ast = ASTs(filename=fake, filter1=filters[0],
                           filter2=filters[1])
                ast.completeness(combined_filters=True, interpolate=True)
                print('Using {0:f} completeness fraction from {1:s}'
                      .format(comp_frac, fake))
                vmax, imax = \
                    ast.get_completeness_fraction(comp_frac,
                                                  bright_lim=bright_lim)
                dove = 'completeness'
            else:
                vmax = np.max(mag1)
                imax = np.max(mag2)
                dove = 'data'

        print('From {}: vmax={} imax={}'.format(dove, vmax, imax))

    template['vimin'] = vimin
    template['vimax'] = vimax
    template['vmin'] = vmin
    template['imin'] = imin
    template['vmax'] = vmax
    template['imax'] = imax
    # if fake color and mag limits within color and mag limits above:
    if fake is not None:
        template = within_limits(template, fake)

    param_files = calcsfh_input_parameter(power_law_imf=power_law_imf,
                                          zinc=zinc, max_tbins=max_tbins,
                                          **template)

    param_files = np.atleast_1d(param_files)
    params = []
    for i, pstr in enumerate(param_files):
        if i == 0:
            fn = param
        else:
            fn = param.replace(PARAMEXT, '_par{0:d}{1:s}'.format(i, PARAMEXT))
        write_param(fn, pstr)
        params.append(fn)
    return params


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description="make calcsfh param file",
                                     fromfile_prefix_chars='@')

    parser.add_argument('--imf', default=None,
                        help='IMF power law (None if using calcsfh flag)')

    parser.add_argument('--bf', type=float, default=0.0,
                        help='Binary fraction [0.0]')

    parser.add_argument('--tbin', type=float, default=0.05,
                        help='age bin width(s) [0.05]')

    parser.add_argument('--vstep', type=float, default=0.15,
                        help='mag step size [0.15]')

    parser.add_argument('--vistep', type=float, default=0.05,
                        help='color step size [0.05]')

    parser.add_argument('--tmin', type=float, default=6.6,
                        help='min log age [6.6]')

    parser.add_argument('--tmax', type=float, default=10.24,
                        help='max log age [10.24]')

    parser.add_argument('--vmax', type=float, default=None,
                        help='faint V limit ')

    parser.add_argument('--imax', type=float, default=None,
                        help='faint I limit')

    parser.add_argument('--zinc', action='store_true',
                        help='use zinc [False]')

    parser.add_argument('--dmod', type=float, nargs=2, default=[10., 10.],
                        help='dmod0, dmod1')

    parser.add_argument('--av', type=float, nargs=2, default=[0.0, 0.0],
                        help='av0, av1')

    parser.add_argument('--dav', type=float, default=0.05,
                        help='Av step -- NOT -dAv flag [0.05]')

    parser.add_argument('--ddmod', type=float, default=0.10,
                        help='dmod step [0.1]')

    parser.add_argument('-s', '--slice', type=float, default=40.,
                        help='cut out mags outside of this value [40.]')

    parser.add_argument('-i', '--interactive', action='store_true',
                        help='find cmd limits interactively')

    parser.add_argument('-f', '--filters', type=str, default=None,
                        help=('comma separated filter names (if filename does '
                              'not follow: PID_TARGET_FILTER1_FILTER2.ext)'))

    parser.add_argument('--fake', type=str, default=None,
                        help=('match fake file to calculate completeness '
                              'and to set faint mag limits of param file'))

    parser.add_argument('-c', '--comp_frac', type=float, default=0.50,
                        help=('completeness fraction as faint mag limit '
                              '(use with --fake=) [0.50]'))

    parser.add_argument('-b', '--bright_lim', type=float, default=20,
                        help=('bright limit to consider for calculating '
                              'completeness (use with --fake=) [20]'))

    parser.add_argument('-p', '--param', type=str,
                        help='template match param file')

    parser.add_argument('--hstflag', type=str,
                        help='if HST filters: uvis or wfc')

    parser.add_argument('--overwrite', action='store_true',
                        help='overwrite')

    parser.add_argument('--max_tbins', type=int, default=100,
                        help=('maximum time bins per param file '
                              '(would create more parameter files)'))

    parser.add_argument('--phot', type=str, help='photometry file match or fits')

    return parser.parse_args(argv)


def main(argv=None):
    """main function for match_param"""
    args = parse_args(argv)

    assert(not isinstance(args.imf, str)), 'Only set IMF if it is a power law.'
    if args.phot.endswith('fits'):
        args.phot = make_phot(args.phot)[0]

    mag1, mag2 = np.loadtxt(args.phot, unpack=True)
    inds, = np.nonzero((np.abs(mag1) < args.slice) &
                       (np.abs(mag2) < args.slice))
    mag1 = mag1[inds]
    mag2 = mag2[inds]

    args.param = args.param or replaceext(args.phot, PARAMEXT)

    if args.filters is None:
        try:
            _, filters = parse_pipeline(args.phot)
        except IndexError:
            print("Could not read filters from filename")
            filters = args.filters
            pass
    else:
        filters = args.filters.split(',')
    if args.hstflag is None:
        if 'uvis' in args.phot.lower():
            args.hstflag == 'uvis'
        elif 'wfc' in args.phot.lower():
            args.hstflag == 'wfc'

    if not os.path.isfile(args.param) or args.overwrite:
        print('Making param file')
        param_kw = {'imf': args.imf,
                    'bf': args.bf,
                    'tbin': args.tbin,
                    'vstep': args.vstep,
                    'vistep': args.vistep,
                    'tmin': args.tmin,
                    'tmax': args.tmax,
                    'dmod0': args.dmod[0],
                    'dmod1': args.dmod[1],
                    'ddmod': args.ddmod,
                    'av0': args.av[0],
                    'av1': args.av[1],
                    'dav': args.dav,
                    'vmax': args.vmax,
                    'imax': args.imax}

        power_law_imf = True
        if args.imf is None:
            power_law_imf = False

        match_param(mag1, mag2, filters, args.param, fake=args.fake,
                    interactive=args.interactive, comp_frac=args.comp_frac,
                    param_kw=param_kw, power_law_imf=power_law_imf,
                    zinc=args.zinc, bright_lim=args.bright_lim,
                    overwrite=args.overwrite, hstflag=args.hstflag,
                    max_tbins=args.max_tbins)

    else:
        print('{} file found, not overwriting'.format(args.param))

    match_diagnostic(args.param, args.phot, fake=args.fake)


if __name__ == "__main__":
    sys.exit(main())
