#!/usr/bin/env python
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from match.scripts import ssp
from match.scripts import cmd
import numpy as np
import glob
import sys
import os
import shutil
import matplotlib.image as mpimg
from matplotlib.patches import Rectangle
import seaborn as sns
# no seaborn gridlines on plot, and white bg:
sns.set_context('paper')
sns.set(font='serif')
sns.set_style("white", {'font.family' : 'serif', 'font.serif': ['Times', 'Palatino', 'serif']})
plt.axes(frameon=False)
from MIST_codes.scripts import read_mist_models as rmm
from match.MISTscripts.read_fake import plt_truth
from match.MISTscripts.fileio import *
from TIGS import photplot as pp

params = {'axes.labelsize': 30,'axes.titlesize':20, 'font.size': 14, 'legend.fontsize': 14, 'xtick.labelsize': 14, 'ytick.labelsize': 14}
mpl.rcParams.update(params)

# Add something to perform diagnostic checks...
# + show that the results are accurate via e.g. plot errors and values as function of star number.
# + show sensitivity to completeness?
# + show effects of alteration of binsize
# + show effects of altering binary fraction? (maybe useful)

# General idea: show that match is accurate (or not) in smallish survey size regime that these clusters possess.


def marginalize(clustername, photobase, jobidstr, outname, params = ['lage', 'vvcrit', 'logZ'], 
                local=False, saveagevsvvc=True,debug=False, interp=True, vvc='all'):

    """
        This function calls certain methods from Philip Rosenfield/Daniel Weisz's MATCH scripts. Particularly,

        + ssp.combine_files(); to compile various v/vcrit MATCH runs into a single run where all best fits will be judged.
        + ssp.SSP(); to create a useable SSP object from the MATCH .scrn files.
            > methods therein of the SSP object are also called.

        + cmd.CMD(); to create a CMD object from the MATCH .cmd files.
            > mehthods therein of the CMD object are called as well.

        Arguments:
        ----------
        + clustername (str): Name of the cluster whose data will be utilized.
        
    """

    #global match_dir
    #global scrn_dir
    #global outpath
    #global plotpath

    # There are two directories here because on my local machine, I keep supercomp
    # output in a separate dir that I DL to, and local runs w/ MATCH are stored in
    # their own location.
    #if not local:
        # downloaded runs (wouldn't be used on remote machine):
    #    match_dir = os.path.join(os.environ['MATCH_ODYOUT'], clustername)
    #else:
        # local runs (should be used on remote machine, or your local machine):
    #    match_dir = get_workdir(data=clustername)#os.path.join(os.environ['MATCH_DIR'], 'MIST', clustername)

    # path to MATCH .scrn files:
    #scrn_dir = os.path.join(get_scrndir(data=clustername), photobase, 'SLURM{:s}'.format(jobidstr))
    # path to MATCH .cmd and .out files:
    #outpath = os.path.join(get_outdir(data=clustername), photobase, 'SLURM{:s}'.format(jobidstr))
    # path to place output plots in:
    #plotpath = os.path.join(outpath, 'plots')

    # list of .scrn files:
    if vvc != 'all':
        flist = glob.glob(os.path.join(scrn_dir, '*vvcrit{:.1f}*.scrn'.format(vvc)))
        params = ['lage', 'logZ']
    else:
        flist = glob.glob(os.path.join(scrn_dir, '*.scrn'))

    #print(scrn_dir)
    #print(flist)

    # print .scrn files found:
    #for afile in flist:
        #print("Reading {:s}...".format(afile))

    # calling MATCH scripts' combine_files() to combine various v/vcrit runs:
    ssp.combine_files(flist, outfile=outname)

    # Create a match SSP object from the combined scrns of the given runs:
    combineddata = ssp.SSP(outname)

    # Corner plots:
    # 'lage', 'vvcrit', 'logZ', 'dmod', 'Av'
    #pdffig, pdfax = combineddata.pdf_plots(params, twod=False, quantile=True, cmap=plt.cm.Reds, interpolateq=interp)
    #plt.close()
    #pdffig, pdfax = combineddata.pdf_plots(params, twod=True, quantile=True, cmap=plt.cm.Reds, interpolateq=interp)
    #plt.savefig(os.path.join(plotpath, 'pdfcornerplot.png'))
    #plt.close()

    # quantile dictionary (the quantiles will be used to draw bounding isochrones for the best-fits at ~+/- 1 sigma parameter values):
    #qdict = {key: value for (key, value) in zip(params, quant)}#{'lage': quant[0], 'vvcrit': quant[1], 'logZ':quant[2]}
    #print(qdict)

    # remove the combined .csv file since it's no longer needed:
    os.remove(os.path.join(os.getcwd(), outname))

    # A plot of the best fit ages found vs. the corresp. v/vcrit values...
    # This is also getting the name of the saved figure for the plot, the matplotlib figure for the plot,
    # the axis object too, the best parameters for each v/vcrit run and the corresp. v/vcrits:
    vvcritagename, bestdict, vvcrits = agevsvvcrit(flist, jobidstr, outpath = plotpath, save=saveagevsvvc) # vbests, vvcrits
    # Move the best fit ages vs v/vcrit plot:
    if saveagevsvvc:
        shutil.move(os.path.join(os.getcwd(), vvcritagename), os.path.join(outpath, vvcritagename))

    # Write a .txt file containing the best fit parameters found, along with their assoc. logP:
    combo_bestdict = {param: getbest(param, combineddata) for param in params}
    with open(os.path.join(plotpath, 'marginalized_params.txt'), 'w+') as f:
        for i, param in enumerate(params):
            #if param == 'lage':
            #    f.write('Best age: {:.2e}, logP = {:.4f}\n'.format(10**bestlist[i][0], bestlist[i][1]))
            #else:
            #    f.write('Best {:s}: {:.2f}, logP = {:.4f}\n'.format(param, bestlist[i][0], bestlist[i][1]))
            if param == 'Av':
                # Converting Av to E(B-V):
                f.write('Best E(B-V): {:.2f}, logP = {:.4f}\n'.format(combo_bestdict[param][0]/3.1, combo_bestdict[param][1]))
            else:
                f.write('Best {:s}: {:.2f}, logP = {:.4f}\n'.format(param, combo_bestdict[param][0], combo_bestdict[param][1]))

    #return vbests, vvcrits, qdict, bestdict#, vvcritagefig
    return vvcrits, bestdict, combo_bestdict

def getbest(param, data):

    """
       This method gets the best-fit (max likelihood) parameter value from an ssp object (SSP class defined 
       in Phil R. and Dan W.'s -- RW's --  match scripts).

       Argument:
       ---------
           + param (str): String for the parameter of interest (would be e.g.: lage, logZ, dmod, see 
              marginalize method of SSP class in RW's scripts)
           + data (SSP): RW script's SSP class instance.

       Returns:
       --------
           + best (float): best-fit value of the parameter corresponding to the string designated by the 
              'param' argument.
           + lnP (float): Associated log probability of the best-fit parameter.

    """
    # Index where the marginalized dist. for the given param is at a maximum probability.
    bestindex = np.where(data.marginalize(param)[1] == np.amax(data.marginalize(param)[1]))[0][0]
    # 0 indexes the parameters values (1 would index the corresponding lnP)
    # This is the max lnP parameter value.
    best, lnP = data.marginalize(param)
    best = best[bestindex]
    lnP = lnP[bestindex]

    return float(best), float(lnP)

# auxillary function for main(). Not majorly useful on its own.
def plotisocmd(ax, vvcrit, best_dict, filters,  extent, 
               txtlbl_kwargs = {"frameon":False,"prop":dict(size=14),"loc":4}, rot=False, texts=['vvc', 'age', 'feh']):

    # blue-red color, red, and blue mag names.
    color_name = '{:s}-{:s}'.format(filters[0], filters[1])
    redmag_name = filters[1]
    bluemag_name = filters[0]
    txtdict = {}
    anctxt=''
    if 'age' in texts:
        txtdict['age'] = 'age = {:.1f} Myrs'.format(10**best_dict[vvcrit]['lage']['val']/1e6)
    if 'feh' in texts:
        txtdict['feh'] = '[Fe/H] = {:.2f}'.format(best_dict[vvcrit]['feh']['val']) 
        if 'age' in texts:
            txtdict['age'] = txtdict['age'] + '\n'

    # best isochrone (for current v/vcrit):
    if rot:
        if 'vvc' in texts:
            txtdict['vvc'] = r'$P\left(\frac{\Omega}{\Omega_c} = 0.3\right)$'
            if len(list(txtdict.keys())) > 1:
                txtdict['vvc'] += '\n'
        #anctxt = (r'$P\left(\frac{\Omega}{\Omega_c} = 0.3\right)$' + '\n' \
        #                           'age = {:.1f} Myrs\n' \
         #                       '[Fe/H] = {:.2f}'.format(10**best_dict[vvcrit]['lage']['val']/1e6, 
            #                                             best_dict[vvcrit]['feh']['val']))
    else:
        if 'vvc' in texts:
            txtdict['vvc'] = r'$\frac{\Omega}{\Omega_c} = $' + '{:.1f}'.format(vvcrit)
            if len(list(txtdict.keys())) > 1:
                txtdict['vvc'] += '\n'
    for k in ['vvc', 'age', 'feh']:
        try:
            anctxt += txtdict[k]
        except KeyError:
            pass

    anctxt = (anctxt)
        #anctxt = (r'$\frac{\Omega}{\Omega_c} = $' + '{:.1f}\n' \
        #                           'age = {:.1f} Myrs\n' \
       #                         '[Fe/H] = {:.2f}'.format(vvcrit, 
        #                                                 10**best_dict[vvcrit]['lage']['val']/1e6, 
         #                                                best_dict[vvcrit]['feh']['val']))

    # the best fitting isochrone:
    iso = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'], vvcrit, ebv=best_dict[vvcrit]['ebv'], exttag='TP')
    iso.set_isodata(best_dict[vvcrit]['lage']['val'], color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'], texts=anctxt)

    # older, redder isochrone (1 sigma away in metallicity and age):
    isou = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'] + best_dict[vvcrit]['feh']['uperr'], vvcrit, 
                      ebv=best_dict[vvcrit]['ebv'], exttag='TP')
    isou.set_isodata(best_dict[vvcrit]['lage']['val'] + best_dict[vvcrit]['lage']['uperr'], 
                     color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])

    # younger, bluer isochrone (1 sigma away in [Fe/H], age but in other direction):
    isol = rmm.ISOCMD(best_dict[vvcrit]['feh']['val'] - best_dict[vvcrit]['feh']['loerr'], vvcrit, 
                                ebv=best_dict[vvcrit]['ebv'], exttag='TP')
    isol.set_isodata(best_dict[vvcrit]['lage']['val'] - best_dict[vvcrit]['lage']['loerr'], 
                     color_name, redmag_name, dmod=best_dict[vvcrit]['dmod'])

    # plot best, 1 sig older age/[Fe/H], 1 sig lower age/[Fe/H] isochrones.
    _ = iso.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, label=txtlbl_kwargs)
    _ = isou.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')
    _ = isol.isoplot(ax, xlims=extent[0:2], ylims=extent[2:], shade = vvcrit, alpha=0.6, ls='--')
                    

    # isol age:
    isolage = best_dict[vvcrit]['lage']['val'] - best_dict[vvcrit]['lage']['loerr']
    isoage = best_dict[vvcrit]['lage']['val']
    isouage = best_dict[vvcrit]['lage']['val'] + best_dict[vvcrit]['lage']['uperr']
    isodmod = best_dict[vvcrit]['dmod']

    # (x, y), i.e., (color, red mag) points of each isochrone in a list:
    mist_pts = [
                isoget_colmags(isol, [color_name, bluemag_name], lage=isolage, dmod=isodmod),
                isoget_colmags(iso, [color_name, bluemag_name], lage=isoage, dmod=isodmod),
                isoget_colmags(isou, [color_name, bluemag_name], lage=isouage, dmod=isodmod)
               ]

    xlims = extent[0:2]
    ylims = extent[2:]

    return ax, iso, isou, isol, mist_pts, (xlims, ylims)

# colmagnames is a list [color_name, magnitude_name]
# helper function for isoplotcmd; cumbersome on its own.
def isoget_colmags(iso, colmagnames, lage, dmod):

    c_dict = iso.get_data([colmagnames[0]], phasemask=[], lage=lage, dmod=dmod)
    m_dict = iso.get_data([colmagnames[1]], phasemask=[], lage=lage, dmod=dmod)

    # converts to values of color and mag: (col_val, mag_val)
    cm_pts = (list(c_dict.values())[0], list(m_dict.values())[0])

    return cm_pts

def main(cluster_name, photbase, jobid, filters, dmod, ebv, local, vvc, 
         outname=None, savebest=True, bestax=None, justbest=False, txtlbl_kwargs = {"frameon":False,"fontsize":14,"loc":4},
         texts=['vvc', 'age', 'feh'], hess_datapts = True, hess_modeliso = True):

    global match_dir
    global scrn_dir
    global outpath
    global plotpath

    # path/to/Cluster
    match_dir = get_workdir(data=cluster_name)
    # path to MATCH .scrn files:
    scrn_dir = os.path.join(get_scrndir(data=cluster_name), photbase, 'SLURM{:s}'.format(jobid))
    # path to MATCH .cmd and .out files:
    outpath = os.path.join(get_outdir(data=cluster_name), photbase, 'SLURM{:s}'.format(jobid))
    # path to place output plots in:
    plotpath = os.path.join(outpath, 'plots')
    #print(outpath)

    rot = False
    if 'rot' in cluster_name:
        rot = True

    if outname != None:
        # Create corner PDF plot:
        if vvc != 'all':
            flist = glob.glob(os.path.join(scrn_dir, '*vvcrit{:.1f}*.scrn'.format(vvc)))
            params = ['lage', 'logZ']
        else:
            flist = glob.glob(os.path.join(scrn_dir, '*.scrn'))

        ssp.combine_files(flist, outfile=outname)

        # Create a match SSP object from the combined scrns of the given runs:
        combineddata = ssp.SSP(outname)

        # Corner plots:
        # 'lage', 'vvcrit', 'logZ', 'dmod', 'Av'
        #pdffig, pdfax = combineddata.pdf_plots(['lage','logZ','vvcrit'], twod=False, quantile=True, cmap=plt.cm.Reds, interpolateq=True)
        #plt.close()
        #pdffig, pdfax = combineddata.pdf_plots(['lage','logZ','vvcrit'], twod=True, quantile=True, cmap=plt.cm.Reds, interpolateq=True)
        #plt.savefig(os.path.join(plotpath, 'pdfcornerplot.png'))
        #plt.close()

    #vbests, vvcrits, qdict, bfdict = marginalize(cluster_name, photbase, sys.argv[3], sys.argv[4], local=local)
    #vvcrits, best_dict, qdict, bfdict = marginalize(cluster_name, photbase, sys.argv[3], sys.argv[4], local=local, incvvc=incvvc)
    # by default, readssp() will return an AVERAGE (over all runs) of the bestfit parameter value and uncertainty.
    # the original, unaveraged values are derived from the sspcombine solutions.
    if vvc == 'all':
        vvcrits = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6]
    else:
        vvcrits = [vvc]
    #print(vvcrits)
    best_dict = {}
#    for vvc in vvcrits:

        # read sspcombine output for age results (subsampled solution being used):
#        sspbest_age, ssp_upageerr, ssp_loageerr = readssp(photname_base=photbase, data=cluster_name, 
#                                                                param='age', vvc=vvc, jobidstr=jobid)
        # read sspcombines for metallicity results:
#        sspbest_logz, ssp_uplogzerr, ssp_lologzerr = readssp(photname_base=photbase, data=cluster_name, 
#                                                                  param='logZ', vvc=vvc, jobidstr=jobid)

        # dictionary whose keys are the corresponding vvc values for the best fit values found above (dmod, Av fixed):
#        best_dict[vvc] = {'lage': {'val': sspbest_age, 'uperr': ssp_upageerr, 'loerr': ssp_loageerr}, 
#                          'feh': {'val': sspbest_logz, 'uperr': ssp_uplogzerr, 'loerr': ssp_lologzerr}, 
#                          'dmod': dmod, 
#                          'ebv': ebv}

    #plt.clf()
    # clearing out pre-existing output files that this script may have produced in prev. runs:
    if not justbest:
        for f in glob.glob(os.path.join(outpath,'cmd*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(plotpath,'cmd*.png')):
            os.remove(f)
        for f in glob.glob(os.path.join(plotpath,'*.txt')):
            os.remove(f)

    # if req., correct 2MASS filter names to those recognizable by MIST.
    for i, afilter in enumerate(filters):
        if afilter == 'J':
            filters[i] = '2MASS_J'
        elif afilter == 'Ks':
            filters[i] = '2MASS_Ks'
        elif afilter == 'B':
            filters[i] = 'Bessell_B'
        elif afilter == 'V':
            filters[i] = 'Bessell_V'


    # Retrieve the magnitude limits used in the fit (from the calcsfh .param file used):
    param_fname = glob.glob(os.path.join(match_dir, 'csfh_param', 'calcsfh_{:s}.param'.format(photbase)))[0]
    with open(param_fname, 'r') as f:
        lines = f.readlines()
        # Color limits:
        color_lims = lines[4].split(' ')
        minvmi = float(color_lims[3])
        maxvmi = float(color_lims[4])
        ibin = float(color_lims[0])
        vibin = float(color_lims[1])
        ilims = lines[6].split(' ')
        imin = float(ilims[1])
        imax = float(ilims[0])
        exclusions = lines[7].split('\n')[0]
        ntbins= int(lines[8])
        tbins=lines[9:8+ntbins]
        for i, aline in enumerate(tbins):
            tbins[i] = np.array(list(map(float, aline.split('\n')[0].split(' ')[-2:])))

    # Manage exclude gates (right now only works for one region):
    # If there are exclude gates...
    if int(exclusions[0]) > 0:
        # Get the exclude gate points:
        expts = np.array(list(map(float, exclusions.split(' ')[1:-1])))
        ex_x = expts[::2]
        ex_y = expts[1::2]
        ex_xy = zip(ex_x, ex_y)

    else:
        expts = None
        ex_xy = None

    # Also get the residuals (???):
    # For overplotting data points, get the data points CMD x, y values from the photometry file:
    photf_path = os.path.join(match_dir, 'photometry', '{:s}.phot'.format(photbase))
    photf_v_raw = np.genfromtxt(photf_path, usecols=(0,))
    photf_i_raw = np.genfromtxt(photf_path, usecols=(1,))
    if any(photf_v_raw >= 99.999) or any(photf_i_raw >= 99.999):
        print("Bad photometric values detected in {:s}. Cutting them out...".format(photf_path))
        photf_v = photf_v_raw[(photf_v_raw < 99.999) & (photf_i_raw < 99.999)]
        photf_i = photf_i_raw[(photf_v_raw < 99.999) & (photf_i_raw < 99.999)]
    else:
        photf_v = photf_v_raw
        photf_i = photf_i_raw

    photf_vmi = photf_v - photf_i
    # Tuple storing the x, y values for plotting the data on a CMD (using color on x, red filter on y):
    photf_pts = (photf_vmi, photf_i)
    # match uses blue filter on y axis.
    photf_pts_alt = (photf_vmi, photf_v)

    bestfit = 99999.9
    # Plot best isochrone for each v/vcrit individually:
    #=================================================================
    allfig = plt.figure()
    allax = allfig.add_subplot(111)

    probfig, probax = plt.subplots(1, 1)

    # iterate through all v/vcrits used for fitting:
    for j, vvcrit in enumerate(vvcrits):

        print(outpath)
        # get the .cmd file for making a MATCH pg style plot:
        if 'flat' in cluster_name:
            fname = glob.glob(os.path.join(outpath,'*vvcflat*.cmd'))[0]
        else:
            try:
                fname = glob.glob(os.path.join(outpath,'*vvc{:.1f}*.cmd'.format(vvcrit)))[0]
            except IndexError:
                fname = glob.glob(os.path.join(outpath,'*vvcrit{:.1f}*.cmd'.format(vvcrit)))[0]
        # use Phil's code to make a CMD object from the .cmd file:
        a_cmd = cmd.CMD(fname, ymag='I')
        # store the axis extents of the CMD data
        extent = a_cmd.extent

        # read sspcombine output for age results (subsampled solution being used):
        sspbest_age, ssp_upageerr, ssp_loageerr = readssp(photname_base=photbase, data=cluster_name, 
                                                                param='age', vvc=vvcrit, jobidstr=jobid)
        # read sspcombines for metallicity results:
        sspbest_logz, ssp_uplogzerr, ssp_lologzerr = readssp(photname_base=photbase, data=cluster_name, 
                                                                  param='logZ', vvc=vvcrit, jobidstr=jobid)

        # dictionary whose keys are the corresponding vvc values for the best fit values found above (dmod, Av fixed):
        best_dict[vvcrit] = {'lage': {'val': sspbest_age, 'uperr': ssp_upageerr, 'loerr': ssp_loageerr}, 
                              'feh': {'val': sspbest_logz, 'uperr': ssp_uplogzerr, 'loerr': ssp_lologzerr}, 
                              'dmod': dmod, 
                              'ebv': ebv,
                              'fit': a_cmd.fit}

        # to plot the fit statistic vs. v/vcrit:
        probax.scatter(vvcrit, a_cmd.fit)

        if True:#not justbest:   

            fig = plt.figure()
            ax = fig.add_subplot(111)

            # scatter plot of stars from data file (observations):
            ax = pp.plotphot(*photf_pts, ax=ax) 
            # X out masked data points.
            if expts is not None:
                for k in range(len(photf_pts[0])):
                    if (min(ex_x) < photf_vmi[k] < max(ex_x)) & (min(ex_y) < photf_v[k] < max(ex_y)):
                        ax.scatter(photf_vmi[k], photf_i[k], s=40, c='k', marker='x')  

            # Plot MIST isochrones corresponding to MATCH best-fit parameters.
            ax, iso, isou, isol, mist_pts, axlims = plotisocmd(ax=ax, vvcrit=vvcrit, best_dict=best_dict, filters=filters, 
                                                               extent=extent, txtlbl_kwargs = txtlbl_kwargs, rot=rot, texts=texts)

            ax.set_title('{:s}'.format(cluster_name.split('rot')[0]))
            # saving the plot of data points w/ isochrones overlaid:
            savename = os.path.join(plotpath, 'cmd_vvc{:.1f}_m2lnP{:.2f}.png'.format(vvcrit, a_cmd.fit))
            fig.savefig(savename, dpi=300)

            # create a MATCH pg style plot using the .cmd file:
            # using photf_pts_alt because MATCH does its CMDs with V-I vs. V & not sure if this is changeable.
            pgcmd_kwargs = {}
            if hess_datapts:
                pgcmd_kwargs['photf_pts'] = photf_pts_alt
            if hess_modeliso:
                pgcmd_kwargs['mist_pts'] = mist_pts
            pgcmd_kwargs['best_list'] = [best_dict[vvcrit]['lage']['val'], best_dict[vvcrit]['feh']['val']]
            pgcmd_kwargs['ymag'] = 'V'           

            # four panel plot of all hess diagrams:
#            a_cmd.pgcmd(**pgcmd_kwargs)
            for m in range(4):
                # each of the 4 plotted separately:
                pgcmd_kwargs['hess_i'] = m
                a_cmd.plthess(**pgcmd_kwargs)
        
            # closing figure before moving on to next plots.
            fig.clf()

        # track which cmd has the overall highest lnP & corresponding v/vcrit.
        if a_cmd.fit < bestfit:

            bestfit = a_cmd.fit
            bestvvc = vvcrit
            bestcmd = a_cmd
  
    # save plot of fit statistic vs. v/vc:
    probax.set_xlabel(r'$\Omega/\Omega_c$')
    probax.set_ylabel(r'$-2\ln$P')
    probfig.savefig(os.path.join(plotpath, 'lnpvsrot.png'))

    # write best fit solutions +/- uncertainties to a summary file:
    # -------------------------------------------------------------
    print('WRITING RESULTS...')
    outlines = ['Hyades\n','_______\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n',
                '2MASS J, Ks\n','-------\n', '\n', '\n', 
                '\n','Praesepe\n','________\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n', 
                '2MASS J, Ks\n','-------\n', '\n', '\n',
                '\n','Pleiades\n_','_______\n', 
                'Tycho B, V\n','--------\n', '\n', '\n', '\n', 
                '2MASS J, Ks\n','-------\n', '\n', '\n', '\n']

    # writing results to a .txt file for reference (don't add fake data results though):
    if 'Fakedata' not in cluster_name:
        with open('/n/home12/sgossage/match_results.txt', 'r+') as sumfile:

            # get current lines of the summary file:
            inlines  = sumfile.readlines()
#            print(inlines)
            #print(len(inlines))
            #print(len(outlines))

            for i, line in enumerate(inlines):
                # maintain lines that already contain best fits
                if ('Best' in line) | ('v/vcrit = 0.0' in line):
                    print(line)
                    outlines[i] = line

            # update best fit line for the current run:
            if 'Hyades' in cluster_name:
                n = 4
            elif 'Praesepe' in cluster_name:
                n = 16
            elif 'Pleiades' in cluster_name:
                n = 28

            if '2MASSJK' in photbase:
                n += 5

            # print these lines when updating...
            # ...if using a distribution of roation rates...
            if ('_rot' in cluster_name) | ('_flat' in cluster_name):
                n += 2
                outlines[n] = 'Best (rotation distribution) v/vcrit = {:.1f}, ' \
                              'age = {:.1f} + {:.1f} - {:.1f} Myrs,' \
                              '[Fe/H] = {:.2f} + {:.2f} - {:.2f}, -2lnP ={:.2f} \n'.format(bestvvc, 
                                                                                (10**best_dict[bestvvc]['lage']['val'])/1e6,
                                                                                (10**(best_dict[bestvvc]['lage']['val'] + \
                                                                                best_dict[bestvvc]['lage']['uperr']) - \
                                                                                10**best_dict[bestvvc]['lage']['val'])/1e6,
                                                                                (10**best_dict[bestvvc]['lage']['val'] - \
                                                                                10**(best_dict[bestvvc]['lage']['val'] - \
                                                                                best_dict[bestvvc]['lage']['loerr']))/1e6,
                                                                                best_dict[bestvvc]['feh']['val'],  
                                                                                best_dict[bestvvc]['feh']['uperr'],
                                                                                best_dict[bestvvc]['feh']['loerr'],
                                                                                best_dict[bestvvc]['fit'])

            # ...or else if not using a distribution of rotation rates.
            else:
                outlines[n] = 'Best v/vcrit = {:.1f}, ' \
                              'age = {:.1f} + {:.1f} - {:.1f} Myrs,' \
                              '[Fe/H] = {:.2f} + {:.2f} - {:.2f}, -2lnP = {:.2f}\n'.format(bestvvc, 
                                                                            (10**best_dict[bestvvc]['lage']['val'])/1e6, 
                                                                            (10**(best_dict[bestvvc]['lage']['val'] + \
                                                                            best_dict[bestvvc]['lage']['uperr']) - \
                                                                            10**best_dict[bestvvc]['lage']['val'])/1e6,
                                                                            (10**best_dict[bestvvc]['lage']['val'] - \
                                                                            10**(best_dict[bestvvc]['lage']['val'] - \
                                                                            best_dict[bestvvc]['lage']['loerr']))/1e6,
                                                                            best_dict[bestvvc]['feh']['val'], 
                                                                            best_dict[bestvvc]['feh']['uperr'],
                                                                            best_dict[bestvvc]['feh']['loerr'],
                                                                            best_dict[bestvvc]['fit'])
                n += 1
                outlines[n] = 'v/vcrit = {:.1f}, ' \
                              'age = {:.1f} + {:.1f} - {:.1f} Myrs,' \
                              '[Fe/H] = {:.2f} + {:.2f} - {:.2f}, -2lnP = {:.2f}\n'.format(0.0, 
                                                                                (10**best_dict[0.0]['lage']['val'])/1e6, 
                                                                                (10**(best_dict[0.0]['lage']['val'] + \
                                                                                best_dict[0.0]['lage']['uperr']) - \
                                                                                10**best_dict[0.0]['lage']['val'])/1e6,
                                                                                (10**best_dict[0.0]['lage']['val'] - \
                                                                                10**(best_dict[0.0]['lage']['val'] - \
                                                                                best_dict[0.0]['lage']['loerr']))/1e6,
                                                                                best_dict[0.0]['feh']['val'], 
                                                                                best_dict[0.0]['feh']['uperr'],
                                                                                best_dict[0.0]['feh']['loerr'],
                                                                                best_dict[0.0]['fit'])

            # return to top of file & then write lines out.
            sumfile.seek(0)
            sumfile.writelines(outlines)
    #=========================================================


    if bestax == None:
        bestfig = plt.figure()
        bestax = bestfig.add_subplot(111)
        bestax.set_title('{:s}'.format(cluster_name.split('rot')[0]))

    #bestax.scatter(*photf_pts, lw=0, s=8, c='r')
    bestax = pp.plotphot(*photf_pts, ax=bestax) 
#    if expts != None:
##        for k in range(len(photf_pts[0])):
#            if (min(ex_x) < photf_vmi[k] < max(ex_x)) & (min(ex_y) < photf_v[k] < max(ex_y)):
#                bestax.scatter(photf_vmi[k], photf_i[k], s=40, c='k', marker='x') 

    bestax, bestiso, bestisou, bestisol, mist_pts, axlims = plotisocmd(ax=bestax, vvcrit=bestvvc, best_dict=best_dict, 
                                                                       filters=filters, extent=extent, txtlbl_kwargs = txtlbl_kwargs, rot=rot, texts=texts)

    if savebest:
        # save the figure.
        bestfig.savefig(os.path.join(plotpath, 'bestcmd_vvcrit{:.1f}_m2lnP{:.2f}.png'.format(vvcrit, a_cmd.fit)), dpi=300)
        bestfig.clf()

    return bestax, allax, (bestiso, bestisou, bestisol), bestcmd, photf_pts, mist_pts, ex_xy, axlims

if __name__ == "__main__":

    """
        This is a script runnable as e.g.

        >> ./plotruns.py Hyades hyades_debruTIGS_TychoBV.phot 12345678 somename.csv Tycho_B Tycho_V 
     
        The intent is to plot the best-fit models found via a MATCH run; this makes use of Phil R./Dan 
        W.'s (RW) match scripts and also Jieun Choi's MIST codes (available on her Github site). The 
        former is mainly used to conglomerate MATCH output across separate runs made at various v/vcrit 
        values (as this parameter is not explored in the usual way by MATCH, RW's scripts concatenate 
        the output of the runs and marginalize to acheive a best-fit over any new parameters not used 
        within MATCH itself, e.g. v/vcrit). The latter is used to plot the best-fit MIST models found 
        by the run(s).

    """

# NOTE:
# sys.argv[x] = cluster name, photometry file, job array slurm id, outputfilename, bluer band filter, redder band filter, dmod, E(B-V), and then optional flags.
    cluster_name = sys.argv[1]
    photbase = os.path.splitext(sys.argv[2])[0]
    jobid = sys.argv[3]
    outname=sys.argv[4]

    # blue and red filters:
    v_filter = sys.argv[5]
    i_filter = sys.argv[6]
    # organize filter inputs from user (blue, red is the expected order):
    filters = [v_filter, i_filter]

    # dmod and E(B-V) used:
    dmod = float(sys.argv[7])
    ebv = float(sys.argv[8])

    # check for script options from user:
    if '-local' in sys.argv[9:]:
        local = True
    else:
        local = False
    if '-truth' in sys.argv[9:]:
        truth = True
    else:
        truth = False

    for arg in sys.argv[9:]:
        print('.')
        if '-vvcin=' in arg:
            vvc = float(arg.split('-vvcin=')[-1])
            break
        else:
            vvc = 'all'

    _ = main(cluster_name=cluster_name, photbase = photbase, jobid = jobid, filters = filters, dmod = dmod, ebv = ebv,
                     local=local, vvc = vvc, outname=outname)
