import os
import glob
import numpy as np

def get_workdir(data='Fakedata'):

    # Environment path to /n/home12/sgossage/match2.6
    matchdir_envkey = 'MATCH_DIR'
    matchdir = os.environ.get(matchdir_envkey)

    if matchdir != None:

        # append '_rot' to a data directory name to tell this function to access MISTrot (rotating models), rather than MIST.
        if '_rot' == data[-4::]:
            mistdir = os.path.join(matchdir, "MISTrot")
            data = data.split('_rot')[0]
        elif '_flat' == data[-5::]:
            mistdir = os.path.join(matchdir, "MISTrot_flat")
            data = data.split('_flat')[0]
        else:
            mistdir = os.path.join(matchdir, "MISTrot_nodist")

        workdir = os.path.join(mistdir, data)
        #print(workdir)
    
        # returns the path to the work directory...e.g. match2.6/MIST/Hyades:
        return  os.path.join(mistdir, data) 
    
    else:
        print("MATCH_DIR environment variable (path to e.g. ../matchx.xx, where x.xx = version #) not set.)")
        return

def get_photodir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'photometry')

def get_scrndir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'scrns')

def get_outdir(data='Fakedata'):

    return os.path.join(get_workdir(data), 'output')

def get_photof(photbase, data='Fakedata'):

    return glob.glob(os.path.join(get_photodir(data), photbase+'.phot'))[0]

def get_jobids(data, photbase):
    
    # searches in the scrn directory for SLURM jobid subdirectories and extracts the jobids present: 
    return [direc.split('SLURM')[-1] for direc in glob.glob(os.path.join(get_scrndir(data), photbase, 'SLURM*'))]

def readssp(photname_base, data='Fakedata', SFR='0.001', param='age', sub=None, vvc='all', jobidstr=None):

    # gets the work dir, e.g., proper/path/to/Fakedata
    workdir = get_workdir(data)

    if param=='lage':
        param = 'age'

    # The ssp output directory is used to extract the solutions found by MATCH.
    sspoutdir = os.path.join(workdir, "sspcombines")

    # If desired, look in a specified subdirectory for the best fits:
    if sub != None:
        sspoutdir = os.path.join(sspoutdir, sub)

    # Now go to the sspcombine output directories and get ages, etc.
    if jobidstr == None:
        # use all SLURM* subdirs if a jobid isn't specified:
        sspoutdir = os.path.join(sspoutdir, photname_base, 'SLURM*')
    else:
        sspoutdir = os.path.join(sspoutdir, photname_base, 'SLURM{:s}'.format(jobidstr))

    print("Reading from {:s}...".format(sspoutdir))


    if 'flat' in data:
        sspfiles = glob.glob(os.path.join(sspoutdir, "*vvcflat*.ssp".format(vvc)))
    elif vvc != 'all':
        # if told to, just get the vvc specified by vvc:
        sspfiles = glob.glob(os.path.join(sspoutdir, "*vvc{:.1f}*.ssp".format(vvc)))
        if len(sspfiles) == 0:
            sspfiles = glob.glob(os.path.join(sspoutdir, "*vvcrit{:.1f}*.ssp".format(vvc)))
    else:
        sspfiles = glob.glob(os.path.join(sspoutdir, "*.ssp"))

    print(sspfiles)

    # model resolution
    dT = 0.02
    dZ = 0.02

    if param == 'age':
        res = dT
    elif param == 'logZ':
        res = dZ

    # age solutions (could be adapted to any parameter):
    bestvals = np.array([]) #[] 
    # associated uncertainties:
    up_uncerts = np.array([])
    lo_uncerts = np.array([])
    # Look through all files found in the sspcombines/ directory:
    for j, f in enumerate(sspfiles):
        # Open the file, read lines, close file:
        currf = open(f)
        lines = currf.readlines()
        currf.close()
        print(f)
        # take the best fit from the 2nd line of the ssp file.
        try:
            bestline = lines[1]
            beststrs = bestline.split(',')
            # loop over param best-fits:
            for astr in beststrs:
                if '{:s}='.format(param) in astr:
                    bestval = float(astr.split('=')[-1])
            print(bestval)
        except IndexError:
            # empty file, skip it.
            continue
        # take the uncertainties from the direct marginalized distribution.        

        solns = {}
        p = {}
        inblock = False
        for l, line in enumerate(lines):

            if 'Direct marginalized distributions' in line:
                inblock = True
                m = 1
                # search the block for solutions:
                while inblock:
                    if ("{:s} = ".format(param) in lines[l+m]) | ("{:s} < ".format(param) in lines[l+m]) | ("{:s} > ".format(param) in lines[l+m]):
                        try:
                            # getting 16, 50, 84th percentiles:
                            if "{:s} = ".format(param) in lines[l+m]:
                                try:
                                    p[16], p[50], p[84] = map(float, ((lines[l+m].split('(')[-1]).split('from sum)')[0]).split(' - '))
                                except ValueError:
                                    lobound = float(((lines[l+m].split('(')[-1]).split(' from sum)')[0]).split('>')[-1])
                                #print(f)
                                print(p)
                            # get lower bound if unbounded (implement one for upper bound):
                            elif "{:s} > ".format(param) in lines[l+m]:
                                print(lines[l+m].split('(')[-1])
                                print((lines[l+m].split('(')[-1]).split('from sum)')[0])
                                print(((lines[l+m].split('(')[-1]).split('from sum)')[0]).split('>')[-1])
                                lobound = float(((lines[l+m].split('(')[-1]).split(' from sum)')[0]).split('>')[-1])
                            break
                        except Exception as emsg:
                            print(emsg)
                            print('Failed to get errors from direct marginalized distribution, will be set to 0.0 uncertainty.')
                            break
                    # here's where the end of the block is, so break:
                    elif (lines[l+m] == '\n') | (lines[l+m] == lines[-1]):
                        break
                    # keep searching until solution is found, or break out.
                    else:
                        m += 1

                #inblock = False
        
        # assign found values to array.
        try:
            bestvals = np.append(bestvals, lobound)
            up_uncerts = np.append(up_uncerts, 0.0)
            lo_uncerts = np.append(lo_uncerts, 0.0)
        except UnboundLocalError:
            bestvals = np.append(bestvals, bestval)
            try:
                up_uncerts = np.append(up_uncerts, abs(bestval - p[16]))
                lo_uncerts = np.append(lo_uncerts, abs(p[84] - bestval))
            except KeyError:
                # no uncertainties found
                up_uncerts = np.append(up_uncerts, 0.0)
                lo_uncerts = np.append(lo_uncerts, 0.0)

    # switch back to the work directory and return the found solutions and uncertainties, along with star number used
    # in finding the respective solutions:

    # MATCH reports @ lower bin edge in best fit solution for age.
    # also add in quadrature the unertainty due to parameter resolution. 
    # From e-mail between D. Weisz and A.E. Dolphin:
    #
    # "For reference, a uniform draw from 0 to 1 has a mean of 0.5 and 
    #  standard deviation of sqrt(1/12) = 0.29.  So, one-sigma resolution 
    # for any variable sampled over uniform spacing is the step size / sqrt(12)." - AED
    #if param == 'age':
    bestvals = bestvals + (res/2.0)

    if len(up_uncerts) > 1:
        up_uncerts = np.sqrt(up_uncerts**2 + (res/np.sqrt(12))**2)
    elif up_uncerts[0] != 0.0:
        up_uncerts = np.sqrt(up_uncerts**2 + (res/np.sqrt(12))**2)
    if len(lo_uncerts) > 1:
        lo_uncerts = np.sqrt(lo_uncerts**2 + (res/np.sqrt(12))**2)
    elif lo_uncerts[0] != 0.0:
        lo_uncerts = np.sqrt(lo_uncerts**2 + (res/np.sqrt(12))**2)

    print(np.mean(bestvals))
    print(np.sqrt(np.sum(up_uncerts**2))/len(up_uncerts))
    print(np.sqrt(np.sum(lo_uncerts**2))/len(lo_uncerts))

    return np.mean(bestvals), np.sqrt(np.sum(up_uncerts**2))/len(up_uncerts), np.sqrt(np.sum(lo_uncerts**2))/len(lo_uncerts)

# Count lines in a file:
def file_len(fname):
    print(fname)
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1
