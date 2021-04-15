#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# ABOUT
# Visualize HDF5 state files calculated by FLEXI. Only particle data for now
#
# USAGE
# python3 visualize [-f path] [-o output] [-H height] [-T simtime]

# Libaries
import os, sys
import argparse
import h5py
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt
import matplotlib
from tikzplotlib import save as tikz_save
from progressbar import *                   # just a simple progress bar

# use LaTeX fonts in the plot
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

widgets = ['Processing: ', Percentage(), ' ', Bar(marker='%',left='[',right=']'),
           ' ', ETA()]                      #see docs for other options

# read user arguments
def read_args(args):
    parser = argparse.ArgumentParser(description='Visualize HDF5 state files calculated by FLEXI')
    parser.add_argument('-f','--filepath', type=str, nargs='?', help = 'Location to be visualized')
    parser.add_argument('-m','--mode',     type=str, nargs='?', help = 'Output mode')
    parser.add_argument('-o','--output',   type=str, nargs='?', help = 'Output folder')
    parser.add_argument('-H','--height',   type=str, nargs='?', help = 'Cube half height')
    parser.add_argument('-T','--simtime',  type=str, nargs='?', help = 'Simulation time')

    try:
        args = parser.parse_known_args(args)[0]
    except:
        parser.print_help()
        exit()

    return args

# main routine
class visualize(object):
    def __init__(self, args):
        # get work folder and output folder
        self.case_dir   = args.filepath
#        self.sgs_name   = os.path.basename(self.case_dir)
        # Path not empty after trailing edge
        if (os.path.split(self.case_dir)[1]):
            self.sgs_name   = os.path.split(self.case_dir)[1]
        else:
            # Try before the trailing slash
            self.sgs_name   = os.path.split(os.path.dirname(self.case_dir))[1]

        # output folder if existing
        if args.output:
            self.output_dir = '{}'.format(args.output)
        else:
            self.output_dir = '{}/output'.format(self.case_dir)

        # channel height if existing
        try:
            if float(args.height) <= 0:
                print('Cube height less or equal zero. Exiting ...')
                exit()
            else:
                self.height = float(args.height)
        except:
            print('Missing cube height. Exiting ..')
            exit()

        self.simtime      = float(args.simtime)

        # if the output folder does not exist, create it
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir, exist_ok=True)

    def case(self):
        # get state and mesh files
        state_files = [kk for kk in os.listdir(self.case_dir) if 'h5' in kk and 'State_0' in kk and not 'ERROR' in kk and not 'TimeAvg' in kk and not 'restart' in kk]
        mesh_file   = [jj for jj in os.listdir(self.case_dir) if 'h5' in jj and 'mesh' in jj]

        # sort files by name
        if len(state_files) <= 0:
            print('Zero state files. Exiting ...')
            exit()

        state_files = sorted(state_files)

        # extract particle data
        with h5py.File(os.path.join(self.case_dir, state_files[0]), 'r') as h5:
            data = {key: val[:] for key, val in h5.items()}
            nSpecies = int(max(data['PartData'][:,6]))
            print('Found FLEXI state file with {} species'.format(nSpecies))
            print('Found datasets for {}'.format(list(data.keys())))

        # Just a progress bar
        pbar = ProgressBar(widgets=widgets, maxval=len(state_files))
        pbar.start()

        # iterate over all state files
        for i,j in enumerate(state_files):
            with h5py.File(os.path.join(self.case_dir, j), 'r') as h5:
#                print('Processing file {} of {}, t={}'.format(i+1,len(state_files),h5.attrs['Time']))
                pbar.update(i)

                data = {key: val[:] for key, val in h5.items()}

                partInt = data['PartData'].shape[0]
                partVar = data['PartData'].shape[1]
                nVarPathStart = partVar - 3

                # run in requested mode
                if args.mode=='time':

                    # just sum it up
                    if i == 0:
                        count  = np.zeros((  nSpecies,len(state_files)))
                        result = np.zeros((8,nSpecies,len(state_files)))
                        deviat = np.zeros((8,nSpecies,len(state_files)))

                    for k in range(partInt):
                        # Overall particle count for concentration
                        # remember that Python indices start at 0
                        species = int(data['PartData'][k][6]-1)
                        count[species,i] += 1

                        partPathX   = data['PartData'][k][nVarPathStart + 0]
                        partPathY   = data['PartData'][k][nVarPathStart + 1]
                        partPathZ   = data['PartData'][k][nVarPathStart + 2]

                        partVeloX   = data['PartData'][k][3]
                        partVeloY   = data['PartData'][k][4]
                        partVeloZ   = data['PartData'][k][5]

                        partPathAbs = (partPathX**2. + partPathY**2. + partPathZ**2.)**0.5
                        delta       = partPathAbs - result[0,species,i]
                        result[0,species,i] += delta/count[species,i]
                        delta2      = partPathAbs - result[0,species,i]
                        deviat[0,species,i] += delta2*delta

                        partPathAbs = (partPathX**2.)**0.5
                        delta       = partPathAbs - result[1,species,i]
                        result[1,species,i] += delta/count[species,i]
                        delta2      = partPathAbs - result[1,species,i]
                        deviat[1,species,i] += delta2*delta

                        partPathAbs = (partPathY**2.)**0.5
                        delta       = partPathAbs - result[2,species,i]
                        result[2,species,i] += delta/count[species,i]
                        delta2      = partPathAbs - result[2,species,i]
                        deviat[2,species,i] += delta2*delta

                        partPathAbs = (partPathZ**2.)**0.5
                        delta       = partPathAbs - result[3,species,i]
                        result[3,species,i] += delta/count[species,i]
                        delta2      = partPathAbs - result[3,species,i]
                        deviat[3,species,i] += delta2*delta

                        VeloMagnAbs = (partVeloX**2. + partVeloY**2. + partVeloZ**2.)**0.5
                        delta       = VeloMagnAbs - result[4,species,i]
                        result[4,species,i] += delta/count[species,i]
                        delta2      = VeloMagnAbs - result[4,species,i]
                        deviat[4,species,i] += delta2*delta

                        VeloMagnAbs = (partVeloX**2.)**0.5
                        delta       = VeloMagnAbs - result[5,species,i]
                        result[5,species,i] += delta/count[species,i]
                        delta2      = VeloMagnAbs - result[5,species,i]
                        deviat[5,species,i] += delta2*delta

                        VeloMagnAbs = (partVeloY**2.)**0.5
                        delta       = VeloMagnAbs - result[6,species,i]
                        result[6,species,i] += delta/count[species,i]
                        delta2      = VeloMagnAbs - result[6,species,i]
                        deviat[6,species,i] += delta2*delta

                        VeloMagnAbs = (partVeloZ**2.)**0.5
                        delta       = VeloMagnAbs - result[7,species,i]
                        result[7,species,i] += delta/count[species,i]
                        delta2      = VeloMagnAbs - result[7,species,i]
                        deviat[7,species,i] += delta2*delta


                    # plot the whole thing
                    if i == len(state_files)-1:

                        # set figure size and label
                        plt.figure(figsize=(8,6),dpi=300)
#                        plt.legend(fontsize=8)
                        plt.xlabel('t [sec]')
                        plt.ylabel('Path length (mean $\pm\sigma^2)$', usetex=True)

                        # Done with the progress bar
                        pbar.finish()

                        print('Generating output ...')
                        x_range = np.linspace(0,self.simtime,len(state_files))
                        for m in range(0,nSpecies):
                            # Remember that delta is the channel HALF height
                            plt.margins(0.)

                            # Absolute magnitude
                            if nSpecies > 1:
                                plt.plot(x_range,result[0][m][:])
                                variance = deviat[0][m][:]/count[m][:]
                                plt.fill_between(x_range,result[0][m][:]-variance,result[0][m][:]+variance,alpha=0.2)

                                plt.gca().legend(['$St^+=10$','$St^+=5$','$St^+=1$','$St^+=0.5$','$St^+=0.1$','$St^+=0.01$'],fontsize=8,loc=2)
                            # Magnitude in Cartesian direction
                            else:
                                plt.plot(x_range,result[1][m][:])
                                variance = deviat[1][m][:]/count[m][:]
                                plt.fill_between(x_range,result[1][m][:]-variance,result[1][m][:]+variance,alpha=0.2)

                                plt.plot(x_range,result[2][m][:])
                                variance = deviat[2][m][:]/count[m][:]
                                plt.fill_between(x_range,result[2][m][:]-variance,result[2][m][:]+variance,alpha=0.2)

                                plt.plot(x_range,result[3][m][:])
                                variance = deviat[3][m][:]/count[m][:]
                                plt.fill_between(x_range,result[3][m][:]-variance,result[3][m][:]+variance,alpha=0.2)

                                plt.gca().legend(['Tracer (x)','Tracer (y)','Tracer (z)'],fontsize=8,loc=2)
                                plt.ylim(0,5)

                            plt.xlim(0,self.simtime)
#                            plt.ylim(0,20)

#                        plt.savefig('Particle_Concentration_unscaled.pdf')
                        pic_namedir = os.path.join(self.output_dir, 'DispTime_SGS_{}'.format(self.sgs_name))
                        print('Writing output to {}'.format(pic_namedir))
#                        plt.savefig('{}.pdf'.format(pic_namedir))
                        tikz_save('{}.tikz'.format(pic_namedir))
                        plt.close()

                        # set figure size and label
                        plt.figure(figsize=(8,6),dpi=300)
#                        plt.legend(fontsize=8)
                        plt.xlabel('t [sec]')
                        plt.ylabel('u^\prime(t)$', usetex=True)

                        for m in range(0,nSpecies):
                            # Remember that delta is the channel HALF height
                            plt.margins(0.)

                            # Absolute magnitude
                            if nSpecies > 1:
                                plt.gca().legend(['$St^+=10$','$St^+=5$','$St^+=1$','$St^+=0.5$','$St^+=0.1$','$St^+=0.01$'],fontsize=8,loc=2)
                                plt.plot(x_range,result[4][m][:])
                                variance = deviat[4][m][:]/count[m][:]
                                plt.fill_between(x_range,result[4][m][:]-variance,result[4][m][:]+variance,alpha=0.2)

                            # Magnitude in Cartesian direction
                            else:
                                plt.plot(x_range,result[5][m][:])
                                variance = deviat[5][m][:]/count[m][:]
                                plt.fill_between(x_range,result[5][m][:]-variance,result[5][m][:]+variance,alpha=0.2)

                                plt.plot(x_range,result[6][m][:])
                                variance = deviat[6][m][:]/count[m][:]
                                plt.fill_between(x_range,result[6][m][:]-variance,result[6][m][:]+variance,alpha=0.2)

                                plt.plot(x_range,result[7][m][:])
                                variance = deviat[7][m][:]/count[m][:]
                                plt.fill_between(x_range,result[7][m][:]-variance,result[7][m][:]+variance,alpha=0.2)

                                plt.gca().legend(['Tracer'],fontsize=8,loc=2)
                                plt.ylim(0,2.0)

                            plt.xlim(0,self.simtime)
#                            plt.ylim(0,20)

#                        plt.savefig('Particle_Concentration_unscaled.pdf')
                        pic_namedir = os.path.join(self.output_dir, 'DispVelo_SGS_{}'.format(self.sgs_name))
                        print('Writing output to {}'.format(pic_namedir))
#                        plt.savefig('{}.pdf'.format(pic_namedir))
                        tikz_save('{}.tikz'.format(pic_namedir))
                        plt.close()

                if args.mode=='concentration':
                    if i == 0:
                        # Allocate array to hold particles in first state file (only works for initial emission)
                        nVar = 8

                        # Find number of particles in each species
                        nPartSpecies = np.zeros(nSpecies,dtype=int)

                        for k in range(partInt):
                            nPartSpecies[int(data['PartData'][k][6]-1)] = nPartSpecies[int(data['PartData'][k][6]-1)] + 1

                        for m in range(nSpecies):
                            print('Found {} particles for species {}'.format(nPartSpecies[m],m+1))

                        # Allocate array to hold particles for each species
                        result  = np.zeros((nSpecies,np.amax(nPartSpecies),nVar))
                        counter = np.zeros(nSpecies,dtype=int)

                    # Operate only on requested sim time
                    time = h5.attrs['Time']
                    # Arbitrary barrier at dt = 10e-8. Might change this later
                    if abs(float(time) - float(self.simtime)) < 10e-8:
                        # Done with the progress bar
                        pbar.finish()

                        print('Found state file for t={} sec. Visualizing ...'.format(self.simtime), end = '')

                        # set figure size and label
                        plt.figure(figsize=(8,6),dpi=300)
                        plt.xlabel('$x/\delta$', usetex=True)
                        plt.ylabel('$<C>/C_0>$', usetex=True)

                        # Sort particles into species
                        counter[:] = 0
                        for k in range(partInt):
                            PartSpecies = int(data['PartData'][k][6]-1)
                            result[PartSpecies,counter[PartSpecies],:] = data['PartData'][k][:]
                            counter[PartSpecies] = counter[PartSpecies] + 1

                        # plot the whole thing
                        # scale y to uniform axis
                        result[:,:,1] = result[:,:,1]/self.height
                        result[:,:,1] = result[:,:,1] + 1

                        # repeat towards both side for kde. coordinates are already scaled with delta. flip the coordinates for mirroring, range [0,2] --> [-2,4]
                        result_kde = np.zeros((nSpecies,3*np.amax(counter)))
                        for k in range(nSpecies):
                            result_kde[k,  0         :  counter[k]-1] = -result[k,0:counter[k]-1,1]
                            result_kde[k,  counter[k]:2*counter[k]-1] =  result[k,0:counter[k]-1,1]
                            result_kde[k,2*counter[k]:3*counter[k]-1] = -result[k,0:counter[k]-1,1] + 4.

                        # plot (normalized) histogram of the data
                        for k in range(nSpecies):
                            plt.hist(result[k,:,1]-1, weights=np.ones_like(result[k,:,1])/np.sum(result[k,:,1])*120,histtype='stepfilled',bins = int(120),label='FLEXI (Hist)',alpha=0.1,color=plt.cm.RdYlGn(k/nSpecies))
                            plt.hist(result[k,:,1]-1, weights=np.ones_like(result[k,:,1])/np.sum(result[k,:,1])*120,histtype='step'      ,bins = int(120),label='FLEXI (Hist)',alpha=0.3,color=plt.cm.RdYlGn(k/nSpecies))

                            # test values for the bw_method option ('None' is the default value)
#                            bw_values =  [3./2.*self.height]
                            bw_values =  [self.height]

                            # generate a list of kde estimators for each bw
                            kde = [scipy.stats.gaussian_kde(result_kde[k,:],bw_method=bw) for bw in bw_values]

                            # plot density estimates
                            plt.plot(np.linspace(-1.,1.,120),6*kde[0](np.linspace(0.,2.,120)),lw=1,color=plt.cm.RdYlGn(k/nSpecies))

                        plt.gca().legend(['$St^+=10$','$St^+=5$','$St^+=1$','$St^+=0.5$','$St^+=0.1$','$St^+=0.01$'],fontsize=8,loc=9)
                        plt.xlim(-1.,1.,0)
                        plt.ylim( 0.,4.,0)
#                        plt.savefig('Particle_Concentration.pdf')
                        pic_namedir = os.path.join(self.output_dir, 'ConvHeight_{}_t_{}'.format(self.pict_dir,time[0]))
                        tikz_save('{}.tikz'.format(pic_namedir))
                        print(' done')
                        sys.exit()

# Run everything
if __name__ == '__main__':
    args = read_args(sys.argv[1:])
    visualize(args).case()
