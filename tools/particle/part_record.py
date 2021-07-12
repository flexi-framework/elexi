import numpy as np
import h5py as h5
import argparse,sys,math,os,csv

########################################################
# Evaluate particle record planes
# Author: Anna Schwarz
# Date: 31.03.2021
#
# Record plane writes
# PartState(1:6), Species, Diameter, Force, Sphericity
########################################################


def user_input(args):
  parser = argparse.ArgumentParser(description='Plot rebound data')
  optional = parser._action_groups.pop()
  required = parser.add_argument_group('required arguments')
#  optional.add_argument('-c','--write_to_csv',        help="Write .h5 data to .csv")
  optional.add_argument('-b','--write_to_binary', default='', help="Write .h5 data to .dat")
  optional.add_argument('-c','--write_to_csv',    default='', help="Write .h5 data to .csv")
  optional.add_argument('-d','--dir',             default='', help="Directory of record points")
  required.add_argument('plotfiles',              default='', nargs='+', help="Record point files")
  parser._action_groups.append(optional)

  try:
    args = parser.parse_known_args(args)[0]
  except:
    parser.print_help()
    exit()

  return args

def write_csv(filename,writedata,VarNames):
  """Write data to .csv"""
  with open(filename,'w', newline='') as c:
    writer=csv.writer(c, delimiter=',')
    writer.writerow(VarNames)
    for i in writedata:
      writer.writerow(i)

def write_binary(filename,writedata,nSpecies,VarNames,SphericityIC,DiameterIC,ForceIC):
  """Write data to binary, .dat"""
  with open(filename,'wb') as f:
    np.array(VarNames.shape[0]).astype('int32').tofile(f,sep="")
    VarNames.tofile(f,sep="")
    nSpecies.tofile(f,sep="")
    SphericityIC.tofile(f,sep="")
    DiameterIC.tofile(f,sep="")
    np.array(ForceIC).astype('U20').tofile(f,sep="")
    np.array(np.array(writedata).shape).astype('int32').tofile(f,sep="")
    np.array(writedata).astype('float64').tofile(f,sep="")

def readin_hdf5(args):
  """ Reading hdf5 data"""
  print('Reading data from h5 files.')
  # Read attributes: VarNames and nSpecies
  # VarNamesPart: PartState(1:6), Species, Diameter, Force, Sphericity
#  plotfiles  = [kk for kk in os.listdir(args.dir) if 'h5' in kk and 'recordpoints_part' in kk]
#  with h5.File(os.path.join(args.dir,plotfiles[0]),'r') as f:
  with h5.File(args.plotfiles[0],'r') as f:
    # Decode bytes
    VarNamesPart = f.attrs['VarNamesPart'].astype('U13')
    nSpecies     = f.attrs['nSpecies'][0].astype('int32')
    if args.write_to_binary != '':
      SphericityIC = f.attrs['SphericityIC'].astype('float32')
      DiameterIC   = f.attrs['DiameterIC'].astype('float32')
      # Read in the forces
      DragForce    = f.attrs['DragForce'].astype('int32')
      LiftForce    = f.attrs['LiftForce'].astype('int32')
      VirtForce    = f.attrs['VirtForce'].astype('int32')
      UndiForce    = f.attrs['UndiForce'].astype('int32')
      BassForce    = f.attrs['BassForce'].astype('int32')
      ForceIC = ['F_D' for i in range(nSpecies)]
      for i in range(nSpecies):
        if LiftForce[i] == 1:
          ForceIC[i] = ForceIC[i] + '+F_L'
        if VirtForce[i] == 1:
          ForceIC[i] = ForceIC[i] + '+F_V'
        if UndiForce[i] == 1:
          ForceIC[i] = ForceIC[i] + '+F_U'
        if BassForce[i] == 1:
          ForceIC[i] = ForceIC[i] + '+F_B'

  recorddata = {i: {k: [] for k in range(len(VarNamesPart))} for i in range(nSpecies)}
  writedata  = []
  for i in args.plotfiles:
    with h5.File(os.path.join(args.dir,i),'r') as f:
      for k,m in enumerate(f['RecordData'][:,6]):
        if m > 0:
          writedata.append(f['RecordData'][k,:])
          for l in range(len(VarNamesPart)):
            recorddata[m-1][l].append(f['RecordData'][k,l])

  # Write data to csv
  if args.write_to_binary != '':
    write_binary(args.write_to_binary, writedata, nSpecies, VarNamesPart, SphericityIC, DiameterIC, ForceIC)

  if args.write_to_csv != '':
    write_csv(args.write_to_csv, writedata, VarNamesPart)

readin_hdf5(user_input(sys.argv[1:]))

