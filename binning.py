from datetime import datetime
startTime = datetime.now()
import tracemalloc
    #from scippnexus import data
#import scippnexus as snx
import numpy as np
import math

    #qimport mantid_args
#from scippnexus import NXdetector
from scippneutron.conversion import graph

import scipp as sc
import scipp.binning
import scippneutron as scn
import scippnexus as snx

#%matplotlib widget


import matplotlib.pyplot as plt
import ipywidgets as ipw
from matplotlib.colors import LogNorm

    #import multiprocessing as mp

import argparse, sys
import h5py
# h5py.enable_ipython_completer()

# starting the momory monitoring
tracemalloc.start()


parser=argparse.ArgumentParser()
#parser.add_argument('--res', help='input res file',required=True)
parser.add_argument('--nex', help='input nexus file',required=True)
parser.add_argument('--bins', help='number of time binns (must be integer)',required=False,default=100)
parser.add_argument('--maxprop', help='nomralisation to maximum probility',required=False,default=1)

args=parser.parse_args()
if args.nex != None:
   filename=args.nex
   print("the used file is: %s" % filename)
if args.bins != None:
   t_step=int(args.bins)
   print("number of time bins:", t_step)
if args.maxprop != None:
   max_prop=int(args.maxprop)
   print("nomalisation of max propailitytie:", max_prop)


#Number op pixels per detecotr dimmension
pix = 1280

#number of detectors
n_det = 3


def CDist2(A,B):
#calculate distance betweenn two points
        dist = len3dvec(twoP_to_vec(A, B))
        return dist



def len3dvec(vec):
## calculates lengh of a 3D vecor
## input as list
        a = np.sqrt(vec[0]**2 + vec[1]**2 + vec[2]**2)
        return a

def twoP_to_vec(A,B):
#creates vector between two points
        vec = np.array([B[0]-A[0], B[1]-A[1], B[2]-A[2]])

        return vec

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])

def read_in_data(filename):

    f = h5py.File(filename)
    a = f['entry1/data']['bank01_events_dat_list_p_x_y_n_id_t']['events'][...]
    #a[0]
    d = np.matrix.transpose(a)
    print("shape of event list (p_x_y_n_id_t)", d.shape)
    print("start extracting data")
    #alocate units to events and create seperate list for each parameter
    t_list = sc.array(dims=['x'], unit='s', values=d[5])
    id_list = sc.array(dims=['x'], unit=None, values=d[4], dtype='int64')
    #print(x_list.shape, y_list.shape, t_list.shape,id_list.shape)

    weights = sc.array(dims=['x'], unit='counts', values=d[0]) #change to integer for measured data
    # normalise neutron with the highest probability to set value
    weights = weights * (max_prop/weights.max()) #delete for actual data
    # weights = sc.ones_like(x_list)
    # weights.unit = 'counts'
    da = sc.DataArray(data=weights, coords={ 't': t_list, 'id': id_list})
    

    #make sure alll IDs are reconised:
    print("id min",id_list.values.min())
    print("id max",id_list.values.max())

    ids1 = sc.arange('id', 1, 1638401, unit=None)
    ids2 = sc.arange('id', 2000001, 3638401, unit=None)
    ids3 = sc.arange('id', 4000001, 5638401, unit=None)
    ids = sc.concat([ids1, ids2, ids3], 'id')
    #grouping by IDs
    return da, ids

print("read in data")

da, ids = read_in_data(filename)


print("beginn grouping")
grouped = da.group(ids)
del da

print("beginn hinstogram")
group_t = grouped.hist(t=t_step)
del grouped
#group_t = da.hist(t=t_step)
#t_edges = sc.linspace('t', t_list.min().value, np.nextafter(t_list.max().value,np.inf), t_step, unit=t_list.unit)
#group_t = sc.binning.make_binned(da, edges=[t_edges], groups=[ids]).hist()
print("endn grouping")
print("extractin geometry informations")
f = h5py.File(filename)
origen = f['entry1/data']['bank01_events_dat_list_p_x_y_n_id_t']['distance'][0].decode()
origen = list(np.float_(origen.split()))

xml  = str(f['entry1/instrument/instrument_xml/data'][...][0]).split('\\n')

comp = False
det = False
source = False
sample = False
sample_pos = [0,0,0]
source_pos = [0,0,0]
d_list = []
rot_l1 = []
fast_l = []
slow_l = []


for i in range(len(xml)):
    ls = xml[i].replace('<t',' ').replace('>',' ').replace('"',' ').replace('<',' ').replace('\\t',' ').split()
    # print(xml[i])
    
    if len(ls) >= 1:
        if ls[0] == 'component':
            det = False
            source = False
            sample = False
            comp = True
            if ls[2].split('-')[0] == 'MonNDtype':

                det = True
                d_list.append([int(ls[2].split('-')[1])])
            elif ls[2] == 'sourceMantid-type':
                source = True
            elif ls[2] == 'sampleMantid-type':
                sample = True
            comp = True
     #   if ls[1].split('-')[0] == 'type="MonNDtype':
            # print("1",ls)
    if len(ls) >= 1:
        if ls[0] == 'type':
            comp = False
    if len(ls) >= 1:
        if ls[0] == 'location' or ls[0] == 'location':
            # print("3",ls)
            if comp == True and det == True:
                # print("2",ls)
                xyz = [float(ls[2]),float(ls[4]),float(ls[6])]
                print("xyz", xyz)
                rot =float(ls[8]) 
                rot_xyz =[ float(ls[8]), float(ls[10]),float(ls[12]),float(ls[14])] 
                print("rotation of detector",rot, rot_xyz) 
                d_list[len(d_list)-1].append(xyz)
                rot_l1.append(rot_xyz)
            elif comp == True and source == True: 
              source_pos = [float(ls[2]),float(ls[4]),float(ls[6])]
            elif comp == True and sample == True: 
               sample_pos = np.array([float(ls[2]),float(ls[4]),float(ls[6])])
# print(len(d_list))
# print("sample and source position",sample_pos,source_pos)
# print("distance between sample and source",CDist2(source_pos, sample_pos))
# print("detector positions, relative to source at 0,0,0:",d_list)
# print("rotation list",rot_l)
#shift from rleative position to source to relative to sample
ds_l = []
sample_pos = sample_pos * [-1,-1,-1]
print("sample_pos",sample_pos)
for i in range(len(d_list)):
    det_pos = np.array([d_list[i][1][0],d_list[i][1][1],d_list[i][1][2]]) * [-1.0, -1.0,-1.0]
    

    rel_xyz= np.round(twoP_to_vec(sample_pos,det_pos),2)
    print("rel_xyz",rel_xyz)
    rel_xyz = rel_xyz # * [-1,-1,-1]
    print("sample to detector dist",len3dvec(rel_xyz))
    # rel_pos= twoP_to_vec(det_pos, sample_pos)
    # print("detector position",sample_pos,det_pos,rel_pos)
    
    # rel_xyz = [ d_list[i][1][0]-sample_pos[0], d_list[i][1][1]-sample_pos[1], (d_list[i][1][2]-sample_pos[2])]
    
    # print("relative position",rel_xyz)
    # rel_xyz = [ sample_pos[0]-d_list[i][1][0],sample_pos[1]-d_list[i][1][1], sample_pos[2]-d_list[i][1][2]]
    # print("position",rel_xyz,rel_pos)
    # print("relative position",rel_xyz)
    ds_l.append(rel_xyz)
vec_f = [-1,0,-0]
vec_s = [0,-1,0]
rot_l = [rot_l1[0],rot_l1[2],rot_l1[1]] #not corect but a temporary work around
for i in range(len(rot_l)):
    # if rot_l[i][2] < 0:
        # theta = np.radians(rot_l[i][0]) 
    # else:
    theta = np.radians(-rot_l[i][0])
    v = [rot_l[i][1],rot_l[i][2],rot_l[i][3]]
    v1 = [0,0,1]
    print("v",i,v,"theta", theta,rot_l[i][0])
    fast_vec = np.round(np.dot(rotation_matrix(v,  theta), vec_f),2)
    fast_v_round = np.array([np.round(fast_vec[0],1), np.round(fast_vec[1],1), np.round(fast_vec[2],1)])

    fast_l.append(fast_vec)
    slow_l.append(np.dot(rotation_matrix(v, theta), vec_s))
    # slow_l.append([0,1,0])

print("relative to sample position",ds_l)
print("fast axis",fast_l)
print("slow axis", slow_l)
print(rot_l)

print((f['entry1']['simulation']['Param'].keys()))
print((f['entry1']['simulation']['Param']['XtalPhiX']))
phix=float(list(str(f['entry1']['simulation']['Param']['XtalPhiX'][...][0]))[2])
phiy=float(list(str(f['entry1']['simulation']['Param']['XtalPhiY'][...][0]))[2])
phiz=float(list(str(f['entry1']['simulation']['Param']['XtalPhiZ'][...][0]))[2])
#str(phix[0])
#int(list(str(phix[0]))[2])
print(phix,phiy,phiz)
cor=[phix,phiy, phiz]
cryst_or = np.array(cor)
cryst_or

no = filename.split('/')
print(no)
name_out= no[-1].split('.')[0]
print(name_out)
file_out = './'+name_out+'_'+str(t_step)+'_out.h5'
print(file_out)

print("start writing to file")

with h5py.File(file_out, 'w') as fo:
## create output nexus file
   fo.attrs['default'] = 'NMX_data'
   nxentry = fo.create_group('NMX_data')
   nxentry.attrs["NX_class"] = 'NXentry'
   nxentry.attrs['default'] = 'data'
   nxentry.attrs['name'] = 'NMX1'
   #nxentry.__setitem__('beamline','NMX')
   nxentry.__setitem__('name','NMX')
   nxentry.__setitem__('definition','TOFRAW')
   nxentry.attrs['name'] = "NMX"
   #nxentry.attrs['name'].__setattr__('name','NMX') 

#SAMPLE
   nx_sample = nxentry.create_group('NXsample')
   nx_sample.__setitem__('name','Single_crystal')


#Instrument
   nx_instrument = nxentry.create_group('NXinstrument')

   nx_detector = nxentry.create_group('NXdetector')
   det_origen = nx_detector.create_dataset('origen', data=ds_l) 
   det_origen.attrs['units'] = 'm'

   fast_axis = nx_detector.create_dataset('fast_axis', data=fast_l) 
   slow_axis = nx_detector.create_dataset('slow_axis', data=slow_l) 

   nx_beam = nxentry.create_group('NXbeam')


   
   proton = nxentry.create_dataset('proton_charge', data=1)    
   
   
   nx_det1 = nxentry.create_group('detector_1')   
   counts = nx_det1.create_dataset('counts', data=[group_t.values], compression="gzip", compression_opts=4)

   t_spectrum = nx_det1.create_dataset('t_bin', data=group_t.coords['t'].values, compression="gzip", compression_opts=4)
   t_spectrum.attrs['units'] = 'ms'
   t_spectrum.attrs['long_name'] = 't_bin TOF (ms)'

   pixel_id = nx_det1.create_dataset('pix_id', data=group_t.coords['id'].values, compression="gzip", compression_opts=4)
   pixel_id.attrs['units'] = ''
   pixel_id.attrs['long_name'] = 'pixel ID'


#SOURCE
   nx_inst = nxentry.create_group('instrument')
   nx_inst.attrs['nr_detector'] = len(d_list)
   nx_source = nxentry.create_group('NXsource')
   nx_source.__setitem__('name','European Spallation Source')
   nx_source.__setitem__('short_name','ESS')
   nx_source.__setitem__('type','Spallation Neutron Source')
   nx_source.__setitem__('distance',-CDist2(source_pos, sample_pos))
   nx_source.__setitem__('probe','neutron')
   nx_source.__setitem__('target_material','W')






#    c_or = nxinst.create_dataset('crystal_orientation',data=cryst_or)
#    c_or.attrs['unit'] = 'degrees'
#    c_or.attrs['long_name'] = 'crystal orientation in Phi (degrees)'
    

   fo.close()
   del fo

current, peak = tracemalloc.get_traced_memory()
print('2 current memory [GB]: {}, peak memory [GB]: {} '.format(round((current/(1024*1024*1024)), 4), round((peak /(1024*1024*1024) ), 4) ))
# stopping the library
tracemalloc.stop()
print("neded time (h:mm:ss): ", datetime.now() - startTime)
