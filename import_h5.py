#!/usr/bin/env python3

# Forked from import_h5_4.py at 05/25/2021 by SJ
# To be publicly released
#
# Difference (wrt. import_h5_4.py):
# - code cleanup
# - hard codes removed
#
# Initial code in 09/25/2019 by Seongsu Jeong
#

import h5py
import glob
import numpy as np
import datetime
#from numpy.lib.nanfunctions import _nanquantile_dispatcher
#from pyproj import Proj, transform
from pyproj import CRS, Transformer
import os
from osgeo import ogr, osr
from scipy import interpolate
import grup
import multiprocessing
import sys
import time


def subsample_xy(vec_x,vec_y,scale_inv=10):
    try:
        nvx=vec_x.size
        nvy=vec_y.size
    except:
        nvx=len(vec_x)
        nvy=len(vec_y)


    if nvx==nvy:
        xsub=vec_x[0::scale_inv]
        ysub=vec_y[0::scale_inv]
        if xsub[-1]!=vec_x[-1]:
            xsub[-1]=vec_x[-1]
            ysub[-1]=vec_y[-1]
        return xsub,ysub
        
    else:
        print('ERROR: Size of x and y vectors do not match')
        return None,None


def xy_2_wktline(vec_x,vec_y):
    #LINESTRING (30 10, 10 30, 40 40)
    form_wkt_linestring='LINESTRING ({})'

    try:
        nvx=vec_x.size
        nvy=vec_y.size
    except:
        nvx=len(vec_x)
        nvy=len(vec_y)

    if nvx==nvy:
        num_node=nvx
    else:
        print('ERROR: Size of x and y vectors do not match')
        num_node=None
        return None
    
    if num_node!=None:
        list_str_xy=[None]*num_node

        for i in range(num_node):
            list_str_xy[i]='{} {}'.format(vec_x[i],vec_y[i])
        
        str_coords=','.join(list_str_xy)

        str_wkt=form_wkt_linestring.format(str_coords)
        return str_wkt


def merge_subtatalog(list_filename_subcat,filename_shp_out):
    #Prepare to load the input subcats
    drvin=ogr.GetDriverByName("ESRI Shapefile")
    datasrc_in=drvin.Open(list_filename_subcat[0],0)
    lyr_in=datasrc_in.GetLayer()
    feat0=lyr_in.GetNextFeature()
    geom0=feat0.GetGeometryRef()
    srs_cat=geom0.GetSpatialReference()
    datasrc_in=None
    
    #Prepare for the out shapefile
    drv=ogr.GetDriverByName("ESRI Shapefile")
    datasrc=drv.CreateDataSource(filename_shp_out)
    lyr=datasrc.CreateLayer('footprint',srs_cat,ogr.wkbLineString)
    
    field_source=ogr.FieldDefn('atl06',ogr.OFTString)
    field_source.SetWidth(64)
    lyr.CreateField(field_source)
    field_source=ogr.FieldDefn('source',ogr.OFTString)
    field_source.SetWidth(128)
    lyr.CreateField(field_source)
    lyr.CreateField(ogr.FieldDefn("track", ogr.OFTInteger))
    field_gt=ogr.FieldDefn('gt',ogr.OFTString)
    field_gt.SetWidth(4)
    lyr.CreateField(field_gt)
    lyr.CreateField(ogr.FieldDefn("cycle", ogr.OFTInteger))
    lyr.CreateField(ogr.FieldDefn("region", ogr.OFTInteger))
    
    #iterate through the sub catalogs
    for i,filename_subcat in enumerate(list_filename_subcat):
        print('Processing {} of {} subcat: {}'.format(i+1,len(list_filename_subcat),os.path.basename(filename_subcat)))
        datasrc_in=drvin.Open(filename_subcat,0)
        lyr_in=datasrc_in.GetLayer()
        for feat_in in lyr_in:
            feat_out=ogr.Feature(lyr.GetLayerDefn())
            feat_out.SetField('atl06',os.path.basename(feat_in.GetField('atl06')))
            feat_out.SetField('source',feat_in.GetField('source'))
            feat_out.SetField('track',feat_in.GetField('track'))
            feat_out.SetField('gt',feat_in.GetField('gt'))
            feat_out.SetField('cycle',feat_in.GetField('cycle'))
            feat_out.SetField('region',feat_in.GetField('region'))
            feat_out.SetGeometry(feat_in.GetGeometryRef())
            lyr.CreateFeature(feat_out)
            feat_out=None
    datasrc=None
    print('Subcat merging completed: {}'.format(filename_shp_out))
    


def hdf5_to_npy(list_hdf5,path_npy,epsg_out=3031,part=0,npart=1):
    scale_gtrk_inv=100
    
    #Deal with projection
    #inproj=Proj(init='epsg:4326')
    #outproj=Proj(init='epsg:{}'.format(epsg_out))
    crs_in=CRS.from_epsg(4326)
    crs_out=CRS.from_epsg(epsg_out)

    proj = Transformer.from_crs(crs_in, crs_out)

    srs=osr.SpatialReference()
    srs.ImportFromEPSG(epsg_out)

    #prepare to write shape file for catalog
    filename_catalog='{}/catalog_{}of{}.shp'.format(path_npy,part,npart)
    drv=ogr.GetDriverByName("ESRI Shapefile")

    datasrc=drv.CreateDataSource(filename_catalog)
    
    lyr=datasrc.CreateLayer('footprint',srs,ogr.wkbLineString)

    field_source=ogr.FieldDefn('atl06',ogr.OFTString)
    field_source.SetWidth(64)
    lyr.CreateField(field_source)
    field_source=ogr.FieldDefn('source',ogr.OFTString)
    field_source.SetWidth(128)
    lyr.CreateField(field_source)
    lyr.CreateField(ogr.FieldDefn("track", ogr.OFTInteger))
    field_gt=ogr.FieldDefn('gt',ogr.OFTString)
    field_gt.SetWidth(4)
    lyr.CreateField(field_gt)
    lyr.CreateField(ogr.FieldDefn("cycle", ogr.OFTInteger))
    lyr.CreateField(ogr.FieldDefn("region", ogr.OFTInteger))


    form_filename_npy='{PARENT}/Track{TRACK}/C{CYCLE}R{REGION}/{GT}.npy'
    list_str_gt=['gt1l','gt1r','gt2l','gt2r','gt3l','gt3r']
    for i,filename_h5 in enumerate(list_hdf5):
        if i%npart==part:
            print(i+1,'/',len(list_hdf5),'-',filename_h5)
            h5name=filename_h5.split('/')[-1]
            seg_h5name=h5name.split('_')
            str_reftrack=seg_h5name[2][0:4]
            str_cycle=seg_h5name[2][4:6]
            str_region=seg_h5name[2][6:8]
            
            with h5py.File(filename_h5,'r') as hin: 
                t0_gps=np.array(hin['ancillary_data']['atlas_sdp_gps_epoch'])[0]
                for str_gt in list_str_gt:
                    if str_gt in hin.keys():
                        filename_npy=form_filename_npy.format(PARENT=path_npy,
                                                            TRACK=str_reftrack,
                                                            GT=str_gt.upper(),
                                                            CYCLE=str_cycle,
                                                            REGION=str_region)
                        #print(filename_npy)
                        if not os.path.exists(filename_npy):
                            gtoi=hin[str_gt]
                            lis=gtoi['land_ice_segments']
                            lat=np.array(lis['latitude'])
                            lon=np.array(lis['longitude'])
                            hgt=np.array(lis['h_li'])
                            hgt_sigma=np.array(lis['h_li_sigma'])
                            qs=np.array(lis['atl06_quality_summary'])
                            t=np.array(lis['delta_time'])
                            geoid=np.array(lis['dem']['geoid_h'])
                            numobs=len(lon)

                            #get rid of values that does not not sense (i.e. too low or too high elev.)
                            flag_valid_obs=(hgt<10000.0) & (hgt>-1000.0)
                            
                            #xmap,ymap=transform(inproj,outproj,lon,lat)
                            #xmap,ymap=proj.transform(lon[flag_valid_obs],lat[flag_valid_obs])
                            xmap,ymap=proj.transform(lat[flag_valid_obs],lon[flag_valid_obs]) #NOTE: New version of proj takes WGS84 coord. as [lat, lon] rather than [lon lat]
                            
                            list_x=xmap
                            list_y=ymap
                            list_z=hgt[flag_valid_obs]
                            list_z_sigma=hgt_sigma[flag_valid_obs]
                            list_t=t[flag_valid_obs]
                            list_qs=qs[flag_valid_obs]
                            list_geoid=geoid[flag_valid_obs]
                            nparr_out=np.array([list_x,list_y,list_z,list_t,list_z_sigma,list_geoid,list_qs]).transpose()
                            if len(list_x)!=0:
                                xsub,ysub=subsample_xy(list_x,list_y,scale_inv=scale_gtrk_inv)
                                feature = ogr.Feature(lyr.GetLayerDefn())
                                # Set the attributes using the values from the delimited text file
                                feature.SetField('source', filename_npy)
                                feature.SetField('atl06', os.path.basename(filename_h5))
                                feature.SetField('track', int(str_reftrack))
                                feature.SetField('gt', str_gt)
                                feature.SetField('cycle', int(str_cycle))
                                feature.SetField('region', int(str_region))

                                str_wkt_footprint=xy_2_wktline(xsub,ysub)
                                feat_geom=ogr.CreateGeometryFromWkt(str_wkt_footprint)
                                feature.SetGeometry(feat_geom)
                                lyr.CreateFeature(feature)
                                feature=None

                                #save the numpy array
                                if not os.path.exists(os.path.dirname(filename_npy)):
                                    os.makedirs(os.path.dirname(filename_npy))
                                np.save(filename_npy,nparr_out)
                    else:
                        print(str_gt,'does not exist in h5 file. Skipping.')

    datasrc=None

    
def hdf5_to_npy_parallel(list_hdf5,path_npy,epsg_out=3031,numworker=4,wait_between_worker=5):
    procs=[]
    for i in range(numworker):
        #hdf5_to_npy(list_hdf5,path_npy,epsg_out=3031,part=0,npart=1):
        proc=multiprocessing.Process(target=hdf5_to_npy,
                                     args=(list_hdf5,path_npy,epsg_out,i,numworker))
        procs.append(proc)
        proc.start()
        time.sleep(wait_between_worker)
    
    for proc in procs:
        proc.join()


if __name__=='__main__':
    str_usage='''
    Usage example:
    import_h5_4.py [hdf5 path] [npy path] [target EPSG] [number of workers]
    '''
    
    if len(sys.argv)==5:
        path_hdf5=sys.argv[1]
        path_import=sys.argv[2]
        epsg_out=int(sys.argv[3])
        ncpu=int(sys.argv[4])
        flag_proceed=True
    else:
        print('Please strictly follow the argument format as below.')
        print(str_usage)
        flag_proceed=False

    if flag_proceed:
        #print out the settings
        print('Input HDF5 path  :',path_hdf5)
        print('Output path      :',path_import)
        print('Output EPSG      :',epsg_out)
        print('# workers        :',ncpu)

        if not os.path.exists(path_import):
            os.makedirs(path_import)
        
        list_atl06=glob.glob('{}/ATL06_*.h5'.format(path_hdf5))
        hdf5_to_npy_parallel(list_atl06,path_import,epsg_out=epsg_out,numworker=ncpu,wait_between_worker=1)
        list_subcat=glob.glob('{}/catalog_*of*.shp'.format(path_import))
        merge_subtatalog(list_subcat,'{}/catalog_all.shp'.format(path_import))
        