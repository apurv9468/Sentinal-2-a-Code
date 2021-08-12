#================================////////////////////////////////===============
##====================Below paramteres are to be changed...(keeping intermediate_functions.py file in same directory)
import rasterio
import rasterio as rio
import pandas as pd
import gdal, ogr, os, osr,shutil,json,gc
import numpy as np
from sentinelhub import opensearch,get_tile_info,AwsTileRequest,WebFeatureService, BBox, CRS, DataCollection, SHConfig
import sentinelhub
import georaster as gr
from shapely import speedups
speedups.disable()
from rasterio.plot import show
from rasterio.mask import mask
from shapely.geometry import mapping
import matplotlib.pyplot as plt
import geopandas as gpd
import earthpy as et
import earthpy.plot as ep
import earthpy.spatial as es
import cartopy as cp
import sys


path=r'D:\Apurv_consolidate'
startdate="2016-01-01"
enddate="2021-04-01"

tiles_list = ['18GYP']  #

#tiles_list=['19HDE']#],'19HCD']#,'19HCC','19HCB','19HDD','19HDC','19HDB']
dataCoveragePercentage_condition=80
cloudyPixelPercentage_condition=50
snowCover_condition=50
import geopandas as gpd
#Himporting  19HCD.shp file(all files should be here)
#crop_extent_soap = gpd.read_file(path+'\\19HCD.shp')
#crop_extent_soap = gpd.read_file(r'C:\Users\Administrator\Desktop\shape file\19HCD.shp')
#================================////////////////////////////////===============


os.chdir(path)
sys.path.append(path)
from intermediate_functions import download_tile_aws,all_tiles_validate,data_tile,jp2_to_gtiff_10,jp2_to_gtiff_60,jp2_to_gtiff_20,tiff2_raster,apply_logic,array2raster
def check_disk(Drive_letter,size):
    import shutil
    if shutil.disk_usage(Drive_letter+':')[2]<size*1024*1024*1024:
        input("Press the <ENTER> key to continue...") 

dates_list=[]
date_list=pd.date_range(startdate,enddate,freq='D')
for date in date_list:
    date=date.strftime("%Y-%m-%d")
    dates_list.append(date)
len(dates_list)

### code for  passing data all at once !!
#path=r'C:\Users\Administrator\AwsData\19HDB'
#ates_list=['2019-02-03']
array_index=dict()
extn = {
    'tileinfo': '.json',
    'jp2': '.jp2',
    'tif_translate': 'gdal_translate.tif',
    'tif_index': 'index.tif',
    'tif_logic': 'logic.tif',
    'csv_index': 'index.csv',
    'png_stats': 'calctime.png'
    }
#tiles_list=['19HDB']#,'19HCC','19HCB','19HDD','19HDC','19HDB']


path_c=path
#dates_list=dates_list[dates_list.index('2020-03-05'):]#12RUU,2017-06-07,0
#     date = dates_list[0]
for tile_name in tiles_list:
    folpath=tile_name+' '+startdate[:7]+'_'+enddate[:7]
    if not os.path.exists(folpath):
        os.makedirs(folpath)
    path=path_c+'\\'+folpath
    for date in dates_list:
        tile_name=tile_name
        aws_index=0
        string_name=tile_name+','+date+','+str(aws_index)
        filtered_data= path+'\\filtered_data'
        check_disk(path[0],5)#checking the disk space to continue
        tile_found,tile_origin_coordinates=download_tile_aws(tile_name,date,path,dataCoveragePercentage_condition,cloudyPixelPercentage_condition,snowCover_condition)
        if(tile_found==True):
            if(all_tiles_validate(filtered_data,string_name,tile_name,date)==True):
                tile_path=filtered_data+'\\'+string_name
                #shutil.copy2(tile_path+'/'+'B02-gdal_translate.tif',tile_path+'/'+'B02-gdal_translate_copy.tif')
                bands=data_tile(filtered_data,string_name)
                #print(tile_dataframe)
                ## Code for converting all 10m jp2 bands itno Gtfiff and resizing it to 20M resoltuion 
                #path of all 10M resolution
                path_10 = bands[bands['res']==10]['path'].values
                #name os all 10M resoltuion bands 
                band_10 = bands[bands['res']==10]['band'].values
                for i in range(len(band_10)):
                     jp2_to_gtiff_10(band_10[i],path_10[i])

                ## code for iteraring over all the bands of 20M resolution and converting it itno Gtiff
                path_20 = bands[bands['res']==20]['path'].values
                band_20 = bands[bands['res']==20]['band'].values

                path_60 = bands[bands['res']==60]['path'].values
                band_60 = bands[bands['res']==60]['band'].values
                for i in range(len(band_20)):
                     jp2_to_gtiff_20(band_20[i], path_20[i])

                for i in range(len(band_60)):
                     jp2_to_gtiff_60(band_60[i],path_60[i])        
                # storing all band names and path names in band_all and path_all respectively 
                path_all = bands['path'].values
                band_all = bands['band'].values

#                print(path_all)
#                print(band_all)

                #creating a dictionary for storing band name as key and its raster matrix as value 
                r=dict()
                meta_data_raster=[]
                for i in range(len(band_all)):
                    # storing band name as key and raster matrix as its value in r dictionary
                    r[band_all[i]],meta_data_raster_=tiff2_raster(band_all[i]+'-gdal_translate',path_all[i])
                    meta_data_raster.append(meta_data_raster_)
                    
                    
                
                logic=apply_logic(r)
                # rams index calculation
                index = np.divide(r['B04'],r['B12'])   # Discriminates for rams relative to hematite-goethite
                print("Index Matrx: ",index)
                # find pixels that pass every condition
                valid_pixel = logic['t1'] & logic['t2'] & logic['t3'] & logic['t4'] & logic['t5'] & logic['t6']
                print( "Valid Pixel matrix:",valid_pixel)


                print("Shape of Index matrix :", index.shape)
                print("Shape of valid pixel matrix", valid_pixel.shape)

                #print(np.sum(index))
                #np.multiply(index,valid_pixel)
                #np.matmul(index,valid_pixel)
                #index
                #type(index)

                index= np.multiply(valid_pixel,index)
                print(np.sum(index ))
                
                
                array_index[string_name]=index
                #array_index.append(index)
                #shutil.copy2(tile_path+'/'+'B02-gdal_translate.tif',tile_path+'/'+'B02-gdal_translate_copy.tif')
                rasterOrigin =tile_origin_coordinates
                pixelWidth = 5490
                pixelHeight = 5490
                
                path_out_index = os.path.join(tile_path,'index_file.tif')
                path_converted_index=os.path.join(tile_path,'index_file_converted.tif')
                meta_data_raster[0]['dtype']='float64' 
                # saving index file in the tif format !!
                with rio.open(path_out_index, 'w', **meta_data_raster[0]) as ff:
                     ff.write(index, 1)
                
                
                command='gdalwarp'+str(' "') +str(path_out_index)+str('" "') +str(path_converted_index)+str('" ')+str('-t_srs EPSG:4326')
                os.system(command)
                gc.collect()
                

            else:
                print("Some bands are missing for the tile",tile_name,"for the date", date)
                print("----Not proceeded--")
        gc.collect()
path=path_c


print("-------Sum before replacing-------")
for key,value in array_index.items():
    print(np.sum(array_index[key]))
    
print("-----sum after rep;acing all the values less than 1 with 0")
for key,value in array_index.items():
    array_index[key][array_index[key]<1]=0
    print(np.sum(array_index[key]))
    

for key,value in array_index.items():
        shape=array_index[key].shape
        break
    
final_matrix=np.zeros(shape)

for key ,value in array_index.items():
    print("Adding the index matrix coressponding to",key)
    final_matrix=final_matrix+array_index[key]

path_out = path+'\\averaged_tif\\averaged_index.tif'

os.chdir(path)
if not os.path.exists('averaged_tif'):
    os.makedirs('averaged_tif')


#path_out = r"C:\Users\Administrator\AwsData\averaged_tif\averaged_index.tif"
with rio.open(path_out, 'w', **meta_data_raster[0]) as ff:
    ff.write(final_matrix.astype('float64'), 1)
    
'''    
with rio.open(path_out, 'w', **meta_data_raster) as ff:
    ff.write(final_matrix.astype('uint16'), 1)
'''
path_converted_index = path+'\\averaged_tif\\averaged_index_converted_epsg_4326.tif'

#path_converted_index=r"C:\Users\Administrator\AwsData\averaged_tif\averaged_index_converted_epsg_4326.tif"
command='gdalwarp'+str(' ') +str(path_out)+str(' ') +str(path_converted_index)+str(' ')+str('-t_srs EPSG:4326')
print(command)
os.system(command)


#soap_chm_path =r'C:\Users\Administrator\AwsData\test_data\19HCD,2019-11-30,0\index_file_converted.tif'
soap_chm_path=os.path.join(tile_path,'index_file_converted.tif')
soap_chm_path=r'C:\Apurv_consolidate\19KBB 2016-01_2021-04\filtered_data\19KBB,2016-03-01,0\B01-gdal_translate.tif'
with rio.open(soap_chm_path) as src:
    print(np.array(src))
    lidar_chm_im = src.read(masked = True)[0] #numpy array and applying mask
    extent = rio.plot.plotting_extent(src) #Returns an extent in the format needed
    soap_profile = src.profile


soap_chm_path1=r'C:\Apurv_consolidate\19KBB 2016-01_2021-04\filtered_data\19KBB,2016-03-01,0\B03-gdal_translate.tif'
soap_chm_path2=r'C:\Apurv_consolidate\19KBB 2016-01_2021-04\filtered_data\19KBB,2016-03-01,0\B02-gdal_translate.tif'
with rasterio.open(soap_chm_path1) as src1:
    with rasterio.open(soap_chm_path2) as src2:
        df1=pd.DataFrame(src1.read()[0])
        df2=pd.DataFrame(src2.read()[0])
        df1=df1.iloc[:1000,:1000]
        df2=df2.iloc[:1000,:1000]
        df3=df1/df2
        df3.to_csv('ratio-gdal.csv')
    
ep.plot_bands(lidar_chm_im,
               cmap='terrain',
               extent=extent,
               #title="Lidar Canopy Height Model (CHM)\n NEON SOAP Field Site",
               cbar=True,
             scale=True)


# In[ ]:
#==========using shape file========
#How this 19HCD.shp file is created?
# plot the data
fig, ax = plt.subplots(figsize = (6, 6))
crop_extent_soap.plot(ax=ax)
ax.set_title("Shape_file_19HCD", 
             fontsize = 16)
ax.set_axis_off();


fig, ax = plt.subplots(figsize=(20, 20))
ep.plot_bands(lidar_chm_im,
              cmap='gist_earth',
              extent=extent,
              ax=ax,
              cbar=False)
crop_extent_soap.plot(ax=ax, alpha=.6, color='r');


# In[ ]:


with rio.open(soap_chm_path) as src:
    lidar_chm_crop, soap_lidar_meta = es.crop_image(src, crop_extent_soap)

# Update the metadata to have the new shape (x and y and affine information)

soap_lidar_meta.update({"driver": "GTiff",
                 "height": lidar_chm_crop.shape[1],
                 "width": lidar_chm_crop.shape[2],
                 "transform": soap_lidar_meta["transform"]})


cr_ext = rio.transform.array_bounds(soap_lidar_meta['height'], 
                                            soap_lidar_meta['width'], 
                                            soap_lidar_meta['transform'])

bound_order = [0,2,1,3]
cr_extent = [cr_ext[b] for b in bound_order]
cr_extent, crop_extent_soap.total_bounds


# In[ ]:
lidar_chm_crop_ma = np.ma.masked_equal(lidar_chm_crop[0], -9999.0) 
ep.plot_bands(lidar_chm_crop_ma, cmap='gist_earth', cbar=False)

# Save to disk so you can use the file later.
path_out = path+"\\b02_cropped\\b02_cropped.tif"

os.chdir(path)
if not os.path.exists('b02_cropped'):
    os.makedirs('b02_cropped')

#
with rio.open(path_out, 'w', **soap_lidar_meta) as ff:
    ff.write(lidar_chm_crop[0], 1)
