
## Code for dowanliding data with tile info and applying filtering criteria to dowanload tiles within a time interval 
import gdal
import os
import json
import rasterio
import pandas as pd
import numpy as np
import shutil
import gdal, ogr, os, osr
import numpy as np
import rasterio as rio
import gc
from sentinelhub import WebFeatureService, BBox, CRS, DataCollection, SHConfig
from sentinelhub import AwsTileRequest
import sentinelhub
from sentinelhub import opensearch
from sentinelhub  import get_tile_info
import os
import json
import pandas as pd
import georaster as gr

def download_tile_aws(tile_name,date,path,dataCoveragePercentage_condition,cloudyPixelPercentage_condition,
                      snowCover_condition):
   # print("=====start here======")
  
    #bands = ['B8A']
    bands = ['B01','B02']
    metafiles = ['tileInfo', 'preview','metadata','TCI']
    metafiles = ['tileInfo', 'preview', 'qi/MSK_CLOUDS_B00']
    metafiles_2=['productInfo', 'tileInfo', 'metadata']

    data_folder = path
    filtered_data= path+'\\filtered_data'


    #for date in dates_list:
    try:
        request = AwsTileRequest(
                               tile=tile_name,
                               time=date,
                               aws_index=None,
                               bands=bands,
                               metafiles=metafiles,
                               data_folder=data_folder,
                               data_collection=DataCollection.SENTINEL2_L1C)
        request.save_data()  # This is where the download is triggered


    ## applying filter creiteria for land and cloud coverage !!
    ## applying filter creiteria for land and cloud coverage !!

        path=path
        aws_index=0
        string_name=tile_name+','+date+','+str(aws_index)
#        string_name=date#+','+str(aws_index)

        json_path=path+'\\'+string_name
    #print(json_path)
        json_file_name=json_path+'\\'+'tileInfo.JSON'
    #print(json_file_name)

    # opening the json file 
        f = open(json_file_name)
    # reading the data in tleiinfo json
        tile_data=json.load(f)

    # data for snow cover over the tile
        tile_data_snow=get_tile_info(tile_name,date, aws_index=None, all_tiles=True) 
    #applying filtering criteria======================================================================================================================================

        if (tile_data['dataCoveragePercentage']> dataCoveragePercentage_condition) and (tile_data['cloudyPixelPercentage' ]<cloudyPixelPercentage_condition) and (tile_data_snow[0]['properties']['snowCover']<snowCover_condition):

        ## if all filters are satisfied then we will downloads all the bands of the data in the filetred date folder !
            request = AwsTileRequest(
                      tile=tile_name,
                      time=date,
                      aws_index=aws_index,
                      bands=None,
                      metafiles=metafiles,
                      data_folder=filtered_data,
                      data_collection=DataCollection.SENTINEL2_L1C)

            request.save_data()  # This is where the download is triggered
            return True,tile_data["tileOrigin"]["coordinates"]
        else:
            print("tile is found but does not clear the filtering criteria on",date)
            return False,[0,0]
    except:
            print("No tile is found on this particular Date:",date)
            return False,[0,0]

        
        
## code for validating the presence  of all the 13 Jp2 bands in the folder
def all_tiles_validate(filtered_data,string_name,tile_name,date):
    
    path_folder=filtered_data+'\\'+string_name
    os.listdir(path_folder)

    bands_list=[
     'B01.jp2',
     'B02.jp2',
     'B03.jp2',
     'B04.jp2',
     'B05.jp2',
     'B06.jp2',
     'B07.jp2',
     'B08.jp2',
     'B09.jp2',
     'B10.jp2',
     'B11.jp2',
     'B12.jp2',
     'B8A.jp2']
    
    flag_band=True
    count=0;
    band_missing=[]

    for i in bands_list:
        if i  not in os.listdir(path_folder):
            band_missing.append(i)
    
    
    metafiles = ['tileInfo', 'preview','metadata','TCI']
    metafiles = ['tileInfo', 'preview', 'qi/MSK_CLOUDS_B00']
 
    request = AwsTileRequest(
                               tile=tile_name,
                               time=date,
                               aws_index=None,
                               bands=band_missing,
                               metafiles=metafiles,
                               data_folder=filtered_data,
                               data_collection=DataCollection.SENTINEL2_L1C,
                                )
    request.save_data()
  
    for i in bands_list:
        if i  not in os.listdir(path_folder):
            flag_band=False
            break
    if flag_band==True:
        print("----All 13  jp2 bands Found-----")
        return True
        pass
    else:
        print("----Some jp2 Bands are missing in the folder----- ")
        return False
    
    
## code for creating the data frame having info about bands

def data_tile(filtered_data,string_name):

##creating a dictionary to store the extension type of eevery file 
    path_folder=filtered_data+'\\'+string_name
    extn = {
    'tileinfo': '.json',
    'jp2': '.jp2',
    'tif_translate': 'gdal_translate.tif',
    'tif_index': 'index.tif',
    'tif_logic': 'logic.tif',
    'csv_index': 'index.csv',
    'png_stats': 'calctime.png'
    }

    # creating a datframe of band with tier resoltuion 
    bands = pd.DataFrame({'band': ['B02','B03','B04','B08'], 'res': 10}).append(pd.DataFrame({'band': ['B05', 'B06','B07','B8A', 'B11', 'B12'], 'res': 20}))
    bands=bands.append(pd.DataFrame({'band': ['B01', 'B09','B10'], 'res': 60}),ignore_index=True)
    bands = bands.reset_index(drop=True)
    bands

    bands['file'] = bands['band'].apply(lambda x: x+extn['jp2'])

    bands['path'] = path_folder
    bands['exists'] = bands.apply(lambda x: os.path.exists(os.path.join(x['path'], x['file'])), axis=1)
    return bands


## Code for converting jp2 files to gtiff 
## and increasing the resolution of 10m band to 20m
## so that they can be divided in the later part for calcautng index !
## 10m resoltuion band is of higher pixel and it needs to be reduced to lower pixel band of 20m
# Convert 10m JP2 files to GTiff
print("gdal_translate jp2 -> tif -outsize 50%")


def jp2_to_gtiff_10(band, path):
    print(" ", band)
    extn = {
    'tileinfo': '.json',
    'jp2': '.jp2',
    'tif_translate': 'gdal_translate.tif',
    'tif_index': 'index.tif',
    'tif_logic': 'logic.tif',
    'csv_index': 'index.csv',
    'png_stats': 'calctime.png'
    }
    jp2_loc = os.path.join(path, band+extn['jp2']+'" ')
    #print(jp2_loc)
    tif_loc = os.path.join(path, band+'-'+extn['tif_translate'])
    #print(tif_loc)
    command='gdal_translate -outsize 50% 50% -of GTiff "'+str(jp2_loc)+str(' "')+str(tif_loc)+'"'#
    # here we outsize the band(Reduce the pixels) by 50% from 10M to 20M resoltuion
    #print(command)
    os.system(command) #The command can also be passed as a string, instead of a variable



# Code for Converting 20m JP2 files to GTiff 
## Note: Here we will not resize the bands
print("gdal_translate jp2 -> tif")

def jp2_to_gtiff_20(band, path):
    print(" ", band)
    extn = {
    'tileinfo': '.json',
    'jp2': '.jp2',
    'tif_translate': 'gdal_translate.tif',
    'tif_index': 'index.tif',
    'tif_logic': 'logic.tif',
    'csv_index': 'index.csv',
    'png_stats': 'calctime.png'
    }
    jp2_loc = os.path.join(path, band+extn['jp2']+'" ')
    #print(jp2_loc)
    tif_loc = os.path.join(path, band+'-'+extn['tif_translate'])
    #print(tif_loc)
    command='gdal_translate -of GTiff "'+str(jp2_loc)+str(' "')+str(tif_loc)+'"'#


    '''jp2_loc = os.path.join(path, band+extn['jp2']+' ')
    print(jp2_loc)
    tif_loc = os.path.join(path, band+'-'+extn['tif_translate'])
    print(tif_loc)
    #Here we do not need to resize the 20M and 60M band 
    command='gdal_translate -of GTiff '+str(jp2_loc)+str( )+str(tif_loc)
    print(command)'''

    os.system(command) #The command can also be passed as a string, instead of a variable
    
    
    
def jp2_to_gtiff_60(band, path):
    print(" ", band)
    extn = {
    'tileinfo': '.json',
    'jp2': '.jp2',
    'tif_translate': 'gdal_translate.tif',
    'tif_index': 'index.tif',
    'tif_logic': 'logic.tif',
    'csv_index': 'index.csv',
    'png_stats': 'calctime.png'
    }
    jp2_loc = os.path.join(path, band+extn['jp2']+'" ')
    #print(jp2_loc)
    tif_loc = os.path.join(path, band+'-'+extn['tif_translate'])
    #print(tif_loc)
    command='gdal_translate -outsize 300% 300% -of GTiff "'+str(jp2_loc)+str(' "')+str(tif_loc)+'"'#
    os.system(command) #The command can also be passed as a string, instead of a variable
      
'''jp2_loc = os.path.join(path, band+extn['jp2']+' ')
print(jp2_loc)
tif_loc = os.path.join(path, band+'-'+extn['tif_translate'])
print(tif_loc)
#Here we do not need to resize the 20M and 60M band 
command='gdal_translate -of GTiff '+str(jp2_loc)+str( )+str(tif_loc)
print(command)'''

## coverting all the tif files into rasters obejcts 

print("gtiff---->>>>raster")

def tiff2_raster(band,path):
    
    path_tif=path+'\\'+band+'.tif'
    # raeding tif file into raster format 
    raster=rasterio.open(path_tif)
    type(raster)
    #reading the only 1 band that is present in the rastor format
    band = raster.read(1)
    raster_matrix=band/10000
    meta_data_raster=raster.meta
    
    print("Shape of the raster_matrix --",band.shape)
    print("Number of Bands Present in the Raster object:",raster.count)
    print("Addtional Meta data of Rastor Object:",raster.meta)
    print("Raster Conversion Successful")
    
    return raster_matrix,meta_data_raster


## Applying logic test and creating logic dictionary

def apply_logic(r):
#Creating a dictionary named logic  having test name as key and its coreesponging boolean(true or false)
#matrix as the value
    logic = dict()

    logic['t1'] = r['B11'] > 0.1 # Filters out water
    logic['t2'] = r['B02'] < 0.4 # Filters out cloud
    logic['t4'] = r['B04'] > 0.1 # Filters out shadow


    logic['t3']=np.divide(r['B04'],r['B8A'])>0.6      ## filters out the vegeation 
    logic['t5'] = np.divide(r['B05'],r['B11'])<0.85   # Filters out topographical effects
    logic['t6'] =np.divide(r['B06'],r['B8A']) > 1     # Filters for the 900 wavelength Fe absorption feature
    return logic

#'B05-gdal_translate.tif', 'B11-gdal_translate.tif', how=align_fun) < 0.85 # Filters out topographical effects
# B06-gdal_translate.tif', 'B8A-gdal_translate.tif', how=align_fun) > 1 # Filters for the 900 wavelength Fe absorption feature
#!gdal_calc.py -A B04-gdal_translate.tif -B B8A-gdal_translate.tif --outfile=logic_3.tif --calc="A/B"
#'B04-gdal_translate.tif', 'B8A-gdal_translate.tif', how=align_fun) > 0.6 # Filters out vegetation


def array2raster(newRasterfn,rasterOrigin,pixelWidth,pixelHeight,array):
    cols = array.shape[1]
    rows = array.shape[0]
    originX = rasterOrigin[0]
    originY = rasterOrigin[1]

    driver = gdal.GetDriverByName('GTiff')
    outRaster = driver.Create(newRasterfn, cols, rows, 1,gdal.GDT_Float32)
    outRaster.SetGeoTransform((originX, pixelWidth, 0, originY, 0, pixelHeight))
    outband = outRaster.GetRasterBand(1)
    outband.WriteArray(array)
    outRasterSRS = osr.SpatialReference()
    outRasterSRS.ImportFromEPSG(4326)
    outRaster.SetProjection(outRasterSRS.ExportToWkt())
    outband.FlushCache()
    output_raster = None
    return False,[0,0]
        
