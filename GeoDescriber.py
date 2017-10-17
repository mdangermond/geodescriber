#-------------------------------------------------------------------------------
#
#   Copyright 2016 Esri
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#
#   you may not use this file except in compliance with the License.
#
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#
#   distributed under the License is distributed on an "AS IS" BASIS,
#
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#
#   See the License for the specific language governing permissions and
#
#   limitations under the License.?
#
#
#
#-------------------------------------------------------------------------------
# Name:         GeoDescriber
#               Describe a polygon based on Landcover, Landform, Bioclimate, Lithology.
#               Analyzes ten landscape services from the Living Atlas and synthesizes text.
# Authors:      Deniz Karagulle and Michael Dangermond
# Created:      10/14/2015, revised 10/17/2017
# Notes:        Please change line 87 to use your own polygon feature class.
#               Change line 65 to a valid scratch/temporary geodatabase,
#               and change line 98-100 to use your own credentials.
#
#               As of 9/20/2017 a bug in arcpy.JoinField_management() in 10.5.1 joins
#               a field to an in_memory table four times. Until this is resolved use the
#               function JoinField_Workaround() instead.
#-------------------------------------------------------------------------------


try:
    # Import modules and start the clock
    import os, sys, math, time, os.path
    import arcpy
    import traceback
    import warnings
    arcpy.CheckOutExtension("Spatial")
    #from arcpy import env
    from arcpy.sa import *
    from random import shuffle
    import tempfile
    import shutil
    import multiprocessing as mp
    import numpy as np
    from uuid import uuid4

    global currentTime
    startTime = time.clock()
    currentTime = 0
    print("Start time is "+str(startTime))
    TEMP= os.getenv("TEMP")

    arcpy.env.overwriteOutput=True
    #output is the geodatabase containing some rasters created by GeoDescriber.
    output= r"C:\gis\GeoDescriber\GeoDescriber.gdb"
    arcpy.env.workspace = "in_memory"
    arcpy.env.scratchWorkspace ="in_memory"
    arcpy.env.compression = "LZW"
    #os.environ["TEMP"] = arcpy.env.scratchWorkspace
    #os.environ["TMP"] = arcpy.env.scratchWorkspace
    tempspace = output
    #tempspace = r"C:\gis\GeoDescriber\notinmem.gdb"
    inmem = "in_memory"

    #output=arcpy.GetParameterAsText(0)
    #set the cellsize to the global default for the population and ELU datasets: 231.9156058
    cellsize = 231.9156058
    #Set the spatial reference to the Mollweide Projection (54009), the equal area projection shown to have
    #the least loss of accuracy in counting up the area of geographic phenomena from a global dataset.
    sr = arcpy.SpatialReference(54009)
    arcpy.env.outputCoordinateSystem = sr
    arcpy.env.overwriteOutput=True
    arcpy.SetLogHistory(False)
    #---------------------------------------------------------------------------
    #---------Provide the path to the polygon you would like to characterize----
    #inFeatureClass = arcpy.GetParameterAsText(0)
    inFeatLyr=r"C:\gis\GeoDescriber\GeoDescriber.gdb\Austria"

except:
    print("A connection to a valid Spatial Analyst license, a valid")
    print("ArcInfo license, and the internet is required to run this script.")
    raise

#-------------------------------------------------------------------------------
#-------------------Use your own Credentials for Landscape Account--------------

# set parameter for username and password
userName = ""
#userName = arcpy.GetParameterAsText(2)
passWord = ""
#passWord = arcpy.GetParameterAsText(3)

#-------------------------------------------------------------------------------
#--------------------------dictionaries-----------------------------------------
# These python dictionaries translate every possible category of the
# four main ecophysiographic criteria into language that can be
# made into better-understood sentences.

# lithology dictionary. Mostly removes capital letters from classes and makes them plural
lithology_dict = {
'Acid Plutonics': 'acid plutonics',
'Acid Volcanic': 'acid volcanics',
'Basic Plutonics': 'basic plutonics',
'Basic Volcanics': 'basic volcanics',
'Carbonate Sedimentary Rock': 'carbonate sedimentary rocks',
'Evaporite': 'areas of evaporite',
'Ice and Glaciers': 'areas of ice and glaciers',
'Intermediate Plutonics': 'intermediate plutonics',
'Intermediate Volcanics': 'intermediate volcanics',
'Metamorphics': 'metamorphics',
'Mixed Sedimentary Rock': 'mixed sedimentary rocks',
'Non-defined': 'areas of undefined lithology',
'Pyroclastics': 'pyroclastics',
'Siliciclastic Sedimentary Rock': 'non-carbonate sedimentary rocks',
'Unconsolidated Sediment': 'areas of unconsolidated sediment',
'None': 'areas of undefined lithology'
}

# bioclimate dictionary. Makes each class lower case and work better in a sentence.
bioclimate_dict = {
'Arctic': 'arctic',
'Cold Dry': 'cold and dry',
'Cold Moist': 'cold and moist',
'Cold Semi-Dry': 'cold and semi-dry',
'Cold Very Dry': 'cold and very dry',
'Cold Very Wet': 'cold and very wet',
'Cold Wet': 'cold and wet',
'Cool Dry': 'cool and dry',
'Cool Moist': 'cool and moist',
'Cool Semi-Dry': 'cool and semi-dry',
'Cool Very Dry': 'cool and very dry',
'Cool Very Wet': 'cool and very wet',
'Cool Wet': 'cool and wet',
'Hot Dry': 'hot and dry',
'Hot Moist': 'hot and moist',
'Hot Semi-Dry': 'hot and semi-dry',
'Hot Very Dry': 'hot and very dry',
'Hot Very Wet': 'hot and very wet',
'Hot Wet': 'hot and wet',
'Very Cold Dry': 'very cold and dry',
'Very Cold Moist': 'very cold and moist',
'Very Cold Semi-Dry': 'very cold and semi-dry',
'Very Cold Very Dry': 'very cold and very dry',
'Very Cold Very Wet': 'very cold and very wet',
'Very Cold Wet': 'very cold and wet',
'Very Hot Dry': 'very hot and dry',
'Very Hot Moist': 'very hot and moist',
'Very Hot Semi-Dry': 'very hot semi-dry',
'Very Hot Very Dry': 'very hot very dry',
'Very Hot Very Wet': 'very hot very wet',
'Very Hot Wet': 'very hot and wet',
'Warm Dry': 'warm and dry',
'Warm Moist': 'warm and moist',
'Warm Semi-Dry': 'warm and semi-dry',
'Warm Very Dry': 'warm and very dry',
'Warm Very Wet': 'warm and very wet',
'Warm Wet': 'warm and wet',
'None': 'undefined climate'
}

# land cover dictionary. Simplifies classes and makes them work better in a sentence.
landcover_dict = {
'Cropland, rainfed': 'rainfed cropland',
'Cropland, rainfed - Herbaceous cover': 'rainfed herbaceous cropland',
'Cropland, rainfed - Tree or shrub cover': 'rainfed tree or shrublike cropland',
'Cropland irrigated or post-flooding': 'irrigated or flooded cropland',
'Mosaic cropland (>50%) / natural vegetation (Tree, shrub, herbaceous cover) (<50%)': 'a mix of cropland with some natural vegetation',
'Mosaic natural vegetation (Tree, shrub, herbaceous cover) (>50%) / cropland (<50%)': 'a mix of natural vegetation with some cropland',
'Tree cover, broadleaved, evergreen, closed to open (>15%)': 'broadleaved evergreen forest',
'Tree cover, broadleaved, deciduous, closed to open (>15%)': 'broadleaved deciduous forest',
'Tree cover, broadleaved, deciduous, closed (>40%)': 'dense broadleaved deciduous forest',
'Tree cover, broadleaved, deciduous, open (15-40%)': 'sparse broadleaved deciduous forest',
'Tree cover, needleleaved, evergreen, closed to open (>15%)': 'needle-leaved evergreen forest',
'Tree cover, needleleaved, evergreen, closed (>40%)': 'dense needle-leaved evergreen forest',
'Tree cover, needleleaved, evergreen, open (15-40%)': 'sparse needle-leaved evergreen forest',
'Tree cover, needleleaved, deciduous, closed to open (>15%)': 'needle-leaved deciduous forest',
'Tree cover, needleleaved, deciduous, closed (>40%)': 'dense needle-leaved deciduous forest',
'Tree cover, needleleaved, deciduous, open (15-40%)': 'sparse needle-leaved deciduous forest',
'Tree cover, mixed leaf type (broadleaved and needleleaved)': 'mixed forest',
'Mosaic Trees and shrub (>50%) / herbaceous cover (<50%)': 'trees or shrubs with some herbaceous cover',
'Mosaic herbaceous cover (>50%) / Trees and shrub (<50%)': 'herbaceous cover with some trees or shrubs',
'Shrubland': 'shrubland',
'Shrubland evergreen': 'evergreen shrubland',
'Shrubland deciduous': 'deciduous shrubland',
'Grassland': 'grassland',
'Lichens and mosses': 'lichens and mosses',
'Sparse vegetation (tree, shrub, herbaceous cover) (<15%)': 'sparse vegetation',
'Sparse shrub (<15%)': 'sparse shrubland',
'Sparse herbaceous cover (<15%)': 'sparse herbaceous cover',
'Tree cover, flooded, fresh or brakish water': 'forest flooded by fresh or brackish water',
'Tree cover, flooded, saline water': 'forest flooded by salt water',
'Shrub or herbaceous cover, flooded, fresh/saline/brakish water': 'flooded herbaceous cover or shrubland',
'Urban areas': 'urban areas',
'Bare areas': 'bare ground',
'Consolidated bare areas': 'consolidated bare ground',
'Unconsolidated bare areas': 'unconsolidated bare ground',
'Water bodies': 'bodies of water',
'Permanent snow and ice': 'permanent ice and snow',
'None': 'undefined land cover'
}

# land cover dictionary. Simplifies classes and makes them work better in a sentence.
short_landcover_dict = {
'Cropland, rainfed': 'cropland',
'Cropland, rainfed - Herbaceous cover': 'cropland',
'Cropland, rainfed - Tree or shrub cover': 'tree or shrub cropland',
'Cropland irrigated or post-flooding': 'cropland',
'Mosaic cropland (>50%) / natural vegetation (Tree, shrub, herbaceous cover) (<50%)': 'cropland and natural vegetation',
'Mosaic natural vegetation (Tree, shrub, herbaceous cover) (>50%) / cropland (<50%)': 'natural vegetation and cropland',
'Tree cover, broadleaved, evergreen, closed to open (>15%)': 'forest',
'Tree cover, broadleaved, deciduous, closed to open (>15%)': 'forest',
'Tree cover, broadleaved, deciduous, closed (>40%)': 'forest',
'Tree cover, broadleaved, deciduous, open (15-40%)': 'forest',
'Tree cover, needleleaved, evergreen, closed to open (>15%)': 'forest',
'Tree cover, needleleaved, evergreen, closed (>40%)': 'forest',
'Tree cover, needleleaved, evergreen, open (15-40%)': 'forest',
'Tree cover, needleleaved, deciduous, closed to open (>15%)': 'forest',
'Tree cover, needleleaved, deciduous, closed (>40%)': 'forest',
'Tree cover, needleleaved, deciduous, open (15-40%)': 'forest',
'Tree cover, mixed leaf type (broadleaved and needleleaved)': 'forest',
'Mosaic Trees and shrub (>50%) / herbaceous cover (<50%)': 'trees or shrubs and herbaceous cover',
'Mosaic herbaceous cover (>50%) / Trees and shrub (<50%)': 'herbaceous cover and trees or shrubs',
'Shrubland': 'shrubs',
'Shrubland evergreen': 'shrubs',
'Shrubland deciduous': 'shrubs',
'Grassland': 'grassland',
'Lichens and mosses': 'lichens and mosses',
'Sparse vegetation (tree, shrub, herbaceous cover) (<15%)': 'sparse vegetation',
'Sparse shrub (<15%)': 'sparse shrubs',
'Sparse herbaceous cover (<15%)': 'sparse cover',
'Tree cover, flooded, fresh or brakish water': 'flooded forest',
'Tree cover, flooded, saline water': 'flooded forest',
'Shrub or herbaceous cover, flooded, fresh/saline/brakish water': 'flooded herbaceous cover or shrubs',
'Urban areas': 'urban areas',
'Bare areas': 'bare ground',
'Consolidated bare areas': 'bare ground',
'Unconsolidated bare areas': 'bare ground',
'Water bodies': 'bodies of water',
'Permanent snow and ice': 'permanent ice and snow',
'None': 'undefined land cover'
}

# landform dictionary. The translations work better in a sentence.
landform_dict = {
'Flat or Nearly Flat Plains': 'flat or nearly flat plains',
'High Hills': 'high hills',
'High Mountains': 'high mountains',
'Irregular Plains with Low Hills': 'irregular plains with low hills',
'Irregular Plains with Moderate Relief': 'irregular plains with moderate relief',
'Low Mountains': 'low mountains',
'Moderate Hills': 'moderate hills',
'Scattered High Hills': 'scattered high hills',
'Scattered High Mountains': 'scattered high mountains',
'Scattered Low Mountains': 'scattered low mountains',
'Scattered Moderate Hills': 'scattered moderate hills',
'Smooth Plains with some local relief': 'smooth plains with some local relief',
'Surface Water': 'bodies of surface water',
'Tablelands with Considerable Relief': 'tablelands with considerable relief',
'Tablelands with High Relief': 'tablelands with high relief',
'Tablelands with Moderate Relief': 'tablelands with moderate relief',
'Tablelands with Very High Relief': 'tablelands with very high relief',
'None': 'undefined landforms'
}

#this matrix attempts to find a shorter, more generic name
#in order to avoid excessive repetition of the proper name
shortname_matrix = [
[" area ","area is","area are"],
["study area","study area is","study area are"],
[" land","land is","land are"],
[" lands","lands are","lands are"],
[" mount ","mountain is","mountain are"],
["mountain","mountain is","mountain are"],
["mountains","mountains are","mountains are"],
["habitat","habitat is","habitat are"],
["neighborhood","neighborhood is","neighborhood are"],
["polder","polder is","polder are"],
[" hill ","hill is","hill are"],
[" hills ","hills are","hills are"],
["polders","polders are","polders are"],
[" site ","site is ","site are"],
[" green ","green is","green are"],
[" greens ","greens is","greens are"],
[" garden ","garden is","garden are"],
[" gardens ","gardens are","gardens are"],
[" quarter","quarter is","quarter are"],
["province","province is","province are"],
["provinces","provinces are","provinces are"],
["municipality","municipality is","municipality are"],
["municipalities","municipalities are","municipalities are"],
["marsh","marsh is","marsh are"],
["marshes","marshes are","marshes are"],
[" flat ","flat is","flat are"],
[" flats ","flats are","flats are"],
[" port ","port is","port are"],
[" ports ","ports are","ports are"],
[" moor ","moor is","moor are"],
["moors","moors are","moors are"],
[" veld","veld is","veld are"],
["veldt","veldt is","veldt are"],
["spring","spring is","spring are"],
["heath","heath is","heath are"],
["icefield","icefield is","icefield are"],
["grassland","grassland is","grassland are"],
["grasslands","grasslands are","grasslands are"],
["canyon","canyon is","canyon are"],
["viewshed","viewshed is","viewshed are"],
[" range ","range is","range are"],
[" plain ","plain is","plain are"],
[" field ","field is","field are"],
[" fort ","fort is","fort are"],
["oblast","oblast is","oblast are"],
["kray","kray is","kray are"],
["krays","krays are","krays are"],
["forest","forest is","forest are"],
["monument","monument is","monument are"],
["property","property is","property are"],
["canton","canton is","canton are"],
["basin","basin is","basin are"],
[" beach ","beach is","beach are"],
[" cape ","cape is","cape are"],
["crater","crater is","crater are"],
["harbor","harbor is","harbor are"],
[" gap ","gap is","gap are"],
["glacier","glacier is","glacier are"],
["island","island is","island are"],
["islands","islands are","islands are"],
["ithsmus","ithsmus is","ithsmus are"],
["plain ","plain is","plain are"],
["plains","plains are","plains are"],
["harbour","harbour is","harbour are"],
[" slope","slope is","slope are"],
["swamp","swamp is","swamp are"],
["valley","valley is","valley are"],
["valleys","valleys are","valleys are"],
[" wood ","wood is","wood are"],
[" woods ","woods are","woods are"],
["tract","tract is","tract are"],
["block","block is","block are"],
[" city ","city is","city are"],
[" town ","town is","town are"],
["village","village is","village are"],
["metropolitan area","metropolitan area is","metropolitan area are"],
["region","region is","region are"],
["district","district is","district are"],
[" erg ","erg is","erg are"],
["delta","delta is","delta are"],
["reef","reef is","reef are"],
["desert","desert is","desert are"],
[" lake ","lake is","lake are"],
["peninsula","peninsula is","peninsula are"],
[" alps","alps are","alps are"],
["volcano","volcano is","volcano are"],
["lagoon","lagoon is","lagoon are"],
[" base ","base is","base are"],
["plateau","plateau is","plateau are"],
["flood plain","flood plain is","flood plain are"],
[" ridge ","ridge is","ridge are"],
["riding","riding is","riding are"],
["constituency","constituency is","constituency are"],
[" camp ","camp is","camp are"],
["pass","pass is","pass are"],
["zone","zone is","zone are"],
["tableland","tableland is","tableland are"],
[" downs","downs are","downs are"],
["state","state is","state are"],
["wilderness","wilderness is","wilderness are"],
["reservoir","reservoir is","reservoir are"],
[" mts","mountains are","mountains are"],
["territory","territory is","territory are"],
["farm","farm is","farm are"],
["farms","farms are","farms are"],
["plantation","plantation is","plantation are"],
["plantations","plantations are","plantations are"],
["group","group is","group are"],
["complex","complex is","complex are"],
["depression","depression is","depression are"],
["republic","republic is","republic are"],
["highland","highland is","highland are"],
["highlands","highlands are","highlands are"],
["prefecture","prefecture is","prefecture are"],
["territories","territories are","territories are"],
["refuge","refuge is","refuge are"],
["enclave","enclave is","enclave are"],
["exclave","exclave is","exclave are"],
["battlefield","battlefield is","battlefield are"],
["oasis","oasis is","oasis are"],
["oases","oases are","oases are"],
["jungle","jungle is","jungle are"],
["jungles","jungles are","jungles are"],
["easement","easement is","easement are"],
["corridor","corridor is","corridor are"],
["sanctuary","sanctuary is","sanctuary are"],
["barrens","barrens are","barrens are"],
["serra","serra are","serra are"],
["sierra","sierra are","sierra are"],
["county","county is","county are"],
["cordillera","cordillera are","cordillera are"],
["bosque","bosque is","bosque are"],
["pampa","pampa is","pampa are"],
["parque","parque is","parque are"],
[" dome","dome is","dome are"],
["gorge","gorge is","gorge are"],
["ravine","ravine is","ravine are"],
[" market","market is","market are"],
["anbaugebiet","anbaugebiet is","anbaugebiet are"],
["reservation","reservation is","reservation are"],
["park","park is","park are"],
["parks","parks are","parks are"],
["watershed","watershed is","watershed are"],
["reservations","reservations are","reservations are"],
["reserve","reserve is","reserve are"],
["reserves","reserves are","reserves are"],
["preserve","preserve is","preserve are"],
["preserves","preserves are","preserves are"],
["conservation area","conservation area is","conservation area are"],
["recreation area","recreation area is","recreation area are"]
]

#-------------------------------------------------------------------------------
#--------------------------functions--------------------------------------------
#
# Calculate the percent of each class in the study area. Put the percents into a table.
# Then join the table back to the original raster so the class text can be used in characterization.
# functions for multi-processing
def p(m):
    global lock
    if lock: lock.acquire()
    #print(m)
    if lock: lock.release()

# multiprocessing pool
def initPool(lk):
    global lock
    lock = lk

# The getResult() function retrieves raster layers from the server asynchronously (using multiprocessing).
# Then it saves each raster layer as a TIF image in a temporary folder. TIF images in a temporary
# folder avoids the worry of locking file geodatabases.
# There are three steps to bringing a raster down off a server for analysis.
# 1. CreateGISServerConnectionFile
# 2. MakeImageServerLayer
# 3. CopyRaster

def getResult(layerInfo):
    """getResult(layerInfo)

    Connects to ArcGIS Online to asynchronously create image server layers.
    Saves each raster as a TIFF image in a temp folder.
    Parameters for the layers are stored in the layerInfo variable.

    """
    try:
        m = ""
        id = str(layerInfo.get('id', str(uuid4().fields[-1])[:5]))
        imageLayer = layerInfo['name']
        outRaster = os.path.join(layerInfo['scratchFolder'], "{0}_{1}.{2}".format(imageLayer, id, "TIF"))
        p(". [{0}]: Begin work on: {1}".format(id, imageLayer))
        p("  . [{0}] Begin: Importing ArcPy".format(id))
        arcpy.env.rasterStatistics = "None 10 10 (0 255)"
        arcpy.env.pyramid = "None -1"
        arcpy.env.overwriteOutput = True
        #arcpy.env.outputCoordinateSystem = sr
        p("  . [{0}] Done:  Importing ArcPy".format(id))
        p("  . [{0}] Begin: Create GIS Server Connection".format(id))
        #print("create gis server connection file...")
        arcpy.mapping.CreateGISServerConnectionFile("USE_GIS_SERVICES", layerInfo['scratchFolder'], layerInfo['service'], layerInfo['serviceURL'],
                                               "ARCGIS_SERVER", '', '', userName, passWord, "SAVE_USERNAME")
        p ( arcpy.GetMessages())
        m += arcpy.GetMessages() + "\n"
        p("  . [{0}] Done:  CreateGISServerConnection".format(id))
        p("  . [{0}] Begin: Make image server layer".format(id))
        #make image server layer and copy raster.
        #print("make image server layer...")
        arcpy.MakeImageServerLayer_management(layerInfo['url'], imageLayer, layerInfo['extentlayer'],"#","#","#","#","#",layerInfo['cellsize'],"#",layerInfo['processingTemplate'])
        p ( arcpy.GetMessages())
        m += arcpy.GetMessages() + "\n"
        p("  . [{0}] Done:  Make image server layer".format(id))
        p("  . [{0}] Begin: Copy raster".format(id))
        #print("copy raster...")
        arcpy.management.CopyRaster(imageLayer, outRaster)
        m += arcpy.GetMessages() + "\n"
        p("  . [{0}] Done:  Copy raster".format(id))
        if not arcpy.Exists(outRaster):
            raise Exception("Result raster '{0}' does not exist.".format(outRaster))
        p("  . [{0}] Begin: Saving results".format(id))
        layerInfo['path'] = outRaster
        layerInfo['messages'] = m
        layerInfo['id'] = id
        p("  . [{0}] Done:  Saving results".format(id))
        p(". [{0}] Done working on: {1}".format(id, imageLayer))
    except Exception as e:
        layerInfo['exception'] = e.message
        return None
    return layerInfo

#This function is nearly identical to getResult() but modified a bit for a second run-through.
#Sometimes not every image returns from the server as requested (3/14/16).
def getResult2(layerInfo2):
    """getResult2(layerInfo2)

    If any datasets are missing after getResult(), getResult2() Connects to
    ArcGIS Online to sequentially create image server layers.
    Each raster is copied to a TIFF image in a temp folder.
    Parameters for the layers are stored in the layerInfo2 variable.

    """
    try:
        m = ""
        id = str(layerInfo2.get('id', str(uuid4().fields[-1])[:5]))
        imageLayer = layerInfo2['name']
        outRaster = os.path.join(layerInfo2['scratchFolder'], "{0}_{1}.{2}".format(imageLayer, id, "TIF"))
        p(". [{0}]: Retrying: {1}".format(id, imageLayer))
        p("  . [{0}] Begin: Importing ArcPy".format(id))
        arcpy.env.rasterStatistics = "None 10 10 (0 255)"
        arcpy.env.pyramid = "None -1"
        arcpy.env.overwriteOutput = True
        #arcpy.env.outputCoordinateSystem = sr
        p("  . [{0}] Done:  Importing ArcPy".format(id))
        p("  . [{0}] Begin: Create GIS Server Connection".format(id))
        #print("create gis server connection file...")
        arcpy.mapping.CreateGISServerConnectionFile("USE_GIS_SERVICES", layerInfo2['scratchFolder'], layerInfo2['service'], layerInfo2['serviceURL'],
                                               "ARCGIS_SERVER", '', '', userName, passWord, "SAVE_USERNAME")
        p ( arcpy.GetMessages())
        m += arcpy.GetMessages() + "\n"
        #print("make image server layer...")
        p("  . [{0}] Done:  CreateGISServerCOnnection".format(id))
        p("  . [{0}] Begin: Make image server layer".format(id))
        arcpy.management.MakeImageServerLayer(layerInfo2['url'], imageLayer, layerInfo2['extentlayer'],"#","#","#","#","#",layerInfo2['cellsize'],"#",layerInfo2['processingTemplate'])
        p ( arcpy.GetMessages())
        m += arcpy.GetMessages() + "\n"
        p("  . [{0}] Done:  Make image server layer".format(id))
        p("  . [{0}] Begin: Copy raster".format(id))
        #print("copy raster...")
        arcpy.management.CopyRaster(imageLayer, outRaster)
        m += arcpy.GetMessages() + "\n"
        p("  . [{0}] Done:  Copy raster".format(id))
        if not arcpy.Exists(outRaster):
            raise Exception("Result raster '{0}' does not exist.".format(outRaster))
        p("  . [{0}] Begin: Saving results".format(id))
        layerInfo2['path'] = outRaster
        layerInfo2['messages'] = m
        layerInfo2['id'] = id
        p("  . [{0}] Done:  Saving results".format(id))
        p(". [{0}] Done working on: {1}".format(id, imageLayer))
    except Exception as e:
        layerInfo2['exception'] = e.message
        return None
    return layerInfo2

#Finds the percentage of classes within a featurelayer, then joins a table of percentages to the studyarea table.
def percent(featurelayer, newlayer, casefield):
    """percent(featurelayer, newlayer, casefield)

    Returns the percent of the study area occupied by every
    significant class found in the study area. These percentages
    are used later to craft sentences in the description.

    """

    if arcversion == '10.5.1':
        ##tempspace = "C:\\gis\\GeoDescriber\\notinmem.gdb"
        print (featurelayer + " is ready for the percent function")
        statFields=[["Count","SUM"]]
        arcpy.Statistics_analysis(featurelayer, newlayer, statFields, casefield)
        arcpy.AddField_management(newlayer, "percent", "DOUBLE")
        sum1=0
        with arcpy.da.SearchCursor(newlayer,["SUM_Count"]) as cursors:
            for row1 in cursors:
                sum1=row1[0]+sum1
        with arcpy.da.UpdateCursor(newlayer,["SUM_Count","percent"]) as cursor:
            for row in cursor:
                row[1]=(row[0]*100)/sum1
                cursor.updateRow(row)
        joinedField=["FREQUENCY","percent"]
        JoinField_Workaround(featurelayer, casefield, newlayer, casefield, joinedField)
    else:
        statFields=[["Count","SUM"]]
        stattable=arcpy.Statistics_analysis (featurelayer, newlayer, statFields, casefield)
        arcpy.AddField_management (stattable, "percent", "DOUBLE")
        sum1=0
        with arcpy.da.SearchCursor(stattable,["SUM_Count"]) as cursors:
            for row1 in cursors:
                sum1=row1[0]+sum1
        with arcpy.da.UpdateCursor(stattable,["SUM_Count","percent"]) as cursor:
            for row in cursor:
                row[1]=(row[0]*100)/sum1
                cursor.updateRow(row)
        joinedField=["FREQUENCY","percent"]
        arcpy.JoinField_management(featurelayer, casefield, stattable, casefield, joinedField)

#find and return the largest value in the percent table.
def largest(layer,field):
    """largest(layer,field)

    Returns the largest area significant class in the study area.

    """
    theitems = []
    rows = arcpy.SearchCursor(layer)
    for row in rows:
        theitems.append(row.getValue(field))
    del rows
    theitems.sort()
    max1= theitems[-1]
    return max1

#find string
def findString(featureCl,largVal):
    """findString(featureCl,largVal)

    Returns the name of the significant class which is the
    largest area in the study area.

    """
    global sideopslist
    with arcpy.da.SearchCursor(featureCl,["ClassName","percent"]) as cursor:
        for row in cursor:
            if row[1]== largVal:
                randomstr = str(row[0])
                sideopslist.append(row[1])
                sideopslist.append(row[0])
                sideopslist.append(i)
                return randomstr

# Find the other land cover, lithology, bioclimate and landform classes that are bigger than 10% but
# smaller than the largest in their classes.
def restofValues (featurelayer, largVal, dictionary):
    """restofValues (featurelayer, largVal, dictionary)

    Returns the names of significant classes making up over 10%
    of the study area but are smaller than the largest class.

    """
    global restopslist
    stringvalues=""
    row_pairs=set()
    z = 0
    with arcpy.da.UpdateCursor(featurelayer,["ClassName","percent"]) as cursor:
            for row in cursor:
                row_pair= tuple(sorted(row))
                if row_pair in row_pairs:
                    cursor.deleteRow()
                else:
                    row_pairs.add(row_pair)
                    if row[1] > int(10) and row[1] != largVal:
                        z = z + 1
                        restopslist.append(row[1])
                        restopslist.append(row[0])
                        restopslist.append(i)
                        stringvalues += str(int(row[1])) +"% is " + dictionary[str(row[0])] + ", "
            restvaluescount = z
            if restvaluescount == 1:
                stringval= stringvalues[:-2]
                restvalues = "but "+stringval+". "
            else:
                stringval=stringvalues[:-2]
                pos = stringval.rfind(',')
                stringv = stringval[:pos] + ', and' + stringval[pos+1:]
                restvalues="while "+stringv+". "
            restvaluescount = 0
    return restvalues


def JoinField_Workaround (indataset,infield,jointable,joinfld,workaroundfields):
    """JoinField_Workaround(indataset,infield,jointable,joinfld,workaroundfields) # joinfield

    Joins a field to a table, but works around the bug in 10.5.1 for joining a
    field to a dataset that is in_memory. Trying to use joinfield_management on
    with the unfixed bug joins the field four times. This gets around that problem.
    """

    try:
        if arcpy.Exists(tempspace+"\\junkg"):
       	    arcpy.Delete_management(tempspace+"\\junkg")
        if arcpy.Exists(tempspace+"\\junktbl"):
       	    arcpy.Delete_management(tempspace+"\\junktbl")
        arcpy.CopyRaster_management(indataset,tempspace+"\\junkg")
        arcpy.CopyRows_management(jointable,tempspace+"\\junktbl")
        if arcpy.Exists(indataset):
       	    arcpy.Delete_management(indataset)
        arcpy.JoinField_management(tempspace+"\\junkg", infield, tempspace+"\\junktbl", joinfld, workaroundfields)
        arcpy.CopyRaster_management(tempspace+"\\junkg",indataset)
        feature = tempspace+"\\junktbl"
        if arcpy.Exists(feature):
       	    arcpy.Delete_management(feature)
        feature = tempspace+"\\junkg"
        if arcpy.Exists(feature):
       	    arcpy.Delete_management(feature)
    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
        # Write Python error messages to log
        err= pymsg + "\n"
        print(err)


#-------------------------------------------------------------------------------
#-----------------retrieve landscape6 and landscape7 rasters--------------------

def CleanUp():
    """CleanUp()

    Cleans up temporary rasters and feature datasets used in
    GeoDescriber, both in_memory and on disk.

    """
    try:
        print("cleaning up...")
        cleanupg = ['aspectg','aspectindexg','Bioclimate_R','Biomass_R','con_extent','cong','Diversity_R','Ecophysdiv_R','Elevation_R','epfcea','extractg','Landcover_R','Landform_R','Lithology_R','mask_extent','mask30','maskext2','mw_zonedg','northupg','Population_R','projcea','Slope_R','thiesmw','thiespts','Water_R','Water30m_R','stattbl','Bioclimates_CR','Bioclimates_ST','Landcover_CR','Landcover_ST','Landform_CR','Landform_ST','Lithology_CR','Lithology_ST','thiesmw','Population_RS','conglf','temp','temppts','cf0']
        for fd in cleanupg:
            feature = os.path.join(output,fd)
            if arcpy.Exists(feature):
                arcpy.Delete_management(feature)
                print("deleting "+feature+"...")
            feature = os.path.join(inmem,fd)
            if arcpy.Exists(feature):
                arcpy.Delete_management(feature)
                print("deleting "+feature+"...")

    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
        # Write Python error messages to log
        err= pymsg + "\n"
        print(err)

    print("Cleanup done.")

def GeoDescriber():
    """

    GeoDescriber creates a field called description on a feature
    class, it analyzes each polygon in the feature class, then it
    writes a landscape description in the field for each polygon.

    """
    global sideopslist
    global lock
    global sr
    global GeoDescriberTries
    global article

    try:
        #The GeoDescriber() function loops through each polygon (using the polygon ID).
        #First it creates a feature layer from one polygon (epf) then it makes a raster
        #for that polygon called mask_extent.
        qstr = "OBJECTID" + " = " + str(intPolyID)
        print("making single poly feature class")
        epf = arcpy.MakeFeatureLayer_management(FL_MollPrj,"EachPolyFeat", qstr)
        thefieldname="OBJECTID"
        arcpy.env.cellSize = cellsize
        thisPolyTime = time.clock()
        currentTime = time.clock()
        print(str(currentTime-startTime)+" seconds. Converting single poly feature to raster...")
        arcpy.CopyFeatures_management(epf,inmem+"\\cf0")
        rasterExt=arcpy.PolygonToRaster_conversion(inmem+"\\cf0", "OBJECTID" ,inmem+"\\mask_extent")
        #rasterExt=arcpy.PolygonToRaster_conversion(epf, "OBJECTID" ,inmem+"\\mask_extent")
        extentRaster=Con(arcpy.Raster(rasterExt),1)
        er=extentRaster.save(inmem+"\\con_extent")
        fieldsdk=["SHAPE@"]
        #obtain the minimum and maximum x and y for the polygon. Used as a template when making image server layers.
        with arcpy.da.SearchCursor(epf, fieldsdk) as cursordk:
            for rowdk in cursordk:
                geom = rowdk[0]
                ext = geom.extent  # or row.Shape.extent
                #print ext
                extentlayer_feature= str(ext.XMin) + " " + str(ext.YMin) + " " + str(ext.XMax) + " " + str(ext.YMax)
                #print extentlayer_feature

        layerInfo2 = []
        lock = None
        mp.freeze_support()

        tempFolder = tempfile.mkdtemp(prefix="agd_")    # Warning: this folder is deleted at the end of this script.

        p(". Scratch workspace is: {0}".format(tempFolder))

        serviceURLLandform = os.path.join(tempFolder, "landscape7", "World_Landforms_Improved_Hammond_Method" + ".ImageServer")
        serviceURLLandcover = os.path.join(tempFolder, "landscape7", "World_Land_Cover_ESA_2010" + ".ImageServer")
        serviceURLLithology = os.path.join(tempFolder, "landscape6", "World_Lithology" + ".ImageServer")
        serviceURLBioclimate = os.path.join(tempFolder, "landscape6", "World_Bioclimates" + ".ImageServer")
        serviceURLElevation = os.path.join(tempFolder, "landscape6", "World_Elevation_GMTED" + ".ImageServer")
        serviceURLPopulation = os.path.join(tempFolder, "landscape6", "World_Population_Estimated" + ".ImageServer")
        serviceURLSlope = os.path.join(tempFolder, "landscape6", "World_Slope_GMTED" + ".ImageServer")
        serviceURLWater= os.path.join(tempFolder,"landscape6", "World_Surface_Water_30m_BaseVue_2013" + ".ImageServer")
        serviceURLDiversity = os.path.join(tempFolder,"landscape7", "World_Ecophysiographic_Diversity_2015" + ".ImageServer")
        serviceURLBiomass = os.path.join(tempFolder,"landscape6", "World_Biomass" + ".ImageServer")

        service7="landscape7.ags"
        serviceURL7="http://landscape7.arcgis.com/arcgis/rest/services"
        service6="landscape6.ags"
        serviceURL6="http://landscape6.arcgis.com/arcgis/rest/services"

        outRaster = "in_memory"

        #Parameters needed by the getResult and processResult methods are stored in this dictionary.
        layerInfo = [
            {'name': "Elevation", 'url': serviceURLElevation, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value", "Count"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Population" , 'url': serviceURLPopulation, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value","Count"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Landform", 'url': serviceURLLandform, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List': ["Value", "ClassName"], 'service': service7, 'serviceURL':serviceURL7, 'cellsize':cellsize, 'processingTemplate':"Ecophysiographic_Facet_Landform_Classes.rft"},
            {'name': "Lithology", 'url': serviceURLLithology, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value", "EF_Litho"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Bioclimate", 'url': serviceURLBioclimate, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["BioClim","Bioclimate"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Landcover", 'url': serviceURLLandcover, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["ELU_ID","ClassName"],'service': service7, 'serviceURL':serviceURL7, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Slope" , 'url': serviceURLSlope, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value","Count"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Water" , 'url': serviceURLWater, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value","Count"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':30,  'processingTemplate':'#'},
            {'name': "Diversity" , 'url': serviceURLDiversity, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value","ecoPhysdiv"],'service': service7, 'serviceURL':serviceURL7, 'cellsize':cellsize, 'processingTemplate':'#'},
            {'name': "Biomass" , 'url': serviceURLBiomass, 'scratchFolder': tempFolder, 'extentlayer': extentlayer_feature, 'List':["Value","Count"],'service': service6, 'serviceURL':serviceURL6, 'cellsize':cellsize, 'processingTemplate':'#'}
            ]

#------------------------------------------------------------------------
#
# Use parallel processing with the getResult function to retrieve images
# for analysis and cut them to the shape of the study area.
#
#------------------------------------------------------------------------

        currentTime = time.clock()
        print(str(currentTime-startTime)+" seconds. Starting multiprocessing...")

        t0 = time.time()
        arcpy.env.outputCoordinateSystem = sr
        lock = mp.Lock()
        pool = mp.Pool(5, initializer=initPool, initargs=(lock,))
        results = [pool.apply_async(getResult, args=(a,)) for a in layerInfo]
        L = []
        if arcversion == '10.5.1':
            for z in results:
                rasterInfo = z.get()
                if rasterInfo is None:
                    continue
                L = L + [rasterInfo]
            for r in L:
                print("fetching " + r['name'] + " from server")
                print(r['path'])
                arcpy.env.outputCoordinateSystem = sr
                outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), r['path'])
                p ( arcpy.GetMessages())
                od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                #Replace nodata from population estimate raster with zero. That way the sums
                #still work. Nodata will will cause GeoDescriber() to fail.
                if r['name'] == "Population":
                    if not arcpy.sa.Raster(inmem+"\\Population_R").maximum > 0:
                        outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), 0)
                        od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                if r['name'] == "Biomass":
                    if not arcpy.sa.Raster(inmem+"\\Biomass_R").maximum > 0:
                        p ( arcpy.GetMessages())
                        saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                        ##arcpy.CopyRaster_management(inmem+"\\"+ r['name'] +"_R", r"C:\gis\GeoDescriber\current.gdb"+"\\"+ r['name'] +"_R")
                    else:
                        p ( arcpy.GetMessages())
                        saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                        JoinField_Workaround(inmem+"\\"+ r['name'] +"_R", "Value", r['path'], "Value", r['List'])
                        p ( arcpy.GetMessages())
                        outputRaster = None
                else:
                    p ( arcpy.GetMessages())
                    saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                    JoinField_Workaround(inmem+"\\"+ r['name'] +"_R", "Value", r['path'], "Value", r['List'])
                    p ( arcpy.GetMessages())
                    outputRaster = None
        else:
            for z in results:
                rasterInfo = z.get()
                if rasterInfo is None:
                    continue
                L = L + [rasterInfo]
            for r in L:
                print("fetching " + r['name'] + " from server")
                arcpy.env.outputCoordinateSystem = sr
                outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), r['path'])
                p ( arcpy.GetMessages())
                od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                #Replace nodata from population estimate raster with zero. That way the sums
                #still work. Nodata will will cause GeoDescriber() to fail.
                if r['name'] == "Population":
                    if not arcpy.sa.Raster(inmem+"\\Population_R").maximum > 0:
                        outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), 0)
                        od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                if r['name'] == "Biomass":
                    if not arcpy.sa.Raster(inmem+"\\Biomass_R").maximum > 0:
                        p ( arcpy.GetMessages())
                        saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                    else:
                        p ( arcpy.GetMessages())
                        saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                        d = arcpy.JoinField_management(saved_raster, "Value", r['path'], "Value", r['List'])
                        p ( arcpy.GetMessages())
                        outputRaster = None
                else:
                    p ( arcpy.GetMessages())
                    saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                    d = arcpy.JoinField_management(saved_raster, "Value", r['path'], "Value", r['List'])
                    p ( arcpy.GetMessages())
                    outputRaster = None
        p(". Shutting down worker pool. Elapsed time: {0:.2f} seconds.".format(time.time()-t0))
        pool.close()
        pool.join()
        #if os.path.exists(tempFolder):
        #    shutil.rmtree(tempFolder)
        p("*** Process complete. {0} jobs in {1:.2f} seconds.".format(len(layerInfo), time.time()-t0))

#
#Check to see if all of the image datasets were retrieved from landscape6 and landscape7.
#If the receipt of all the datasets is not confirmed, copy data for the missing rasters
#from layerInfo to layerInfo2.

        try:
            for d in layerInfo:
                print("checking for " + d['name'])
                dataset = inmem +"\\"+ d['name'] +"_R"
                if not arcpy.Exists(dataset):
                    print(dataset+" dataset was not received from server. Trying again...")
                    layerInfo2.append(d)

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

#----------------------------------------------------------------------------
#
# If any of the ten rasters needed for analysis are missing, retrieve them again,
# this time one at a time. Otherwise the script will skip this step.
#
#----------------------------------------------------------------------------

        try:
            print layerInfo2
            if layerInfo2 != []:
            #a second round of parallel processing if any rasters are missing
                lock = None
                mp.freeze_support()
                p(". Scratch workspace is: {0}".format(tempFolder))
                print("Retrieving lost datasets again from the server...")
                t0 = time.time()
                arcpy.env.outputCoordinateSystem = sr
                lock = mp.Lock()
                #note that this time there is only one process in the pool at a time.
                pool = mp.Pool(1, initializer=initPool, initargs=(lock,))
                results = [pool.apply_async(getResult2, args=(a,)) for a in layerInfo2]
                L2 = []
                for z in results:
                    rasterInfo = z.get()
                    if rasterInfo is None:
                        continue
                    L2 = L2 + [rasterInfo]
                    print(L2)
                for r in L2:
                    if arcversion == '10.5.1':
                        print("fetching " + r['name'] + " from server")
                        ##arcpy.CopyRaster_management(r['path'], "C:\\gis\\GeoDescriber\\junk\\" + r['name'] + ".TIF")
                        arcpy.env.outputCoordinateSystem = sr
                        outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), r['path'])
                        p ( arcpy.GetMessages())
                        od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                        if r['name'] == "Population":
                            if not arcpy.sa.Raster(inmem+"\\Population_R").maximum > 0:
                                outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), 0)
                                od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                        if r['name'] == "Biomass":
                            if not arcpy.sa.Raster(inmem+"\\Biomass_R").maximum > 0:
                                p ( arcpy.GetMessages())
                                saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                                arcpy.BuildRasterAttributeTable_management(inmem+"\\"+ r['name'] +"_R")
                                ##arcpy.CopyRaster_management(inmem+"\\"+ r['name'] +"_R", r"C:\gis\GeoDescriber\current.gdb"+"\\"+ r['name'] +"_R")
                            else:
                                p ( arcpy.GetMessages())
                                saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                                JoinField_Workaround(inmem+"\\"+ r['name'] +"_R", "Value", r['path'], "Value", r['List'])
                                p ( arcpy.GetMessages())
                                outputRaster = None
                        else:
                            p ( arcpy.GetMessages())
                            saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                            JoinField_Workaround(inmem+"\\"+ r['name'] +"_R", "Value", r['path'], "Value", r['List'])
                            p ( arcpy.GetMessages())
                            outputRaster = None
                    else:
                        print("fetching " + r['name'] + " from server")
                        arcpy.env.outputCoordinateSystem = sr
                        outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), r['path'])
                        p ( arcpy.GetMessages())
                        od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                        if r['name'] == "Population":
                            if not arcpy.sa.Raster(inmem+"\\Population_R").maximum > 0:
                                outputRaster = Con(arcpy.Raster(inmem+"\\con_extent"), 0)
                                od = outputRaster.save(inmem +"\\"+ r['name'] +"_R")
                        if r['name'] == "Biomass":
                            if not arcpy.sa.Raster(inmem+"\\Biomass_R").maximum > 0:
                                p ( arcpy.GetMessages())
                                saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                            else:
                                p ( arcpy.GetMessages())
                                saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                                d = arcpy.JoinField_management(saved_raster, "Value", r['path'], "Value", r['List'])
                                p ( arcpy.GetMessages())
                                outputRaster = None
                        else:
                            p ( arcpy.GetMessages())
                            saved_raster=arcpy.Raster(inmem+"\\"+ r['name'] +"_R")
                            d = arcpy.JoinField_management(saved_raster, "Value", r['path'], "Value", r['List'])
                            p ( arcpy.GetMessages())
                            outputRaster = None
                p(". Shutting down worker pool. Elapsed time: {0:.2f} seconds.".format(time.time()-t0))
                pool.close()
                pool.join()
                if os.path.exists(tempFolder):
                    shutil.rmtree(tempFolder)
                p("*** Round 2 process complete. {0} jobs in {1:.2f} seconds.".format(len(layerInfo2), time.time()-t0))

        except arcpy.ExecuteError:
            print()
            print("Did you provide a valid username and password?")
            print("Do you have enough RAM for a polygon this large?")
            print("Also, shapefiles are not yet supported by this script...")
            print("Be sure your extent layer is a polygon feature class.")
            print()
            raise
        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            #Delete the temporary folder used to store TIF rasters retrieved from the server.
            if os.path.exists(tempFolder):
                shutil.rmtree(tempFolder)

            # Derive the aspect of the study area terrain. Reclass the aspect into 8 directions. Later the script will
            # add 180 degrees together facing all 8 directions to find a general aspect trend.
            webmerc = arcpy.SpatialReference(3857)
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Projecting elevation to webmerc to derive aspect...")
            arcpy.ProjectRaster_management(inmem+"\\Elevation_R",inmem+"\\northupg",webmerc,"NEAREST",cellsize)
            arcpy.BuildRasterAttributeTable_management(inmem+"\\northupg", "Overwrite")
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Deriving aspect from webmerc elevation raster...")
            #xbz = arcpy.sa.Aspect(inmem+"\\northupg")
            xbz = Aspect(inmem+"\\northupg")
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Projecting aspect raster back to Mollweide...")
            arcpy.ProjectRaster_management(xbz,inmem+"\\aspectg",sr,"NEAREST",cellsize)
            #arcpy.BuildRasterAttributeTable_management(inmem+"\\aspectg", "Overwrite")
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Generating remap range for aspect raster...")
            aspectRemapRange = RemapRange([[0,45,1],[45,90,2],[90,135,3],[135,180,4],[180,225,5],[225,270,6],[270,315,7],[315,360,8]])
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Reclassifying aspect raster into aspect index raster...")
            aspectindexg = arcpy.sa.Reclassify(arcpy.Raster(inmem+"\\aspectg"),"Value",aspectRemapRange)
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds. Saving and cleaning up after aspect job...")
            aspectindexg.save(inmem+"\\aspectindexg")
            #clean up and free up memory
            cleanupg = ['northupg','aspectg']
            for fd in cleanupg:
                feature = os.path.join(inmem,fd)
                if arcpy.Exists(feature):
                    arcpy.Delete_management(feature)
                    print("deleting "+feature+"...")

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

    #-----------------------------------------------------------------------------------------
    #-----------------------------------------------------------------------------------------

        try:
            # Calculate the percentage of every class for the four main ecophysiographic
            # criteria in the study area, bioclimate, landform, lithology, and land cover.
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Calculating percentages...")
            list_FC=[inmem+"\\Bioclimate_R", inmem+"\\Landform_R", inmem+"\\Lithology_R", inmem+"\\Landcover_R"]
            if arcversion == '10.5.1':
                for fc in list_FC:
                    if fc.endswith("Bioclimate_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName":
                                arcpy.DeleteField_management(fc,"ClassName")
                        arcpy.AddField_management(inmem+"\\Bioclimate_R","ClassName", "TEXT", "", "", "250")
                        with arcpy.da.UpdateCursor(inmem+"\\Bioclimate_R", ['Bioclimate','ClassName']) as cursor:
                            for row in cursor:
                                row[1] = row[0]
                                cursor.updateRow(row)
                        lookmeupl = arcpy.sa.Lookup(inmem+"\\Bioclimate_R","ClassName")
                        lookmeupl.save(inmem+"\\Bioclimate_R")
                        percent(inmem+"\\Bioclimate_R", inmem+"\\Bioclimates_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Bioclimate_R",inmem+"\\Bioclimates_CR")
                    elif fc.endswith("Landform_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName_1":
                                for field in fields:
                                    if field.name == "ClassName":
                                        arcpy.DeleteField_management(fc,"ClassName")
                                arcpy.AddField_management(inmem+"\\Landform_R","ClassName", "TEXT", "", "", "250")
                                with arcpy.da.UpdateCursor(inmem+"\\Landform_R", ['ClassName_1','ClassName']) as cursor:
                                    for row in cursor:
                                        row[1] = row[0]
                                        cursor.updateRow(row)
                                fields = arcpy.ListFields(fc)
                                for field in fields:
                                    if field.name == "ClassName_1":
                                        arcpy.DeleteField_management(fc,"ClassName_1")
                        lookmeup2 = arcpy.sa.Lookup(inmem+"\\Landform_R","ClassName")
                        lookmeup2.save(inmem+"\\Landform_R")
                        percent(inmem+"\\Landform_R",inmem+"\\Landform_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Landform_R",inmem+"\\Landform_CR")
                    elif fc.endswith("Lithology_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName":
                                arcpy.DeleteField_management(fc,"ClassName")
                        arcpy.AddField_management(inmem+"\\Lithology_R","ClassName", "TEXT", "", "", "250" )
                        with arcpy.da.UpdateCursor(inmem+"\\Lithology_R", ("EF_Litho","ClassName")) as cursor:
                            for row in cursor:
                                row[1]=row[0]
                                cursor.updateRow(row)
                        lookmeup3 = arcpy.sa.Lookup(inmem+"\\Lithology_R","ClassName")
                        lookmeup3.save(inmem+"\\Lithology_R")
                        percent(inmem+"\\Lithology_R",inmem+"\\Lithology_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Lithology_R",inmem+"\\Lithology_CR")
                    elif fc.endswith("Landcover_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName_1":
                                for field in fields:
                                    if field.name == "ClassName":
                                        arcpy.DeleteField_management(fc,"ClassName")
                                arcpy.AddField_management(inmem+"\\Landcover_R","ClassName", "TEXT", "", "", "250")
                                with arcpy.da.UpdateCursor(inmem+"\\Landcover_R", ['ClassName_1','ClassName']) as cursor:
                                    for row in cursor:
                                        row[1] = row[0]
                                        cursor.updateRow(row)
                                fields = arcpy.ListFields(fc)
                                for field in fields:
                                    if field.name == "ClassName_1":
                                        arcpy.DeleteField_management(fc,"ClassName_1")
                        lookmeup4 = arcpy.sa.Lookup(inmem+"\\Landcover_R","ClassName")
                        lookmeup4.save(inmem+"\\Landcover_R")
                        percent(inmem+"\\Landcover_R",inmem+"\\Landcover_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Landcover_R",inmem+"\\Landcover_CR")
                    else:
                        p("*** dictionaries starting ")
            else:
                for fc in list_FC:
                    if fc.endswith("Bioclimate_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName":
                                arcpy.DeleteField_management(fc,"ClassName")
                                print("Deleting ClassName...")
                        arcpy.AddField_management(inmem+"\\Bioclimate_R","ClassName", "TEXT", "", "", "250" )
                        with arcpy.da.UpdateCursor(inmem+"\\Bioclimate_R", ("Bioclimate","ClassName")) as cursor:
                            for row in cursor:
                                row[1]=  row[0]
                                cursor.updateRow(row)
                        lookmeupl = arcpy.sa.Lookup(inmem+"\\Bioclimate_R","ClassName")
                        lookmeupl.save(inmem+"\\Bioclimate_R")
                        percent(inmem+"\\Bioclimate_R", inmem+"\\Bioclimates_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Bioclimate_R",inmem+"\\Bioclimates_CR")
                    elif fc.endswith("Landform_R") == True:
                        lookmeupl = arcpy.sa.Lookup(inmem+"\\Landform_R","ClassName")
                        lookmeupl.save(inmem+"\\Landform_R")
                        percent(inmem+"\\Landform_R",inmem+"\\Landform_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Landform_R",inmem+"\\Landform_CR")
                    elif fc.endswith("Lithology_R") == True:
                        fields = arcpy.ListFields(fc)
                        for field in fields:
                            if field.name == "ClassName":
                                arcpy.DeleteField_management(fc,"ClassName")
                                print("Deleting ClassName...")
                        arcpy.AddField_management(inmem+"\\Lithology_R","ClassName", "TEXT", "", "", "250" )
                        with arcpy.da.UpdateCursor(inmem+"\\Lithology_R", ("EF_Litho","ClassName")) as cursor:
                            for row in cursor:
                                row[1]=  row[0]
                                cursor.updateRow(row)
                        lookmeupl = arcpy.sa.Lookup(inmem+"\\Lithology_R","ClassName")
                        lookmeupl.save(inmem+"\\Lithology_R")
                        percent(inmem+"\\Lithology_R",inmem+"\\Lithology_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Lithology_R",inmem+"\\Lithology_CR")
                    elif fc.endswith("Landcover_R") == True:
                        lookmeupl = arcpy.sa.Lookup(inmem+"\\Landcover_R","ClassName")
                        lookmeupl.save(inmem+"\\Landcover_R")
                        percent(inmem+"\\Landcover_R",inmem+"\\Landcover_ST","ClassName")
                        arcpy.CopyRows_management(inmem+"\\Landcover_R",inmem+"\\Landcover_CR")
                    else:
                        p("*** dictionaries starting ")

            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Done calculating percentages.")


        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        #--------------------------------------------------------------------------------
        #
        #                  study area facts (studyarealist)
        #
        #--------------------------------------------------------------------------------

        #assemble statistics for entire study area into a single list called studyarealist.

        try:
            #assemble statistics for entire study area into a single list
            dummyvalue = -9999
            studyarealist = []

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)
            #What's this area named?

        try:
            #search for a name field. If there isn't one, the script will use
            #the generic term "study area"
            #studyarealist[0] = name
            polyname = "study area"
            field_names = [f.name for f in arcpy.ListFields(inFeatLyr)]
            for field in reversed(field_names):
                fieldsupper = field.upper()
                if "NAM" in fieldsupper:
                    fields = []
                    fields.append(field)
                    fields.append("OBJECTID")
                    with arcpy.da.SearchCursor(inFeatLyr,fields) as cursor:
                        for row in cursor:
                            if row[1] == intPolyID:
                                polyname = row[0]
                elif "NOM" in fieldsupper:
                    fields = []
                    fields.append(field)
                    fields.append("OBJECTID")
                    with arcpy.da.SearchCursor(inFeatLyr,fields) as cursor:
                        for row in cursor:
                            if row[1] == intPolyID:
                                polyname = row[0]
                elif "NAAM" in fieldsupper:
                    fields = []
                    fields.append(field)
                    fields.append("OBJECTID")
                    with arcpy.da.SearchCursor(inFeatLyr,fields) as cursor:
                        for row in cursor:
                            if row[1] == intPolyID:
                                polyname = row[0]
            studyarealist.append(polyname)

        except IOError:
            print("An error has occurred. Please provide a valid pathname to your extent layer and try again.")
            raise

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            #What's this area's short name? studyarealist[1]
            #What's this area's short name singular version? studyarealist[2]
            #What's this area's short name plural version? studyarealist[3]
            shortname = "study area"
            shortname_sing = "study area is"
            shortname_plur = "study area are"
            for row in shortname_matrix:
                if row[0] in studyarealist[0].lower():
                    shortname = row[0]
                    shortname_sing = row[1]
                    shortname_plur = row[2]
            studyarealist.append(shortname)
            studyarealist.append(shortname_sing)
            studyarealist.append(shortname_plur)
            del shortname
            del shortname_sing
            del shortname_plur

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            #Generate a grid with north, south, east, and west sides, coded, to evaluate the study area.
            #First change the projection to cylindrical equal area so directions are true.
            cea = arcpy.SpatialReference(54034)
            arcpy.env.outputCoordinateSystem = cea
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Generating N/S/E/W. Spatial reference is cylindrical equal area.")
            arcpy.env.overwriteOutput=True
            arcpy.env.cellSize = cellsize
            arcpy.env.extent = "MAXOF"
            qstr = "OBJECTID" + " = " + str(intPolyID)
            epfcea0 = arcpy.Project_management(epf,output+"\\epfcea",cea)
            epf2 = arcpy.MakeFeatureLayer_management(epfcea0,"EachPolyFeat", qstr)

            #Obtain the envelope for the polygon, then make its limits the limits of the
            #grid which divides the shape into north, south, east, west, northeast, northwest,
            #southeast, southwest, and center. Create points at the center of each of the nine zones,
            #then generate thiessen polygons from the points to generate the zones.
            fieldsdk = ["SHAPE@"]
            with arcpy.da.SearchCursor(epf2, fieldsdk) as cursordk:
                for rowdk in cursordk:
                    #if rowdk[0] == intPolyID:
                    #geom = rowdk[1]
                    geom = rowdk[0]
                    ext = geom.extent  # or row.Shape.extent
                    ceaxmin = ext.XMin
                    ceaymin = ext.YMin
                    ceaxmax = ext.XMax
                    ceaymax = ext.YMax

            ceay1 = ((ceaymax - ceaymin) * 0.4) + ceaymin
            ceay2 = ((ceaymax - ceaymin) * 0.6) + ceaymin
            ceax1 = ((ceaxmax - ceaxmin) * 0.4) + ceaxmin
            ceax2 = ((ceaxmax - ceaxmin) * 0.6) + ceaxmin

            point11 = arcpy.Point(((ceax1-ceaxmin)/2)+ceaxmin,((ceay1-ceaymin)/2)+ceaymin)
            point12 = arcpy.Point(((ceax1-ceaxmin)/2)+ceaxmin,((ceay2-ceay1)/2)+ceay1)
            point13 = arcpy.Point(((ceax1-ceaxmin)/2)+ceaxmin,((ceaymax-ceay2)/2)+ceay2)
            point21 = arcpy.Point(((ceax2-ceax1)/2)+ceax1,((ceay1-ceaymin)/2)+ceaymin)
            point22 = arcpy.Point(((ceax2-ceax1)/2)+ceax1,((ceay2-ceay1)/2)+ceay1)
            point23 = arcpy.Point(((ceax2-ceax1)/2)+ceax1,((ceaymax-ceay2)/2)+ceay2)
            point31 = arcpy.Point(((ceaxmax-ceax2)/2)+ceax2,((ceay1-ceaymin)/2)+ceaymin)
            point32 = arcpy.Point(((ceaxmax-ceax2)/2)+ceax2,((ceay2-ceay1)/2)+ceay1)
            point33 = arcpy.Point(((ceaxmax-ceax2)/2)+ceax2,((ceaymax-ceay2)/2)+ceay2)

            point10x = ((ceax1-ceaxmin)/2)+ceaxmin
            point20x = ((ceax2-ceax1)/2)+ceax1
            point30x = ((ceaxmax-ceax2)/2)+ceax2
            point1y = ((ceay1-ceaymin)/2)+ceaymin
            point2y = ((ceay2-ceay1)/2)+ceay1
            point3y = ((ceaymax-ceay2)/2)+ceay2

            row_values = [(11,(point10x,point1y)),(12,(point10x,point2y)),(13,(point10x,point3y)),(21,(point20x,point1y)),(22,(point20x,point2y)),(23,(point20x,point3y)),(31,(point30x,point1y)),(32,(point30x,point2y)),(33,(point30x,point3y))]

            feature_class = "temppts"
            arcpy.CreateFeatureclass_management(inmem,feature_class,"POINT")
            arcpy.AddField_management(inmem+"\\"+feature_class,"ZONE","SHORT")
            cursor = arcpy.da.InsertCursor(inmem+"\\"+feature_class,["ZONE","SHAPE@XY"])
            for row in row_values:
                cursor.insertRow(row)
            del cursor
            arcpy.env.extent = ext

            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Creating Thiessen Polygons for N/S/E/W")
            arcpy.CreateThiessenPolygons_analysis(inmem+"\\"+feature_class, inmem+"\\thiespts", "ALL")

            #return spatial reference to mollweide for the most accurate calculations (equal area projection)
            sr = arcpy.SpatialReference(54009)
            arcpy.env.outputCoordinateSystem = sr
            arcpy.env.cellSize = cellsize
            arcpy.env.extent = "MAXOF"

            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Projecting and rasterizing Thiessen polygons")

            #project thiessen polygons into mollweide then create a raster of the thiessen zones.
            epfcea = arcpy.Project_management(inmem+"\\thiespts",output+"\\thiesmw",sr)
            zonedg = arcpy.PolygonToRaster_conversion(output+"\\thiesmw","ZONE")
            mw_zonedg = arcpy.sa.Con(rasterExt,zonedg)
            mw_zonedg.save(inmem+"\\mw_zonedg")
            feature = os.path.join(output,"thiesmw")
            if arcpy.Exists(feature):
                arcpy.Delete_management(feature)
                print("deleting "+feature+"...")

            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("spatial reference is back to MW")
            arcpy.env.overwriteOutput=True

            #zero out the variables that keep a tally of the spatial position
            #of geographic phenomena
            zonesw= 0
            zonew = 0
            zonenw = 0
            zones = 0
            zonec = 0
            zonen = 0
            zonese = 0
            zonee = 0
            zonene = 0
            neq_n = -9999
            seq_n = -9999
            swq_n = -9999
            nwq_n = -9999
            zes_n = -9999
            zws_n = -9999
            zss_n = -9999
            zns_n = -9999
            zonesw_n = -9999
            zonew_n = -9999
            zonenw_n = -9999
            zones_n = -9999
            zonec_n = -9999
            zonen_n = -9999
            zonese_n = -9999
            zonee_n = -9999
            zonene_n = -9999

            #count the number of cells in each of the nine zones in the study area.
            #Transfer them to nine variables to be used later for comparison with
            #counts of parts of the study area.
            print("populating variables from the nine parts of the study area")
            fieldsmwz = ["Value","Count"]
            with arcpy.da.SearchCursor(inmem+"\\mw_zonedg",fieldsmwz) as cursormwz:
                for rowmwz in cursormwz:
                    if rowmwz[0] == 11:
                        zonesw = rowmwz[1]
                    if rowmwz[0] == 12:
                        zonew = rowmwz[1]
                    if rowmwz[0] == 13:
                        zonenw = rowmwz[1]
                    if rowmwz[0] == 21:
                        zones = rowmwz[1]
                    if rowmwz[0] == 22:
                        zonec = rowmwz[1]
                    if rowmwz[0] == 23:
                        zonen = rowmwz[1]
                    if rowmwz[0] == 31:
                        zonese = rowmwz[1]
                    if rowmwz[0] == 32:
                        zonee = rowmwz[1]
                    if rowmwz[0] == 33:
                        zonene = rowmwz[1]

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            print("appending elevation statistics to the study area...")
            #Obtain the mean and standard deviation of the elevation for the whole study area,
            #and store the values in meanelevation and stdelevation. Used later in comparison.
            a=arcpy.sa.ZonalStatisticsAsTable(extentRaster,"Value",inmem+"\\Elevation_R",inmem+"\\stattbl","DATA")
            with arcpy.da.SearchCursor(inmem+"\\stattbl",["COUNT","MIN","MAX","MEAN","MEDIAN","STD"]) as cursor:
                for row in cursor:
                    #What is its area? studyarealist[4]
                    studyarealist.append(row[0])
                    #What is its lowest elevation? studyarealist[5]
                    studyarealist.append(row[1])
                    #What is its highest elevation? studyarealist[6]
                    studyarealist.append(row[2])
                    #What is its mean elevation? studyarealist[7]
                    studyarealist.append(row[3])
                    #What is its median elevation? studyarealist[8]
                    studyarealist.append(row[4])
                    #What is one standard deviation below the mean elevation? studyarealist[9]
                    studyarealist.append(row[3]-row[5])
                    #What is one standard deviation above the mean elevation? studyarealist[10]
                    studyarealist.append(row[3]+row[5])
            print("done appending elevation statistics to the study area...")

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

            #If bodies of water make up a significant class, use the 30m resolution raster
            #to estimate its area with more accuracy. This improves estimates in areas where
            #there are a lot of small lakes, such as Minnesota

        try:
            print("getting statistics on bodies of water...")
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            divzero = 0
            bowno = -9999
            bowyes = -9999

            xyy=Con(arcpy.Raster(inmem+"\\Water_R") == 11, 1, 0)
            xyz=Con(rasterExt,xyy)
            xyz.save(inmem+"\\Water30m_R")
            arcpy.BuildRasterAttributeTable_management(inmem+"\\Water30m_R", "Overwrite")
            with arcpy.da.SearchCursor(inmem+"\\Water30m_R",["Value","Count"]) as cursor:
                for row in cursor:
                    if row[0] == 0 and row[1] == 0:
                        divzero == 1
            if divzero == 0:
                with arcpy.da.SearchCursor(inmem+"\\Water30m_R",["Value","Count"]) as cursor:
                    for row in cursor:
                        if row[0] == 0:
                            bowno = row[1]
                        if row[0] == 1:
                            bowyes = row[1]
                    #What percentage of the study area is bodies of water, at 30m resolution? studyarealist[11]
                studyarealist.append((bowyes/bowno) * 100)
            else:
                studyarealist.append(-9999)
            if arcpy.Exists(xyy):
            	arcpy.Delete_management(xyy)
            if arcpy.Exists(xyz):
            	arcpy.Delete_management(xyz)

##            arcpy.env.cellSize = 30
##
##            print("studyarealist[4] (cell count) is "+str(studyarealist[4]))
##
##            # If a raster has 100 cells or less, you can calculate the area of the
##            # bodies of water by making datasets in memory. Otherwise you may run out of memory.
##            # In memory datasets save a lot of time, however. Experiment with this.
##            # if you have enough RAM you can safely set this threshold higher than 10,000 cells.
##            if studyarealist[4] > 100:
##                # find the percentage bodies of water, and write temporary datasets to disk.
##                # This is slower, but with less risk of running out of memory.
##                #xyy=Con(Con(arcpy.PolygonToRaster_conversion(output+"\\proj", thefieldname, output+"\\mask30"), inmem+"\\Water_R") == 11, 1, 0)
##                print("making a 30m mask for bodies of water...")
##                xyw=arcpy.PolygonToRaster_conversion(output+"\\proj", thefieldname, output+"\\mask30")
##                print("isolating bodies of water...")
##                xyx=Con(xyw, inmem+"\\Water_R")
##                print("converting bodies of water to binary raster...")
##                xyy=Con(xyx == 11, 1, 0)
##                xyy.save(output+"\\Water30m_R")
##                feature = os.path.join(output,'mask30')
##                if arcpy.Exists(feature):
##            	   arcpy.Delete_management(feature)
##                arcpy.BuildRasterAttributeTable_management(output+"\\Water30m_R", "Overwrite")
##                with arcpy.da.SearchCursor(output+"\\Water30m_R",["Value","Count"]) as cursor:
##                    for row in cursor:
##                        if row[0] == 0 and row[1] == 0:
##                            divzero == 1
##                if divzero == 0:
##                    with arcpy.da.SearchCursor(output+"\\Water30m_R",["Value","Count"]) as cursor:
##                        for row in cursor:
##                            if row[0] == 0:
##                                bowno = row[1]
##                            if row[0] == 1:
##                                bowyes = row[1]
##                        #What percentage of the study area is bodies of water, at 30m resolution? studyarealist[11]
##                    studyarealist.append((bowyes/bowno) * 100)
##                else:
##                    studyarealist.append(-9999)
##                arcpy.env.cellSize = cellsize
##                if arcpy.Exists(xyy):
##                	arcpy.Delete_management(xyy)
##                if arcpy.Exists(xyx):
##                	arcpy.Delete_management(xyx)
##                if arcpy.Exists(xyw):
##                	arcpy.Delete_management(xyx)
##                feature = os.path.join(output,'Water30m_R')
##                if arcpy.Exists(feature):
##                	arcpy.Delete_management(feature)
##            else:
##                # find the percentage bodies of water, but in_memory (fastest way)
##                xyy=Con(Con(arcpy.PolygonToRaster_conversion(output+"\\proj", thefieldname, inmem+"\\mask30"), inmem+"\\Water_R") == 11, 1, 0)
##                xyy.save(inmem+"\\Water30m_R")
##                feature = os.path.join(inmem,'mask30')
##                if arcpy.Exists(feature):
##            	   arcpy.Delete_management(feature)
##                arcpy.BuildRasterAttributeTable_management(inmem+"\\Water30m_R", "Overwrite")
##                with arcpy.da.SearchCursor(inmem+"\\Water30m_R",["Value","Count"]) as cursor:
##                    for row in cursor:
##                        if row[0] == 0 and row[1] == 0:
##                            divzero == 1
##                if divzero == 0:
##                    with arcpy.da.SearchCursor(inmem+"\\Water30m_R",["Value","Count"]) as cursor:
##                        for row in cursor:
##                            if row[0] == 0:
##                                bowno = row[1]
##                            if row[0] == 1:
##                                bowyes = row[1]
##                        #What percentage of the study area is bodies of water, at 30m resolution? studyarealist[11]
##                    studyarealist.append((bowyes/bowno) * 100)
##                else:
##                    studyarealist.append(-9999)
##                arcpy.env.cellSize = cellsize
##                feature = os.path.join(inmem,'Water30m_R')
##                if arcpy.Exists(feature):
##                	arcpy.Delete_management(feature)
##                arcpy.env.cellSize = cellsize
##                feature = os.path.join(output,'Water30m_R')
##                if arcpy.Exists(feature):
##                	arcpy.Delete_management(feature)
##                feature = xyy
##                if arcpy.Exists(feature):
##                	arcpy.Delete_management(feature)

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            currentTime = time.clock()
            print(str(currentTime-startTime)+" seconds have elapsed so far")
            print("Calculating mean ELU diversity")
            print("looking up diversity")
            ecophysdivg=arcpy.sa.Lookup(inmem+"\\Diversity_R","ecoPhysdiv")
            print("saving diversity")
            ecophysdivg.save(inmem+"\\Ecophysdiv_R")
            print("populating list with the diversity")
            a=arcpy.sa.ZonalStatisticsAsTable(extentRaster,"Value",inmem+"\\Ecophysdiv_R",inmem+"\\stattbl","DATA")
            with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN","STD"]) as cursor:
                for row in cursor:
                    #What is the mean ELU diversity in the study area? studyarealist[12]
                    studyarealist.append(row[0])
                    #What is the median ELU diversity in the study area? studyarealist[13]
                    #median is not available from floating point raster, using mean as dummy value!
                    studyarealist.append(row[0])
                    #What is one standard deviation below the mean ELU diversity? studyarealist[14]
                    studyarealist.append(row[0]-row[1])
                    #What is one standard deviation above the mean ELU diversity? studyarealist[15]
                    studyarealist.append(row[0]+row[1])
                    #dummy value studyarealist[16]
                    studyarealist.append(dummyvalue)
                    #dummy value studyarealist[17]
                    studyarealist.append(dummyvalue)
                    #dummy value studyarealist[18]
                    studyarealist.append(dummyvalue)
                    #dummy value studyarealist[19]
                    studyarealist.append(dummyvalue)
                    #dummy value studyarealist[20]
                    studyarealist.append(dummyvalue)
                    #dummy value studyarealist[21]
                    studyarealist.append(dummyvalue)

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        try:
            #number of cells in the southwest. studyarealist[22]
            studyarealist.append(zonesw)
            #number of cells in the west. studyarealist[23]
            studyarealist.append(zonew)
            #number of cells in the northwest. studyarealist[24]
            studyarealist.append(zonenw)
            #number of cells in the south. studyarealist[25]
            studyarealist.append(zones)
            #number of cells in the center. studyarealist[26]
            studyarealist.append(zonec)
            #number of cells in the north. studyarealist[27]
            studyarealist.append(zonen)
            #number of cells in the southeast. studyarealist[28]
            studyarealist.append(zonese)
            #number of cells in the east. studyarealist[29]
            studyarealist.append(zonee)
            #number of cells in the northeast. studyarealist[30]
            studyarealist.append(zonene)

            neq_n = studyarealist[26]+studyarealist[27]+studyarealist[29]+studyarealist[30]
            seq_n = studyarealist[26]+studyarealist[25]+studyarealist[29]+studyarealist[28]
            swq_n = studyarealist[22]+studyarealist[23]+studyarealist[25]+studyarealist[26]
            nwq_n = studyarealist[23]+studyarealist[24]+studyarealist[26]+studyarealist[27]
            zes_n = studyarealist[28]+studyarealist[29]+studyarealist[30]
            zws_n = studyarealist[22]+studyarealist[23]+studyarealist[24]
            zss_n = studyarealist[22]+studyarealist[25]+studyarealist[28]
            zns_n = studyarealist[24]+studyarealist[27]+studyarealist[30]
            zonesw_n = studyarealist[22]
            zonew_n = studyarealist[23]
            zonenw_n = studyarealist[24]
            zones_n = studyarealist[25]
            zonec_n = studyarealist[26]
            zonen_n = studyarealist[27]
            zonese_n = studyarealist[28]
            zonee_n = studyarealist[29]
            zonene_n = studyarealist[30]

            article = ""
            if studyarealist[0][0].isupper():
                article = ""
            else:
                article = "the "
            print(article)

            print(studyarealist)
            #print("The study area list is "+str(len(studyarealist))+" tokens in length. (should be 31)")

            if len(studyarealist) != 31:
                GeoDescriberTries +=1
                return

        except ValueError:
            print("value error")

        except IndexError:
            print("index error")
            neq_n = 0
            seq_n = 0
            swq_n = 0
            nwq_n = 0
            zes_n = 0
            zws_n = 0
            zss_n = 0
            zns_n = 0
            zonesw_n = 0
            zonew_n = 0
            zonenw_n = 0
            zones_n = 0
            zonec_n = 0
            zonen_n = 0
            zonese_n = 0
            zonee_n = 0
            zonene_n = 0

        except:
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)

        #---------------------------------------------------------------------------------
        #
        #                        look up significant classes
        #
        #
        #--------------------------------------------------------------------------------
        if studyarealist[4] > 15:
            try:
                global i
                global restopslist
                p("*** look up significant classes...")
                sideopslist = []
                restopslist = []
                alllist = []

                #Assemble characterization text.
                currentTime = time.clock()
                print(str(currentTime-startTime)+" seconds have elapsed so far")
                print("Assembling characterization text...")
                AllPoly=[inmem+"\\Bioclimate_R", inmem+"\\Landform_R", inmem+"\\Lithology_R", inmem+"\\Landcover_R"]
                for i in AllPoly:
                    if i.endswith("Bioclimate_R") == True:
                        Bio = i
                        #Bio
                        largVal_bio=largest(Bio,"percent")
                        key_bioclimate= findString(Bio, largVal_bio)
                        bioclimate_string = bioclimate_dict[key_bioclimate]
                        #find other values bigger than 10% smaller than the largest value
                        bioclimate_rest= restofValues(Bio, largVal_bio, bioclimate_dict)
                    elif i.endswith("Landform_R") == True:
                        Landform=i
                        #Landform
                        largVal_landform=largest(Landform,"percent")
                        key_landform=findString(Landform,largVal_landform)
                        landform_string = landform_dict[key_landform]
                        #find other values bigger than 10% smaller than the largest value
                        landform_rest= restofValues(Landform, largVal_landform, landform_dict)
                    elif i.endswith("Lithology_R") == True:
                        Lithology=i
                        #Lithology
                        largVal_lit=largest(Lithology,"percent")
                        ##tr_lit=int(round((largVal_lit/10)-.5))
                        key_lithology=findString(Lithology,largVal_lit)
                        lithology_string = lithology_dict[key_lithology]
                        #find rest of the values bigger than 10% smaller than the biggest value
                        lithology_rest= restofValues(Lithology, largVal_lit, lithology_dict)
                    elif i.endswith("Landcover_R") == True:
                        Landcover= i
                        largVal_landcover=largest(Landcover,"percent")
                        key_landcover=findString(Landcover, largVal_landcover)
                        landcover_string = landcover_dict[key_landcover]
                        #find rest of the values bigger than 10%  smaller than the biggest value
                        landcover_rest= restofValues(Landcover, largVal_landcover, landcover_dict)

                        polygonfc=output+"\\proj"

                        #The area of the polygonfc determines some text, if it's called a place, area, or region.
                        largVal_polygonfc=largest(polygonfc,"SHAPE_Area")

                        if largVal_polygonfc < warnthreshold:
                            warnings.warn("Warning: "+polyname+" is less than 1000 pixels. This is less than the minimum designed study area size.", UserWarning, stacklevel=2)

            except IOError:
                print("This does not appear to be a valid polygon feature class. A valid")
                print("feature class contains at least some land area and can not be")
                print("completely over ocean.")
                raise

            except:
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)


            #----------------------------------side operations--------------------------------
            #
            # Order a list so the four predominant ecophysiographic criteria appear in order from highest
            # to lowest by percentage. If lithology appears first, demote it to second place.
            # Follow up listing each of the other criteria that are over 10% with the predominant criteria,
            # so all the lithology are together, all the bioclimate, landform and land cover are together
            # in their own paragraphs from predominant value, in descending order to the least significant value.
            #
            #--------------------------------------------------------------------------------

        p("*** side operations..")
        dummyvalue = -9999
        alllist = []
        descendingorderlist = []
        descendingrestlist = []
        junklist = []
        connectors = [['In addition, '],['Furthermore, '],['Also, '],['Moreover, ']]

        if studyarealist[4] > 15:
            try:
                while len(sideopslist) != 0:
                    sideoppct = sideopslist.pop(0)
                    sideopclass = sideopslist.pop(0)
                    sideoppath = sideopslist.pop(0)
                    sideoppctrnd = int(sideoppct + .5)
                    junklist.append(sideoppct)
                    junklist.append(sideopclass)
                    junklist.append(sideoppath)
                    descendingorderlist.append(junklist)
                    junklist = []
                descendingorderlist.sort(reverse=True)
                x = 0
                # move lithology paragraph so that it is the last paragraph.
                for triplet in descendingorderlist:
                    if triplet[2].endswith("Lithology_R") == True and x <= 2:
                        lithofirst = descendingorderlist.pop(x)
                        descendingorderlist.insert(3,lithofirst)
                    x += 1
                x = 0
                # now if landcover appears first or second, move it to third place. The characterization will now lead with either landform or bioclimate.
                for triplet in descendingorderlist:
                    if triplet[2].endswith("Landcover_R") == True and x <= 1:
                        lcfirst = descendingorderlist.pop(x)
                        descendingorderlist.insert(2,lcfirst)
                    x += 1
                for triplet in descendingorderlist:
                    if triplet[2].endswith("Bioclimate_R") == True:
                        topbioclimatepercent = triplet[0]
                    if triplet[2].endswith("Landform_R") == True:
                        toplandformpercent = triplet[0]
                    if triplet[2].endswith("Landcover_R") == True:
                        toplandcoverpercent = triplet[0]
                    if triplet[2].endswith("Lithology_R") == True:
                        toplithopercent = triplet[0]
                while len(restopslist) != 0:
                    sideoppct = restopslist.pop(0)
                    sideopclass = restopslist.pop(0)
                    sideoppath = restopslist.pop(0)
                    sideoppctrnd = int(sideoppct + .5)
                    junklist.append(sideoppct)
                    junklist.append(sideopclass)
                    junklist.append(sideoppath)
                    descendingrestlist.append(junklist)
                    junklist = []
                descendingrestlist.sort(reverse=True)
                for triplet in descendingorderlist:
                    alllist.append(triplet)
                    for resttriplet in descendingrestlist:
                        if triplet[2] == resttriplet[2]:
                            alllist.append(resttriplet)

            except:
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)

            #---------------------------------------------------------------------------------
            #
            #                            significant class facts
            #
            #  This section analyzes every lithology, landcover, landform, and bioclimate that
            #  makes up over 10% of the study area and populates classlist with the analysis facts.
            #
            #---------------------------------------------------------------------------------

        if studyarealist[4] > 15:
            try:
                #p("Analyzing significant ecophysiographic phenomena in detail...")
                currentTime = time.clock()
                print(str(currentTime-startTime)+" seconds have elapsed so far")
                print("Analyzing significant ecophysiographic phenomena in detail...")

                # clear counters which count off how many times a loop has been executed on each ecophysiographic class.
                landcovertimes = 0
                lithologytimes = 0
                bioclimatetimes = 0
                landformtimes = 0
                classlist = []
                count0 = 0
                classname0 = []

                # clear variables which store the text that is used to create the output text fields.
                # description = this variable collects the text for the bulleted fact sheet, with percentages written out.
                # description = the text for the four paragraph characterization, using adjectives to approximate quantity.
                preface = ""
                description = ""
                description = ""


                # Now that alllist is ordered in the order that sentences will be written, perform some analysis
                # to find out what to say about each class. The analysis results will populate the classlist list object.
                # bioclimatelist,
                for analclass in alllist:
                    thepct = analclass[0]
                    theclass = analclass[1]

                    # BIOCLIMATES side operations
                    # Extract the bioclimate class from the bioclimates raster. Output the minimum, mean, and maximum elevation
                    # for that bioclimates raster into variables. Use those variables in a sentence which gives the minimum and
                    # maximum elevation as context for that bioclimate. If the mean of the bioclimate is greater than the mean
                    # plus the standard deviation of the con_extent, add "at the higher elevations" to the sentence. If the mean
                    # of the bioclimate is less than the mean minus the standard deviation of con_extent, add "at the lower
                    # elevations" to the sentence.
                    #
                    # If any bioclimate classes are over 10% of the study area, find out more details about those classes.
                    #if analclass[2] == inmem+"\\Bioclimate_R":
                    if analclass[2].endswith("Bioclimate_R") == True:
                        bioclimatelist = []
                        #What is the pathname to this bioclimate subclass? bioclimatelist[0]
                        bioclimatelist.append(analclass[2])
                        #What's the bioclimate called? bioclimatelist[1]
                        bioclimatelist.append(analclass[1])
                        currentTime = time.clock()
                        print(str(currentTime-startTime)+" seconds have elapsed so far")
                        print("-----performing side operations on bioclimate class " + analclass[1]+"-----")
                        bioclimatetimes += 1
                        attExtract = ExtractByAttributes(analclass[2], "ClassName = '"+analclass[1]+"'")

                        feature = os.path.join(inmem,'extractg')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)

                        attExtract.save(inmem+"\\extractg")
                        # For each significant bioclimate, find out if its mean elevation is above or below
                        # one standard deviation of the study area mean elevation.
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Elevation_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["Count","MIN","MAX","MEAN","MEDIAN"]) as cursor:
                            for row in cursor:
                                #What is this bioclimate subclass' area? bioclimatelist[2]
                                bioclimatelist.append(row[0])
                                #What is this bioclimate subclass' percentage of the overall study area? bioclimatelist[3]
                                bioclimatelist.append(analclass[0])
                                #What is this bioclimate subclass' rank in order of greatest to least percent of the study area? bioclimatelist[4]
                                bioclimatelist.append(bioclimatetimes)
                                #What is this bioclimate subclass' lowest elevation? bioclimatelist[5]
                                bioclimatelist.append(row[1])
                                #What is this bioclimate subclass' highest elevation? bioclimatelist[6]
                                bioclimatelist.append(row[2])
                                #What is this bioclimate subclass' mean elevation? bioclimatelist[7]
                                bioclimatelist.append(row[3])
                                #What is this bioclimate subclass' median elevation? bioclimatelist[8]
                                bioclimatelist.append(row[4])
                        #dummy bioclimatelist[9]
                        bioclimatelist.append(dummyvalue)
                        #dummy bioclimatelist[10]
                        bioclimatelist.append(dummyvalue)
                        # Take the extracted bioclimate class from the bioclimates raster and use it as a conditional raster to
                        # output the landcover on just that bioclimate class. If a particular land cover is over half of the
                        # bioclimate class, write a sentence that says most of this bioclimate zone is covered by a particular land
                        # cover class. If a particular land cover class is over 89%, change 'over half' to 'almost all'. And if
                        # that land cover class is over 99%, change 'over half' to 'all'.
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landform_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            ##arcpy.CopyRaster_management(inmem+"\\cong", r"C:\gis\GeoDescriber\current.gdb\cong")
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top landform in the area covered by this subclass? bioclimatelist[11]
                        bioclimatelist.append(classname0)
                        #What percentage of the area covered by this subclass is the top landform? bioclimatelist[12]
                        bioclimatelist.append((count0/bioclimatelist[2])*100)
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Lithology_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top rock type in the area covered by this subclass? bioclimatelist[13]
                        bioclimatelist.append(classname0)
                        #What percentage of the area covered by this subclass is the top rock type? bioclimatelist[14]
                        bioclimatelist.append((count0/bioclimatelist[2])*100)
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landcover_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top land cover in the area covered by this subclass? bioclimatelist[15]
                        bioclimatelist.append(classname0)
                        #What percentage of the area covered by this subclass is the top land cover? bioclimatelist[16]
                        bioclimatelist.append((count0/bioclimatelist[2])*100)
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Ecophysdiv_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN"]) as cursor:
                            for row in cursor:
                                #What is the mean ELU diversity in this climate zone? bioclimatelist[17]
                                bioclimatelist.append(row[0])
                                #What is the median ELU diversity in this climate zone? bioclimatelist[18]
                                #this is a dummy value because there is no median in floating point
                                bioclimatelist.append(row[0])
                        #dummy value bioclimatelist[19]
                        bioclimatelist.append(dummyvalue)
                        #dummy value bioclimatelist[20]
                        bioclimatelist.append(dummyvalue)
                        #dummy value bioclimatelist[21]
                        bioclimatelist.append(dummyvalue)

                        neq = 0
                        seq = 0
                        swq = 0
                        nwq = 0
                        zes = 0
                        zws = 0
                        zss = 0
                        zns = 0
                        zonesw = 0
                        zonew = 0
                        zonenw = 0
                        zones = 0
                        zonec = 0
                        zonen = 0
                        zonese = 0
                        zonee = 0
                        zonene = 0

                        try:
                            feature = os.path.join(inmem,'cong')
                            if arcpy.Exists(feature):
                            	arcpy.Delete_management(feature)
                            cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\mw_zonedg")
                            cong.save(inmem+"\\cong")
                            with arcpy.da.SearchCursor(inmem+"\\cong",["Value","Count"]) as cursor:
                                for rowmwz in cursor:
                                    if rowmwz[0] == 11:
                                        zonesw = rowmwz[1]
                                    if rowmwz[0] == 12:
                                        zonew = rowmwz[1]
                                    if rowmwz[0] == 13:
                                        zonenw = rowmwz[1]
                                    if rowmwz[0] == 21:
                                        zones = rowmwz[1]
                                    if rowmwz[0] == 22:
                                        zonec = rowmwz[1]
                                    if rowmwz[0] == 23:
                                        zonen = rowmwz[1]
                                    if rowmwz[0] == 31:
                                        zonese = rowmwz[1]
                                    if rowmwz[0] == 32:
                                        zonee = rowmwz[1]
                                    if rowmwz[0] == 33:
                                        zonene = rowmwz[1]

                            del cong
                            #number of cells in the southwest. bioclimatelist[22]
                            bioclimatelist.append(zonesw)
                            #number of cells in the west. bioclimatelist[23]
                            bioclimatelist.append(zonew)
                            #number of cells in the northwest. bioclimatelist[24]
                            bioclimatelist.append(zonenw)
                            #number of cells in the south. bioclimatelist[25]
                            bioclimatelist.append(zones)
                            #number of cells in the center. bioclimatelist[26]
                            bioclimatelist.append(zonec)
                            #number of cells in the north. bioclimatelist[27]
                            bioclimatelist.append(zonen)
                            #number of cells in the southeast. bioclimatelist[28]
                            bioclimatelist.append(zonese)
                            #number of cells in the east. bioclimatelist[29]
                            bioclimatelist.append(zonee)
                            #number of cells in the northeast. bioclimatelist[30]
                            bioclimatelist.append(zonene)

                            neq = bioclimatelist[26]+bioclimatelist[27]+bioclimatelist[29]+bioclimatelist[30]
                            seq = bioclimatelist[26]+bioclimatelist[25]+bioclimatelist[29]+bioclimatelist[28]
                            swq = bioclimatelist[22]+bioclimatelist[23]+bioclimatelist[25]+bioclimatelist[26]
                            nwq = bioclimatelist[23]+bioclimatelist[24]+bioclimatelist[26]+bioclimatelist[27]
                            zes = bioclimatelist[28]+bioclimatelist[29]+bioclimatelist[30]
                            zws = bioclimatelist[22]+bioclimatelist[23]+bioclimatelist[24]
                            zss = bioclimatelist[22]+bioclimatelist[25]+bioclimatelist[28]
                            zns = bioclimatelist[24]+bioclimatelist[27]+bioclimatelist[30]

                            maxzone = max([neq,'northeastern part',neq_n],[seq,'southeastern part',seq_n],[swq,'southwestern part',swq_n],[nwq,'northwestern part',nwq_n],[zes,'east side',zes_n],[zws,'west side',zws_n],[zss,'south side',zss_n],[zns,'north side',zns_n],[zonesw,'southwesternmost portion',zonesw_n],[zonew,'westernmost portion',zonew_n],[zonenw,'northwesternmost portion',zonenw_n],[zones,'southernmost portion',zones_n],[zonec,'most central portion',zonec_n],[zonen,'northernmost portion',zonen_n],[zonese,'southeasternmost portion',zonese_n],[zonee,'easternmost portion',zonee_n],[zonene,'northeasternmost portion',zonene_n])

                            #How many cells are in the most likely part of the bioclimate class?[31]
                            bioclimatelist.append(maxzone[0])
                            #What phrase describes the most likely part of this bioclimate class?[32]
                            bioclimatelist.append(maxzone[1])
                            #What percentage of the most likely part of the bioclimate class is this bioclimate class?[33]
                            if maxzone[2] != 0:
                                bioclimatelist.append((maxzone[0]/maxzone[2])*100)
                            else:
                                bioclimatelist.append(0)

                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0

                        except ValueError:
                            print("value error")

                        except IndexError:
                            print("index error")
                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0

                        classlist.append(bioclimatelist)
                        #print bioclimatelist
                        print("length of bioclimatelist is "+str(len(bioclimatelist)))

                        if len(bioclimatelist) != 34:
                            GeoDescriberTries +=1
                            return

                    # LANDFORMS side operations
                    # Extract the landform class from the landforms raster and
                    # if any landform classes are over 10% of the study area, find out more details about those classes.
                    elif analclass[2].endswith("Landform_R") == True:
                    #elif analclass[2] == output+"\\Landform_R":
                        landformlist = []
                        #What is the pathname to this landform subclass? landformlist[0]
                        landformlist.append(analclass[2])
                        #What's the landform called? landformlist[1]
                        landformlist.append(analclass[1])
                        currentTime = time.clock()
                        print(str(currentTime-startTime)+" seconds have elapsed so far")
                        print("-----performing side operations on landform class " + analclass[1]+"-----")
                        landformtimes += 1
                        attExtract = ExtractByAttributes(analclass[2], "ClassName = '"+analclass[1]+"'")

                        feature = os.path.join(inmem,'extractg')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)

                        attExtract.save(inmem+"\\extractg")
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Elevation_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["Count","MIN","MAX","MEAN","MEDIAN"]) as cursor:
                            for row in cursor:
                                #What is this landform subclass' area? landformlist[2]
                                landformlist.append(row[0])
                                #What is this landform subclass' percentage of the overall study area? landformlist[3]
                                landformlist.append(analclass[0])
                                #What is this landform subclass' rank in order of greatest to least percent of the study area? landformlist[4]
                                landformlist.append(landformtimes)
                                #What is this landform subclass' lowest elevation? landformlist[5]
                                landformlist.append(row[1])
                                #What is this landform subclass' highest elevation? landformlist[6]
                                landformlist.append(row[2])
                                #What is this landform subclass' mean elevation? landformlist[7]
                                landformlist.append(row[3])
                                #What is this landform subclass' median elevation? landformlist[8]
                                landformlist.append(row[4])
                        # Take the extracted landform class from the landforms raster and use it as a conditional raster to
                        # output the landcover on just that landform class. If a particular land cover is over half of the
                        # landform class, write a sentence that says most of this landform zone is covered by a particular land
                        # cover class. If a particular land cover class is over 89%, change 'over half' to 'almost all'. And if
                        # that land cover class is over 99%, change 'over half' to 'all'.
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Bioclimate_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top bioclimate in the area covered by this subclass? landformlist[9]
                        landformlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top bioclimate? landformlist[10]
                        landformlist.append((count0/landformlist[2])*100)
            			#dummy values landformlist[11] and landformlist[12]
                        landformlist.append(dummyvalue)
                        landformlist.append(dummyvalue)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Lithology_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top rock type in the area covered by this subclass? landformlist[13]
                        landformlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top rock type? landformlist[14]
                        landformlist.append((count0/landformlist[2])*100)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landcover_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top land cover in the area covered by this subclass? landformlist[15]
                        landformlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top land cover? landformlist[16]
                        landformlist.append((count0/landformlist[2])*100)
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Ecophysdiv_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN"]) as cursor:
                            for row in cursor:
                                #What is the mean ELU diversity in this landform zone? landformlist[17]
                                landformlist.append(row[0])
                                #What is the median ELU diversity in this landform zone? landformlist[18]
                                #this is a dummy value because there is no median in floating point
                                landformlist.append(row[0])
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Slope_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN","MEDIAN"]) as cursor:
                            for row in cursor:
                                #What is the mean slope percentage of this landform subclass? landformlist[19]
                                landformlist.append(row[0])
                                #What is the median slope percentage of this landform subclass? landformlist[20]
                                landformlist.append(row[1])
                        # Take the extracted landform class raster and use it as a conditional raster, this time with the aspect
                        # raster aspectindexg. Next, sum up aspect counts for the 180 degrees which face all 8 directions. Of all
                        # 8 directions, take the largest number. If this number does not equal half of the cell count for the whole
                        # area of the extracted landform class, populate the variable aspectstatement with the phrase "not facing
                        # any direction, by majority. " Otherwise populate the variable aspectstatement with the phrase "generally
                        # facing ____. ", for example "north", filling in the blank with the direction with the largest cell count
                        # in its 180 degree face.
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        xbx=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\aspectindexg")
                        xbx.save(inmem+"\\cong")
                        aspecttuples = []
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Value","Count"]) as cursor:
                            for row in cursor:
                                aspecttuples.append([row[0],row[1]])
                            aspectdict0 = dict(aspecttuples)
                            aspectkeys = [-1,1,2,3,4,5,6,7,8]
                            for key in aspectkeys:
                                try:
                                    aspectdict0[key]
                                except:
                                    aspecttuples.append([key,0])
                            aspectdict = dict(aspecttuples)
                            face_n = int(aspectdict[7] + aspectdict[8] + aspectdict[1] + aspectdict[2])
                            face_ne = int(aspectdict[8] + aspectdict[1] + aspectdict[2] + aspectdict[3])
                            face_e = int(aspectdict[1] + aspectdict[2] + aspectdict[3] + aspectdict[4])
                            face_se = int(aspectdict[2] + aspectdict[3] + aspectdict[4] + aspectdict[5])
                            face_s = int(aspectdict[3] + aspectdict[4] + aspectdict[5] + aspectdict[6])
                            face_sw = int(aspectdict[4] + aspectdict[5] + aspectdict[6] + aspectdict[7])
                            face_w = int(aspectdict[5] + aspectdict[6] + aspectdict[7] + aspectdict[8])
                            face_nw = int(aspectdict[6] + aspectdict[7] + aspectdict[8] + aspectdict[1])
                            all_faces = int(aspectdict[1] + aspectdict[2] + aspectdict[3] + aspectdict[4] + aspectdict[5] + aspectdict[6] + aspectdict[7] + aspectdict[8] + aspectdict[-1])
                            largest_face = max(face_n, face_nw, face_w, face_sw, face_s, face_se, face_e, face_ne)
                            aspectstatement = ""
                            #aspectstatement = " "
                            #aspectstatement = "not facing any direction, by majority"
                            if float(largest_face)/float(all_faces) > .583:
                                if largest_face == face_n:
                                    aspectstatement = ", generally facing north"
                                elif largest_face == face_nw:
                                    aspectstatement = ", generally facing northwest"
                                elif largest_face == face_ne:
                                    aspectstatement = ", generally facing northeast"
                                elif largest_face == face_se:
                                    aspectstatement = ", generally facing southeast"
                                elif largest_face == face_s:
                                    aspectstatement = ", generally facing south"
                                elif largest_face == face_sw:
                                    aspectstatement = ", generally facing southwest"
                                elif largest_face == face_w:
                                    aspectstatement = ", generally facing west"
                                elif largest_face == face_e:
                                    aspectstatement = ", generally facing east"
                            #What aspect does this landform face, by majority?  landformlist[21]
                            landformlist.append(aspectstatement)

                        neq = 0
                        seq = 0
                        swq = 0
                        nwq = 0
                        zes = 0
                        zws = 0
                        zss = 0
                        zns = 0
                        zonesw = 0
                        zonew = 0
                        zonenw = 0
                        zones = 0
                        zonec = 0
                        zonen = 0
                        zonese = 0
                        zonee = 0
                        zonene = 0

                        try:
                            feature = os.path.join(inmem,'cong')
                            if arcpy.Exists(feature):
                            	arcpy.Delete_management(feature)
                            cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\mw_zonedg")
                            cong.save(inmem+"\\cong")
                            with arcpy.da.SearchCursor(inmem+"\\cong",["Value","Count"]) as cursor:
                                for rowmwz in cursor:
                                    if rowmwz[0] == 11:
                                        zonesw = rowmwz[1]
                                    if rowmwz[0] == 12:
                                        zonew = rowmwz[1]
                                    if rowmwz[0] == 13:
                                        zonenw = rowmwz[1]
                                    if rowmwz[0] == 21:
                                        zones = rowmwz[1]
                                    if rowmwz[0] == 22:
                                        zonec = rowmwz[1]
                                    if rowmwz[0] == 23:
                                        zonen = rowmwz[1]
                                    if rowmwz[0] == 31:
                                        zonese = rowmwz[1]
                                    if rowmwz[0] == 32:
                                        zonee = rowmwz[1]
                                    if rowmwz[0] == 33:
                                        zonene = rowmwz[1]

                            del cong
                            #number of cells in the southwest. landformlist[22]
                            landformlist.append(zonesw)
                            #number of cells in the west. landformlist[23]
                            landformlist.append(zonew)
                            #number of cells in the northwest. landformlist[24]
                            landformlist.append(zonenw)
                            #number of cells in the south. landformlist[25]
                            landformlist.append(zones)
                            #number of cells in the center. landformlist[26]
                            landformlist.append(zonec)
                            #number of cells in the north. landformlist[27]
                            landformlist.append(zonen)
                            #number of cells in the southeast. landformlist[28]
                            landformlist.append(zonese)
                            #number of cells in the east. landformlist[29]
                            landformlist.append(zonee)
                            #number of cells in the northeast. landformlist[30]
                            landformlist.append(zonene)
                            neq = landformlist[26]+landformlist[27]+landformlist[29]+landformlist[30]
                            seq = landformlist[26]+landformlist[25]+landformlist[29]+landformlist[28]
                            swq = landformlist[22]+landformlist[23]+landformlist[25]+landformlist[26]
                            nwq = landformlist[23]+landformlist[24]+landformlist[26]+landformlist[27]
                            zes = landformlist[28]+landformlist[29]+landformlist[30]
                            zws = landformlist[22]+landformlist[23]+landformlist[24]
                            zss = landformlist[22]+landformlist[25]+landformlist[28]
                            zns = landformlist[24]+landformlist[27]+landformlist[30]

                            maxzone = max([neq,'northeastern part',neq_n],[seq,'southeastern part',seq_n],[swq,'southwestern part',swq_n],[nwq,'northwestern part',nwq_n],[zes,'east side',zes_n],[zws,'west side',zws_n],[zss,'south side',zss_n],[zns,'north side',zns_n],[zonesw,'southwesternmost portion',zonesw_n],[zonew,'westernmost portion',zonew_n],[zonenw,'northwesternmost portion',zonenw_n],[zones,'southernmost portion',zones_n],[zonec,'most central portion',zonec_n],[zonen,'northernmost portion',zonen_n],[zonese,'southeasternmost portion',zonese_n],[zonee,'easternmost portion',zonee_n],[zonene,'northeasternmost portion',zonene_n])

                            #How many cells are in the most likely part of the landform class?[31]
                            landformlist.append(maxzone[0])
                            #What phrase describes the most likely part of this landform class?[32]
                            landformlist.append(maxzone[1])
                            #What percentage of the most likely part of the landform class is this landform class?[33]
                            if maxzone[2] != 0:
                                landformlist.append((maxzone[0]/maxzone[2])*100)
                            else:
                                landformlist.append(0)

                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0

                        except ValueError:
                            print("value error")
                            print("starting over")

                        except IndexError:
                            print("index error")
                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0
                            print("starting over")

                        classlist.append(landformlist)
                        #print landformlist
                        print("length of landformlist is "+str(len(landformlist)))

                        if len(landformlist) != 34:
                            GeoDescriberTries +=1
                            return

                    # LITHOLOGY side operations
                    # If any lithology classes are over 10% of the study area, find out more details about those classes.
                    elif analclass[2].endswith("Lithology_R") == True:
                    #elif analclass[2] == output+"\\Lithology_R":
                        lithologylist = []
                        #What is the pathname to this lithology subclass? lithologylist[0]
                        lithologylist.append(analclass[2])
                        #What's the rock type called?  lithologylist[1]
                        lithologylist.append(analclass[1])
                        currentTime = time.clock()
                        print(str(currentTime-startTime)+" seconds have elapsed so far")
                        print("-----performing side operations on lithology class " + analclass[1]+"-----")
                        lithologytimes += 1
                        attExtract = ExtractByAttributes(analclass[2], "ClassName = '"+analclass[1]+"'")

                        feature = os.path.join(inmem,'extractg')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)

                        attExtract.save(inmem+"\\extractg")
                        # Extract the lithology class from the lithology raster. Output the minimum and maximum elevation
                        # for that lithology raster into variables. Use those variables in a phrase which gives the minimum and
                        # maximum elevation as context for that bioclimate. " found at surface elevations between ______ and
                        # _______m. "
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Elevation_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["Count","MIN","MAX","MEAN","MEDIAN"]) as cursor:
                            for row in cursor:
                                #What is this rock type subclass' area? lithologylist[2]
                                lithologylist.append(row[0])
                                #What is this rock type subclass' percentage of the overall study area? lithologylist[3]
                                lithologylist.append(analclass[0])
                                #What is this rock type subclass' rank in order of greatest to least percent of the study area? lithologylist[4]
                                lithologylist.append(lithologytimes)
                                #What is this rock type subclass' lowest elevation? lithologylist[5]
                                lithologylist.append(row[1])
                                #What is this rock type subclass' highest elevation? lithologylist[6]
                                lithologylist.append(row[2])
                                #What is this rock type subclass' mean elevation? lithologylist[7]
                                lithologylist.append(row[3])
                                #What is this rock type subclass' median elevation? lithologylist[8]
                                lithologylist.append(row[4])
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Bioclimate_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top bioclimate in the area covered by this subclass? lithologylist[9]
                        lithologylist.append(classname0)
                        #What percentage of the area covered by this subclass is the top bioclimate? lithologylist[10]
                        lithologylist.append((count0/lithologylist[2])*100)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landform_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top landform in the area covered by this subclass? lithologylist[11]
                        lithologylist.append(classname0)
                        #What percentage of the area covered by this subclass is the top landform? lithologylist[12]
                        lithologylist.append((count0/lithologylist[2])*100)
                        #dummy lithologylist[13]
                        lithologylist.append(dummyvalue)
                        #dummy lithologylist[14]
                        lithologylist.append(dummyvalue)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landcover_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landcover_CR", "Value", ["Value","ClassName"])
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top land cover in the area covered by this subclass? lithologylist[15]
                        lithologylist.append(classname0)
                        #What percentage of the area covered by this subclass is the top land cover? lithologylist[16]
                        lithologylist.append((count0/lithologylist[2]) * 100)
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Ecophysdiv_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN"]) as cursor:
                            for row in cursor:
                                #What is the mean ELU diversity in this lithology zone? landformlist[17]
                                lithologylist.append(row[0])
                                #What is the median ELU diversity in this lithology zone? landformlist[18]
                                #this is a dummy value because there is no median in floating point
                                lithologylist.append(row[0])

                        #dummy value lithologylist[19]
                        lithologylist.append(dummyvalue)
                        #dummy value lithologylist[20]
                        lithologylist.append(dummyvalue)
                        #dummy value lithologylist[21]
                        lithologylist.append(dummyvalue)

                        neq = 0
                        seq = 0
                        swq = 0
                        nwq = 0
                        zes = 0
                        zws = 0
                        zss = 0
                        zns = 0
                        zonesw = 0
                        zonew = 0
                        zonenw = 0
                        zones = 0
                        zonec = 0
                        zonen = 0
                        zonese = 0
                        zonee = 0
                        zonene = 0

                        try:
                            feature = os.path.join(inmem,'cong')
                            if arcpy.Exists(feature):
                            	arcpy.Delete_management(feature)
                            cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\mw_zonedg")
                            cong.save(inmem+"\\cong")
                            with arcpy.da.SearchCursor(inmem+"\\cong",["Value","Count"]) as cursor:
                                for rowmwz in cursor:
                                    if rowmwz[0] == 11:
                                        zonesw = rowmwz[1]
                                    if rowmwz[0] == 12:
                                        zonew = rowmwz[1]
                                    if rowmwz[0] == 13:
                                        zonenw = rowmwz[1]
                                    if rowmwz[0] == 21:
                                        zones = rowmwz[1]
                                    if rowmwz[0] == 22:
                                        zonec = rowmwz[1]
                                    if rowmwz[0] == 23:
                                        zonen = rowmwz[1]
                                    if rowmwz[0] == 31:
                                        zonese = rowmwz[1]
                                    if rowmwz[0] == 32:
                                        zonee = rowmwz[1]
                                    if rowmwz[0] == 33:
                                        zonene = rowmwz[1]
                            del cong
                            #number of cells in the southwest. lithologylist[22]
                            lithologylist.append(zonesw)
                            #number of cells in the west. lithologylist[23]
                            lithologylist.append(zonew)
                            #number of cells in the northwest. lithologylist[24]
                            lithologylist.append(zonenw)
                            #number of cells in the south. lithologylist[25]
                            lithologylist.append(zones)
                            #number of cells in the center. lithologylist[26]
                            lithologylist.append(zonec)
                            #number of cells in the north. lithologylist[27]
                            lithologylist.append(zonen)
                            #number of cells in the southeast. lithologylist[28]
                            lithologylist.append(zonese)
                            #number of cells in the east. lithologylist[29]
                            lithologylist.append(zonee)
                            #number of cells in the northeast. lithologylist[30]
                            lithologylist.append(zonene)

                            neq = lithologylist[26]+lithologylist[27]+lithologylist[29]+lithologylist[30]
                            seq = lithologylist[26]+lithologylist[25]+lithologylist[29]+lithologylist[28]
                            swq = lithologylist[22]+lithologylist[23]+lithologylist[25]+lithologylist[26]
                            nwq = lithologylist[23]+lithologylist[24]+lithologylist[26]+lithologylist[27]
                            zes = lithologylist[28]+lithologylist[29]+lithologylist[30]
                            zws = lithologylist[22]+lithologylist[23]+lithologylist[24]
                            zss = lithologylist[22]+lithologylist[25]+lithologylist[28]
                            zns = lithologylist[24]+lithologylist[27]+lithologylist[30]

                            maxzone = max([neq,'northeastern part',neq_n],[seq,'southeastern part',seq_n],[swq,'southwestern part',swq_n],[nwq,'northwestern part',nwq_n],[zes,'east side',zes_n],[zws,'west side',zws_n],[zss,'south side',zss_n],[zns,'north side',zns_n],[zonesw,'southwesternmost portion',zonesw_n],[zonew,'westernmost portion',zonew_n],[zonenw,'northwesternmost portion',zonenw_n],[zones,'southernmost portion',zones_n],[zonec,'most central portion',zonec_n],[zonen,'northernmost portion',zonen_n],[zonese,'southeasternmost portion',zonese_n],[zonee,'easternmost portion',zonee_n],[zonene,'northeasternmost portion',zonene_n])

                            #How many cells are in the most likely part of the lithology class?[31]
                            lithologylist.append(maxzone[0])
                            #What phrase describes the most likely part of this lithology class?[32]
                            lithologylist.append(maxzone[1])
                            #What percentage of the most likely part of the lithology class is this lithology class?[33]
                            if maxzone[2] != 0:
                                lithologylist.append((maxzone[0]/maxzone[2])*100)
                            else:
                                lithologylist.append(0)

                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0

                        except ValueError:
                            print("value error")
                            print("starting over")

                        except IndexError:
                            print("index error")
                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0
                            print("starting over")

                        classlist.append(lithologylist)
                        #print lithologylist
                        print("length of lithologylist is "+str(len(lithologylist)))

                        if len(lithologylist) != 34:
                            GeoDescriberTries +=1
                            return

                    # LAND COVER side operations
                    # Extract the land cover class from the landcover raster. Use it as a conditional raster to cut out a piece of elevation data
                    # that matches the area of the landcover class raster. Use zonal statistics to find the mean elevation for the land cover.
                    # If the mean of the landcover is greater than the mean plus the standard deviation of the con_extent, add
                    # ", found at the higher elevations" to landcoverelevstr. If the mean of the bioclimate is less than the mean minus
                    # the standard deviation of con_extent, add ", found at the lower elevations" to the string landcoverelevstr.
                    elif analclass[2].endswith("Landcover_R") == True:
                    #elif analclass[2] == output+"\\Landcover_R":
                        landcoverlist = []
                        #What is the pathname to this land cover subclass? landcoverlist[0]
                        landcoverlist.append(analclass[2])
                        #What's this land cover called? landcoverlist[1]
                        landcoverlist.append(analclass[1])
                        currentTime = time.clock()
                        print(str(currentTime-startTime)+" seconds have elapsed so far")
                        print("-----performing side operations on land cover class " + analclass[1]+"-----")
                        landcovertimes += 1
                        # Compare the mean elevation of the significant land cover class with the mean elevation of the study area.
                        # If the class elevation is below or above one standard deviation from the mean elevation of the study area,
                        # add a clause to landcoverelevstr saying it's found at the higher or lower elevations. Otherwise leave it blank.
                        attExtract = ExtractByAttributes(analclass[2], "ClassName = '"+analclass[1]+"'")

                        feature = os.path.join(inmem,'extractg')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)

                        attExtract.save(inmem+"\\extractg")
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Elevation_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["Count","MIN","MAX","MEAN","MEDIAN"]) as cursor:
                            for row in cursor:
                                #What is this land cover class' area? landcoverlist[2]
                                landcoverlist.append(row[0])
                                #What is this land cover class' percentage of the overall study area? landcoverlist[3]
                                landcoverlist.append(analclass[0])
                                #What is this land cover subclass' rank in order of greatest to least percent of the study area? landcoverlist[4]
                                landcoverlist.append(landcovertimes)
                                #What is this land cover subclass' lowest elevation? landcoverlist[5]
                                landcoverlist.append(row[1])
                                #What is this land cover type subclass' highest elevation? landcoverlist[6]
                                landcoverlist.append(row[2])
                                #What is this land cover type subclass' mean elevation? landcoverlist[7]
                                landcoverlist.append(row[3])
                                #What is this land cover type subclass' median elevation? landcoverlist[8]
                                landcoverlist.append(row[4])
                        feature = os.path.join(inmem,'cong')
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Bioclimate_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Bioclimates_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top bioclimate in the area covered by this subclass? landcoverlist[9]
                        landcoverlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top bioclimate? landcoverlist[10]
                        landcoverlist.append((count0/landcoverlist[2])*100)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Landform_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Landform_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top landform in the area covered by this subclass? landcoverlist[11]
                        landcoverlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top landform? landcoverlist[12]
                        landcoverlist.append((count0/landcoverlist[2])*100)
                        if arcpy.Exists(feature):
                        	arcpy.Delete_management(feature)
                        cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\Lithology_R")
                        cong.save(inmem+"\\cong")
                        if arcversion == '10.5.1':
                            JoinField_Workaround(inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        else:
                            d=arcpy.JoinField_management (inmem+"\\cong", "Value", inmem+"\\Lithology_CR", "Value", ["Value","ClassName"])
                        sumLandCover = 0
                        with arcpy.da.SearchCursor(inmem+"\\cong",["Count","Value","ClassName"]) as cursor:
                            for row in sorted(cursor):
                                count0 = row[0]
                                classname0 = row[2]
                        #What is the top rock type in the area covered by this subclass? landcoverlist[13]
                        landcoverlist.append(classname0)
                        #What percentage of the area covered by this subclass is the top rock type? landcoverlist[14]
                        landcoverlist.append((count0/landcoverlist[2])*100)
                        landcoverlist.append(dummyvalue) #dummy landcoverlist[15]
                        landcoverlist.append(dummyvalue) #dummy landcoverlist[16]
                        a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Ecophysdiv_R",inmem+"\\stattbl","DATA")
                        with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN"]) as cursor:
                            for row in cursor:
                                #What is the mean ELU diversity in this land cover zone? landcoverlist[17]
                                landcoverlist.append(row[0])
                                #What is the median ELU diversity in this land cover zone? landcoverlist[18]
                                #this is a dummy value because there is no median in floating point
                                landcoverlist.append(row[0])

                        if landcoverlist[1] in ['Urban areas']:
                            feature = os.path.join(inmem,'Population_RS')
                            if arcpy.Exists(feature):
                            	arcpy.Delete_management(feature)
                            if arcversion == '10.5.1':
                                outputRaster0 = Con(arcpy.Raster(inmem+"\\extractg"),arcpy.Raster(inmem+"\\Population_R"), 0)
                                outputRaster = Con(outputRaster0 != 65535, outputRaster0)
                            else:
                                outputRaster = Con(arcpy.Raster(inmem+"\\extractg"),arcpy.Raster(inmem+"\\Population_R"), 0)
                            od = outputRaster.save(inmem +"\\Population_RS")
                            a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Population_RS",inmem+"\\stattbl","DATA")
                            if not arcpy.sa.Raster(inmem+"\\Population_RS").maximum > 0:
                                landcoverlist.append(0)
                            else:
                                with arcpy.da.SearchCursor(inmem+"\\stattbl",["SUM"]) as cursor:
                                    for row in cursor:
                                        #What is the population of this subclass? landcoverlist[19]
                                        landcoverlist.append(row[0])
                        else:
                            #Population will not be in text description, then just make landcoverlist[19] equal to -9999
                            landcoverlist.append(-9999)

                        if arcpy.sa.Raster(inmem+"\\Biomass_R").maximum > 0:
                            a=arcpy.sa.ZonalStatisticsAsTable(inmem+"\\extractg","Value",inmem+"\\Biomass_R",inmem+"\\stattbl","DATA")
                            with arcpy.da.SearchCursor(inmem+"\\stattbl",["MEAN"]) as cursor:
                                for row in cursor:
                                    #What is the mean biomass per acre of this subclass? landcoverlist[20]
                                    landcoverlist.append(row[0])
                        else:
                            #What is the mean biomass per acre of this subclass? landcoverlist[20]
                            #put dummy value in landcoverlist[20]
                            landcoverlist.append(dummyvalue)
                        landcoverlist.append(dummyvalue) #dummy landcoverlist[21]

                        neq = 0
                        seq = 0
                        swq = 0
                        nwq = 0
                        zes = 0
                        zws = 0
                        zss = 0
                        zns = 0
                        zonesw = 0
                        zonew = 0
                        zonenw = 0
                        zones = 0
                        zonec = 0
                        zonen = 0
                        zonese = 0
                        zonee = 0
                        zonene = 0

                        try:
                            feature = os.path.join(inmem,'cong')
                            if arcpy.Exists(feature):
                            	arcpy.Delete_management(feature)
                            cong=arcpy.sa.Con(inmem+"\\extractg",inmem+"\\mw_zonedg")
                            cong.save(inmem+"\\cong")
                            with arcpy.da.SearchCursor(inmem+"\\cong",["Value","Count"]) as cursor:
                                for rowmwz in cursor:
                                    if rowmwz[0] == 11:
                                        zonesw = rowmwz[1]
                                    if rowmwz[0] == 12:
                                        zonew = rowmwz[1]
                                    if rowmwz[0] == 13:
                                        zonenw = rowmwz[1]
                                    if rowmwz[0] == 21:
                                        zones = rowmwz[1]
                                    if rowmwz[0] == 22:
                                        zonec = rowmwz[1]
                                    if rowmwz[0] == 23:
                                        zonen = rowmwz[1]
                                    if rowmwz[0] == 31:
                                        zonese = rowmwz[1]
                                    if rowmwz[0] == 32:
                                        zonee = rowmwz[1]
                                    if rowmwz[0] == 33:
                                        zonene = rowmwz[1]
                            del cong
                            #number of cells in the southwest. landcoverlist[22]
                            landcoverlist.append(zonesw)
                            #number of cells in the west. landcoverlist[23]
                            landcoverlist.append(zonew)
                            #number of cells in the northwest. landcoverlist[24]
                            landcoverlist.append(zonenw)
                            #number of cells in the south. landcoverlist[25]
                            landcoverlist.append(zones)
                            #number of cells in the center. landcoverlist[26]
                            landcoverlist.append(zonec)
                            #number of cells in the north. landcoverlist[27]
                            landcoverlist.append(zonen)
                            #number of cells in the southeast. landcoverlist[28]
                            landcoverlist.append(zonese)
                            #number of cells in the east. landcoverlist[29]
                            landcoverlist.append(zonee)
                            #number of cells in the northeast. landcoverlist[30]
                            landcoverlist.append(zonene)

##                            if len(landcoverlist) == 30 and landcoverlist[1] == 'Water bodies':
##                                landcoverlist.insert(19,-9999)
##                            if len(landcoverlist) == 30 and landcoverlist[1] == 'Urban areas':
##                                landcoverlist.insert(19,-9999)

                            if len(landcoverlist) == 30:
                                landcoverlist.insert(19,-9999)

                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            neq = landcoverlist[26]+landcoverlist[27]+landcoverlist[29]+landcoverlist[30]
                            seq = landcoverlist[26]+landcoverlist[25]+landcoverlist[29]+landcoverlist[28]
                            swq = landcoverlist[22]+landcoverlist[23]+landcoverlist[25]+landcoverlist[26]
                            nwq = landcoverlist[23]+landcoverlist[24]+landcoverlist[26]+landcoverlist[27]
                            zes = landcoverlist[28]+landcoverlist[29]+landcoverlist[30]
                            zws = landcoverlist[22]+landcoverlist[23]+landcoverlist[24]
                            zss = landcoverlist[22]+landcoverlist[25]+landcoverlist[28]
                            zns = landcoverlist[24]+landcoverlist[27]+landcoverlist[30]

                            maxzone = max([neq,'northeastern part',neq_n],[seq,'southeastern part',seq_n],[swq,'southwestern part',swq_n],[nwq,'northwestern part',nwq_n],[zes,'east side',zes_n],[zws,'west side',zws_n],[zss,'south side',zss_n],[zns,'north side',zns_n],[zonesw,'southwesternmost portion',zonesw_n],[zonew,'westernmost portion',zonew_n],[zonenw,'northwesternmost portion',zonenw_n],[zones,'southernmost portion',zones_n],[zonec,'most central portion',zonec_n],[zonen,'northernmost portion',zonen_n],[zonese,'southeasternmost portion',zonese_n],[zonee,'easternmost portion',zonee_n],[zonene,'northeasternmost portion',zonene_n])
                            print "maxzone is "+ str(maxzone)

                            #How many cells are in the most likely part of the landcover class?[31]
                            landcoverlist.append(maxzone[0])
                            #What phrase describes the most likely part of this landcover class?[32]
                            landcoverlist.append(maxzone[1])
                            #What percentage of the most likely part of the landcover class is this landcover class?[33]
                            if maxzone[2] != 0:
                                landcoverlist.append((maxzone[0]/maxzone[2])*100)
                            else:
                                landcoverlist.append(0)

                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0

                        except ValueError:
                            print("value error")
                            print("starting over")

                        except IndexError:
                            print("index error")
                            neq = 0
                            seq = 0
                            swq = 0
                            nwq = 0
                            zes = 0
                            zws = 0
                            zss = 0
                            zns = 0
                            zonesw = 0
                            zonew = 0
                            zonenw = 0
                            zones = 0
                            zonec = 0
                            zonen = 0
                            zonese = 0
                            zonee = 0
                            zonene = 0
                            print("starting over")

                        classlist.append(landcoverlist)
                        print landcoverlist
                        print("length of landcoverlist is "+str(len(landcoverlist)))

                        if len(landcoverlist) != 34:
                            GeoDescriberTries +=1
                            return

                feature = os.path.join(inmem,'extractg')
                if arcpy.Exists(feature):
                	arcpy.Delete_management(feature)

            except:
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)


                #---------------------------------------------------------------------------------
                #
                #                               synthesis
                #
                #               Synthesize human readable text from classlist.
                #               Use it to populate the description variable.
                #
                #---------------------------------------------------------------------------------

        if studyarealist[4] > 15:

            try:
                divfragment = ""
                divstatement = ""
                studyareasqmeters = studyarealist[4] * cellsize
                if studyareasqmeters < 53000:
                	if studyarealist[12] < 0.1:
                		divfragment = "An area with very uniform landscape, "
                	elif studyarealist[12] >= 0.1 and studyarealist[12] < 0.5:
                		divfragment = "An area with uniform landscape, "
                	elif studyarealist[12] >= 0.5 and studyarealist[12] < 1:
                		print("")
                	elif studyarealist[12] >= 1 and studyarealist[12] < 2:
                		divstatement = studyarealist[0] + " has rather high landscape diversity. "
                	elif studyarealist[12] >= 2 and studyarealist[12] < 10:
                		divstatement = "For such a small area, "+studyarealist[0]+" has very high landscape diversity. "
                	else: #studyarealist[12] > 10
                		divstatement = "For such a small area, "+studyarealist[0]+" has extraordinarily high landscape diversity. "
                elif studyareasqmeters >= 53000 and studyareasqmeters < 1000000:
                	if studyarealist[12] < 0.1:
                		divfragment = "A very uniform landscape, "
                	elif studyarealist[12] >= 0.1 and studyarealist[12] < 0.5:
                		divfragment = "A uniform landscape, "
                	elif studyarealist[12] >= 0.5 and studyarealist[12] < 1:
                		print("")
                	elif studyarealist[12] >= 1 and studyarealist[12] < 2:
                		print("")
                	elif studyarealist[12] >= 2 and studyarealist[12] < 10:
                		divstatement = "For a small area, "+studyarealist[0]+" has very high landscape diversity. "
                	else: #studyarealist[12] > 10
                		divstatement = "For a small area, "+studyarealist[0]+" has extraordinarily high landscape diversity. "
                elif studyareasqmeters >= 1000000 and studyareasqmeters < 10000000:
                	if studyarealist[12] < 0.1:
                		divstatement = studyarealist[0]+" has extraordinarily low landscape diversity. "
                	elif studyarealist[12] >= 0.1 and studyarealist[12] < 0.5:
                		divstatement = studyarealist[0]+" has very low landscape diversity. "
                	elif studyarealist[12] >= 0.5 and studyarealist[12] < 1:
                		print("")
                	elif studyarealist[12] >= 1 and studyarealist[12] < 2:
                		print("")
                	elif studyarealist[12] >= 2 and studyarealist[12] < 10:
                		divstatement = studyarealist[0]+" is a highly diverse landscape. "
                	else: #studyarealist[12] > 10
                		divstatement = studyarealist[0]+" has extraordinarily high diversity in its landscape. "
                elif studyareasqmeters >= 10000000 and studyareasqmeters < 100000000:
                	if studyarealist[12] < 0.1:
                		divstatement = "For such a large area, "+studyarealist[0]+" has extraordinarily low diversity in its landscapes. "
                	elif studyarealist[12] >= 0.1 and studyarealist[12] < 0.5:
                		divstatement = "For such a large area, "+studyarealist[0]+" has very low diversity in its landscapes. "
                	elif studyarealist[12] >= 0.5 and studyarealist[12] < 1:
                		print("")
                	elif studyarealist[12] >= 1 and studyarealist[12] < 2:
                		print("")
                	elif studyarealist[12] >= 2 and studyarealist[12] < 10:
                		divfragment = "A rather large area with very diverse landscapes, "
                	else: #studyarealist[12] > 10
                		divfragment = "A rather large area with extraordinarily diverse landscapes, "
                else: #if studyareasqmeters >= 100000000:
                	if studyarealist[12] < 0.1:
                		divstatement = "For such a large area, "+studyarealist[0]+" has extraordinarily little diversity in its landscapes. "
                	elif studyarealist[12] >= 0.1 and studyarealist[12] < 0.5:
                		divstatement = "For such a large area, "+studyarealist[0]+" has very low landscape diversity. "
                	elif studyarealist[12] >= 0.5 and studyarealist[12] < 1:
                		divstatement = studyarealist[0]+" has rather low landscape diversity. "
                	elif studyarealist[12] >= 1 and studyarealist[12] < 2:
                		print("")
                	elif studyarealist[12] >= 2 and studyarealist[12] < 10:
                		divfragment = "A very diverse area, "
                	else: #studyarealist[12] > 10
                		divfragment = "An extraordinarily diverse area, "

            except:
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)

        if studyarealist[4] > 15:
            try:
                #Define and clear the variables used to test whether
                #elevation is a factor in temperature change.
                elev = -9999
                junk = ""
                jjunk = ""
                temptext = ""
                elevationalready = 0
                tempelevlist=[]
                tempelevlist0=[]
                temprangelist=[]
                temprangelist0=[]

                #Check to see if elevation is a factor in temperature change.
                #Find the range of temperatures (without moisture) in the study area.
                #Put them into a list called tempelevlist0.
                templist = ['very hot','hot','warm','cool','cold','very cold']
                for row in classlist:
                    if row[0].endswith("Bioclimate_R") == True and row[3] > 10:
                    #if row[0]== output+"\\Bioclimate_R" and row[3] > 10:
                        junk = bioclimate_dict[row[1]].rfind(' and')
                        jjunk =bioclimate_dict[row[1]][:junk]
                        tempelevlist0.append(row[8])
                        tempelevlist0.append(jjunk)
                        tempelevlist.append(tempelevlist0)
                    tempelevlist0 = []
                for temp in templist:
                    elev = -9999
                    x = 0
                    elevsum = 0
                    for row in sorted(tempelevlist):
                        if temp == row[1]:
                            x = x + 1
                            elevsum = elevsum + row[0]
                            elev = float(elevsum / x)
                    if elev > -1000:
                        temprangelist0.append(elev)
                        temprangelist0.append(temp)
                        temprangelist.append(temprangelist0)
                    temprangelist0 = []
                    elev = -9999
                n = 1

                #Elevation can be considered a factor in temperature change when there is a
                #change in the median elevation of 300m. When this occurs, synthesize text
                #where elevation is found to be a factor in temperature change. For example:
                #"As the elevation increases, bioclimate temperatures drop from cool to cold."
                #Since this may be when we bring up the subject of elevation for the first time,
                #This may also be the time to give an elevation range for the whole study area.
                while n < len(temprangelist):
                    if temprangelist[n][0]-temprangelist[n-1][0] > 300:
                        if n == 1:
                            temptext0000 = ""
                            temptext000 = ""
                            temptext0000 = divfragment+studyarealist[0]+" rises in elevation from {:,}".format(studyarealist[5])+" to {:,}".format(studyarealist[6])+" meters. "
                            temptext000 = temptext0000[0].capitalize()+temptext0000[1:-1]+" "
                            if studyarealist[5] == -4:
                                temptext = temptext000.replace(" rises in elevation from -4 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == -3:
                                temptext = temptext000.replace(" rises in elevation from -3 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == -2:
                                temptext = temptext000.replace(" rises in elevation from -2 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == -1:
                                temptext = temptext000.replace(" rises in elevation from -1 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == 0:
                                temptext = temptext000.replace(" rises in elevation from 0 to "," rises in elevation from sea level to ")
                            elif studyarealist[5] == 1:
                                temptext = temptext000.replace(" rises in elevation from 1 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == 2:
                                temptext = temptext000.replace(" rises in elevation from 2 to "," rises in elevation from about sea level to ")
                            elif studyarealist[5] == 3:
                                temptext = temptext000.replace(" rises in elevation from 3 to "," rises in elevation from about sea level to ")
                            else:
                                temptext = temptext000
                            temptext += "As the elevation increases, temperatures drop from "+temprangelist[n-1][1]+" to "+temprangelist[n][1]
                        else:
                            temptext += " and from "+temprangelist[n-1][1]+" to "+temprangelist[n][1]
                        elevationalready = 1
                    n = n + 1
                if temptext != "":
                    temptext0 = ""
                    temptext00 = ""
                    i = temptext.rfind(" and from")
                    temptext0 = temptext[:i].replace(" and from ", ", ")
                    temptext00 = temptext0 + temptext[i:].replace(" and from",", and from")
                    temptext = temptext00+". "
                tempsearch = [[', and from very hot to hot. ','As the elevation increases, temperatures drop from very hot to hot. '],[', and from hot to warm. ','As the elevation increases, temperatures drop from hot to warm. '],[', and from warm to cool. ','As the elevation increases, temperatures drop from warm to cool. '],[', and from cool to cold. ','As the elevation increases, temperatures drop from cool to cold. '],[', and from cold to very cold. ','As the elevation increases, temperatures drop from cold to very cold. ']]
                for pairs in tempsearch:
                    if temptext == pairs[0]:
                        temptext = pairs[1]

            except:
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)


        if studyarealist[4] > 15:
            try:
                #Define and clear variables used in language synthesis.
                n = 0
                nolandformsonlitho = 1
                nolandformsonlithostr0 = ""
                nolandformsonlithostr00 = ""
                nolandformsonlithostr000 = ""
                alllitholist = []
                bioclimatefractions = []
                landformfractions = []
                landcoverfractions = []
                lithologyfractions = []
                bioclimatefractiontext = ""
                landformfractiontext = ""
                landcoverfractiontext = ""
                lithologyfractiontext = ""
                bioclimatestr0 = ""
                landformstr0 = ""
                landcoverstr0 = ""
                lithologystr0 = ""
                elevrangestr = ""
                elevrangestr0 = ""
                landcoverredundency = 0
                lithologyredundency = 0
                bodyofwaterredundency = 0
                predescription = ""
                description = ""
                opener = ""
                biolocstr = ""
                litholocstr = ""
                landformlocstr = ""
                landcoverlocstr = ""
                if divstatement == "":
                    description = ""
                else:
                    opener = divstatement[0].capitalize()+divstatement[1:]
                    description = opener.replace('Study area ', 'The study area ')
                #conjunctions = ['','about a','another','yet another','still another','and another','another','another','another','another','another','another','another']
                conjunctions = ['','about a','about a','about a','about a','about a','about a','about a','about a','about a','about a','about a','about a']
                conjunctions_simple = ['','a','a','a','a','a','a','a','a','a','a','a','a']

                #Loop through the classlist, passing data assembled in the "significant class facts" phase
                #to conditional statements to synthesize the text.
                currentTime = time.clock()
                print(str(currentTime-startTime)+" seconds have elapsed so far")
                print("Synthesizing text.")
                for row in classlist:
                    #surface water occurs in both land cover and lithology classes. Flag the redundency so that
                    #it is only mentioned once in the synthesized text.
                    if row[11] == 'Surface Water':
                        bodyofwaterredundency = 1

                    #############################################################################
                    #if this is a bioclimate class, assemble the text in the following section:
                    if row[0].endswith("Bioclimate_R") == True:
                    #if row[0] == output+"\\Bioclimate_R":
                        bioclimatefractiontext = ""
                        bioclimatefractiontext0 = ""
                        landcoveronbioclim = ""
                        bioclimateelevationtext = ""

                        #If a bioclimate class is mostly on one side of a watershed, assemble text for that into the variable biolocstr.
                        #Otherwise the variable biolocstr is left blank.
                        biolocstr = ""
                        if row[31]/row[2] > .66:
                        	biolocstr = "Most of this "+bioclimate_dict[row[1]]+" bioclimate zone is in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .9:
                        	biolocstr = "Almost all of this "+bioclimate_dict[row[1]]+" bioclimate zone is in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .99:
                        	biolocstr = "This "+bioclimate_dict[row[1]]+" bioclimate zone is in the "+row[32]+" of the "+studyarealist[1]+". "

                        if row[4] == 1:
                            description += temptext

                        #If the median elevation of the bioclimate class occurs one half of a standard deviation above the mean,
                        #add text that it occurs at the higher elevations. If it's one half of a standard deviation below
                        #the mean, add text that it occurs at the lower elevations.
                        if row[7] > (studyarealist[7] + ((studyarealist[10]-studyarealist[7])/2)):
                            bioclimateelevationtext = "at the higher elevations, "
                        if row[7] < (studyarealist[7] - ((studyarealist[7]-studyarealist[9])/2)):
                            bioclimateelevationtext = "at the lower elevations, "
                        if row[4] > 1:
                            n = 1
                            bioclimatearticle = "the "
                        else:
                            n = 0
                            bioclimatearticle = article

                        #If the bioclimate class has mostly one land cover type, include an adjective phrase which
                        #gives detail about what land cover class is found there.
                        if row[16] < 95.000001 and landcoverredundency == 0:
                            if row[16] >= 50:
                                landcoveronbioclim = ", which is mostly covered by "+short_landcover_dict[row[15]]
                            if row[16] >= 89:
                                landcoveronbioclim = ", which is almost completely covered by "+short_landcover_dict[row[15]]
                            if row[16] >= 99:
                                landcoveronbioclim = ", which is completely covered by "+short_landcover_dict[row[15]]

                        #assemble descriptive text for the bioclimate class, capitalize it, and append it to the description variable.
                        if row[3] > 99.0:
                            bioclimatestr0 = bioclimateelevationtext+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 80:
                            bioclimatestr0 = bioclimateelevationtext+"most of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 60:
                            bioclimatestr0 = bioclimateelevationtext+"much of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 51:
                            bioclimatestr0 = bioclimateelevationtext+"just over half of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 49:
                            bioclimatestr0 = bioclimateelevationtext+"half of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 40:
                            bioclimatestr0 = bioclimateelevationtext+"just under half of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 30:
                            bioclimatestr0 = bioclimateelevationtext+conjunctions[row[4]]+" third of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]
                        elif row[3] > 20:
                            bioclimatestr0 = bioclimateelevationtext+conjunctions[row[4]]+" quarter of "+bioclimatearticle+studyarealist[n]+" has a "+ bioclimate_dict[row[1]] +" bioclimate"+landcoveronbioclim+ ". "+biolocstr
                            description += bioclimatestr0[0].capitalize()+bioclimatestr0[1:]

                        #If there is more than one fraction of a class, group them together and make all the fractions into one sentence.
                        else:
                            bioclimatefractions.append(row)
                            if row[4] == bioclimatetimes:
                                if len(bioclimatefractions) < 1:
                                    bioclimatefractiontext = ""
                                if len(bioclimatefractions) == 1:
                                    bioclimatefractiontext = bioclimateelevationtext+conjunctions_simple[row[4]]+ " fraction of the "+studyarealist[n]+" has a "+bioclimate_dict[row[1]]+" bioclimate"+landcoveronbioclim+". "
                                else:
                                    bioclimatefractiontext += "fractions of the "+studyarealist[1]+" have "
                                    for row in bioclimatefractions:
                                        bioclimatefractiontext += bioclimate_dict[row[1]] + ", "
                                    bioclimatefractiontext0 = bioclimatefractiontext[:-2]+" bioclimates. "
                                    k = bioclimatefractiontext0.rfind(",")
                                    if len(bioclimatefractions) == 2:
                                        bioclimatefractiontext = bioclimatefractiontext0[:k] + " and" + bioclimatefractiontext0[k+1:]
                                    else:
                                        bioclimatefractiontext = bioclimatefractiontext0[:k] + ", and" + bioclimatefractiontext0[k+1:]
                        if row[4] == bioclimatetimes:
                            if len(bioclimatefractiontext) > 2:
                                description += bioclimatefractiontext[0].capitalize()+bioclimatefractiontext[1:-1]+"</p><p>"
                            else:
                                description += "</p><p>"

                    #############################################################################
                    #if this is a landform class, assemble the text in the following section:
                    if row[0].endswith("Landform_R") == True:
                    #if row[0] == output+"\\Landform_R":
                        landformfractiontext = ""
                        landformfractiontext0 = ""
                        landformslopetext = ""
                        landformaspecttext = ""
                        landformelevationtext = ""
                        landformelevationtext0 = ""
                        landformlithologytext = ""
                        if row[4] > 1:
                            n = 1
                            landformarticle = "the "
                        else:
                            n = 0
                            landformarticle = article


                        #If a landform class is mostly on one side of a watershed, assemble text for that into the variable landformlocstr.
                        #Otherwise the variable landformlocstr is left blank.
                        landformlocstr = ""
                        if row[31]/row[2] > .66:
                        	landformlocstr = "Most of these "+landform_dict[row[1]]+" are on the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .9:
                        	landformlocstr = "Almost all of these "+landform_dict[row[1]]+" are on the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .99:
                        	landformlocstr = "These "+landform_dict[row[1]]+" are on the "+row[32]+" of the "+studyarealist[1]+". "

                        #If landforms are hills or mountains, and the mean percent slope is > 40%,
                        #append the adjective "extremely steep" to the description.
                        if row[1] != 'Surface Water':
                            if 'ills' in landform_dict[row[1]]:
                                landformaspectext = row[21]
            ##                    if row[19] > 8:
            ##                        landformslopetext = "moderately sloped "
            ##                    if row[19] > 20:
            ##                        landformslopetext = "steep "
                                if row[19] > 40:
                                    landformslopetext = "extremely steep "
                            if 'ountains' in landform_dict[row[1]]:
                                landformaspecttext = row[21]
            ##                    if row[19] > 8:
            ##                        landformslopetext = "moderately sloped "
            ##                    if row[19] > 20:
            ##                        landformslopetext = "steep "
                                if row[19] > 40:
                                    landformslopetext = "extremely steep "

                            #If the median elevation of the bioclimate class occurs a standard deviation above the mean,
                            #add text that it occurs at the higher elevations. If it's a standard deviation below
                            #the mean, add text that it occurs at the lower elevations.
                            if row[7] > studyarealist[10]:
                                landformelevationtext = "at higher elevations, "
                            if row[7] < studyarealist[9]:
                                landformelevationtext = "at lower elevations, "

                            #If the landform class has mostly one rock type, include an adjective phrase which
                            #gives detail about what lithology is found there.
                            if row[14] < 95.000001 and lithologyredundency == 0:
                                if row[14] >= 50:
                                    landformlithologytext = ", mostly composed of "+lithology_dict[row[13]]
                                if row[14] >= 89:
                                    landformlithologytext = ", almost all composed of "+lithology_dict[row[13]]
                                if row[14] >= 99:
                                    landformlithologytext = ", completely composed of "+lithology_dict[row[13]]

                            #assemble descriptive text for the landform class, capitalize it, and append it to the description variable.
                            if row[3] > 99.0:
                                landformstr0 = landformelevationtext+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 80:
                                landformstr0 = landformelevationtext+"most of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 60:
                                landformstr0 = landformelevationtext+"much of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 51:
                                landformstr0 = landformelevationtext+"just over half of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 49:
                                landformstr0 = landformelevationtext+"half of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 40:
                                landformstr0 = landformelevationtext+"just under half of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 30:
                                landformstr0 = landformelevationtext+conjunctions[row[4]]+" third of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]
                            elif row[3] > 20:
                                landformstr0 = landformelevationtext+conjunctions[row[4]]+" quarter of "+landformarticle+studyarealist[n]+" is "+landformslopetext+ landform_dict[row[1]] +landformaspecttext+landformlithologytext+". "+landformlocstr
                                description += landformstr0[0].capitalize()+landformstr0[1:]

                            #If there is more than one fraction of a class, group them together, making the fractions into one sentence.
                            else:
                                landformfractions.append(row)
                                if row[4] == landformtimes:
                                    if len(landformfractions) < 1:
                                        landformfractiontext = ""
                                    if len(landformfractions) == 1:
                                        landformfractiontext = landformelevationtext+conjunctions_simple[row[4]]+" fraction of the "+studyarealist[2]+" "+landformslopetext+landform_dict[row[1]]+landformaspecttext+landformlithologytext+". "
                                    else:
                                        landformfractiontext += "fractions of the "+studyarealist[1]+" form "
                                        for row in landformfractions:
                                            landformfractiontext += landform_dict[row[1]] + ", "
                                        landformfractiontext0 = landformfractiontext[:-2]+". "
                                        k = landformfractiontext0.rfind(",")
                                        if len(landformfractions) == 2:
                                            landformfractiontext = landformfractiontext0[:k] + " and" + landformfractiontext0[k+1:]
                                        else:
                                            landformfractiontext = landformfractiontext0[:k] + ", and" + landformfractiontext0[k+1:]
                            if row[4] == landformtimes:
                                if len(landformfractiontext) > 2:
                                    description += landformfractiontext[0].capitalize()+landformfractiontext[1:-1]+"</p><p>"
                                else:
                                    description += "</p><p>"

                    #############################################################################
                    #if this is a land cover class, assemble the text in the following section:
                    if row[0].endswith("Landcover_R") == True:
                    #if row[0] == output+"\\Landcover_R":
                        landcoverfractiontext = ""
                        landcoverfractiontext0 = ""
                        landcoverpoptext = ""
                        landcoverelevstr = ""
                        densitytext = ""
                        waterbodypoptext = ""
                        covertext = " covers "
                        thelandcover = landcover_dict[row[1]]

                        #If a land cover class is mostly on one side of a watershed, assemble text for that into the variable landcoverlocstr.
                        #Otherwise the variable landcoverlocstr is left blank.
                        landcoverlocstr = ""
                        landcoverlocstr0 = ""
                        landcoverlocstr00 = ""
                        if row[31]/row[2] > .66:
                        	landcoverlocstr0 = "Most of this "+landcover_dict[row[1]]+" is in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .9:
                        	landcoverlocstr0 = "Almost all of this "+landcover_dict[row[1]]+" is in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .99:
                        	landcoverlocstr0 = "This "+landcover_dict[row[1]]+" is in the "+row[32]+" of the "+studyarealist[1]+". "
                        landcoverlocstr00 = landcoverlocstr0.replace('this urban areas is','these urban areas are')
                        landcoverlocstr = landcoverlocstr00.replace('This urban areas is','These urban areas are')

                        if row[4] > 1:
                            n = 1
                            landcoverarticle = "the "
                        else:
                            n = 0
                            landcoverarticle = article
                        if len(row) > 20:
                            if 'Tree cover' in row[1] and row[20] > 7250:
                                densitytext = "dense, "
                        ##if row[1] in ['Cropland, rainfed','Cropland, rainfed - Herbaceous cover','Cropland, rainfed - Tree or shrub cover','Cropland irrigated or post-flooding','Mosaic cropland (>50%) / natural vegetation (Tree, shrub, herbaceous cover) (<50%)','Mosaic natural vegetation (Tree, shrub, herbaceous cover) (>50%) / cropland (<50%)']:
                            ##landcoverpoptext = "with an estimated population of "+"{:,}".format(int(round(row[19],-2)))+", "
                        if row[1] in ['Urban areas']:
                            landcoverpoptext = "with an estimated population of "+"{:,}".format(int(round(row[19],-3)))+", "
                            covertext = " cover "
            ##            if row[1] in ['Water bodies']:
            ##                waterbodypoptext = " (approximately "+str(round(studyarealist[11],1))+"% of the "+studyarealist[1]+")"
            ##                covertext = " cover "
                        if row[1] in ['Lichens and mosses']:
                            covertext = " cover "

                        #If the land cover class is mostly bare ground or sparsely vegetated,
                        #give detail about type of bare ground is found there (from the lithology).
                        if row[1] in ['Bare areas','Consolidated bare areas','Unconsolidated bare areas']:
                            if row[14] > 50:
                                thelandcover = landcover_dict[row[1]]+", mostly "+lithology_dict[row[13]]+","
                            if row[14] > 90:
                                thelandcover = landcover_dict[row[1]]+", almost all "+lithology_dict[row[13]]+","
                            if row[14] > 99:
                                thelandcover = landcover_dict[row[1]]+", composed of "+lithology_dict[row[13]]+","
                        if row[1] in ['Sparse shrub (<15%)','Sparse herbaceous cover (<15%)','Sparse vegetation (tree, shrub, herbaceous cover) (<15%)','Lichens and mosses']:
                            if row[14] > 50:
                                thelandcover = landcover_dict[row[1]]+", mostly on "+lithology_dict[row[13]]+","
                            if row[14] > 90:
                                thelandcover = landcover_dict[row[1]]+", almost all of it on "+lithology_dict[row[13]]+","
                            if row[14] > 99:
                                thelandcover = landcover_dict[row[1]]+" on "+lithology_dict[row[13]]

                        if row[7] > studyarealist[10]:
                            landcoverelevstr = "at the higher elevations, "
                        if row[7] < studyarealist[9]:
                            landcoverelevstr = "at the lower elevations, "

                        #assemble descriptive text for the land cover class, capitalize it, and append it to the description variable.
                        #If one of the significant classes is bodies of water, specify the quantity from the 30m water bodies raster.
                        if row[1] in ['Water bodies'] and studyarealist[11] > 12 and studyarealist[11] < 100:
                            #landcoverstr0 = "Bodies of water cover "+str(round(studyarealist[11],1))+"% of the "+studyarealist[1]+". "
                            landcoverstr0 = "Bodies of water cover about "+str(round(studyarealist[11],0))[:-2]+"% of the "+studyarealist[1]+". "
                            description += landcoverstr0
                        else:
                            if row[3] > 99.0:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+" "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 80:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+"most of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 60:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+"much of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 51:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+"just over half of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 49:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+"half of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 40:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+"just under half of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 30:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+conjunctions[row[4]]+" third of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]
                            elif row[3] > 20:
                                landcoverstr0 = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+conjunctions[row[4]]+" quarter of "+landcoverarticle+studyarealist[n]+". "+landcoverlocstr
                                description += landcoverstr0[0].capitalize()+landcoverstr0[1:]

                            #If there is more than one fraction of a class, group them together, making the fractions into one sentence.
                            else:
                                landcoverfractions.append(row)
                                if row[4] == landcovertimes:
                                    if len(landcoverfractions) < 1:
                                        landcoverfractiontext = ""
                                    if len(landcoverfractions) == 1:
                                        landcoverfractiontext = landcoverpoptext+landcoverelevstr+densitytext+thelandcover+waterbodypoptext+covertext+conjunctions_simple[row[4]]+" fraction of the "+studyarealist[n]+". "
                                    else:
                                        landcoverfractiontext += "fractions of the "+studyarealist[1]+" are covered by "
                                        for row in landcoverfractions:
                                            landcoverfractiontext += landcover_dict[row[1]] + ", "
                                        landcoverfractiontext0 = landcoverfractiontext[:-2]+". "
                                        k = landcoverfractiontext0.rfind(",")
                                        if len(landcoverfractions) == 2:
                                        	landcoverfractiontext = landcoverfractiontext0[:k] + " and" + landcoverfractiontext0[k+1:]
                                        else:
                                        	landcoverfractiontext = landcoverfractiontext0[:k] + ", and" + landcoverfractiontext0[k+1:]
                        if row[4] == landcovertimes:
                            if len(landcoverfractiontext) > 2:
                                description += landcoverfractiontext[0].capitalize()+landcoverfractiontext[1:-1]+"</p><p>"
                            else:
                                description += "</p><p>"

                    #############################################################################
                    #if this is a lithology class, assemble the text in the following section:
                    if row[0].endswith("Lithology_R") == True:
                    #if row[0] == output+"\\Lithology_R":
                        alllitholist.append(row)
                        lithologyfractiontext = ""
                        lithologyfractiontext0 = ""
                        lithologylandformtext = ""
                        if row[4] > 1:
                            n = 1
                            lithologyarticle = "the "
                        else:
                            n = 0
                            lithologyarticle = article

                        #If a lithology class is mostly on one side of a watershed, assemble text for that into the variable litholocstr.
                        #Otherwise the variable litholocstr is left blank.
                        litholocstr = ""
                        litholocstr0 = ""
                        litholocstr00 = ""
                        if row[31]/row[2] > .66:
                        	litholocstr0 = "Most of these "+lithology_dict[row[1]]+" are in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .9:
                        	litholocstr0 = "Almost all of these "+lithology_dict[row[1]]+" are in the "+row[32]+" of the "+studyarealist[1]+". "
                        if row[31]/row[2] > .99:
                        	litholocstr0 = "These "+lithology_dict[row[1]]+" are in the "+row[32]+" of the "+studyarealist[1]+". "
                        litholocstr00 = litholocstr00.replace('these unconsolidated sediment are','this unconsolidated sediment is')
                        litholocstr = litholocstr0.replace('These unconsolidated sediment are','This unconsolidated sediment is')

                        #If the lithology class is mostly one type of landform, include an adjective phrase which
                        #gives detail about what landform is found there. #most popular landform on the grassland is landformd_dict[row[11]]
                        if row[12] < 95.000001 and lithologyredundency == 0:
                            if row[12] >= 50:
                                lithologylandformtext = ", mostly forming "+landform_dict[row[11]]+","
                                if row[3] > 20:
                                    nolandformsonlitho = 0
                            if row[12] >= 89:
                                lithologylandformtext = ", almost all of which form "+landform_dict[row[11]]+","
                                if row[3] > 20:
                                    nolandformsonlitho = 0
                            if row[12] >= 99:
                                lithologylandformtext = ", forming "+landform_dict[row[11]]+","
                                if row[3] > 20:
                                    nolandformsonlitho = 0
                        if row[3] > 99.0:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 80:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie most of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 60:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie much of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 51:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie just over half of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 49:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie half of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 40:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie just under half of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 30:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie "+conjunctions[row[4]]+" third of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]
                        elif row[3] > 20:
                            lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie "+conjunctions[row[4]]+" quarter of "+lithologyarticle+studyarealist[n]+". "+litholocstr
                            predescription += lithologystr0[0].capitalize()+lithologystr0[1:]

                        #If there is more than one fraction of a class, group them together, making the fractions into one sentence.
                        else:
                            lithologyfractions.append(row)
                            if row[4] == lithologytimes:
                                if len(lithologyfractions) < 1:
                                    lithologyfractiontext = ""
                                if len(lithologyfractions) == 1:
                                    lithologystr0 = lithology_dict[row[1]] +lithologylandformtext+ " underlie "+conjunctions_simple[row[4]]+" fraction of the "+studyarealist[n]+". "
                                else:
                                    lithologyfractiontext += "fractions of the "+studyarealist[1]+" are "
                                    for row in lithologyfractions:
                                        lithologyfractiontext += lithology_dict[row[1]] + ", "
                                    lithologyfractiontext0 = lithologyfractiontext[:-2]+". "
                                    k = lithologyfractiontext0.rfind(",")
                                    if len(lithologyfractions) == 2:
                                        lithologyfractiontext = lithologyfractiontext0[:k] + " and" + lithologyfractiontext0[k+1:]
                                    else:
                                        lithologyfractiontext = lithologyfractiontext0[:k] + ", and" + lithologyfractiontext0[k+1:]
                        if row[4] == lithologytimes:
                            if nolandformsonlitho == 1:
                                nolandformsonlithostr = "A mix of "
                                for row in alllitholist:
                                    if row[3] > 20:
                                        nolandformsonlithostr += lithology_dict[row[1]]+" and "
                                nolandformsonlithostr0 = ""
                                nolandformsonlithostr += " predominate throughout the landforms of "+studyarealist[0]+", "+lithologyfractiontext.replace('fractions of the watershed are areas of ','along with small amounts of ')
                                nolandformsonlithostr0 = nolandformsonlithostr.replace('A mix of  predominate throughout','A mix with no dominant lithology characterizes')
                                nolandformsonlithostr00 = nolandformsonlithostr0.replace('. . ','. ')
                                nolandformsonlithostr000 = nolandformsonlithostr00.replace('and  predominate','predominate')
                                description += nolandformsonlithostr000[:-3]+nolandformsonlithostr000[-3:].replace(', ','. ') +"</p><p>"
                            else:
                                if len(lithologyfractiontext) > 2:
                                    description += predescription+lithologyfractiontext[0].capitalize()+lithologyfractiontext[1:-1]+"</p><p>"
                                else:
                                    description += predescription+"</p><p>"

            except:
                print("if you have a name field in your study area polygon, does the feature")
                print("have a name? Check to make sure the field is populated.")
                # Get the traceback object
                tb = sys.exc_info()[2]
                tbinfo = traceback.format_tb(tb)[0]
                # Concatenate information together concerning the error into a message string
                pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
                # Write Python error messages to log
                err= pymsg + "\n"
                print(err)

        if studyarealist[4] < 16:
            description = "This area is too small to meaningfully describe."

# Add a field to the study area called Description, populate it with four paragraphs (separated by the </p><p> tags).
# The Description paragraphs use adjectives that depict ranges of significant classes, such as "most of this area".
# It is populated by the variable description. inFeatLyr FL_MollPrj
        with arcpy.da.UpdateCursor(inFeatLyr,["OBJECTID","Description"]) as cursor:
            for row in cursor:
                if row[0] == intPolyID:
                    #print(description)
                    if description[-7:] == '</p><p>':
                        row[1]= '<p>'+description[:-7]+'</p>'
                    else:
                        row[1]= '<p>'+description+'</p>'
                    print(row[1])
                cursor.updateRow(row)
                GeoDescriberTries = 1

        del cursor, row
        currentTime = time.clock()
        print("This polygon took "+str(currentTime-thisPolyTime)+" seconds.")
        print("Feature class so far has taken "+str(currentTime-startTime)+" seconds.")

    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
        # Write Python error messages to log
        err= pymsg + "\n"
        print(err)



if __name__ == "__main__":

    # Geodescriber will count how many times it tried to describe a polygon.
    # This count is stored in the variable GeoDescriberTries.
    global GeoDescriberTries
    GeoDescriberTries = 1

    # Make sure all temporary datasets are removed from both in_memory and disk.
    CleanUp()

    ##global currentTime
    currentTime = 0
    arcversion = arcpy.GetInstallInfo()['Version']
    if arcversion == '10.5.1':
        print("running 10.5.1")
    else:
        print("not running 10.5.1")
    print("0 seconds have elapsed so far")
    print("Adding field...")
    arcpy.AddField_management (inFeatLyr, "Description", "TEXT", "", "", "50000")
    #print(str(currentTime-startTime)+" seconds have elapsed so far")
    sr = arcpy.SpatialReference(54009)
    arcpy.env.outputCoordinateSystem = sr
    arcpy.env.overwriteOutput = True
    warnthreshold = (cellsize * cellsize) * 1000
    arcpy.env.cellSize = cellsize

    try:
        print("deleting proj...")
        cleanupg = ['proj']
        for fd in cleanupg:
            feature = os.path.join(output,fd)
            if arcpy.Exists(feature):
                arcpy.Delete_management(feature)
                print("deleting "+feature+"...")
            feature = os.path.join(inmem,fd)
            if arcpy.Exists(feature):
                arcpy.Delete_management(feature)
                print("deleting "+feature+"...")

    except:
        # Get the traceback object
        tb = sys.exc_info()[2]
        tbinfo = traceback.format_tb(tb)[0]
        # Concatenate information together concerning the error into a message string
        pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
        # Write Python error messages to log
        err= pymsg + "\n"
        print(err)

    print("Projecting proj into mollweide...")
    FL_MollPrj =arcpy.Project_management(inFeatLyr,output+"\\proj",sr)
    listFeatureIDs=[]
    with arcpy.da.SearchCursor(FL_MollPrj,["OID@"]) as cursorid:
        for rowid in cursorid:
            rowidk=rowid[0]
            listFeatureIDs = listFeatureIDs + [rowidk]

    #try three times to describe each polygon.
    for intPolyID in listFeatureIDs:
        try:
            if GeoDescriberTries == 1:
                CleanUp()
                GeoDescriber()
            if GeoDescriberTries == 2:
                print("second try to run GeoDescriber()")
                CleanUp()
                GeoDescriber()
            if GeoDescriberTries == 3:
                print("third and final try to run GeoDescriber()")
                CleanUp()
                GeoDescriber()
            if GeoDescriberTries > 3:
                print("tried three times to run GeoDescriber()")

        except:
            print("GeoDescriber() did not work")
            # Get the traceback object
            tb = sys.exc_info()[2]
            tbinfo = traceback.format_tb(tb)[0]
            # Concatenate information together concerning the error into a message string
            pymsg = tbinfo + "\n" + str(sys.exc_type)+ ": " + str(sys.exc_value)
            # Write Python error messages to log
            err= pymsg + "\n"
            print(err)