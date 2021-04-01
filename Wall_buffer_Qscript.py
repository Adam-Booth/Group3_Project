## CREATES Buffer
# Output mannually uploaded to become an ee object

from PyQt5.QtGui import *

#file path for hadrians wall shape file
w_fp = r"C:\Users\lgxsv2\OneDrive - The University of Nottingham\PhD\yr_1\group_project\NDVI\hadrians_wall.shp"


#Add this layer to the map
w_layer = iface.addVectorLayer(w_fp, 'wall', 'ogr')


#Output name variable for buffer
outfn = 'buffer_350.shp'


#Buffer distance == distance between vallum and roman road at Magna
bufferDist = 350

#processing buffer
#disolved to create just one polygon
processing.run('native:buffer', {'INPUT':w_fp, "DISTANCE":bufferDist, \
'DISSOLVE':True, 'OUTPUT':outfn})

#Add to map 
iface.addVectorLayer(outfn, "", 'ogr')
