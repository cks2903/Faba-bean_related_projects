#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 09:37:11 2020

@author: CathrineKiel
"""

import pandas as pd
import numpy as np
pd.options.display.precision = 20
pd.set_option('display.precision', 20)
#!pip install pyproj # to install using pip
#!pip install geopandas
from pyproj import Proj
import matplotlib.pyplot as plt
#!pip install geopy
import geopy.distance
import gmplot
from termcolor import colored



# Define functions


#fix this one so that it subtract easting when it has to and add it when it has to
def calculate_y_fromfixedx(start_east,start_north):
    '''This function takes in GPS coordinates as 
    UTM (easting and northing) and predicts how 
    GPS coordinates of plots change as we move
    from west to east, but remain almost at the
    same northing position
    '''
    if Latitude_axis==1:
        plot_numbers=number_plots_y
        meters=plotsize_y
        path_btw_plots=plot_path_y
        plots_until_big_path=no_of_plots_until_big_path_y_axis
        lengthofbigpath=rep_path
        # check if going down the x-axis is going to east or west
        East=df1.Meters_East[2]-df1.Meters_East[0] #if <0 we are going west
        
    indexwithshift=list()
    for i in range(1,blocks_along_y_axis):
        indexwithshift.append(i*plots_until_big_path+1)     
    y_coordinates_north=list()
    y_coordinates_east=list()
    y_coordinates=list()
    x_coordinates=list()

    current_easting=start_east
    current_northing=start_north


    for i in range(1,(plot_numbers+1)):
        
        if i==1:
            current_northing=start_north
            current_easting=start_east
            
        
    
        if i in indexwithshift:
            current_northing+=slope_Northing_avg*(meters+lengthofbigpath)
            
            if East>0:
                current_easting+=(meters+lengthofbigpath)
                
            if East<0:
                current_easting-=(meters+lengthofbigpath) #men nu går vi ikke mod øst men mod vest så minus?
            
        if i%2==0 and rowaxis=='y': #equal numbers
            current_northing+=slope_Northing_avg*(meters)
            if East>0:
                current_easting+=(meters)
                
            if East<0:
                current_easting-=(meters)
                
        else: #all unequal numbers that are not the first in a block
            current_northing+=slope_Northing_avg*(meters+path_btw_plots)
            if East>0:
                current_easting+=(meters+path_btw_plots)
                
            if East<0:
                current_easting-=(meters+path_btw_plots)
                
            
            

    
        y_coordinates_east.append(current_easting)
        y_coordinates_north.append(current_northing)
        y_coordinates.append(i)
        x_coordinates.append(1)
        x0=pd.DataFrame(np.c_[x_coordinates,y_coordinates,y_coordinates_east,y_coordinates_north],columns=['x','y','Meters_East','Meters_North'])
    
    return(x0)


#not fixed yet
def calculate_x_fromfixedy(start_east,start_north,y):
    '''This function takes in GPS coordinates as 
    UTM (easting and northing) and predicts how 
    GPS coordinates of plots change as we move
    from north to south, but remain almost at the
    same easting position (we stay in the same 
    coordinate along that axis).
    Inputs are start_position_as_easting,
    start_position_as_northing, and position as
    integer on the axis where we are not moving
    '''
    if Latitude_axis==1:
        plot_numbers=number_plots_x
        meters=plotsize_x
        plots_until_big_path=no_of_plots_until_big_path_x_axis
        lengthofbigpath=xpath
        path_btw_plots=Betweenparcel_path_x
        North=df1.Meters_South[1]-df1.Meters_South[0] #if <0 we are going south
   
    
    x_coordinates_north=list()
    x_coordinates_east=list()
    y_coordinates=list()
    x_coordinates=list()

    current_easting=start_east
    current_northing=start_north
    
    indexwithshift=list()
    for i in range(1,blocks_along_x_axis):
        indexwithshift.append(i*plots_until_big_path+1) 

    for i in range(1,(plot_numbers+1)):
        
        if i==1:
            current_northing=start_north
            current_easting=start_east
    
        if i in indexwithshift:
            current_easting+=slope_Easting_avg*(meters+lengthofbigpath)
            if North<0:
                current_northing+=-(meters+lengthofbigpath)
            if North>0:
                current_northing+=(meters+lengthofbigpath)
                
        if i%2==0 and rowaxis=='x': #equal numbers
            current_easting+=slope_Easting_avg*(meters)
            if North>0:
                current_northing+=(meters)
                
            if North<0:
                current_northing-=(meters)
                
                
        else:
            current_easting+=slope_Easting_avg*(meters+path_btw_plots)
            if North<0:
                current_northing+=-(meters+path_btw_plots)
            if North>0:
                current_northing+=(meters+path_btw_plots)
            
        x_coordinates_east.append(current_easting)
        x_coordinates_north.append(current_northing)
        y_coordinates.append(y+1)
        x_coordinates.append(i)
        y0=pd.DataFrame(np.c_[x_coordinates,y_coordinates,x_coordinates_east,x_coordinates_north],columns=['x','y','Meters_East','Meters_North'])
    
    return(y0)


# introduction 
input(colored("Welcome to the GPS prediction function. Please read the documentation file to get an overview of input data needed. Press ENTER to proceed","green"))
answer=input(colored("First, you should convert the axis so that the x-axis is the one displaying changes in latitude (Going from North to South or vice versa. Have you done this? (y/n)","green"))
if answer=="n":
    exit()    

# Set parameters for the model

Name = input(colored("Give a name to your project: ","green"))
#Name="Nordic Seed 2019_Core"

#x=1,y=1
Corner1 = input(colored("UTM coordinates of x=1,y=1 corner (ex. 578779.86,6201856.26 : ","green"))
#ex Corner1=578779.86,6201856.26
Corner1 = Corner1.strip('][').split(',')
Corner1=[float(i) for i in Corner1]
#x=endpoint,y=1
Corner2 = input(colored("UTM coordinates of x=endpoint,y=1 corner: ",'green')) 
#Corner2=578783.17,6201829.74
Corner2 = Corner2.strip('][').split(',')
Corner2=[float(i) for i in Corner2]


# x=1, y=endpoint
Corner3 = input(colored("UTM coordinates of x=1,y=endpoint corner: ",'green')) 
#Corner3=578848.44,6201864.77
Corner3 = Corner3.strip('][').split(',')
Corner3=[float(i) for i in Corner3]

# x=endpoint, y=endpoint
Corner4 = input(colored("UTM coordinates of x=endpoint,y=endpoint corner: ",'green')) 
#Corner4=578847.11,6201837.62
Corner4 = Corner4.strip('][').split(',')
Corner4=[float(i) for i in Corner4]

UTM = input(colored("Give UTM area code as number only: ",'green')) 
UTM=str(UTM)+'K'
#UTM='32K'

plotsize_x = input(colored("Plot size along x-axis (in meters): ",'green'))
#plotsize_x=0.75
plotsize_x=float(plotsize_x)
plotsize_y = input(colored("Plot size along y-axis (in meters): ",'green')) 
#plotsize_y=1.78
plotsize_y=float(plotsize_y)
plot_path_y = input(colored("Size of paths between plots along y-axis (in meters): ",'green')) 
#plot_path_y=0.80
plot_path_y=float(plot_path_y)
Betweenparcel_path_x = input(colored("Size of paths between parcels (not plots) along x-axis (in meters): ",'green')) 
#Betweenparcel_path_x=0
Betweenparcel_path_x=float(Betweenparcel_path_x)
plotsinparcel = input(colored("Number of plots in one parcel: ",'green')) 
#plotsinparcel=2
plotsinparcel=int(plotsinparcel)

if plotsinparcel!=1:
    rowaxis = input(colored("Along which axis is a parcel divided into two plots (y/x): ",'green')) 

number_plots_x = input(colored("Number of plots along x axis: ",'green')) 
#number_plots_x=28
number_plots_x=int(number_plots_x)
number_plots_y = input(colored("Number of plots along y axis: ",'green')) 
#number_plots_y=30    
number_plots_y=int(number_plots_y)
blocks_along_y_axis = input(colored("Number of blocks along y-axis: ",'green')) 
#blocks_along_y_axis=3
blocks_along_y_axis=int(blocks_along_y_axis)
if blocks_along_y_axis!=0:
    blocks_along_y_axis=int(blocks_along_y_axis)
    no_of_plots_until_big_path_y_axis = input(colored("Number of plots until big path along y-axis: ",'green')) 
    #no_of_plots_until_big_path_y_axis=10
    no_of_plots_until_big_path_y_axis=int(no_of_plots_until_big_path_y_axis)
    rep_path = input(colored("Length of path between blocks in y-axis: ",'green')) 
    if rep_path!='NA':
        rep_path=int(rep_path)    

if blocks_along_y_axis==0:
    no_of_plots_until_big_path_y_axis=100000000000000000 #just so we never reach that limit
    rep_path=0 #true though

blocks_along_x_axis = input(colored("Number of blocks along x-axis: ",'green')) 
#blocks_along_x_axis=2
if blocks_along_x_axis!=0:
    blocks_along_x_axis=int(blocks_along_x_axis)
    no_of_plots_until_big_path_x_axis = input(colored("Number of plots until big path along x-axis: ",'green')) 
    #no_of_plots_until_big_path_x_axis=14
    no_of_plots_until_big_path_x_axis=int(no_of_plots_until_big_path_x_axis)

    xpath = input(colored("Size of path between blocks along x-axis. (in meters and NA if not known: ",'green')) 
    #xpath='NA'
    if xpath!='NA':
        xpath=int(xpath)

if blocks_along_x_axis==0:
    no_of_plots_until_big_path_x_axis=100000000000000000 #just so we never reach that limit

# Make necessary calculations to estimate big path on x-axis
Easting=(Corner1[0],Corner2[0],Corner3[0],Corner4[0])
Westing=(Corner1[1],Corner2[1],Corner3[1],Corner4[1])
df=pd.DataFrame(np.c_[Easting,Westing],columns=['Meters_East','Meters_South'])
df
myProj = Proj("+proj=utm +zone={}, +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs".format(UTM))
lon, lat = myProj(df['Meters_East'].values, df['Meters_South'].values, inverse=True)
df1=pd.DataFrame(np.c_[df, lon, lat],columns=['Meters_East','Meters_South','Longitude','Latitude'])

if xpath=='NA':
    xpath=(geopy.distance.vincenty((df1.Latitude[2],df1.Longitude[2]), (df1.Latitude[3],df1.Longitude[3])).m)-(plotsize_x*28)


# Identify Latitude axis
dist_corner1to2=np.abs(df1.Latitude[0]-df1.Latitude[1])
dist_corner1to3=np.abs(df1.Latitude[0]-df1.Latitude[2])
if dist_corner1to2>dist_corner1to3:
    Latitude_axis = 1 #x is latitude axis
    print(colored("x-axis is estimated to be the axis with changing latitude",'yellow'))
else:
    Latitude_axis = 2 #y is latitude axis
    print(colored("y-axis is estimated to be the axis with changing latitude",'yellow'))

    
# Calculate slopes of fields
if rowaxis=='x':
    #Across north-south (Latitude). Longitude should be similar
    slope_Northing_1 = np.abs(df1.Meters_East[1]-df1.Meters_East[0])/(plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*(number_plots_x-blocks_along_x_axis)*(1/plotsinparcel))
    slope_Northing_2 = np.abs(df1.Meters_East[2]-df1.Meters_East[3])/(plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*(number_plots_x-blocks_along_x_axis)*(1/plotsinparcel))
    slope_Northing_avg = (slope_Northing_1+slope_Northing_2)/2 #Pr. m

    #Across east-west (Longitude). Latitude should be similar
    slope_Easting_1 = np.abs(df1.Meters_South[0]-df1.Meters_South[2])/(plotsize_y*(number_plots_y-1)+plot_path_y*(number_plots_y-blocks_along_y_axis)+(blocks_along_y_axis-1)*(blocks_along_y_axis-1))
    slope_Easting_2 = np.abs(df1.Meters_South[1]-df1.Meters_South[3])/(plotsize_y*(number_plots_y-1)+plot_path_y*(number_plots_y-blocks_along_y_axis)+(blocks_along_y_axis-1)*(blocks_along_y_axis-1))
    slope_Easting_avg = (slope_Easting_2+slope_Easting_1)/2 #Pr. m
    
if rowaxis=='y':
    #Across north-south (Latitude). Longitude should be similar
    slope_Northing_1 = np.abs(df1.Meters_East[1]-df1.Meters_East[0])/(plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*(number_plots_x-blocks_along_x_axis))
    slope_Northing_2 = np.abs(df1.Meters_East[2]-df1.Meters_East[3])/(plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*(number_plots_x-blocks_along_x_axis))
    slope_Northing_avg = (slope_Northing_1+slope_Northing_2)/2 #Pr. m

    #Across east-west (Longitude). Latitude should be similar
    slope_Easting_1 = np.abs(df1.Meters_South[0]-df1.Meters_South[2])/(plotsize_y*(number_plots_y-1)+plot_path_y*(number_plots_y-blocks_along_y_axis)+(blocks_along_y_axis-1)*(blocks_along_y_axis-1)*(1/plotsinparcel))
    slope_Easting_2 = np.abs(df1.Meters_South[1]-df1.Meters_South[3])/(plotsize_y*(number_plots_y-1)+plot_path_y*(number_plots_y-blocks_along_y_axis)+(blocks_along_y_axis-1)*(blocks_along_y_axis-1)*(1/plotsinparcel))
    slope_Easting_avg = (slope_Easting_2+slope_Easting_1)/2 #Pr. m

# Check if calculated length and width of field fit with what is expected from GPS coordinates

# x-axis
if rowaxis=='x':
    length_x_axis_cal=plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*((number_plots_x/plotsinparcel)-blocks_along_x_axis)

if rowaxis=='y':
    length_x_axis_cal=plotsize_x*(number_plots_x-1)+xpath*(blocks_along_x_axis-1)+Betweenparcel_path_x*((number_plots_x)-blocks_along_x_axis)
  
    
xdist_1=geopy.distance.vincenty((df1.Latitude[0],df1.Longitude[0]), (df1.Latitude[1],df1.Longitude[1])).m
xdist_2=geopy.distance.vincenty((df1.Latitude[2],df1.Longitude[2]), (df1.Latitude[3],df1.Longitude[3])).m
xdist_avg=(xdist_1+xdist_2)/2

difference_x=np.abs(length_x_axis_cal-xdist_avg)
if difference_x>plotsize_x*2:
    print(colored("length of x-axis do calculated from GPS coordinates and from plot size and paths differ more than two plot sizes",'yellow'))
    print(colored("The distance calculated from plot sizes and paths will be used",'yellow'))


# y-axis

if blocks_along_y_axis!=0:
    length_y_axis_cal=plotsize_y*(number_plots_y-1*plotsinparcel)+plot_path_y*(number_plots_y-1*plotsinparcel-blocks_along_y_axis)+(blocks_along_y_axis-1)*rep_path

if blocks_along_y_axis==0:
    length_y_axis_cal=plotsize_y*(number_plots_y-1*plotsinparcel)+plot_path_y*((number_plots_y-1*plotsinparcel)/2-blocks_along_y_axis)
    
ydist_1=geopy.distance.vincenty((df1.Latitude[0],df1.Longitude[0]), (df1.Latitude[2],df1.Longitude[2])).m
ydist_2=geopy.distance.vincenty((df1.Latitude[1],df1.Longitude[1]), (df1.Latitude[3],df1.Longitude[3])).m
ydist_avg=(ydist_1+ydist_2)/2
difference_y=np.abs(length_y_axis_cal-ydist_avg)
if difference_y>plotsize_y*2:
    print(colored("length of y-axis do calculated from GPS coordinates and from plot size and paths differ more than two plot sizes",'yellow'))
    print(colored("The distance calculated from plot sizes and paths will be used",'yellow'))
    


        
# now calculate Latitudes for all 'fixed' Longitudes

if Latitude_axis==1:
    plot_numbers=number_plots_y


All_gpscoordinates=list()
y=calculate_y_fromfixedx(df1.Meters_East[0],df1.Meters_South[0])
for i in range(0,(plot_numbers)): 
    x_values=calculate_x_fromfixedy(y.Meters_East[i],y.Meters_North[i],i)
    All_gpscoordinates.append(x_values) #decimals rounded away

Total = pd.concat(All_gpscoordinates)



# Convert UTM coordinates to longitude and latitude, and make csv
myProj = Proj(("+proj=utm +zone={}, +south +ellps=WGS84 +datum=WGS84 +units=m +no_defs").format(UTM))
lon, lat = myProj(Total['Meters_East'].values, Total['Meters_North'].values, inverse=True)
Total2=pd.DataFrame(np.c_[lon, lat],columns=['Longitude','Latitude'])
Total3=pd.DataFrame(np.c_[Total, lon, lat],columns=['x','y','Meters_East','Meters_North','Longitude','Latitude'])
Total3.to_csv(('{}.csv').format(Name),index=False)


# Last plot all gps coordinates in an axis and perhaps on a map if possible
plt.plot( 'Longitude', 'Latitude', data=Total3, linestyle='none', marker='o',color='green',markersize=2.5)
plt.xlabel('Longitude')
plt.ylabel('Latitude')
plt.text(Total3.Longitude[0],Total3.Latitude[0],'start (1,1)')
plt.text(Total3.Longitude[len(Total3)-1],Total3.Latitude[len(Total3)-1],'end')
filename=('FieldPlan' + '_' + Name)
plt.savefig(filename+'.png',dpi=500)

# save googlemaps based version of field plan
gmap4 = gmplot.GoogleMapPlotter(Total3.Latitude[0],Total3.Longitude[0],18) 
# points on the Google map 
gmap4.heatmap( Total3.Latitude , Total3.Longitude ) 
gmap4.draw( filename+ ".html" )          








