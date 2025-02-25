
library(sf)
library(stringr)
library(matrixStats)
library(dplyr)
library(units)
library(lubridate)
library(terra)
library(tidyterra)
source("C:/Users/Eliza.Carter/Documents/Projects/California/ORI/code/ais_processing_functions.R") # make sure path is correct
fname <- basename(rstudioapi::getSourceEditorContext()$path) #get the name of the current file
start0 <- Sys.time()

###########
# years to loop through
years = c(2023)

###########
## settings used in all iterations
# dsn of study area polygons
fp = "C:/Users/Eliza.Carter/Documents/Projects/California/ORI/data/ORI_Area_of_Interest.shp"
study_areas = st_read(fp, quiet = TRUE) %>%
  st_transform(., crs="EPSG:26911")

# grid resolution of output raster
output_resolution = 100 # meters

###
# output directory
out_dir = "C:/Users/Eliza.Carter/Documents/Projects/California/ORI/data/source_data"  # make sure this is right
#
dir.create(file.path(str_remove(out_dir, basename(out_dir))), showWarnings = FALSE) # create extracted_rasters
dir.create(file.path(out_dir), showWarnings = FALSE) # create monthly_multiband_rasters
if(!(file.exists(paste0(out_dir,"log_file1.txt")))){ # copy log file template if it doesn't exist
  file.copy("log_file_template.txt",paste0(out_dir,"log_file.txt"))
}

###########

# dsn of merged lines for the relevant year
dsn_lines = "C:/Users/Eliza.Carter/Documents/Projects/California/ORI/data/source_data/ori_ais_transects.shp"   # make sure this is right

study_area_name = "ori"
study_area = study_areas
# create buffered (by 2 resolutions) bounding box around the current study area
bbox = st_bbox(study_area) + c(-2, -2, 2, 2)*output_resolution

      
########
# pull the relevant lines datasets
lines_month = list() # empty spatial vector collection
        # read in lines filtered to current time range
line_sf = st_read(dsn_lines, quiet = TRUE)
lines_month = c(lines_month, list(line_sf)) # add to list
      
      # combine the line features and crop to the study area
lines_month = st_crop(do.call(rbind, lines_month), bbox) 
      # add the ship type category as a field
lines_month$ShipTypeCat <- categorize_ship_types(lines_month$vessel_typ) # from ais_processing_functions.R
########
      

#######
# initialize empty raster for bounding box and resolution
empty_ras <- rast(ext(bbox), resolution=output_resolution)
crs(empty_ras) <- crs(lines_month) # set the crs
      
# read in the hex grid
fishnet_polygon <- st_read("C:/Users/Eliza.Carter/Documents/Projects/California/ORI/data/ORI_Area_of_Interest_grid.shp")
fishnet_polygon$fn_id = seq(nrow(fishnet_polygon)) # add an id field
# create a version of the fishnet that is a line for each grid instead of polygons
fishnet_boundary_lines = st_cast(fishnet_polygon, "LINESTRING")

#######
      
# Extract the intersecting points between the lines and each grid cell in the fishnet 
# Do this seperately for each ship type
ship_types = c("Cargo", "Tug_Tow", "Fishing",        
                     "Passenger", "Pleasure", "Other")
      
      for(j in seq(length(ship_types))){
        # get just the lines for this ship (vessel) type
        vessel_lines <- lines_month %>%
          filter(ShipTypeCat == ship_types[j])
        
        ###
        # generate points for all of the intersections between ais lines and the grid cell boundarys
        line_intersect_points = st_intersection(fishnet_boundary_lines, vessel_lines) %>%
          st_cast(., "MULTIPOINT") %>%
          st_cast(., "POINT")
        
        # count the number of point intersections for each fishnet grid (using the fishnet id)
        intsct_data = st_drop_geometry(line_intersect_points) %>%
          group_by(fn_id) %>%
          summarise(!!ship_types[j] := n()) # set column name to ship type
        
        # join the counts for this ship type back to the fishnet polygon
        fishnet_polygon = fishnet_polygon %>%
          left_join(intsct_data, by=c("fn_id"))
        
        ###
        
        #####
        ### count the number of line endpoints in each grid cell
        # 1st step: get points for transect ends, crop out ones that aren't true 
        #  transect ends on the study area boundary
        if(nrow(vessel_lines)>0){ # check there is a vessel line
          line_end_pts <- st_coordinates(st_cast(vessel_lines, "MULTILINESTRING")) %>% # get coordinates of line vertices
            as.data.frame(.) %>% # convert from matrix to data frame
            group_by(c(L2)) %>%  # group by the origninal line it came from
            slice(c(1,n())) %>%  # get the first and last point of the linestring
            st_as_sf(., coords=c("X", "Y"), crs=st_crs(vessel_lines)) %>% # convert to points
            st_crop(., st_bbox(study_area))  # crop to the study area (cuts out edge touching buffered bbox)
          
          
          # get a raster with the count of each endpoints in each cell in the same
          #  multiband format for each ship type
          if(j==1){ trnsct_end_raster <- rast(rasterize(line_end_pts, empty_ras, fun="count"), names=ship_types[j])
          }else{ trnsct_end_raster <- c(trnsct_end_raster,
                                        rast(rasterize(line_end_pts, empty_ras, fun="count"), names=ship_types[j])) }
          
        }else{ # if no vessel line
          if(j==1){ trnsct_end_raster <- empty_ras
          }else{ trnsct_end_raster <- c(trnsct_end_raster, empty_ras) }
        }
        #####
        
      }
######

#write fishnet_polygon table to shapefile
st_write(fishnet_polygon, "C:/Users/Eliza.Carter/Documents/Projects/California/ORI/data/input_data/ORI_Industry.gpkg", layer = "ORI_AIS_Grid_Intersect", delete_layer = TRUE)

######
# convert the polygon to raster object, applied to each ship type
rast_output = lapply(ship_types,
          function(x){rasterize(fishnet_polygon,
                     empty_ras,
                     field=x)})
# convert the list of separate raster layers into one multiband layer
rast_output = rast(rast_output)
      
# add the end point counts and replace NAs
rast_output = subst(rast_output, NA, 0) + subst(trnsct_end_raster, NA, 0) 
      
######
# output file name
oname = paste(out_dir, "ori_ais.tif", sep="/")
# write raster to the output folder
writeRaster(rast_output, oname, overwrite=TRUE)
######

#############
      
end <- Sys.time()
print(end - start) # print runtime for this iteration
      

print("full loop ended")
print(Sys.time() - start0)


