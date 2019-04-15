# load relevant packages

library(intersectr)
library(dplyr)
library(daymetr)
library(rgdal)
library(sf)
library(RNetCDF)

# https://usgs-r.github.io/intersectr/dev/articles/intersectr.html

# output <- tile_outlines

# writeOGR(obj=output,dsn = "D:\\WBEEP", layer = "daymettiles", driver = "ESRI Shapefile")

# data <- download_daymet_ncss(location = c(65.261,-147.185,64.188,-142.928), start = 2000, end = 2018, param = "prcp" )

# set working directory

setwd("D:\\WBEEP")

# Set an analysis projection to do intersection in, projection chosen is NAD 83 Alaska from http://spatialreference.org/ref/epsg/
analysis_projection <- "+init=epsg:3338"

# Use sample daymet precipitation file for 2005 using daymetr package above

nc <- open.nc(con = "prcp_daily_2005_ncss.nc")

# Use Salcha HUC12 geometry 
geom <- sf::read_sf(dsn = "D:\\WBEEP", layer = "SalchaHUC12") %>%
  st_transform(analysis_projection)

# Look at what variables are available.
(nc_var <- ncmeta::nc_vars(nc))
#> # A tibble: 5 x 5
#>      id name                    type      ndims natts
#>   <dbl> <chr>                   <chr>     <dbl> <dbl>
#> 1     0 prcp                    NC_FLOAT      3     9
#> 2     1 time                    NC_DOUBLE     1     6
#> 3     2 y                       NC_FLOAT      1     2
#> 4     3 x                       NC_FLOAT      1     2
#> 5     4 lambert_conformal_conic NC_SHORT      0    11

# Choose a variable and get coordinate variables.
variable_name <- "prcp"
(nc_coord_vars <- ncmeta::nc_coord_var(nc, variable_name))
#> Warning in nc_coord_var_finder(dim, var, att, axe, variable): missing
#> coordinate variables names in coordinates attribute trying to find non-
#> auxiliary coordinate variables.
#> # A tibble: 1 x 6
#>   variable X     Y     Z     T     bounds
#>   <chr>    <chr> <chr> <chr> <chr> <chr> 
#> 1 prcp     x     y     <NA>  time  <NA>

# Pull out the coordinate variables and get them in projection units.
meters_per_km <- 1000
x <- var.get.nc(nc, nc_coord_vars$X, unpack = TRUE) * meters_per_km
y <- var.get.nc(nc, nc_coord_vars$Y, unpack = TRUE) * meters_per_km

# Create cell geometry for this geom
cell_geometry <- create_cell_geometry(X_coords = x, 
                                      Y_coords = y, 
                                      prj = ncmeta::nc_gm_to_prj(ncmeta::nc_grid_mapping_atts(nc)), 
                                      geom = geom, 
                                      buffer_dist = 1000)

# Be explicit about the area intersection geometry and ids.
data_source_cells <- st_sf(select(cell_geometry, grid_ids))
target_polygons <- st_sf(select(geom, HUC12))
st_agr(data_source_cells) <- "constant"
st_agr(target_polygons) <- "constant"

# Create the areal weights
area_weights <- 
  calculate_area_intersection_weights(data_source_cells, target_polygons)

# Run the intersection
intersected <- execute_intersection(nc_file = "prcp_daily_2005_ncss.nc",
                                    variable_name = variable_name,
                                    intersection_weights = area_weights,
                                    cell_geometry = cell_geometry,
                                    x_var = nc_coord_vars$X, 
                                    y_var = nc_coord_vars$Y, 
                                    t_var = nc_coord_vars$T)

# Get data ready to create plots.change third number to get layer of p...
X_inds <- seq(min(cell_geometry$X_ind), max(cell_geometry$X_ind), 1)
Y_inds <- seq(min(cell_geometry$Y_ind), max(cell_geometry$Y_ind), 1)

ids <- intersectr:::get_ids(length(X_inds), length(Y_inds))

grid_data <- var.get.nc(nc, variable_name,
                        start = c(min(X_inds), min(Y_inds), 3),
                        count = c(length(X_inds), length(Y_inds), 1), 
                        unpack = TRUE)

grid_data <- data.frame(grid_data = matrix(grid_data,
                                           ncol = 1,
                                           byrow = TRUE),
                        grid_ids = matrix(ids, ncol = 1))

grid_data$grid_data[grid_data$grid_data < 0] <- NA

grid_data <- left_join(cell_geometry, grid_data, by = "grid_ids")

geom$HUC12 <- as.numeric(geom$HUC12)


geom_data <- select(geom, HUC12) %>%
  left_join(data.frame(HUC12 = as.numeric(names(intersected)[2:ncol(intersected)]),
                       poly_data = as.numeric(intersected[170, 2:ncol(intersected)])), by = "HUC12")

geom_data <- st_transform(geom_data, st_crs(grid_data))

breaks <- seq(0,40,1)
plot(grid_data["grid_data"], border = NA, breaks = breaks)
plot(geom_data["poly_data"], breaks = breaks)


plot(grid_data$geometry)
plot(geom_data["poly_data"], breaks = breaks, add = TRUE)
