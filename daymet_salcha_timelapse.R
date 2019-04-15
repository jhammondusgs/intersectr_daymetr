# timelapse daymet HUC 12 Salcha

# load relevant packages

library(intersectr)
library(dplyr)
library(daymetr)
library(rgdal)
library(sf)
library(RNetCDF)
library(stringr)

setwd("D:\\WBEEP")

analysis_projection <- "+init=epsg:3338"

for(year in 2000:2018){
  # year = 2000
nc <- open.nc(con = paste("prcp_daily_",year,"_ncss.nc",sep = ""))

geom <- sf::read_sf(dsn = "D:\\WBEEP", layer = "SalchaHUC12") %>%
  st_transform(analysis_projection)

(nc_var <- ncmeta::nc_vars(nc))

variable_name <- "prcp"
(nc_coord_vars <- ncmeta::nc_coord_var(nc, variable_name))

meters_per_km <- 1000
x <- var.get.nc(nc, nc_coord_vars$X, unpack = TRUE) * meters_per_km
y <- var.get.nc(nc, nc_coord_vars$Y, unpack = TRUE) * meters_per_km

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
intersected <- execute_intersection(nc_file = paste("prcp_daily_",year,"_ncss.nc",sep = ""),
                                    variable_name = variable_name,
                                    intersection_weights = area_weights,
                                    cell_geometry = cell_geometry,
                                    x_var = nc_coord_vars$X, 
                                    y_var = nc_coord_vars$Y, 
                                    t_var = nc_coord_vars$T)

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

breaks <- seq(0,20,1)

    for(day in 1:365){
      # day = 20
geom_data <- select(geom, HUC12) %>%
  left_join(data.frame(HUC12 = as.numeric(names(intersected)[2:ncol(intersected)]),
                       poly_data = as.numeric(intersected[day, 2:ncol(intersected)])), by = "HUC12")

doy <- str_pad(day, 3, pad = "0")
date <- intersected$time_stamp[day]
date <- substr(date,1,10)

jpeg(filename = paste(year,doy,"daymetP.jpeg",sep = ""), width = 8, height = 6, units = "in",res = 100)
plot(geom_data["poly_data"], breaks = breaks, main = date, pal = rev(topo.colors(20)), key.pos = NULL)
dev.off()

}
}

