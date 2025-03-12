#Groundfish Community Analysis

#library(devtools)
#install_github("vast-lib/tinyVAST", dependencies = TRUE)
library(dplyr)
library(sp)
library(tinyVAST)
library(fmesher)
library(ggplot2)
library(sf)
library(viridisLite)
library(viridis)
library(fastcluster)
library("rnaturalearth")
library(rnaturalearthhires)
#devtools::install_github("ropensci/rnaturalearthhires")
#devtools::install_github("ropensci/rnaturalearthdata")
#devtools::install_github("ropensci/rnaturalearthhires")

#require(rnaturalearthhires)
#remotes::install_github("ropensci/rnaturalearthhires")
#devtools::install_github("ropensci/rnaturalearthhires")

Data = readRDS("groundfish_community_study/dfa_data_origplus_fullcoast_Feb2025.RDS")
#Data = readRDS("dfa_data_origplus_fullcoast_Feb2025.RDS")

# Subset data

Data = subset( Data, year >= 2003 )
Data = subset( Data, year <= 2024 )
min(Data$date)
max(Data$date)


unique(Data$common_name)
Data = subset( Data, common_name == "pacific ocean perch"| common_name == "arrowtooth flounder"| common_name == "pacific halibut"  | common_name == "sablefish" | common_name == "pacific cod" | common_name == "english sole" | common_name == "flathead sole"| common_name == "rex sole"| common_name == "dover sole" | common_name == "lingcod" )#common_name == "yelloweye rockfish" | common_name == "shortraker rockfish" | )


#
tapply( Data$catch_weight, INDEX=list(Data$year, Data$common_name), FUN=\(x) mean(x>0) )
tapply( Data$catch_weight, INDEX=list(Data$year, Data$survey_name), FUN=length )

#
Data0 = na.omit( as.data.frame(Data)[,c('year','common_name','catch_weight','lat_start','lon_start')])
# Replace spaces (not allowed in SEM notation)
Data0$common_name = sapply(Data0$common_name, FUN=gsub, pattern=" ", replace="_", fixed=TRUE )

#CREATE map
  sf_locs = st_as_sf(Data0[,c('lon_start','lat_start')], coords=1:2)
  sf_polys = st_buffer( sf_locs, dist=0.2 )
  domain = st_union( sf_polys )
  saveRDS( domain, file="groundfish_community_study/domain.RDS" )

#
Land = ne_countries( scale=10, continent = "North America")

#
sf_grid = st_make_grid( domain, cellsize=0.1 )
sf_grid = st_intersection( sf_grid, domain )
grid = st_coordinates(st_centroid( sf_grid ))

#
graph = fm_mesh_2d( loc=domain, cutoff=0.5, refine=TRUE )
#graph = fm_mesh_2d( loc=Data0[,c('lon_start','lat_start')], cutoff=0.5, refine=TRUE )

#
DF = expand.grid( "grid"=1:nrow(grid), "common_name"=unique(Data0$common_name) )
DF = cbind( DF, 'lon_start'=grid[DF$grid,'X'], 'lat_start'=grid[DF$grid,'Y'] )

family = list()
for(i in seq_along(unique(Data0$common_name)) ){
  family[[i]] = tweedie()
}

names(family) = unique(Data0$common_name)

  fit = tinyVAST(
    formula = catch_weight ~ 0 + common_name, # interaction(common_name,year),
    data = Data0,
    control = tinyVASTcontrol( profile = "alpha_j",
                               trace = 1,
                               getsd = FALSE ),

      family = list(
        "pacific_ocean_perch" = tweedie(),
       #  "walleye_pollock" = tweedie(),
            "arrowtooth_flounder" = tweedie(),
       "pacific_cod"= tweedie(),
            "pacific_halibut" = tweedie(),
            "rex_sole"= tweedie(),
           "flathead_sole"= tweedie(),
       #    "northern_rock_sole"= tweedie(),
          "english_sole"= tweedie(),
         "sablefish" = tweedie(),
       "dover_sole" = tweedie(),
       "lingcod" = tweedie()
    ),
    sem = "",
    spatial_graph = graph,
    variable_column = "common_name",
    space_columns = c('lon_start','lat_start'),
    distribution_column = "common_name"
   )

  # Ordination
  dens = predict( fit, newdata=DF, what="p_g" )
 dens = matrix(dens, dimnames=list(NULL,"name")) #added "name" as placeholder to replace early and late time period names

  DF = cbind( DF, dens )

  Y_gc = matrix( DF$name, nrow=nrow(grid), ncol=length(unique(Data0$common_name)),
                 dimnames = list(NULL,unique(Data0$common_name)) )

#
Hclust = hclust.vector( Y_gc, method="ward" )

# Plot
k=3
Class_g = cutree( Hclust, k = k )
Class_gz = model.matrix( ~ 0 + factor(Class_g) )
pmean_gz = (t(Y_gc) %*% Class_gz) / nrow(Y_gc)

#name for plot
rnames=c("Arrowtooth fl.", "P.halibut", "Flathead sole", "P.cod", "Rex sole", "P.ocean perch",  "Sablefish", "Dover sole", "English sole","Lingcod")

#png( "Orig.modplusK3.Feb25.png", width=4, height=5, res=900, units="in" )
#   par(mfrow=c(2,1), oma=c(4,2,0,0), mar=c(3,3,4,1), mgp=c(3,0.5,0), tck=-0.02) #mar B,L,T,R
#
#   plotgrid = st_sf( sf_grid, Class=Class_g[1:length(sf_grid)] )
#   plot( plotgrid, reset=FALSE, key.pos=NULL, border=NA, pal=viridis, main="" )
#   plot( Land, add=TRUE, col="grey" )
#   axis(side=2, cex.axis=0.5, font.axis=1)
#   axis(side=1, cex.axis=0.5, font.axis=1)
#   box()
#   cex.lab=3
#   mtext( side=2, text="Latitude (⁰N)", line=2, cex=0.5 )
#   mtext( side=1, text="Longitude (⁰W)", line=2, cex=0.5)
#
#   ##
#   par(mar=c(3,1,3.5,1)) #mar B,L,T,R
#
#   matplot( y=pmean_gz, xlab="", xaxt="n",yaxt="n", type="l", lwd=2, col=viridis(k), lty="solid" )
#
#   axis(1, at=1:nrow(pmean_gz), labels=rnames, las=3, cex.axis=0.5,font.axis=1 )
#   axis(2,cex.axis=0.5,font.axis=1 )
#   mtext( side=2, text="Average log-density", line=2,cex=0.5 )
# #dev.off()
rnames=c("Arrowtooth fl.", "P.halibut", "Flathead sole", "P.cod", "Rex sole", "P.ocean perch",  "Sablefish", "Dover sole", "English sole","Lingcod")


sf_grid <- readRDS("sf_grid.rds")
plotgrid <- readRDS("plotgrid.rds")
Class_g <- readRDS("Class_g.rds")
  # This first plot is the map
  # make sure sf_grid has a CRS
  if (is.na(st_crs(sf_grid))) st_crs(sf_grid) <- 4326  # Assign EPSG:4326 if it's not set

  # Transform sf_grid to EPSG:3338
  sf_grid <- st_transform(sf_grid, crs = 3338)

  # Ensure Land has CRS
  #if (is.na(st_crs(Land))) {
  #  st_crs(Land) <- 4326  # if missing
  #} else {
  #  Land <- st_transform(Land, 4326)  # If it's already set, just transform it to EPSG:4326
  #}

  # transform Land to EPSG:3338
  #Land <- st_transform(Land, crs = 3338)

  # is assigned to sf_grid as before
  plotgrid <- st_sf(sf_grid, Class = Class_g[1:length(sf_grid)])

  # Define the place labels (you can adjust the lat/lon coordinates as needed)
  place_labels <- data.frame(
    type = c("mainland", "mainland", "mainland", "mainland", "survey",
             "peninsula", "survey", "survey", "survey", "survey", "survey"),
    lab = c("Alaska", "Russia", "Canada", "USA", "West Coast",
            "Alaska Peninsula", "Aleutian Islands", "Gulf of Alaska",
            "Bering\nSea\nSlope", "Eastern\nBering Sea", "Northern\nBering Sea"),
    angle = c(0, 0, 0, 0, -45, 45, 0, 30, 0, 0, 0),
    lat = c(63, 61.798276, 58, 42, 40, 56.352495, 53.25, 54.720787,
            57, 57.456912, 62.25),
    lon = c(-154, 173.205231, -122, -120, -125.2,
            -159.029430, -173, -154.794131, -176, -162, -170.5)
  ) %>%
    dplyr::filter(type != "peninsula") %>%
    sf::st_as_sf(coords = c("lon", "lat"),
                 crs = "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0") %>%
    sf::st_transform(crs = 3338)  # Transform labels to the same CRS (EPSG:3338)

  # I don't know we need the bbox_geo -- this was pulled in from some other plots
  bbox_geo <- c(xmin = -175, ymin = 32, xmax = -115, ymax = 60)
  # Create a bounding box in geographic coordinates (EPSG:4326)
  bbox_sf <- st_sfc(
    st_polygon(list(matrix(c(bbox_geo["xmin"], bbox_geo["ymin"],
                             bbox_geo["xmin"], bbox_geo["ymax"],
                             bbox_geo["xmax"], bbox_geo["ymax"],
                             bbox_geo["xmax"], bbox_geo["ymin"],
                             bbox_geo["xmin"], bbox_geo["ymin"]),
                           ncol = 2, byrow = TRUE))),
    crs = 4326
  )

  # Transform bounding box to EPSG:3338
  bbox_sf_projected <- st_transform(bbox_sf, crs = 3338)

  # Get the transformed bounding box extents in EPSG:3338
  bbox_projected <- st_bbox(bbox_sf_projected)
  bbox_projected$ymax <- bbox_projected$ymax + 250000 # this was in orig plot

  # Convert to factor
  plotgrid$Community <- as.factor(plotgrid$Class)

  map_data <- rnaturalearth::ne_countries(
    scale = "medium",
    returnclass = "sf", country = c("united states of america","canada","mexico"))
  coast <- suppressWarnings(suppressMessages(
    sf::st_crop(map_data,
                c(xmin = -179, ymin = 15, xmax = -80, ymax = 70))))
  coast_proj <- sf::st_transform(coast, crs = "EPSG:3338")

  p1 <- ggplot(coast_proj) +
    geom_sf(fill="grey10", col = "grey20") +
    geom_sf(data = plotgrid, aes(fill = Community), color = NA) +
    scale_fill_viridis_d(name = "Community") +
    theme_bw() +
    labs(x = "Longitude (\u00B0W)", y = "Latitude (\u00B0N)") +
    scale_x_continuous(breaks = seq(-180, 180, by = 10)) +  # Longitude gridlines every 10 degrees
    scale_y_continuous(breaks = seq(-90, 90, by = 5)) +  # Latitude gridlines every 5 degrees
    coord_sf(xlim = c(bbox_projected$xmin, bbox_projected$xmax),  # Apply the bounding box
             ylim = c(bbox_projected$ymin, bbox_projected$ymax))  # Apply the bounding box

  # Add place labels to p1 (same as before)
  p1 <- p1 + ggplot2::geom_sf_text(
    data = place_labels %>% dplyr::filter(type == "mainland", lab != "Russia"),
    mapping = aes(label = lab, angle = angle),
    color = "grey60",
    size = 3,
    show.legend = FALSE
  ) +
    coord_sf(
      xlim = c(bbox_projected$xmin, bbox_projected$xmax),
      ylim = c(bbox_projected$ymin, bbox_projected$ymax),
      expand = FALSE  #no additional padding around the plot
    ) +
    scale_x_continuous(breaks = seq(-180, 180, by = 10)) +  # Longitude gridlines every 10 degrees
    scale_y_continuous(breaks = seq(-90, 90, by = 5))  # Latitude gridlines every 5 degrees


  # Second plot here,
  # Convert pmean_gz into a data frame with 'x' as the row index and 'y' as the data
  df <- data.frame(
    x = 1:nrow(pmean_gz),
    y = as.vector(pmean_gz),  # Flatten pmean_gz for y-axis
    group = rep(1:ncol(pmean_gz), each = nrow(pmean_gz))  # Create a group for each column of pmean_gz
  )

  # Create a ggplot line plot
  p2 <- ggplot(df, aes(x = x, y = y, group = group, color = factor(group))) +
    geom_line(lwd = 1.3) +  # Line width and color
    scale_color_viridis(discrete = TRUE, name = "Community") +
    theme_bw() +
    labs(x = NULL,
      y = "Average log-density"
    ) +
    scale_x_continuous(
      breaks = 1:nrow(pmean_gz),  # x-axis ticks at each row of pmean_gz
      labels = rnames,  # row names as x-axis labels
      expand = c(0, 0)  # remove extra space on the left and right
    ) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1)  # Rotate x-axis labels
    )


  library(cowplot)
  p3 <- plot_grid(p1,p2,
                  ncol = 1,
                  align = "v",
                  axis = "lr",
                  rel_heights = c(2,1))

  ggsave(p3, filename="groundfish_community_study/combined_plot.png",height=7,width=6, bg = "white")
