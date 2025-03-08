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
library(fastcluster)
library("rnaturalearth")
library(rnaturalearthhires)
#devtools::install_github("ropensci/rnaturalearthhires")
#devtools::install_github("ropensci/rnaturalearthdata")
#devtools::install_github("ropensci/rnaturalearthhires")

#require(rnaturalearthhires)
#remotes::install_github("ropensci/rnaturalearthhires")
#devtools::install_github("ropensci/rnaturalearthhires")


Data = readRDS("dfa_data_origplus_fullcoast_Feb2025.rds") 

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
  saveRDS( domain, file="domain.RDS" )

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
  par(mfrow=c(2,1), oma=c(4,2,0,0), mar=c(3,3,4,1), mgp=c(3,0.5,0), tck=-0.02) #mar B,L,T,R

  plotgrid = st_sf( sf_grid, Class=Class_g[1:length(sf_grid)] )
  plot( plotgrid, reset=FALSE, key.pos=NULL, border=NA, pal=viridis, main="" )
  plot( Land, add=TRUE, col="grey" )
  axis(side=2, cex.axis=0.5, font.axis=1)
  axis(side=1, cex.axis=0.5, font.axis=1)
  box()
  cex.lab=3
  mtext( side=2, text="Latitude (⁰N)", line=2, cex=0.5 )
  mtext( side=1, text="Longitude (⁰W)", line=2, cex=0.5)

  ##
  par(mar=c(3,1,3.5,1)) #mar B,L,T,R
  
  matplot( y=pmean_gz, xlab="", xaxt="n",yaxt="n", type="l", lwd=2, col=viridis(k), lty="solid" )

  axis(1, at=1:nrow(pmean_gz), labels=rnames, las=3, cex.axis=0.5,font.axis=1 )
  axis(2,cex.axis=0.5,font.axis=1 )
  mtext( side=2, text="Average log-density", line=2,cex=0.5 )
#dev.off()
