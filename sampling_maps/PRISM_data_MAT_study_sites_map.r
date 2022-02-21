require(devtools)
require(reshape2)
require(dplyr)
require(raster)

require(ggplot2)
require(ggmap)

#install_github(repo = "ropensci/prism")
require(prism)

source("~/ggplot_theme.txt")

prism_set_dl_dir("/Users/ericmorrison/PRISM_data")
#options(prism.path = "/Users/ericmorrison/PRISM_data")
#get_prism_normals(type = 'ppt', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmax', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmean', resolution = '4km', annual = T, keepZip = TRUE)
#get_prism_normals(type = 'tmin', resolution = '4km', annual = T, keepZip = TRUE)


#Convert raster to point data

new_file<-3#this number corresponds to the row of the file of interest
RS <- pd_stack(prism_archive_ls()[3]) ##raster file of data
proj4string(RS)<-CRS("+proj=longlat +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs") ##assign projection info

##convert raster to point data frame
df <- data.frame(rasterToPoints(RS))
m.df <- melt(df, c("x", "y"))
names(m.df)[1:2] <- c("lon", "lat") #rename columns


site_coords = read.table("~/GARNAS_neonectria_genome_reseq_10072020/sample_metadata/site_coords.txt", header = T)

minLat=min(site_coords$lat)-2
maxLat=max(site_coords$lat)+2
minLon=min(site_coords$lon)-2
maxLon=max(site_coords$lon)+2


m.df.study_area<-m.df%>%filter(minLat < lat, lat < maxLat, minLon < lon, lon <maxLon)%>%
mutate(MAT = value)%>%
dplyr::select(-value)

min_max = m.df.study_area$MAT %>% range
#min_max[2] = min_max[2]+2
mid_point = (min_max[1]+min_max[2])/2

p1 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT))+
#geom_point(data=site_coords, aes(x=lon, y = lat), color = "black") +
labs(x = "longitude", y = "latitude") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2(expression("Mean\nannual\ntemp. ("*degree*C*")"),
    low='darkslateblue',
    mid="#d1e5f0",
    high = 'red',
    midpoint=mid_point,
    limits = c(min_max[1], min_max[2])
)+
theme(
legend.title = element_text(size = 25),
legend.text = element_text(size = 25),
axis.text = element_blank(),
axis.ticks = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
legend.position = c(0.825,0.225)
)

p2 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT+2))+
#geom_point(data=site_coords, aes(x=lon, y = lat), color = "black") +
labs(x = "longitude", y = "latitude") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2(expression("Mean\nannual\ntemp. ("*degree*C*")"),
low='darkslateblue',
mid="#d1e5f0",
high = 'red',
midpoint=mid_point,
limits = c(min_max[1], min_max[2])
)+
theme(
legend.title = element_text(size = 25),
legend.text = element_text(size = 25),
axis.text = element_blank(),
axis.ticks = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
legend.position = c(0.825,0.225)
)

p3 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT))+
geom_point(data=site_coords, aes(x=lon, y = lat, color = core), size = 4) +
labs(x = "longitude", y = "latitude") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2(expression("Mean\nannual\ntemp. ("*degree*C*")"),
low='darkslateblue',
mid="#d1e5f0",
high = 'red',
midpoint=mid_point,
limits = c(min_max[1], min_max[2])
)+
theme(
legend.title = element_text(size = 25),
legend.text = element_text(size = 25),
axis.text = element_blank(),
axis.ticks = element_blank(),
axis.title.x = element_blank(),
axis.title.y = element_blank(),
legend.position = c(0.825,0.225)
) +
scale_color_manual(values = c("y" = "black", "n" = "grey45"), guide = F)

min_max = m.df.study_area$MAT %>% range
min_max[2] = 15
mid_point = (min_max[1]+min_max[2])/2

p4 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT))+
geom_point(data=site_coords %>% filter(state.name != "WI" & state.name != "NH.HUBB"), aes(x=lon, y = lat, color = sample), size = 3) +
labs(x = "longitude", y = "latitude") +
my_gg_theme +
theme(
legend.title = element_text(size = 20)
) +
scale_fill_gradient2(expression("MAT ("*degree*C*")"),
    low="#053061",
    mid="#d1e5f0",
    high = "#b2182b",
    midpoint=mid_point,
    limits = c(min_max[1], min_max[2])
)+
scale_color_manual("Sample status",
    values = c("y" = "black", "c" = "white"),
    labels = c("y" = "sequenced", "c" = "collected"),
) +
theme(
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 25),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    legend.position = c(0.825,0.225),
    legend.key = element_rect(fill = "grey")
)

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

pdf("figures/MAT_site_map.pdf", height = plot_height, width = plot_width)
print(p3)
dev.off()

pdf("figures/MAT_site_map.potential_sites.pdf", height = plot_height*1.2, width = plot_width*1.2)
print(p4)
dev.off()

#big text

p1 = ggplot()+
geom_raster(data=m.df.study_area, aes(x=lon, y=lat, fill=MAT))+
geom_point(data=site_coords, aes(x=lon, y = lat, color = sample), size = 5) +
my_gg_theme +
labs(y = NULL, x = NULL) +
theme(
legend.title = element_text(size = 35),
legend.text = element_text(size = 30),
axis.text = element_blank(),
axis.ticks = element_blank(),
legend.position = c(0.825,0.225)
) +
scale_fill_gradient2(expression("MAT ("*degree*C*")"), low='darkslateblue',mid='lightblue',high = 'red',midpoint=10) +
scale_color_manual(values = c("y" = "black", "n" = "dark grey"), guide = F)

plot_height = (maxLat-minLat)/2
plot_width = (maxLon-minLon)/2

pdf("PRISM_maps/site_MAT_map.big_text.pdf", height = plot_height, width = plot_width)
print(p1)
print(p2)
dev.off()

