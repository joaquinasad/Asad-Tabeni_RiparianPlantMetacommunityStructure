#===============================================================================
# ESTIMATION OF HABITAT FRAGMENTATION THROUGH THE LANDSCAPE DIVISON INDEX (LDI)
#===============================================================================

setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/FiltrosAmbientalesUrbanos/Fragmentacion")

library(terra)
citation("terra")
library(landscapemetrics)
citation("landscapemetrics")
library(tools)
library(dplyr)
library(sf)

citation()

#Three layers
r1 = rast("Antropic_CoverMapBiomas.tif")  # binary: 0 = natural, 1 = antrópc
r2 = rast("Area_Impermeabilizada_Raster.tif")  # 1 = antropic, NA = natural
r3 = rast("Open_Buildings_Raster.tif") # 1 = antropic, NA = natural

# align resolution and extention
r2 = resample(r2, r1)
r3 = resample(r3, r1)

# reeplace NAs for 0 to make all layers binarys
r_comb = (ifel(is.na(r1), 0, r1) +
             ifel(is.na(r2), 0, r2) +
             ifel(is.na(r3), 0, r3))

# convert in binary layers
r_final = r_comb > 0
r_final = as.numeric(r_final)  

# save results
writeRaster(r_final, "fragmentation_layer.tif", overwrite = TRUE)

# Raster binario de fragmentación: 0 = natural, 1 = antrópico
landcover = rast("fragmentation_layer.tif")


# add buffers (500 m arround sampling sites)

buffer_dir=setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/FiltrosAmbientalesURbanos/Fragmentacion/Buffers")
buffer_files = list.files(buffer_dir, pattern = "\\.gpkg$", full.names = TRUE)

# Function for LDI wihth landscape metrics

calcular_division_natural = function(buffer_path) {
  # Read buffer
  buf = vect(buffer_path)
  
# Match CRS with raster
  buf = project(buf, crs(landcover))
  
# Calculate Landscape Division Index for class = 0 (natural cover)
  div = sample_lsm(landcover,
                    y = buf,
                    what = "lsm_c_division",  # class-level division
                    classes = 0,              # only natural pixels
                    progress = FALSE)
  
 # Add buffer name (file name)
  div$buffer_name = tools::file_path_sans_ext(basename(buffer_path))
  
  return(div)
}

# Loop through all buffers

resultados = lapply(buffer_files, calcular_division_natural)

# Combine into one data frame
resultados_df = bind_rows(resultados)

View(resultados_df)

# save results
setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/Fragmentacion")

write.csv(resultados_df, file = "LDI.csv", row.names = FALSE)


# ----------------------------------------------------------
# OTPUT
# Each row = one buffer
# Columns:
#   layer | level | class | id | metric | value | buffer_name
# ----------------------------------------------------------


#Estimation of imprevious area in the buffers

setwd("C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/FiltrosAmbientalesUrbanos/Fragmentacion")
r1 = rast("Area_Impermeabilizada_Raster.tif")  
r2 = rast("Open_Buildings_Raster.tif") 

r2 = project(r2, r1)
r2 = resample(r2, r1, method = "near")  # 'near' porque es binario

ext_total = union(ext(r1), ext(r2))

#combine both layers
r_comb = (ifel(is.na(r1), 0, r1)) +
             ifel(is.na(r2), 0, r2)

# Convert to binary
r_final = r_comb > 0
r_final = as.numeric(r_final)            


buffer_dir = "C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/FiltrosAmbientalesUrbanos/Fragmentacion/Buffers"
buffers_files = list.files(buffer_dir, pattern = "\\.gpkg$", full.names = TRUE)

buffers_list = lapply(buffers_files, st_read)
buffers = do.call(rbind, buffers_list)
          
buffers = st_transform(buffers, crs(r_final))             
buffers_v = vect(buffers)
vals = terra::extract(r_final, buffers_v)

# Estimate percentaje of pixels with value 1 (imprevious area)
result = aggregate(vals[,2], by = list(vals$ID), FUN = function(x) {
  mean(x == 1, na.rm = TRUE) * 100
})
names(result) = c("ID", "perc_impermeable")

plot_ids = data.frame(ID = seq_len(nrow(buffers)), Plot = buffers$Plot)
result = merge(plot_ids, result, by = "ID", all.x = TRUE)

#Save
write.csv(result[, c("Plot", "perc_impermeable")],
          "C:/Users/ASUS/OneDrive - Facultad de Agronomía - Universidad de Buenos Aires/Escritorio/Doctorado - IADIZA/FiltrosAmbientalesUrbanos/Fragmentacion/porcentaje_impermeable_por_plot.csv",
          row.names = FALSE)
