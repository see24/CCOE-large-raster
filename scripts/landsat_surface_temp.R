# Process Landsat data to get surface temperature in Ottawa area

library(dplyr)
library(purrr)
library(stringr)
library(terra)
library(sf)
library(ggplot2)
library(tidyr)
library(httr) # for dealing with kml sites
library(rvest) # for dealing with kml sites

raw_dat <- "data/data-raw"
raw_dat_pth <- "data/data-raw/Landsat2021SummerOttawa_DL20220401"
int_dat_pth <- "data/data-int"
out_dat_pth <- "data/data-out"

# multiplicative scale factor
msf <- 0.00341802

# additive scale factor
asf <- 149

# Load files --------------------------------------------------------------

# get file names
all_fls <- list.files(raw_dat_pth)

# get the tile IDs
tiles <- str_extract(all_fls, "(?<=SP_).{6}") %>% unique()

# Make a list of list of rasters one for each tile for the ST_B10 band
ST_tile_rast_ls <- map(tiles, 
                      ~list.files(raw_dat_pth, paste0(.x, ".*ST_B10"), 
                                       full.names = TRUE)) %>% 
  map(~map(.x, rast))

names(ST_tile_rast_ls) <- tiles
ST_tile_rast_ls <- map(ST_tile_rast_ls, ~set_names(.x, map_chr(.x, names)))

# Make a list of list of rasters one for each tile for the QA_PIXEL band
QA_tile_rast_ls <- map(tiles, 
                       ~list.files(raw_dat_pth, paste0(.x, ".*QA_PIXEL"), 
                                   full.names = TRUE)) %>% 
  map(~map(.x, rast))

names(QA_tile_rast_ls) <- tiles
QA_tile_rast_ls <- map(QA_tile_rast_ls, ~set_names(.x, map_chr(.x, names)))

# # get a minimum extent for each tile to crop rasters to same extent
# get_min_ext <- function(r_ls){
#   xmn <- map_dbl(r_ls, ~xmin(.x)) %>% max()
#   xmx <- map_dbl(r_ls, ~xmax(.x)) %>% min()
#   ymn <- map_dbl(r_ls, ~ymin(.x)) %>% max()
#   ymx <- map_dbl(r_ls, ~ymax(.x)) %>% min()
#   
#   ext(xmn, xmx, ymn, ymx)
# }
# 
# min_ext_ls <- map(ST_tile_rast_ls, get_min_ext)
# 
# ST_tile_rast_ls_crop <- map2(ST_tile_rast_ls, min_ext_ls, 
#                              ~map(.x, crop, y = .y, 
#                                   filename = file.path(int_dat_pth, names(.x))))
# 
# # combine the lists of rasters into one raster stack per tile
# ST_tile_stk_ls <- map(ST_tile_rast_ls_crop, rast)
# 
# rm(ST_tile_rast_ls_crop)

# Do corrections ----------------------------------------------------------

# This tool works for getting bit values. You read it backwards but 2 bit parts you read forwards?
RStoolbox::decodeQA(21824)

# try to find one "good" scene per tile

freq_ls <- map(QA_tile_rast_ls, ~map_dfr(.x, freq, usenames = TRUE))

decode_qa <- function(x){
  x %>% rowwise() %>%
    mutate(bits = RStoolbox::decodeQA(value), 
           fill = substr(bits, 16, 16),
           dcloud = substr(bits, 15, 15),
           cirrus = substr(bits, 14, 14),
           cloud = substr(bits, 13, 13),
           cshadow = substr(bits, 12, 12),
           snow = substr(bits, 11, 11),
           clear = substr(bits, 10, 10),
           water = substr(bits, 9, 9),
           cloud_conf = substr(bits, 7, 8),
           cshadow_conf = substr(bits, 5, 6),
           snow_conf = substr(bits, 3, 4),
           cirrus_conf = substr(bits, 1, 2)) %>% 
    mutate(across(-layer, as.numeric)) %>% 
    mutate(all_clear = sum(cirrus, cloud , cshadow, dcloud, water) == 0) %>% 
    mutate(all_clear = ifelse(fill == 1, NA, all_clear))
}

freq_df <- map_dfr(freq_ls, decode_qa, .id = "tile")

best_days <- filter(freq_df, value != 1, all_clear) %>% 
  mutate(layer = str_extract(layer, "(?<=SP_\\d{6}_).{8}")) %>% 
  group_by(tile, layer) %>% 
  summarise(totalcount = sum(count, na.rm = TRUE)) %>% 
  group_by(tile) %>% 
  filter(totalcount == max(totalcount, na.rm = TRUE))

ST_tile_best <- map2(ST_tile_rast_ls, best_days$layer, 
                       ~.x[[which(str_detect(names(.x), .y))]])

QA_tile_best <- map2(QA_tile_rast_ls, best_days$layer, 
                     ~.x[[which(str_detect(names(.x), .y))]])

par(mfrow = c(4,2))
walk2(ST_tile_best, QA_tile_best, ~(function(x, y){
  plot(x)
  plot(y)
  })(.x, .y))

#freq_df %>% filter(layer == names(QA_tile_best[[1]]))

make_mask <- function(qa_r, k_r, df){
  rcl <- df %>% filter(layer == names(qa_r)) %>% 
    select(value, all_clear)
  msk <- classify(qa_r, rcl)
  mask(k_r, msk, maskvalue = FALSE, overwrite = TRUE,
       filename = file.path(int_dat_pth, paste0(names(k_r), "masked.tif")))
}

ST_tile_best_mskd <- map2(QA_tile_best, ST_tile_best, make_mask, 
                            df = freq_df)

# Set the pattern to choose data from the same day
pat <- "20210529"

vrt(list.files(int_dat_pth, pattern = pat, full.names = TRUE),
    file.path(int_dat_pth, "vrt_masked.vrt"), overwrite = TRUE)

vrt_masked <- rast(file.path(int_dat_pth, "vrt_masked.vrt"))


conv_to_C <- function(x, msf, asf){(x*msf+asf) - 273.15}

# convert to Celsius
ST_tile_best_C <- app(vrt_masked, fun = conv_to_C, msf = msf, asf = asf,
                      overwrite = TRUE,
                      filename = file.path(out_dat_pth, "ST_landsat_best.tif")) 



