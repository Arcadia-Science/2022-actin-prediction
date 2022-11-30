library(readr)
library(dplyr)
library(purrr)

files <- unlist(snakemake@input[['tsv_summaries']])
# files <- Sys.glob("outputs/shared_feature_residues/2_shared_residue_summaries/*tsv")
atp_files <- files[grepl(pattern = "atp_binding.tsv", x = files)]
lat_files <- files[grepl(pattern = "lateral_actin_contact.tsv", x = files)]
lon_files <- files[grepl(pattern = "longitudinal_actin_contact.tsv", x = files)]

lon <- lon_files %>%
  purrr::set_names() %>%
  purrr::map_dfr(readr::read_tsv, show_col_types = F, .id = "protein") %>%
  dplyr::mutate(protein = gsub("-longitudinal_actin_contact.tsv", "", basename(protein))) %>%
  dplyr::select(protein, lon_feature_count = feature_count, 
                lon_num_matching = num_matching, 
                lon_fraction_matching = fraction_matching)

lat <- lat_files %>%
  purrr::set_names() %>%
  purrr::map_dfr(readr::read_tsv, show_col_types = F, .id = "protein") %>%
  dplyr::mutate(protein = gsub("-lateral_actin_contact.tsv", "", basename(protein))) %>%
  dplyr::select(protein, lat_feature_count = feature_count, 
                lat_num_matching = num_matching, 
                lat_fraction_matching = fraction_matching)

atp <- atp_files %>%
  purrr::set_names() %>%
  purrr::map_dfr(readr::read_tsv, show_col_types = F, .id = "protein") %>%
  dplyr::mutate(protein = gsub("-atp_binding.tsv", "", basename(protein))) %>%
  dplyr::select(protein, atp_feature_count = feature_count, 
                atp_num_matching = num_matching, 
                atp_fraction_matching = fraction_matching)

all_features <- lon %>%
  dplyr::left_join(lat, by = "protein") %>%
  dplyr::left_join(atp, by = "protein") %>%
  dplyr::mutate(w_avg_contacts =  (lat_num_matching + lon_num_matching)/(lat_feature_count + lon_feature_count))

readr::write_tsv(all_features, snakemake@output[['all_features']])
#write_tsv(all_features, "20221129-all-features.tsv")
