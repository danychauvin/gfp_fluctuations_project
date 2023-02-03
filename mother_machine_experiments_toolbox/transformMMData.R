# Setting strains
myframes_complete <- myframes_complete %>% 
  mutate(strain='MG1655')
# Converting GFP units
myframes_complete <- myframes_complete %>%
  # append relevant conversion parameters
  mutate(autofluo_predict=NA, fp_per_dn=NA,
         # case of GFPmut2 (from zaslaver library)
         autofluo_predict = ifelse(strain!='asc662', 133.6, autofluo_predict),
         fp_per_dn = ifelse(strain!='asc662', 0.198, fp_per_dn),
         # case of asc662
         autofluo_predict = ifelse(strain=='asc662', 422.8, autofluo_predict),
         fp_per_dn = ifelse(strain=='asc662', 0.0361 * 4, fp_per_dn)) %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um_oriented_bbox,
         gfp_nb = fluogfp_amplitude * fp_per_dn ) %>% 
  group_by(cell) %>% 
  arrange(time_sec) %>% 
  mutate(g_birth=first(gfp_nb)) %>% 
  mutate(g_div=last(gfp_nb))

myframes <- myframes %>% 
  mutate(strain='MG1655')
# Converting GFP units
myframes <- myframes %>%
  # append relevant conversion parameters
  mutate(autofluo_predict=NA, fp_per_dn=NA,
         # case of GFPmut2 (from zaslaver library)
         autofluo_predict = ifelse(strain!='asc662', 133.6, autofluo_predict),
         fp_per_dn = ifelse(strain!='asc662', 0.198, fp_per_dn),
         # case of asc662
         autofluo_predict = ifelse(strain=='asc662', 422.8, autofluo_predict),
         fp_per_dn = ifelse(strain=='asc662', 0.0361 * 4, fp_per_dn)) %>% 
  # convert to gfp units (after subtracting autofluorescence)
  mutate(fluogfp_amplitude = fluo_amplitude - autofluo_predict * length_um_oriented_bbox,
         gfp_nb = fluogfp_amplitude * fp_per_dn )
