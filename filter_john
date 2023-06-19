# Identify valid crossovers (those that have realistic and contemporaneous values):
nan = np.nan
df_valid = df.query('(`/reference/canopy_ht/canopy_ht_q100` != @nan) and \
                    (`/reference/canopy_ht/canopy_ht_q100` < 100) and \
                    (`/reference/canopy_ht/canopy_ht_q100` > -5) and \
                    (`/reference/canopy_ht/canopy_ht_q100` > -5) and \
                    (`/simulation/als_zt` != @nan) and \
                    (`/simulation/als_zg` != @nan) and \
                    (`/simulation/als_rh100` < 100) and \
                    (`/simulation/als_rh098` > 2.5) and \
                    (`/land_cover_data/urban_proportion` == 0)')
df_valid = df_valid[df_valid['/simulation/als_project'].isin(valid_sites)]
#Identify orbits affected by noise:
def get_tile_id(longitude, latitude, tilesize=72):
    ease2_origin = -17367530.445161499083042, 7314540.830638599582016
    ease2_nbins = int(34704 / tilesize), int(14616 / tilesize)
    ease2_binsize = 1000.895023349556141*tilesize, 1000.895023349562052*tilesize 
    transformer = Transformer.from_crs('epsg:4326', 'epsg:6933', always_xy=True)
    x,y = transformer.transform(longitude, latitude)
    xidx = np.uint16( (x - ease2_origin[0]) / ease2_binsize[0]) + 1
    yidx = np.uint16( (ease2_origin[1] - y) / ease2_binsize[1]) + 1
    tile_id = np.array([f'X{x:03d}Y{y:03d}' for x,y in zip(xidx,yidx)])
    return tile_id

fn = '/gpfs/data1/vclgp/armstonj/l4b_granule_filtering/issgedi_l4b_excluded_granules_r002_20230315a_plusManual.json'
with open(fn, 'r') as f:
    excluded_granules = json.load(f)

tile_id = get_tile_id(df_valid['/geolocation/longitude'].to_numpy(), df_valid['/geolocation/latitude'].to_numpy())
orbit = (df_valid.index // 10000000000000)
granule = (df_valid.index % 100000000000) // 100000000
granule_id = orbit * 10 + granule

granule_mask = np.array([(True if g in excluded_granules[t] else False) if t in excluded_granules else False
                         for g,t in zip(granule_id,tile_id)])

df_valid = df_valid.assign(granule_mask=granule_mask)





c = df.loc[valid,'/gedi/gedi_rh_a2/sensitivity'].to_numpy()
pft = df.loc[valid,'/land_cover_data/modis_pft_1km'].to_numpy()
region = df.loc[valid,'/land_cover_data/continental_region_1km'].to_numpy()
tropics = (pft == 2) & np.isin(region, [4,5,6])
s = np.full(c.shape, 0.95, dtype=np.float32)
s[tropics] = 0.98
q4 = (c > s) & (c < 1)

df_quality = df_valid.query('(`/reference/ground_elev/ground_elev_count` > 10) and \
                            (`/gedi/gedi_rh_a1/sensitivity` > 0.7) and \
                            (`/gedi/gedi_rh_a1/sensitivity` < 1) and \
                            (`/geolocation/dz_pearson` >= 0.9) and \
                            (`/geolocation/offset_pearson` >= 0.9) and \
                            (`/geolocation/dz_count` >= 10) and \
                            (`/geolocation/offset_count` > 40) and \
                            (`/geolocation/pearson` > 0.95) and \
                            (`loss_sum_3x3` == 0) and \
                            (`granule_mask` != True) and \
                            (`mw` <= 202)')
