import geopandas as gpd
import pandas as pd
import numpy as np
import os

################################################################################
# FUNCTIONS
################################################################################
def clean_shp_columns(df):
	'''
	Keeps the good columns in polygon dataframe.
	'''
	return df[['GEOID','geometry']]

def clean_cdc_columns(df):
	'''
	Keeps the good columns in cdc table.
	'''
	return df[['CountyName','TotalPopulation','OBESITY_CrudePrev','TractFIPS']]


def join_dotted_polygons(df):
	'''
	Take a geopandas data frame and join the polygons listed as XXXX.YY (bad) to
	create a single XXXX (good) polygon. 
	'''
	mask_bad  = np.array([_[-2:] != '00' for _ in df['GEOID']])#True 914 of 1654
	mask_good = ~mask_bad #True 740, size 1654
	bad_df    = df[mask_bad] #size 914

	new_ids    = np.array([_[:-2] + '00' for _ in  bad_df['GEOID']]) #914
	unique_ids = np.unique(new_ids) #323

	new_df = {'GEOID':[],'geometry':[]}

	for i,u in enumerate(unique_ids):
		subset = bad_df[u == new_ids]
		union  = subset['geometry'].unary_union
		new_df['geometry'].append(union)
		new_df['GEOID'].append(u)

	new_df = pd.concat([gpd.GeoDataFrame(new_df),df[mask_good]],ignore_index=True)
	return new_df
	

def join_dotted_cdc_tracts(df):
	'''
	Joins the tracts in the cdc with trailing digits different from zero (corres
	ponding to subdivision that can be joined to form a single tract). Tracts 
	are joined by (i) setting last two digits to zero, (ii) adding (sum) the
	population  of individual tracts, and (iii) setting the obesity prevalance 
	to a weighted sum: population * obesity / (sum populations)
	'''	
	mask_bad  = np.array([_[-2:] != '00' for _ in df['TractFIPS']]) #506/1397 T
	mask_good = ~mask_bad #881 True, shape 1397
	bad_df = df[mask_bad] #shape (506, 4)

	new_ids        = np.array([_[:-2] + '00' for _ in  bad_df['TractFIPS']])#506
	new_ids_unique = np.unique(new_ids) #178

	new_df = {'CountyName':[],'TotalPopulation':[],
		'OBESITY_CrudePrev':[],'TractFIPS':[]}
	
	#calculate average obesity
	for i,u in enumerate(new_ids_unique):
		subset = bad_df[u == new_ids] # subdivisions matching the joined area

		#population is weights for obesity weighted averages
		weight = subset['TotalPopulation'].to_numpy().astype(int)
		means  = subset['OBESITY_CrudePrev'].to_numpy().astype(float)
		new_avg = (weight * means).sum() / weight.sum()

		#append all to a new dataframe
		new_df['CountyName'].append(subset['CountyName'].iloc[0])
		new_df['TotalPopulation'].append(weight.sum())
		new_df['OBESITY_CrudePrev'].append(new_avg)
		new_df['TractFIPS'].append(u)

	new_df = pd.concat([pd.DataFrame(new_df),df[mask_good]],ignore_index=True)
	return new_df

def save_shapefile(gdf, name_str):
	'''
	Saves a geopandas dataframe to .shp format
	'''
	dirpath = './%s/' % name_str 

	#if subdir doesn't exist
	if os.path.isdir(dirpath) is False: 
		# make a subdir 
		os.mkdir(dirpath)

	#file not in dir?
	if os.path.isfile(dirpath + '%s.shp' % name_str) is False:
		#save
		print("Saving shapefile to " + dirpath + '%s.shp' % name_str)
		gdf.to_file(dirpath + '%s.shp' % name_str) #<------ the actual save	

def intersect_shp_cdc(cdc_df, shp_df):
	'''
	Take the CDC dataframe and polygon dataframe and return the intersection.
	'''
	# some bool variables to keep track of indices in obth arrays
	keep_shp = np.zeros(len(shp_df)).astype(bool) #true/false size 1654
	keep_cdc = np.zeros(len(cdc_df)).astype(bool) #true/false size 1387

	#for each tract in CDC, match GEOID and add True to mask
	for index, tract_id in cdc_df['TractFIPS'].items():
		mask = np.array(tract_id == shp_df['GEOID'])
		if mask.sum() > 0:
			keep_cdc[index] = True
		keep_shp += mask

	# Select rows
	shp_intersect = shp_df[keep_shp].reset_index(drop=True) 
	cdc_intersect = cdc_df[keep_cdc].reset_index(drop=True)

	#sort
	cdc_intersect = cdc_intersect.sort_values(by=['TractFIPS'],
		ignore_index=True) # Sort by tract id
	shp_intersect = shp_intersect.sort_values(by=['GEOID'],ignore_index=True)

	#merge
	intersect_df  = pd.concat([cdc_intersect,shp_intersect],axis=1) 

	#back to geopandas
	intersect_gdf = gpd.GeoDataFrame(intersect_df)

	return intersect_gdf

################################################################################
# MAIN
################################################################################
if __name__ == '__main__': # <--- execution starts here

	#I.READ FILES
	shp_df = gpd.read_file('./tl_rd22_29_tract/tl_rd22_29_tract.shp')
	cdc_df = pd.read_csv('./CDC_tracts.csv')

	#II.CLEAN COLUMNS
	shp_df = clean_shp_columns(shp_df)
	cdc_df = clean_cdc_columns(cdc_df).astype(str)

	#III.INTERSECTION
	intersect_gdf = intersect_shp_cdc(cdc_df,shp_df)
	save_shapefile(intersect_gdf,'intersection')
	
	#IV. JOIN POLYGONS
	join_gdf = join_dotted_polygons(shp_df)
	save_shapefile(join_gdf,'joined_polygons')

	#V. JOIN CDC TRACTS
	join_cdc = join_dotted_cdc_tracts(cdc_df)

	#VI. INTERSECTION OF TABLES WITH WITH JOINED/MERGED CENSUS TRACTS
	intersect_joined_polygons = intersect_shp_cdc(join_cdc,join_gdf)
	intersect_joined_polygons.to_csv('./final.csv')

