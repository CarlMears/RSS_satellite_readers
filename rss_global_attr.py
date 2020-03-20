def rss_global_attr(year=year,month=month,day=day,platform='F08',sensor='SSM/I',version='V7.0')

    import datetime

    dt = datetime.datetime.today
    date_str = f"{dt.year}{dt.month:02d}{dt.day:02d}T{dt.hour:02d}{dt.minute:02d}{dt.second:02d}Z"
    default_global_attr = {
    'Conventions' : 'CF-1.6',
    'title' : 'daily oceanic retreivals on a 2.5 degree grid',
    'source' : 'imaging microwave radiometer',
    'references' : 'doi10.1175/JCLI-D-15-0744.1, doi10.1175/JCLI-D-16-0768.1',
    'history' : 'netcdf4 file converted from binary file Remote Sensing Systems',
    'Metadata_Conventions' : 'CF-1.6, Unidata Dataset Discovery v1.0, NOAA CDR v1.0, GDS v2.0',
    'standard_name_vocabulary' : 'CF Standard Name Table (v19, 22 March 2012)',
    'date_created' : date_str,
    'license' : 'No constraints on data access or use',
    'cdm_data_type' : 'grid',
    'project' : 'RSS Ocean Retrievals',
    'processing_level' : 'NASA Level 3',
    'creator_name' : 'Remote Sensing Systems',
    'creator_url' : 'http//www.remss.com/',
    'creator_email' : 'support@remss.com',
    'institution' : 'Remote Sensing Systems',
    'geospatial_lat_min : -90.0f, 
    'geospatial_lat_max : 90.0f, 
    'geospatial_lon_min : -180.0f, 
    'geospatial_lon_max : 180.0f, 
    'geospatial_lat_units' : 'degrees_north',
    'geospatial_lon_units' : 'degrees_east',
    'geospatial_lat_resolution : 0.25, 
    'geospatial_lon_resolution : 0.25, 
    'time_coverage_start' : f"{year}-{month:02d}-{day:02d}T000000Z",
    'time_coverage_end' : f"{year}-{month:02d}-{day:02d}T235959Z",
    'time_coverage_duration' : '1D',
    'product_version' : version,
    'platform' : platform,
    'sensor' : sensor,
    'spatial_resolution' : '0.25 degree'
    }

    return default_global_attr