
import gzip
import numpy as np
import xarray as xr
import copy

def read_radiometer_daily_xr(satname='amsre', year=2005, month=1, day=1):

    valid_sats = ['f08', 'f10', 'f11', 'f13', 'f14', 'f15', 'f16', 'f17', 'f18',
                  'amsre', 'amsr2', 'windsat', 'tmi', 'gmi']
    valid_ssmi = ['f08', 'f10', 'f11', 'f13', 'f14', 'f15', 'f16', 'f17', 'f18']

    ssmi_path =    '//athena/public/ftp/ssmi/'
    amsre_path =   '//athena/public/ftp/amsre/bmaps_v07/'
    amsr2_path =   '//athena/public/ftp/amsr2/bmaps_v08/'
    windsat_path = '//athena/public/ftp/windsat/bmaps_v07.0.1/'
    tmi_path =     '//athena/public/ftp/tmi/bmaps_v07.1/'
    gmi_path =     '//athena/public/ftp/gmi/bmaps_v08.2/'




    sat_name = satname.lower()  #I think this copies the string too, so no side effects

    if sat_name in valid_sats:
        if sat_name in valid_ssmi:
            radiometer_data = read_ssmi_daily_xr(sat_name, year=year, month=month, day=day, path=ssmi_path)
            glb_attr = rss_global_attr(year=year,month=month,day=day,platform=satname.upper(),sensor='SSM/I',version='V7.0')
        elif sat_name == 'amsre':
            radiometer_data = read_amsre_daily_xr(year=year, month=month, day=day, path=amsre_path)
            glb_attr=rss_global_attr(year=year,month=month,day=day,platform='AQUA',sensor='AMSRe',version='V7.0') 
        elif sat_name == 'amsr2':
            radiometer_data = read_amsr2_daily_xr(year=year,month=month, day=day,path=amsr2_path)
            glb_attr=rss_global_attr(year=year,month=month,day=day,platform='GCOM-W',sensor='AMSR2',version='V8.0')
        elif sat_name == 'windsat':
            radiometer_data = read_windsat_daily_xr(year=year,month=month,day=day,path=windsat_path)
            glb_attr=rss_global_attr(year=year,month=month,day=day,platform='Coriolis',sensor='WindSat',version='V7.0')
        elif sat_name == 'tmi':
            radiometer_data = read_tmi_daily_xr(year=year,month=month,day=day,path=tmi_path)
            glb_attr=rss_global_attr(year=year,month=month,day=day,platform='TRMM',sensor='TMI',version='V7.1')
        elif sat_name == 'gmi':
            radiometer_data = read_gmi_daily_xr(year=year,month=month,day=day,path=gmi_path)
            glb_attr=rss_global_attr(year=year,month=month,day=day,platform='GPM',sensor='GMI',version='V8.2')
 
        for key in glb_attr.keys():
            print(key, glb_attr[key])
            radiometer_data.attrs[key] = glb_attr[key]
        return radiometer_data
    else:
        raise(ValueError(f"satname {satname} not valid"))


def read_amsre_daily_xr(year = 2005,month=1,day=1,path='//athena/public/ftp/amsre/bmaps_v07/'):

    #construct the file name
    filename = (path + f"y{year:04d}/m{month:02d}/f32_{year:04d}{month:02d}{day:02d}v7.gz")
    f = gzip.open(filename,'rb')
    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    sst       = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_lf   = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (14, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*7
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        sst[node, :, :]       = (byte_data[map_offset + 1, :, :]*0.15) - 3.0
        wspd_lf[node, :, :]   = byte_data[map_offset + 2, :, :]*0.2
        wspd_mf[node, :, :]   = byte_data[map_offset + 3, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 4, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 5, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 6, :, :]*0.1

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'sst':         (('node', 'latitude', 'longitude'),sst),
                    'wspd_lf':     (('node', 'latitude', 'longitude'),wspd_lf),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )
    ds['sst'].attrs['standard_name']='sea_surface_temperature'
    ds['sst'].attrs['units']='degrees K'
    print()
    return ds

def read_amsr2_daily_xr(year = 2017,month=1,day=1,path='//athena/public/ftp/amsr2/bmaps_v08/'):

    #construct the file name
    filename = (path + f"y{year:04d}/m{month:02d}/f34_{year:04d}{month:02d}{day:02d}v8.gz")
    f = gzip.open(filename,'rb')
    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    sst       = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_lf   = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (14, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*7
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        sst[node, :, :]       = (byte_data[map_offset + 1, :, :]*0.15) - 3.0
        wspd_lf[node, :, :]   = byte_data[map_offset + 2, :, :]*0.2
        wspd_mf[node, :, :]   = byte_data[map_offset + 3, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 4, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 5, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 6, :, :]*0.1

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'sst':         (('node', 'latitude', 'longitude'),sst),
                    'wspd_lf':     (('node', 'latitude', 'longitude'),wspd_lf),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )

    print()
    return ds

def read_windsat_daily_xr(year = 2017,month=1,day=1,path=None):

    #construct the file name
    filename = (path + f"y{year:04d}/m{month:02d}/wsat_{year:04d}{month:02d}{day:02d}v7.0.1.gz")
    f = gzip.open(filename,'rb')

    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    sst       = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_lf   = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_aw   = np.zeros((2, 720, 1440), dtype=np.float32)
    wdir      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (18, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*9
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        sst[node, :, :]       = (byte_data[map_offset + 1, :, :]*0.15) - 3.0
        wspd_lf[node, :, :]   = byte_data[map_offset + 2, :, :]*0.2
        wspd_mf[node, :, :]   = byte_data[map_offset + 3, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 4, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 5, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 6, :, :]*0.1
        wspd_aw[node, :, :]   = byte_data[map_offset + 7, :, :]*0.2
        wdir[node, :, :]      = byte_data[map_offset + 8, :, :]*1.5

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'sst':         (('node', 'latitude', 'longitude'),sst),
                    'wspd_lf':     (('node', 'latitude', 'longitude'),wspd_lf),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain),
                    'wspd_aw':     (('node', 'latitude', 'longitude'),wspd_aw),
                    'wdir':        (('node', 'latitude', 'longitude'),wdir)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )

    print()
    return ds

def read_ssmi_daily_xr(satname,year = 2017,month=1,day=1,path='//athena/public/ftp/ssmi/'):

    #construct the file name
    filename = (path + f"{satname}/bmaps_v07/y{year:04d}/m{month:02d}/{satname}_{year:04d}{month:02d}{day:02d}v7.gz")
    f = gzip.open(filename,'rb')
    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (10, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*5
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        wspd_mf[node, :, :]   = byte_data[map_offset + 1, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 2, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 3, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 4, :, :]*0.1

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )

    print()
    return ds

def read_tmi_daily_xr(year = 2007,month=1,day=1,path='//athena/public/ftp/tmi/'):

    #construct the file name
    filename = (path + f"y{year:04d}/m{month:02d}/f12_{year:04d}{month:02d}{day:02d}v7.1.gz")
    f = gzip.open(filename,'rb')
    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    sst       = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_lf   = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (14, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*7
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        sst[node, :, :]       = byte_data[map_offset + 1, :, :]*0.2
        wspd_lf[node, :, :]   = byte_data[map_offset + 2, :, :]*0.2
        wspd_mf[node, :, :]   = byte_data[map_offset + 3, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 4, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 5, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 6, :, :]*0.1

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'sst':         (('node', 'latitude', 'longitude'),sst),
                    'wspd_lf':     (('node', 'latitude', 'longitude'),wspd_lf),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )
    ds['sst'].attrs['standard_name']='sea_surface_temperature'
    ds['sst'].attrs['units']='degrees K'

    print()
    return ds

def read_gmi_daily_xr(year = 2007,month=1,day=1,path='//athena/public/ftp/gmi/'):

    #construct the file name
    filename = (path + f"y{year:04d}/m{month:02d}/f35_{year:04d}{month:02d}{day:02d}v8.2.gz")
    f = gzip.open(filename,'rb')
    #need to use frombuffer because np.fromfile doesn't use gzip
    z = np.frombuffer(f.read(),dtype=np.uint8)
    f.close()

    time_min  = np.zeros((2, 720, 1440), dtype=np.float32)
    sst       = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_lf   = np.zeros((2, 720, 1440), dtype=np.float32)
    wspd_mf   = np.zeros((2, 720, 1440), dtype=np.float32)
    vapor     = np.zeros((2, 720, 1440), dtype=np.float32)
    cloud     = np.zeros((2, 720, 1440), dtype=np.float32)
    rain      = np.zeros((2, 720, 1440), dtype=np.float32)

    byte_data = np.reshape(z, (14, 720, 1440)).astype(np.float32)
    byte_data[byte_data > 250] = np.nan
    
    for node in (0,1):
        map_offset = node*7
        time_min[node, :, :]  = byte_data[map_offset + 0, :, :]*6.0
        sst[node, :, :]       = byte_data[map_offset + 1, :, :]*0.2
        wspd_lf[node, :, :]   = byte_data[map_offset + 2, :, :]*0.2
        wspd_mf[node, :, :]   = byte_data[map_offset + 3, :, :]*0.2
        vapor[node, :, :]     = byte_data[map_offset + 4, :, :]*0.3
        cloud[node, :, :]     = (byte_data[map_offset + 5, :, :]*0.01) - 0.05
        rain[node, :, :]      = byte_data[map_offset + 6, :, :]*0.1

    latitude = -90.0 + 0.125 + np.arange(0,720)*0.25
    longitude = 0.125 + np.arange(0,1440)*0.25
    node = np.arange(0,2)
    
    ds = xr.Dataset(
        data_vars= {'time':        (('node', 'latitude', 'longitude'),time_min),
                    'sst':         (('node', 'latitude', 'longitude'),sst),
                    'wspd_lf':     (('node', 'latitude', 'longitude'),wspd_lf),
                    'wspd_mf':     (('node', 'latitude', 'longitude'),wspd_mf),
                    'vapor':       (('node', 'latitude', 'longitude'),vapor),
                    'cloud':       (('node', 'latitude', 'longitude'),cloud),
                    'rain':        (('node', 'latitude', 'longitude'),rain)},
        coords=    {'node': node,
                    'latitude':latitude,
                    'longitude':longitude}
                    )
    ds['sst'].attrs['standard_name']='sea_surface_temperature'
    ds['sst'].attrs['units']='degrees K'

    print()
    return ds


def rss_global_attr(year=1988,month=1,day=1,platform='F08',sensor='SSM/I',version='V7.0'):

    import datetime
    dt = datetime.datetime.today()
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
    'geospatial_lat_min' : -90.0, 
    'geospatial_lat_max' : 90.0, 
    'geospatial_lon_min' : -180.0, 
    'geospatial_lon_max' : 180.0, 
    'geospatial_lat_units' : 'degrees_north',
    'geospatial_lon_units' : 'degrees_east',
    'geospatial_lat_resolution' : 0.25, 
    'geospatial_lon_resolution' : 0.25, 
    'time_coverage_start' : f"{year}-{month:02d}-{day:02d}T000000Z",
    'time_coverage_end' : f"{year}-{month:02d}-{day:02d}T235959Z",
    'time_coverage_duration' : '1D',
    'product_version' : version,
    'platform' : platform,
    'sensor' : sensor,
    'spatial_resolution' : '0.25 degree'
    }
    return default_global_attr


if __name__ == '__main__':
    import sys
    sys.path.append('C:/job_CCMP/python/')
    import matplotlib.pyplot as plt

    def global_map(a, vmin=0.0, vmax=30.0, cmap=None, plt_colorbar=False,title=''):
        import numpy as np
        import matplotlib.pyplot as plt
        import matplotlib as mpl
        import cartopy.crs as ccrs

        img_extent = [-180.0, 180.0, -90.0, 90.0]
        fig = plt.figure(figsize=(10, 5))  
        ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(),title=title)
        for item in ([ax.title]):
            item.set_fontsize(16)
        map = ax.imshow(np.flipud(np.roll(a, 720, axis=1)), cmap=cmap, origin='upper', transform=ccrs.PlateCarree(),
                        norm=mpl.colors.Normalize(vmin=vmin, vmax=vmax), extent=img_extent)
        if plt_colorbar:
            cbar = fig.colorbar(map, shrink=0.7, orientation='horizontal')
            cbar.ax.tick_params(labelsize=14)
        ax.coastlines()
        ax.set_global()
        return fig, ax


    

    year = 2018
    month = 1
    day = 15
    satname='tmi'
    date_str = f"{month:02d}/{day:02d}/{year:04d}"
    ds = read_radiometer_daily_xr(satname=satname,year=year, month=month, day=day)
    global_map(ds['wspd_mf'].values[0, :, :],plt_colorbar=True,title=f"{satname.upper()}, Wind Speed (LF), {date_str}",vmax=25.0)
    global_map(ds['vapor'].values[0, :, :],plt_colorbar=True,title=f"{satname.upper()}, Vapor, {date_str}",vmax=60.0)
    global_map(ds['cloud'].values[1, :, :],plt_colorbar=True,title=f"{satname.upper()}, Cloud, {date_str}",vmax=0.5)
    global_map(ds['rain'].values[1, :, :],plt_colorbar=True,title=f"{satname.upper()}, Rain Rate, {date_str}",vmax=2.0)

    plt.show()

    test_file = f"C:/job_CCMP/python/ssmi/test_nc/{satname}_{year:04d}{month:02d}{day:02d}.nc"
    ds.to_netcdf(test_file)
    print()

