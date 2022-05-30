import glob
import os
import re


"""
Create symlink files of the form

source_model/SalishSeaCast-VNR002_1h_grid_T_2D_y2016m02d03.nc
source_model/SalishSeaCast-VNR002_1h_grid_T_2D_y2016m02d04.nc
source_model/SalishSeaCast-VNR002_1h_grid_T_2D_y2016m02d05.nc
source_model/SalishSeaCast-VNR002_1h_grid_T_2D_y2016m02d06.nc
source_model/SalishSeaCast-VNR002_1h_grid_T_2D_y2016m02d07.nc

that reference the original files

model/SalishSeaCast-VNR002_1h_grid_T_20160203-20160203.nc
model/SalishSeaCast-VNR002_1h_grid_T_20160204-20160204.nc
model/SalishSeaCast-VNR002_1h_grid_T_20160205-20160205.nc
model/SalishSeaCast-VNR002_1h_grid_T_20160206-20160206.nc
model/SalishSeaCast-VNR002_1h_grid_T_20160207-20160207.nc
"""

p = re.compile(r'(.*)_(\d{4})(\d{2})(\d{2})-\2\3\4.nc')

ncfilelist = glob.glob('model/*.nc')
for ncfile in ncfilelist:
    basename = os.path.basename(ncfile)

    m = re.search(p, basename)
    new_basename = f"{m[1]}_y{m[2]}m{m[3]}d{m[4]}.nc"

    new_ncfile = os.path.join('source_model', new_basename)

    if not os.path.exists(new_ncfile):
        os.symlink( os.path.join('..', ncfile), new_ncfile)
