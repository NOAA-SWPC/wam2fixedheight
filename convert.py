from datetime import datetime
from pathlib import Path
from glob import glob
from subprocess import Popen, PIPE, STDOUT
from netCDF4 import Dataset, date2num
from drivers import get_inputs, INPUT_FILE
from multiprocessing import Pool
from os.path import basename

PARALLEL_PROCESSES = 16
NML_TEMPLATE = '''
&FIXED_HEIGHT
  IDOY        = {:d}
  IHOUR       = {:.3f}
  INPUT_FILE  = '{}'
  OUTPUT_FILE = '{}'
  F107        = {:.2f}
  F107A       = {:.2f}
  AP          = {}
/
'''
INITIAL_DRIVER_DT = datetime(2021,7,14)

OUTPUT_FILE_FMT = 'wam_fixed_height.{}'
TIME_UNITS = 'days since 1970-01-01'

def update_time_field(file, nc_dts):
    ncfid = Dataset(file, 'r+')
    ncfid.variables['time'][:] = nc_dts
    ncfid.close()

def convert(file, driver_file=INPUT_FILE):
    output_file = OUTPUT_FILE_FMT.format(basename(file))

    ncfid = Dataset(file)

    dt = datetime.strptime(ncfid.fcst_date, '%Y%m%d_%H%M%S')
    doy = dt.timetuple().tm_yday
    hour = dt.hour + dt.minute / 60

    nc_dts = [date2num(dt, TIME_UNITS)]

    start_idx = int((dt - INITIAL_DRIVER_DT).total_seconds() // 60)

    f107, f107a, ap = get_inputs(start_idx, input_file=driver_file)
    p = Popen(['./int_driver'], stdout=PIPE, stdin=PIPE, stderr=PIPE)
    p.communicate(input=NML_TEMPLATE.format(doy,
                                            hour,
                                            file,
                                            output_file,
                                            f107,
                                            f107a,
                                            ', '.join(['{:.2f}'.format(a) for a in ap])).encode())

    update_time_field(output_file, nc_dts)

def main():
    files = sorted(glob('gsm10.*'))
    with Pool(processes=PARALLEL_PROCESSES) as p:
        p.map(convert, files)

if __name__ == '__main__':
    main()
