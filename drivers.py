from netCDF4 import Dataset
import numpy as np

INPUT_FILE = 'input_parameters.nc'
START = 2160
DAY = 60*24

def maxindex(i, hours=0):
    return max([0, i - hours*60])

def process(f107, f107a, ap, apa, i):
    f107_o  = np.mean(f107[maxindex(i-DAY):i])
    f107a_o = f107a[i]
    ap_o = [apa[i], ap[i], ap[maxindex(i,3)],
            ap[maxindex(i,6)], ap[maxindex(i,9)],
            apa[maxindex(i,12)], apa[maxindex(i,36)]]
    return f107_o, f107a_o, ap_o

def get_inputs(i, input_file=INPUT_FILE):
    ncfid = Dataset(input_file)

    f107  = ncfid.variables['f107'][:]
    f107a = ncfid.variables['f107d'][:]
    ap    = ncfid.variables['ap'][:]
    apa   = ncfid.variables['apa'][:]

    return process(f107, f107a, ap, apa, i)

def fwrite(fn):
    from scipy.io import FortranFile
    with FortranFile(fn, 'w') as f:
        f.write_record(f107)
        f.write_record(f107a)
        f.write_record(ap)

def main():
    f107, f107a, ap = get_inputs(START)
    fwrite('driver_file', f107, f107a, ap)

if __name__ == '__main__':
    main()

