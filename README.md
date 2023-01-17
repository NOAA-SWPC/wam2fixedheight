# Instructions

The following has been tested on Linux Mint with anaconda3-2021.11, cmake 3.16.3, and the NetCDF dev libraries installed.

```
git clone --branch main https://github.com/NOAA-SWPC/wam2fixedheight
cd wam2fixedheight

./compile.sh
```

`cp` in input_parameters.nc and some gsm10.*.nc files, and then:

`python convert.py`

The above will create equivalent wam_fixed_height.gsm10.*.nc files using multiprocessing, 16 parallel processes as defined by default in PARALLEL_PROCESSES at the top of convert.py

Alternatively, you can write your own wrapper using Python; here’s a very simple example:

```
from convert import convert
convert(‘/path/to/gsm10.YYYYMMDD_HHMMSS.nc’, driver_file=’/path/to/input_parameters.nc’)
```

and a single output would appear in the current working directory (assuming int_driver and msis21.parm are also in the CWD).

