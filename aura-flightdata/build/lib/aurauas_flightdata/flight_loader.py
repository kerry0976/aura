import h5py
import os
import pandas as pd

from . import aura_csv
from . import aura_hdf5
from . import px4_sdlog2
from . import px4_ulog
from . import sentera
from . import sentera2
from . import umn1_mat2
#from . import umn1_mat3 # For Grove data Set 
from . import umn3_hdf5

def load(path):
    flight_data = {}
    flight_format = None

    (root, ext) = os.path.splitext(path)
    aura_hdf5_path = os.path.join(path, 'flight.h5')
    aura_csv_path = os.path.join(path, 'imu-0.csv')
    ulog_path = path + '_sensor_combined_0.csv'
    sentera_path = os.path.join(path, 'imu.csv')
    
    # determine the data log format and call the corresponding loader code

    if ext == '.h5':
        # quick peek
        data = h5py.File(path, 'r')
        if "metadata" in data:
            md = data["/metadata"]
            if md.attrs.get("format", "") == "AuraUAS":
                print("Detected AuraUAS hdf5 format.")
                flight_data = aura_hdf5.load(path)
                flight_format = 'aura_hdf5'
        else:
            print('Detected UMN3 (hdf5) format.')
            flight_data = umn3_hdf5.load(path)
            flight_format = 'umn3'
    elif os.path.exists(aura_hdf5_path):
        # aura hdf5 format
        print("Detected AuraUAS hdf5 format.")
        flight_data = aura_hdf5.load(aura_hdf5_path)
        flight_format = 'aura_hdf5'
    elif os.path.exists(aura_csv_path):
        # aura csv format
        print('Detected aura csv format.')
        flight_data = aura_csv.load(path)
        flight_format = 'aura_csv'
    elif ext == '.mat':
        # umn1
        print('Detected umn1 format.')
        print('Notice: assuming umn1 .mat format')
        flight_data = umn1_mat2.load(path) # modified to adapt new data format (K.S.)
        flight_format = 'umn1'
    elif os.path.exists(ulog_path):
        # px4_ulog
        print('Detected px4 ulog (csv family of files) format.')
        print('Support needs code updates')
        quit()
        flight_data = px4_ulog.load(path)
        flight_format = 'px4_ulog'
    elif ext == '.csv':
        # px4 sdlog2
        print('Detected px4 ulog (single csv file) format.')
        print('Support needs code updates')
        quit()
        flight_data = px4_sdlog2.load(path)
        flight_format = 'px4_sdlog2'
    elif os.path.exists(sentera_path):
        # sentera1 or sentera2
        print('Detected sentera format.')
        print('Notice: assuming original sentera camera format')
        print('Support needs code updates')
        quit()
        flight_data = sentera.load(path)
        # imu_data, gps_data, air_data, filter_data = sentera2.load(path)
        flight_format = 'sentera'
    else:
        print('Unable to determine data log format (or path not valid):', path)

    return flight_data, flight_format

def as_pandas(flight_data):
    result = {}
    # convert to pandas DataFrame's
    for key in flight_data:
        result[key] = pd.DataFrame(flight_data[key])
        result[key].set_index('time', inplace=True, drop=False)
    return result
       
def save(filename, data):
    aura_csv.save_filter_result(filename, data)

