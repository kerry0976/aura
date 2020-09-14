# build linear interpolaters for the 'standard' flight data fields.

import math
import numpy as np
from scipy import interpolate # strait up linear interpolation, nothing fancy

# helpful constants
d2r = math.pi / 180.0

class FlightInterpolate():
    def __init__(self, data):
        # df is a pd.DataFrame indexed by time (in seconds)
        columns = {}
        if not len(data):
            return
        for key in data[0]:
            if type(data[0][key]) != str:
                #print("  field:", key, type(data[0][key]))
                columns[key] = []
        for record in data:
            for key in columns:
                columns[key].append(record[key])
        for key in columns:
            columns[key] = np.array(columns[key])
        self.interp = {}
        for key in columns:
            self.interp[key] = interpolate.interp1d(columns['time'],
                                                    columns[key],
                                                    bounds_error=False,
                                                    fill_value=0.0)

    def query(self, t):
        result = {}
        result['time'] = t
        for key in self.interp:
            result[key] = self.interp[key](t).item()
        return result

class pdFlightInterpolate():
    def __init__(self, df):
        # df is a pd.DataFrame indexed by time (in seconds)
        self.interp = {}
        for column in df.columns:
            self.interp[column] = interpolate.interp1d(df.index, df[column],
                                                       bounds_error=False,
                                                       fill_value=0.0)

    def query(self, t):
        result = {}
        result['time'] = t
        for key in self.interp:
            result[key] = self.interp[key](t).item()
        return result

class InterpolationGroup():
    def __init__(self, data):
        self.group = {}
        for key in data:
            #print("group:", key)
            self.group[key] = FlightInterpolate(data[key])

    def query(self, t, key):
        return self.group[key].query(t)

# emulate realtime linear processing of a data set
class IterateGroup():
    def __init__(self, data):
        self.data = data
        self.counter = {}
        for key in data:
            self.counter[key] = 0

    def size(self):
        if 'imu' in self.data:
            return len(self.data['imu'])
        else:
            return 0
    
    def next(self):
        result = {}
        # next imu record
        i = self.counter['imu']
        if i < len(self.data['imu']):
            t = self.data['imu'][i]['time']
            #print("t = ", t)
            result['imu'] = self.data['imu'][i]
            self.counter['imu'] += 1
            
            # look for new records of other types
            for key in self.data:
                if key != 'imu':
                    #print(" ", key)
                    i = self.counter[key]
                    while i < len(self.data[key]):
                        d = self.data[key][i]
                        if d['time'] <= t:
                            #print("  ", d['time'])
                            result[key] = d
                            self.counter[key] += 1
                            i = self.counter[key]
                            if d['time'] == t:
                                break
                        else:
                            break
        return result
            
