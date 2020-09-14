# nav_wrapper2.py
# 1. simplify calling the C++ ekf filters
# 2. optional gps lag support

# import sys
# sys.path.append('/Users/kerry/Desktop/aura/aura-navigation/src/nav_common/')
# sys.path.append('/Users/kerry/Desktop/aura/aura-navigation/')
import pandas as pd

from aurauas_navigation.structs import IMUdata, GPSdata, NAVconfig, GNSS_measurement
#from aurauas_navigation.ekf15 import EKF15
# from aurauas_navigation.ekf15_mag import EKF15_mag
#from aurauas_navigation.uNavINS import uNavINS
#from aurauas_navigation.uNavINS_BFS import uNavINS_BFS
from aurauas_navigation.openloop import OpenLoop
from aurauas_navigation.ekf17 import EKF17

class filter():
    def __init__(self, nav, gps_lag_sec=0.0, imu_dt=0.02):
        self.name = nav
        if nav == 'EKF15':
            self.filter = EKF15()
        elif nav == 'EKF15_mag':
            self.filter = EKF15_mag()
        elif nav == 'EKF17':
            self.filter = EKF17()
        elif nav == 'uNavINS':
            self.filter = uNavINS()
        elif nav == 'uNavINS_BFS':
            self.filter = uNavINS_BFS()
        else:
            print("Unknown nav filter specified aborting:", nav)
            quit()

        self.openloop = OpenLoop()
        self.gps_lag_frames = int(round(gps_lag_sec / imu_dt))
        print("gps lag frame:", self.gps_lag_frames)
        self.imu_queue = []
        
    def set_config(self, config):
        Cconfig = NAVconfig()
        Cconfig.from_dict(config)
        self.filter.set_config(Cconfig)
        
    def update(self, imu, gnss):
        Cimu = IMUdata()
        Cimu.from_dict(imu)

        # queue delay
        self.imu_queue.insert(0, Cimu)
        Cimu = self.imu_queue.pop()
        #self.Cgps = GPSdata()
        #self.Cgps.from_dict( gps )
        #self.filter.update(Cimu, self.Cgps
        #assert "time" in gnss
        # print(gnss['time'])
        # print(gnss.get('gnss_measument'))
        
        self.Cgnss = GNSS_measurement(gnss['time'], gnss.get('gnss_measurement'))
        # self.Cgnss = GNSS_measurement()
        # self.Cgnss.create(gnss["gnss_measument.time"], gnss["gnss_measument"])
        # self.Cgnss.from_dict(gnss) 
        # self.Cgnss.from_dict(gnss['time'], gnss.get('gnss_measurement'))
        self.filter.update(Cimu, self.Cgnss)

        nav = self.filter.get_nav()

        if len(self.imu_queue):
            # forward propagate from the lagged solution to new
            self.openloop.init_by_nav(nav)
            for imu in reversed(self.imu_queue):
                nav = self.openloop.update(imu)
            
        return nav.as_dict()

    def close(self):
        pass
