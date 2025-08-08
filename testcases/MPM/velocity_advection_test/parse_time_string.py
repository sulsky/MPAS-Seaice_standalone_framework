from datetime import datetime
from netCDF4 import date2num

#--------------------------------------------------------------------

def parse_time_string(xtime):

    date = datetime.strptime(b''.join(xtime.compressed()).decode('utf-8'),"%Y-%m-%d_%H:%M:%S")
    secs = date2num(date, units='secs since 2000-01-01',
                                         calendar='noleap')
    return secs

#--------------------------------------------------------------------

if __name__ == "__main__":

    parse_time_string(xtime)
