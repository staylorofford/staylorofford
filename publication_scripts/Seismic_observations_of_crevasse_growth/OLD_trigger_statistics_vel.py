# -*- coding: utf-8 -*-
"""
Created on Thu Oct 13 10:27:27 2016

@author: taylorsa
""" 

# KEY ASPECTS OF THESE FUNCTIONS:

# HOURLY DATA IS OBSERVED IN THE HOUR PRIOR TO THE HOUR OF OBSERVATION, I.E.
# 1100 ON JULIAN DAY 2 IS MEASURED FROM 1000-1100 ON JULIAN DAY 2

# DAILY DATA IS OBSERVED IN THE DAY PRIOR TO THE DAY OF OBSERVATION, I.E.
# JULIAN DAY 3 IS MEASURED FROM 0000 TO 2400 ON JULIAN DAY 2
# NOTE: THE ABOVE IS NON-ARBITRARY. IT FOLLOWS FROM THE WAY HOURLY DATA IS BINNED   

# standardise plotting segment

from matplotlib import rc
from matplotlib import rcParams
import brewer2mpl

# adjust generic plot params

params = {
   'axes.labelsize': 10,
   'text.fontsize': 10,
   'legend.fontsize': 10,
   'xtick.labelsize': 10,
   'ytick.labelsize': 10,
   'text.usetex': False,
   'figure.figsize': [3.54331, 3.54331] # default is 90 mm wide and high
   }
rcParams.update(params)

# update font to arial, size 10

rc('text', usetex=False)
rc('font', **{'family': 'Arial'})
rcParams.update({'font.size': 10})

# adjust default colour range for line plots

bmap = brewer2mpl.get_map('Set1', 'qualitative', 5)
colors = bmap.mpl_colors
rcParams['axes.color_cycle'] = colors
std_linewidth=2

# tasteful outline for scatterplots

std_edgecolors='grey'
std_linewidths=0.1

# plot only variable saturation when using a colourbar scale

brewer2mpl.get_map('Greys', 'sequential', 8).mpl_colormap
 
def runmean(interval, window_size):
    import numpy as np
    # HJH
    # INPUTS    interval        time series
    #           window_size     number of samples to smooth over
    window= np.ones(int(window_size))/float(window_size)
    return np.convolve(interval, window, 'same')




def interp1d(x,y,xi):
    import numpy as np
    # HJH
    yi=np.interp(xi,x,y)
    return xi, yi





def matlab2datetime(matlab_datenum):
    import datetime as dt
    day = dt.datetime.fromordinal(int(matlab_datenum))
    dayfrac = dt.timedelta(days=matlab_datenum%1) - dt.timedelta(days = 366)
    return day + dayfrac




#def GPS_parse(GPS_in, GPS_vel, GPS_vel_time, start_date, end_date):
def GPS_parse(GPS_in, GPS_vdat, start_date, end_date):
    # parse the geodetic data for input into trigger statistics
    
    import datetime
    import math
    import csv
    import pylab as pl
    import numpy as np
    
    GPS_out=[[] for GPS in range(len(GPS_in))]
    GPS_time=[[] for GPS in range(len(GPS_in))]
    g=-1
    for GPS in GPS_in:
        if GPS==[]: continue
        g+=1
        with open(GPS) as ssv:
            l=-1
            for line in ssv:
                n=0
                l+=1
                lineout=[]
                flineout=[[],[],[]] # extracts only position data
                for i in range(len(line)/25):
                    if i==0:
                        time=float(line[3:25])
                        [odtime, stime]=[int(math.floor(time/86400.0))-366, time-86400*(int(math.floor(time/86400.0)))]
                        jdate=((datetime.datetime.fromordinal(odtime)).timetuple()).tm_yday
                        timefull=str(jdate)+'.'+str(stime/86400.0)[2:] #dec julian days
                        if jdate<start_date: 
                            n=1
                            break
                        if jdate>end_date: 
                            n=1
                            break
                    else:
                        pass
                    i*=25
                    lineout.append(float(line[3+i:25+i]))
                if n==1: continue
                timeout=timefull
                [flineout[0], flineout[1], flineout[2]]=[lineout[1], lineout[2], lineout[3]] # just take positions
                GPS_out[g].append(flineout)
                GPS_time[g].append(timeout)

    # get velocity value and time
    
    GPS_vel_out=[[] for GPS in range(len(GPS_in))]
    GPS_vel_unc=[[] for GPS in range(len(GPS_in))]
    GPS_rel_unc=[[] for GPS in range(len(GPS_in))]
    GPS_vel_time_out=[[] for GPS in range(len(GPS_in))]
    g=-1
    for GPS in GPS_vdat:
        g+=1
        # from HJH
        v_data=np.genfromtxt(GPS,dtype=None,names=['time','v','vmin','vmax'])
        vel_mdatenum=v_data['time']/86400.0
        vel_pdatenum_obj=[matlab2datetime(ii) for ii in vel_mdatenum]
        vel_pdatenum=[datetime.datetime.strftime(ii,'%Y-%m-%d %H:%M:%S') for ii in vel_pdatenum_obj]
        vdatenum=(pl.datestr2num(vel_pdatenum))
        v = v_data['v']*86400.0
#        vmin = v_data['vmin']*86400.0
#        vmax = v_data['vmax']*86400.0
        # now shape it into my format...
        for day in range(len(vdatenum)):
            [odtime, stime]=[int(vdatenum[day]), (vdatenum[day]-(int(vdatenum[day])))]#*86400.0]
            jdate=((datetime.datetime.fromordinal(odtime)).timetuple()).tm_yday
            if jdate<start_date: continue
            if jdate>end_date: continue
            GPS_vel_time_out[g].append(str(jdate)+'.'+str(stime)[2:])
            GPS_vel_out[g].append(v[day])
            GPS_vel_unc[g].append(86400.0*abs(v_data['vmax'][day]-v_data['vmin'][day])/2)
            GPS_rel_unc[g].append(100*GPS_vel_unc[g][-1]/float(GPS_vel_out[g][-1]))
                
    # GPS vel is horizontal velocity to ~ 90% confidence
        
    return GPS_out, GPS_time, GPS_vel_out, GPS_vel_time_out, GPS_rel_unc




def lake_level(lake_level_data, year, start_date, end_date):

    import numpy as np
    import datetime
    import csv
    
    # get lake level data into format and timeseries window
    
    lakedata=[]
    laketime=[]
    with open(lake_level_data) as csvfile:
        for line in csv.reader(csvfile):
            [time, level] = [line[0], line[3]]
            time=datetime.datetime.strptime(time,'%Y-%m-%d %H:%M:%S')
            jday=(time.timetuple()).tm_yday
            if time.timetuple()[0]!=year: continue
            if jday<start_date: continue
            if jday>end_date: continue
            hour=(time.timetuple()[3]+time.timetuple()[4]/60.0+time.timetuple()[5]/3600.0) # add in /24.0 for not binning
            lakedata.append(float(level))
            laketime.append(str(jday)+'.'+str(hour)) # add in [2:] for hour if not binning
            
    # bin into hour
    hlast=0 # assume it begins at 0 hour
    lake_level=0
    ns=0
    lake_levels=[np.nan] # adjust for end-hour binning
    lake_time=[np.nan]
    for t in range(len(laketime)):
        jday = int((laketime[t].split('.')[0]))
        hour = float(laketime[t].split('.')[1])//1
        if hour!=hlast:
            if ns>0:
                lake_level/=ns
                lake_levels.append(lake_level)
                lake_time.append(str(jday)+'.'+str(int(hour)/24.0)[2:])
                ns=0
                hlast=hour
        elif t==len(laketime): #finalise
            lake_level/=ns
            ns=0
            lake_levels.append(lake_level)
            lake_time.append(str(jday)+'.'+str(int(hour)/24.0)[2:])
        else:
            lake_level+=float(lakedata[t])
            ns+=1
            hlast=hour

    # find rate of change with a forward difference
    
    # not using binned data is too noisy
    
#    dy=[]
#    c=0
#    delx=float(laketime[50])-float(laketime[49]) # reliable measure of x step
#    for value in lakedata:
#        dy.append((lakedata[c+1]-lakedata[c])/delx) # forward difference
#        c+=1
#        if c==len(lakedata)-1:
#            dy.append(np.nan) # no derivative at the last value
#            break

    # use bin data
    
    dy=[]
    c=0
    delx=float(lake_time[50])-float(lake_time[49])
    for value in lake_levels:
        dy.append((lake_levels[c+1]-lake_levels[c])/delx)
        c+=1
        if c==len(lake_levels)-1:
            dy.append(np.nan)
            break
#        
    dy[1]=np.nan # correct for edge effect
    dy[2]=np.nan
    lake_levels[1]=np.nan
    lake_levels[2]=np.nan
      
    return laketime, dy, lake_levels




def strain_rate(GPS_pos,GPS_vel):
    
    import math
    import numpy as np
    
    # calculate representative positions of GNSS sensors to 1 m precision

    [GPS1_E,GPS1_N,GPS2_E,GPS2_N]=[[],[],[],[]]
    for i in range(len(GPS_pos[0])):
        try:    GPS1_E.append(int(round(GPS_pos[0][i][0])))
        except: pass
        try:    GPS1_N.append(int(round(GPS_pos[0][i][1])))
        except: pass
        try:    GPS2_E.append(int(round(GPS_pos[1][i][0])))
        except: pass
        try:    GPS2_N.append(int(round(GPS_pos[1][i][1])))
        except: pass
    
    GPS1_Epos=sum(GPS1_E)/len(GPS1_E)
    GPS1_Npos=sum(GPS1_N)/len(GPS1_N)
    GPS2_Epos=sum(GPS2_E)/len(GPS2_E)
    GPS2_Npos=sum(GPS2_N)/len(GPS2_N)
    
    GPS_loc=[[GPS1_Epos, GPS1_Npos], [GPS2_Epos, GPS2_Npos]]
    
    # calculate strain rate over time using fixed sensor positions
    
    [rep_GPS1, rep_GPS2]=[[],[]]
    for a in range(len(GPS_vel[0])): 
        rep_GPS1.append([GPS1_Epos, GPS1_Npos])
        rep_GPS2.append([GPS2_Epos, GPS2_Npos])
    
    strains=[0 for x in range(len(GPS_vel[0]))]
    i=0
    for lower, upper, vlower, vupper in zip(rep_GPS1, rep_GPS2, GPS_vel[0], GPS_vel[1]):
        try:
            hsep=math.sqrt((lower[0]-upper[0])**2+(lower[1]-upper[1])**2)
            eh=(vupper-vlower)/hsep # measure how the upper glacier is approaching the lower
            strains[i]=eh
        except:
            strains[i]=np.nan
        i+=1
        
    return strains




def backwards_difference(y,delx):
    ### Takes a single function value and finds its first derivative using the backwards difference
    ### Specify the function (y) and the spacing between its indendent variable (delx)
    
    dy=[]
    c=1
    dy.append(0) #there is no derivative at the first value
    for value in y:
        dy.append((y[c]-y[c-1])/delx) #calculates the backwards difference at every x
        c+=1
        if c==len(y):
            break
    return dy
        
    
    
    
def chosen_rainfall(start_date,end_date,year,rainfall_times,rainfall_data):
    ### Takes two .csv files (tab delimited) of hours in days (ordinal, here taken as the MATLAB format, and in UTC time)
    ### -> save('temp_data','temp','-ascii','-double','-tabs') is the MATLAB command <-
    ### and hourly mm precipitation of rainfall and converts it into two variables: julian dates.hours for the window and the data within the window
    ### NOTE: the .hours is trimmed to the 7th index for funcitonality with other scripts - minimal data is lost

    import csv
    import datetime
    import math
    import sys
    import numpy as np
    import matplotlib.pyplot as plt
    
    csv.field_size_limit(sys.maxsize)
    
    # parse rainfall data
    
    data_all=[]
    with open(rainfall_data) as tsv:
        for line in csv.reader(tsv,dialect='excel-tab'):
            for i in range(len(line)):
                if i==len(line)-1: break
                data_all.append(float(line[i][3:]))
                
    # extra loop to re-float the data? Not sure this is needed
                
    rainfall_data=[]        
    for data in data_all:
        rainfall_data.append(float(data))
        
    # parse rainfall times

    rain_times_all=[]
    rain_times_day=[]
    rain_times_hour=[]
    drain=[0 for x in range(int(math.ceil(len(rainfall_data)/24.0))+400)] # add in extra days for partial start and end days
    odrain=[0 for x in range(int(math.ceil(len(rainfall_data)/24.0))+400)]
    rain_day_last=0
    rind=-1
    
    with open(rainfall_times) as tsv:
        for line in csv.reader(tsv,dialect='excel-tab'):
            for i in range(len(line)):
                if i==len(line)-1: break
                rain_times_all.append(float(line[i][3:]))
                
    # calculate decimal ordinal days for hourly rainfall and ordinary days for daily rain
    
    hrind=-1
    for time in rain_times_all:
            hrind+=1
            rain_times_day.append(int(math.floor((float(time))-366))) #takes ordinal day and applies shift from MATLAB format to python format
            hrain=str((float(time))-(math.floor(float(time))))
            rain_times_hour.append(hrain)#[:8]) #takes decimal ordinal day (hour, min, sec)
    
            # create daily rain sums
            
            if int(math.ceil((float(time))-366))==rain_day_last:
                rain_day_last=int(math.ceil((float(time))-366))
                drain[rind]+=rainfall_data[hrind]
                continue
            else:
                rind+=1
                odrain[rind]=int(math.ceil(float(time))-366) # keep an index of the ordinal days for the daily rainfall (records rainfall for the day preceeding)
                rain_day_last=int(math.ceil((float(time))-366))
                drain[rind]+=rainfall_data[hrind]
    
    #converts chosen start and end dates into ordinal dates
    
    das=JulianDate_to_MMDDYYY(year,start_date)
    dae=JulianDate_to_MMDDYYY(year,end_date)
    dns=datetime.date.toordinal(datetime.date(das[0],das[1],das[2]));
    dne=datetime.date.toordinal(datetime.date(dae[0],dae[1],dae[2]));
    
    #finds index in date from start of start and end dates for interval of interest
    
    gen = (i for i,x in enumerate(rain_times_day) if x==dns)
    for i in gen: #takes first index of the given value
        s=i  
        break
    
    gen = (i for i,x in enumerate(rain_times_day) if x==dne)
    for i in gen: #takes last index of the given value
        e=i
        
    # get indexes for daily rainfall
    
    sdfs=odrain.index(dns)
    edfs=odrain.index(dne) 
    
    odrain=odrain[sdfs:edfs+1]
    drain=drain[sdfs:edfs+1]
    
    #converts ordinal dates to julian days and then trims the date list to the chosen window (hourly)
    
    julian_rain_days_complete=[]
    for d in rain_times_day: 
        julian_rain_dates=datetime.date.fromordinal(d)
        julian_rain_days_tt=julian_rain_dates.timetuple()
        julian_rain_days_complete.append(julian_rain_days_tt.tm_yday)
    #julian_rain_days=np.unique(julian_rain_days_complete[s:e+1]) #extra variable of just julian days without hour component
    
    #add in an hour componenet to the days given
    #note: days start at 0 hours and end at 23 hours

    julian_rain_dec_days=[]    
    for i in range(e+1-s):     
        julian_rain_dec_days.append(float("{0:.8f}".format((julian_rain_days_complete[i+s])+float(rain_times_hour[i+s]))))
    
    #calculates total mm rainfall during interest period and trims the data to the chosen date window
    
    trimmed_rainfall_data=rainfall_data[s:e+1]
    total_mm=sum(rainfall_data[s:e+1])
    print('Total mm of rainfall during the chosen period is '+str(total_mm)+' mm.')
    
    # APPENDED: create julian day timing for daily rainfall data
    
    odrain_jdays=[]
    for d in odrain: 
        julian_rain_dates=datetime.date.fromordinal(d)
        julian_rain_days_tt=julian_rain_dates.timetuple()
        odrain_jdays.append(julian_rain_days_tt.tm_yday)

    return julian_rain_dec_days,trimmed_rainfall_data, odrain_jdays, drain
    
    


def diurnal_event_statistics(start_date,end_date,event_directory,delay_time):

    import matplotlib.pyplot as plt        
    import os
    from obspy import UTCDateTime
    import numpy as np
    import math
    import datetime
    
    #prepares empty lists for later appending
    
    files=[]
            
    # load in files
    
    files=os.listdir(event_directory)
    files.sort()

    # only allow .MSEED files to be in the event list    
    
    i=0
    for e in range(len(files)):
        try:
            if files[i][-6:]=='.MSEED':
                i+=1
            else:
                files.pop(files.index(files[i]))
        except:
            files.pop(files.index(files[i]))
    
    #removes successive events (stops statistical skewing due to many triggers for one long event)
    s=1 #successor counter
    n=0 #index counter
    r=0 #count number of entries removed
    for i,j in zip(files,files[1:]): #for current and next entry in zip
        try: 
            if s==1:
                first_event_time=UTCDateTime(i[:-7])
                second_event_time=UTCDateTime(j[:-7])
            elif s>1:
                second_event_time=UTCDateTime(j[:-7])
    
            if (first_event_time+delay_time)>=second_event_time:
                files.pop(n+1) #remove the successor
                first_event_time=second_event_time #adapt the remover event time to make it the same as the first
                s+=1
            else:
                s=1
                n+=1
            r+=1
        except:
            pass
        
    popped=0
    flast=0
    thours=[]
    for f in range(len(files)): #strips file names of their days and hours
       if popped==1: f=flast # backshift due to a lost index
       if f==len(files): break # this loop is complete - force exit
       popped=0
       event_file=files[f]
              
       try:
            date=datetime.datetime.strptime(event_file[:10],'%Y-%m-%d')
            datett=date.timetuple()
            jdate=datett.tm_yday
            
            # only load data up to the end date
            if jdate<start_date:
                files.pop(files.index(event_file))
                popped=1
                flast=f # necessary reset to override the for loop indexing ( a differently named variable would probably be better! )
                continue
            elif jdate>end_date:
                files.pop(files.index(event_file))
                popped=1
                flast=f
                continue 
       except:
            pass
            
       try: #pull the closest hour to the event from the file and append it to a list
           hour=int(math.ceil(float(event_file[11:13])+(float(event_file[14:16]))/60+(float(event_file[17:26]))/3600)) #bins events to the closest hour  
           thours.append(hour)
       except:
           print('I do not think that is a event file.')
           continue #otherwise ignore it  
            
    hour_sum=[0]*24 #list of 0 values for the 24 hours in a day, starting at 12am 
    hour_index=range(24)
    ihours=range(25) #include the 24th hour
    
    for hour in thours:
        c=0    
        for ihour in ihours:
            if hour==ihour:
                hour_NZT=hour+12 #get it in NZT (standard time)
                if hour_NZT>=24: c-=24 #shift it back
                cloc=hour_index.index(c+12)
                hour_sum[cloc]+=1 #increase the number of counts for that hour (set to shift)
            c+=1
            
    # 0 is all events occuring in the 23rd hour (23pm-12am, "observation" at hour 24 AKA 0)
    # 1 is all events occuring in the 1st hour (12am-1am, "observation" at hour 1)

    
    plt.bar(range(24),hour_sum,width=1,color='red',lw=4)
    plt.xlim(0,23)
    plt.tick_params(axis='both',which='major',labelsize=12)
    plt.xlabel('Hour (NZST)',fontsize=12, labelpad=15)
    plt.ylabel('Icequake Detections',fontsize=12, labelpad=15)
    #plt.title('Diurnal Distribution of Impulsive Seismic Events', fontsize=18, fontweight='bold',y=1.03)
    ax=plt.gca()
    ax.set_xticks(np.arange(3,24,3))
    ax.set_yticks(np.arange(0,25,5))
#    plt.savefig('IT5_diurnal_event_trend.png',format='eps',dpi=1000)    
    plt.show()
    
 

   
def isclose(a, b, rel_tol=0.001, abs_tol=0.0):
    #Python 2.7 version of math.isclose
    return abs(a-b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)  



    
def JulianDate_to_MMDDYYY(y,jd):
    ###Converts year and julian day into year, month, day
    ###Taken from http://www.science-emergence.com/Articles/MODIS-Convert-Julian-Date-to-MMDDYYYY-Date-format-using-python/ on 14/10/16
    import calendar
    
    month = 1
    while jd - calendar.monthrange(y,month)[1] > 0 and month <= 12:
        jd = jd - calendar.monthrange(y,month)[1]
        month = month + 1
    return y,month,jd
    
    
  
  
def roundtohalf(x):
    #round coordinate to nearest half hour, adapted from Alok Singhai on stackoverflow.com on 18/11/16
    import math
    cellsize=0.5
    rounded = cellsize * math.floor(float(x)/cellsize)
    if (rounded/cellsize) %2==0: rounded+=0.5
    return rounded
    
    
    
def time_reconciler(hourly_event_sum,hourly_event_hours,dec_julian_days,start_date,end_date):
    #converts disjointed hourly event record to an "all-time" record for 
    #compatability with rainfall data
    #NOTE: first input is from event counter, second input is from rainfall
    #(or just all decimal julian days for the window length)
    
    hourly_event_sum_all_hours=[0]*(((end_date-start_date)*24)+24)
    h_in_d=0
    i=0 #index counter
    for h in hourly_event_hours:
        if float(h)>end_date: break
        h=float("{0:.8f}".format(float(h)))#[:9] #trim string
        h_in_d=dec_julian_days.index(h)
        hourly_event_sum_all_hours[h_in_d]=hourly_event_sum[i]
        i+=1

    #output should function with decimal julian day record (i.e. rainfall timing)
    return hourly_event_sum_all_hours
        
        
        
        
def event_counter(start_date,end_date,event_directory,delay_time):
    ### Takes event names in a directory and counts the total number on a given day (for a single year)
    
    import datetime
    import os
    from obspy import UTCDateTime
    import numpy as np
    import matplotlib.pyplot as plt
    import math
    
    #import and organise event "meta"data and strip event dates and hours from it    
    
    #prepares empty lists for later appending
    
    files=[]
    ttimes=[]
            
    # load in files
    
    files=os.listdir(event_directory)
    files.sort()

    # only allow .MSEED files to be in the event list    
    
    i=0
    for e in range(len(files)):
        try:
            if files[i][-6:]=='.MSEED':
                i+=1
            else:
                files.pop(files.index(files[i]))
        except:
            files.pop(files.index(files[i]))
    
    #removes successive events (stops statistical skewing due to many triggers for one long event)
    s=1 #successor counter
    n=0 #index counter
    r=0 #count number of entries removed
    for i,j in zip(files,files[1:]): #for current and next entry in zip
        try: 
            if s==1:
                first_event_time=UTCDateTime(i[:-7])
                second_event_time=UTCDateTime(j[:-7])
            elif s>1:
                second_event_time=UTCDateTime(j[:-7])
    
            if (first_event_time+delay_time)>=second_event_time:
                files.pop(n+1) #remove the successor
                first_event_time=second_event_time #adapt the remover event time to make it the same as the first
                s+=1
            else:
                s=1
                n+=1
            r+=1
        except:
            pass
        
    popped=0
    flast=0
    for f in range(len(files)): #strips file names of their days and hours
       if popped==1: f=flast # backshift due to a lost index
       if f==len(files): break # this loop is complete - force exit
       popped=0
       event_file=files[f]
              
       try:
            date=datetime.datetime.strptime(event_file[:10],'%Y-%m-%d')
            datett=date.timetuple()
            jdate=datett.tm_yday
            
            # only load data up to the end date
            if jdate<start_date:
                files.pop(files.index(event_file))
                popped=1
                flast=f # necessary reset to override the for loop indexing ( a differently named variable would probably be better! )
                continue
            elif jdate>end_date:
                files.pop(files.index(event_file))
                popped=1
                flast=f
                continue 

            #events binned by the hour following...
            hour=math.ceil(float(event_file[11:13])+(float(event_file[14:16]))/60+(float(event_file[17:26]))/3600)
            hour=(float(hour)/24) #convert hours into decimal days
            ttimes=np.append(ttimes,str((float(jdate)+hour))) #special event time format for hourly events (format: julian_date.hour)
       except:
            #print('I do not think that is a event file.')
            continue #otherwise ignore it
    ttimes=np.sort(ttimes)  
    ttimes_all=ttimes # all event times ...
    ttimes_unique=np.unique(ttimes) #collapses array into days.hours 

    #put all event times through an indexing routine to count the number of events in a given hour       
       
    idates=range((end_date-start_date+1)) #establishes index days
    ihours=range(24)
    numtrigdhr=0 #initialise
    events_per_hour=np.zeros(len(ttimes_unique)) #a list of zeros corresponding to the index days and their hours

    for entry in ttimes: #looks at each day.hour in the list of julian dates and hours taken from files
        if numtrigdhr>1: #skip the events that have been removed (but haven't in the iterator)
            numtrigdhr-=1
            continue
            
        m=0 #event match off
        date=int(entry[:3]) #reiterates date every time the inner loop(s) finishes, same for hour
        hour=int(round((float(entry[3:]))*24)) #reinflate hour (0-23)
        d=0 #clear day count
        h=0 #clear hour count
        
        for idate in idates: #look through all index dates
            idate+=start_date #set each index date into the date window of interest
            if date==idate: #if the date from the event file matches the index date
                for ihour in ihours: #then look at all index hours
                
                    if hour==ihour: #if the hour from the event file matches the index hour
                        dh=np.where(ttimes_unique==(str(date)+entry[3:]))[0]
                        dhr=np.where(ttimes==(str(date)+entry[3:]))[0]
                        numtrigdhr=len(dhr)
                        events_per_hour[dh]+=len(dhr) #then increment the number of events at the index position representing the given date and hour from the event file
                        ttimes=np.delete(ttimes,dhr) #then remove that value - it has been found
                        m=1 #event match on
                        break #match found - now look at the next entry
                    else:
                        h+=1 #otherwise increment the hour counter and check the next hour                    
                if m==1:
                    break #continue the break
            else:
                d+=1 #otherwise increment the day counter and check the next day
        if m==1:
            continue #finish the break
            
    #sub-routine to count the total number of events in the events_per_hour file (QC)        
    
    cumulative_events=accumulate_data(events_per_hour)
    print('The total number of events is: '+str(cumulative_events[len(cumulative_events)-1]))
    
    return events_per_hour,ttimes_unique, files, ttimes_all




def smooth_data(data,symmetric_smooth_range):
    # smooth data by taking a running mean from a data point out to the data 1/2 the smooth range away
    # NOTE: the symmetric_smooth_range must be even - it is the number of data points to sample on either side of the central data
    import numpy as np
    smooth_data=[0 for x in range(len(data))] 
    for smooth_cen in range(len(data)):
        missing=0
        smooth_cen+=symmetric_smooth_range # skip to the first viable option (and so on)
        if smooth_cen==len(data)-symmetric_smooth_range: break # leave before the smoothing fails
        if str(data[smooth_cen])!='nan':
            s=data[smooth_cen] # data sum (starting with central data point)
        else:
            s=0
            missing=1
        for r in range(symmetric_smooth_range): # look at data on either side of the central data point
            # add this data to the sum
            if str(data[smooth_cen-r])!='nan':
                s+=data[smooth_cen-r]
            else:
                missing+=1
            if str(data[smooth_cen+r])!='nan':
                s+=data[smooth_cen+r]
            else:
                missing+=1
        # append the mean to the smoothed data list
        try:
            smooth_data[smooth_cen]=(s/(2*symmetric_smooth_range+1-missing))
        except:
            smooth_data[smooth_cen]=np.nan
    
    return smooth_data
            
            
        
        
def central_difference(u,dx):
    # find the first and second derivatives using the central difference
    
    du=[]
    d2u=[]
    for j in range(len(u)):
        if j==0 or j==(len(u)-1): #central difference is not defined for start or end point
            du.append(0) #no data value
            d2u.append(0) #no data value
            continue
        else:
            du.append(((u[j+1]-u[j-1])/(2*dx))) #defined as the absolute value to make later computation easier
            d2u.append((u[j+1]-(2*u[j])+u[j-1])/(dx**2))
            
    return du,d2u



def accumulate_data(events_per_hour):
    
    # calculate cumulture event # at hourly resolution
    
    c=1
    count=1
    cumulative_events=[]
    cumulative_events.append(events_per_hour[0]) #pre-set the first entry so subsequent entries can add it to their value
    for entry in events_per_hour: #look at every entry in the list of event counts per hour
        cumulative_events.append(events_per_hour[c]+cumulative_events[c-1]) #add the event count for a given hour to all of those preceeding
        c+=1 #increment the position in the list
        if c==len(events_per_hour): #when the index reaches the end of the list terminate the loop
            break
        count+=events_per_hour[c] #increment the number of events by the number in that hour
    
    return cumulative_events
    
    
    
    
def accumulate_events(filenames):
    
    # calculate cumulate event # at event resolution    
    
    import datetime    
    
    count=0
    cumulative_events=[]
    data_times=[]
    for i in range(len(filenames)):
            # get event time
            timein=filenames[i][:19]
            timein=datetime.datetime.strptime(timein,'%Y-%m-%dT%H:%M:%S')
            hfrac=(timein.timetuple()[3]+timein.timetuple()[4]/60.0+timein.timetuple()[5]/3600.0)/24.0
            jday=datetime.datetime.timetuple(timein).tm_yday
            event_time=float(jday+hfrac) 
            # calculate cumulate count
            count+=1
            cumulative_events.append(count)
            data_times.append(event_time)
        
    return cumulative_events, data_times




def roundtocell(x, cellsize):
    #round coordinate to nearest cell, adapted from Alok Singhai on stackoverflow.com on 18/11/16
    ans=cellsize * round(float(x)/float(cellsize))
    str_ans=str(ans)
    ans=float(str_ans[:6]) # only allow 1 dp for km NZTM
    return ans
    
    
    
    
def xcorrlag(dat1,dat2,spacing):
    #HJH code
    import numpy as np
    # time series at equal spacing
    # INPUTS    dat1 time series 1
    #           dat2 time series 2
    #           spacing dt (must be equal)

    cc = np.correlate(dat1,dat2,"full")
    
    breaker=0
    for entry in cc:
        if str(entry)=='nan':
            breaker=1
            break
    if breaker==1: return np.nan, np.nan
    
#    print cc
#    print len(cc)
    maxcorr = max(cc)
    mi = [i for i, j in enumerate(cc) if j == maxcorr] # max index
    # zero lag at len(dat1) (if both vectors same length)
    # len(cc) = len(dat1)+len(dat2)-1 when cc "full"
    # max correlation
    
#    lagtime=mi[0] # number of samples from the start of the window
    
#    if ccmode == "rain":
    lagtime=(len(dat1)-mi[0]-1) # in hours
#    elif ccmode == "lake":
#        lagtime=spacing*(len(lakei)-mi[0]-1)

    return lagtime, maxcorr




def drainage_model(outdir, start_date,end_date,rainfall_times, rainfall_data, crevassing_times, crevassing_rate):
    
    import matplotlib.pyplot as plt
    import numpy as np
    
    # model that tries to match rising crevassing peaks with parameters
    # simple version assumes infinite storage maxima over time
    # and that the glacial drainage system cannot overflow into the crevasse
    
    # drainage rate is in mm/hr:
        
    Do_S_all=[[] for x in np.linspace(0.1,2,20)]
    b_S_all=[[] for x in np.linspace(0.01,0.25,25)]
    a_S_all=[[] for x in np.linspace(0.01,0.25,25)]
    
    Do_D_all=[[] for x in np.linspace(0.1,2,20)]
    b_D_all=[[] for x in np.linspace(0.01,0.25,25)]
    a_D_all=[[] for x in np.linspace(0.01,0.25,25)]
    
    di=-1
    for Do in np.linspace(0.1,2,20):
        di+=1
        ai=-1
        for a in np.linspace(0.01,0.25,25):
            ai+=1
            bi=-1
            for b in np.linspace(0.01,0.25,25):
                bi+=1
                print('Model parameters are: '+str(Do)+' , '+str(a)+' , '+str(b))
                # initialise model:
                I=0 # current input
                D=0 # current drainage
                S=0 # current storage
                Sl=0 # previous storage
                Dmax=Do # current drainage maxima
                
                # calculate storage over all times:
                S_all=[0 for x in range(len(rainfall_times))]
                D_all=[0 for x in range(len(rainfall_times))]
                for t in range(len(rainfall_times)):
                    I=rainfall_data[t] # input is rainfall
                    # calculate drainage
                    if I<Dmax: # rain rate is less than drainage rate
                        if Sl==0: # there is no storage
                            D=I
                        elif Sl>=(Dmax-I): # storage is sufficient to fill to maximum
                            Ds=Dmax-I
                            D=I+Ds
                        elif Sl<(Dmax-I): # storage is insufficient "
                            Ds=Sl
                            D=I+Ds
                    elif I>=Dmax: # rain rate is more than drainage rate
                        D=Dmax
                    # calculate storage
                    S=Sl+I-D
                    S_all[t]=S # register this to a time
                    D_all[t]=Dmax
                    # now save the storage and drainage maxima to last time step
                    Sl=S
                    Dmaxl=Dmax
                    # adjust the drainage maxima
                    if (Dmaxl+a*D-b*Dmaxl)>Do:
                        Dmax=Dmaxl+a*D-b*Dmaxl
                    else:
                        Dmax=Do
        
                b_S_all[bi]=S_all
                b_D_all[bi]=D_all
            a_S_all[ai]=b_S_all
            a_D_all[ai]=b_D_all
        Do_S_all[di]=a_S_all
        Do_D_all[di]=a_D_all

    # now compare all storage to the crevassing rate
    
    # find storage maximum for consistent plotting
    all_storage=[]
    for s1 in Do_S_all:
        for s2 in s1:
            for s3 in s2:
                all_storage=all_storage+s3
    s_t_max=max(all_storage)
    
    all_drainage=[]
    for d1 in Do_D_all:
        for d2 in d1:
            for d3 in d2:
                all_drainage=all_drainage+d3
    d_t_max=max(all_drainage)
    
    # reference lists for titles etc.
        
    a_all=np.linspace(0.01,0.25,25)
    b_all=np.linspace(0.01,0.25,25)
    d_all=np.linspace(0.1,2,20)
                
    di=-1
    for Do_S_t in Do_S_all:
        di+=1
        Do_D_t=Do_D_all[di]
        ai=-1
        for a_S_t in Do_S_t:
            ai+=1
            a_D_t=Do_D_t[ai]
            bi=-1
            for S_t in a_S_t:
                bi+=1
                b_D_t=a_D_t[bi]
                print('Model parameters for this plot are: '+str(d_all[di])+' , '+str(a_all[ai])+' , '+str(b_all[bi]))
                
                # this is super crude, but literally just plot them for manual comparison
                rain_major_ticks=np.arange(0,roundtocell(int(max(rainfall_data)),5),5)
                rain_minor_ticks=np.arange(0,roundtocell(int(max(rainfall_data)),5),1)
                drainage_major_ticks=np.arange(0,roundtocell(int(d_t_max),20),20)
                drainage_minor_ticks=np.arange(0,roundtocell(int(d_t_max),20),10)
                storage_major_ticks=np.arange(0,roundtocell(int(s_t_max),20),20)
                storage_minor_ticks=np.arange(0,roundtocell(int(s_t_max),20),10)
                events_major_ticks=np.arange(0,roundtocell(max(crevassing_rate),5)+5,5)
                events_minor_ticks=np.arange(0,roundtocell(max(crevassing_rate),5)+5,1)
                decimal_days_major_ticks=np.arange(start_date,end_date,1)
                decimal_days_minor_ticks=np.arange(start_date,end_date,24)
            
                #plotting...
        #        bar_width=0.2 # use this for all-time plot
                bar_width=0.02 # use this for zooming in
                line_width=1
                
                plt.figure(figsize=(11.69,8.27))
                plt.subplot(411)
                ax1=plt.gca() 
                
                ax1.bar(rainfall_times,rainfall_data,color='b',width=bar_width,lw=line_width)
                ax1.set_ylabel('Rain Rate\n(mm $hr^-1$)',fontsize=10, labelpad=25)
                ax1.set_xticks(decimal_days_major_ticks)
                ax1.set_xticks(decimal_days_minor_ticks,minor='True')
                ax1.set_xlim(start_date,end_date)
                ax1.set_yticks(rain_major_ticks)
                ax1.set_yticks(rain_minor_ticks,minor='True')
                ax1.set_ylim(0,roundtocell(max(rainfall_data),5))
                ax1.tick_params(axis='both',which='major',labelsize=10)
                ax1.tick_params(axis='x',which='minor',labelsize=None)
                ax1.get_yaxis().set_label_coords(-0.06,0.5)
                ax1.grid(which='major',alpha=0.5)
                ax1.grid(which='minor',alpha=0.25)
                
                plt.subplot(412, sharex=ax1)
                ax2=plt.gca()
                
                ax2.bar(rainfall_times,b_D_t,color='green',width=bar_width,lw=line_width)
                ax2.set_ylabel('Drainage Rate\n(mm $hr^-1$)',fontsize=10, labelpad=25)
                ax2.set_xticks(decimal_days_major_ticks)
                ax2.set_xticks(decimal_days_minor_ticks,minor='True')
                ax2.set_xlim(start_date,end_date)
                ax2.set_yticks(drainage_major_ticks)
                ax2.set_yticks(drainage_minor_ticks,minor='True')
                ax2.set_ylim(0,roundtocell(int(d_t_max),20))
                ax2.tick_params(axis='both',which='major',labelsize=10)
                ax2.tick_params(axis='x',which='minor',labelsize=None)
                ax2.get_yaxis().set_label_coords(-0.06,0.5)
                ax2.grid(which='major',alpha=0.5)
                ax2.grid(which='minor',alpha=0.25)
                
                plt.subplot(413, sharex=ax1)
                ax3=plt.gca()
                
                ax3.bar(rainfall_times,S_t,color='purple',width=bar_width,lw=line_width)
                ax3.set_ylabel('Storage\n(mm)',fontsize=10, labelpad=25)
                ax3.set_xticks(decimal_days_major_ticks)
                ax3.set_xticks(decimal_days_minor_ticks,minor='True')
                ax3.set_xlim(start_date,end_date)
                ax3.set_yticks(storage_major_ticks)
                ax3.set_yticks(storage_minor_ticks,minor='True')
                ax3.set_ylim(0,roundtocell(s_t_max,20))
                ax3.tick_params(axis='both',which='major',labelsize=10)
                ax3.tick_params(axis='x',which='minor',labelsize=None)
                ax3.get_yaxis().set_label_coords(-0.06,0.5)
                ax3.grid(which='major',alpha=0.5)
                ax3.grid(which='minor',alpha=0.25)
                
                plt.subplot(414, sharex=ax1)
                ax4=plt.gca()
                
#                print crevassing_rate
#                print rainfall_times
#                print len(crevassing_rate)
#                print len(rainfall_times)
                
                ax4.bar(rainfall_times,crevassing_rate,color='r',width=bar_width,lw=line_width)
                ax4.set_xlabel('Julian Day (2016)', fontsize=10, labelpad=25)
                ax4.set_ylabel('Crevassing Event Count\n(number $hr^-1$)',fontsize=10, labelpad=25)
                ax4.set_xticks(decimal_days_major_ticks)
                ax4.set_xticks(decimal_days_minor_ticks,minor='True')
                ax4.set_xlim(start_date,end_date)
                ax4.set_yticks(events_major_ticks)
                ax4.set_yticks(events_minor_ticks,minor='True')
                ax4.set_ylim(0,roundtocell(max(crevassing_rate),5)+5)
                ax4.tick_params(axis='both',which='major',labelsize=10)
                ax4.tick_params(axis='x',which='minor',labelsize=None)
                ax4.get_xaxis().set_label_coords(0.5,-0.11)
                ax4.get_yaxis().set_label_coords(-0.06,0.5)
                ax4.grid(which='major',alpha=0.5)
                ax4.grid(which='minor',alpha=0.25)      
                
                ax1.set_title('Crevasse Storage Model Results\n Do = '+str(d_all[di])+ ', a = '+str(a_all[ai])+', b = '+str(b_all[bi]), fontsize=16, y=1.03)
                
                plt.savefig(outdir+'/STORAGE/'+str(d_all[di])+'_'+str(a_all[ai])+'_'+str(b_all[bi])+'_mmPhr.pdf',format='pdf',dpi=300)
                plt.close('all')
                
#from scipy import signal
#signals=[signal.gaussian(24, 6, sym=True)[i] for i in range(len(signal.gaussian(24, 6, sym=True)))]
#drainage_model('/Volumes/arc_02/taylorsa/gaussianmodel', 0,3,range(0,72),[0,0,0,0,0,0]+signals+[0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0], range(0,72), 24*[0]+[2,4,8,16,4,2]+42*[0])

def run_event_statistics():        
    
    # this version made for smoothed data
    
    #this runs the stand-alone event statistics script. Otherwise the script
    #acts as a number of functions to be called for more advanced processes
               
    import csv 
    import matplotlib.pyplot as plt
    import numpy as np
    np.set_printoptions(threshold=np.nan)
    import math
    import scipy
    import datetime
    import pylab as pl
    
    start_date=112
    end_date=139
    year=2016
    hour_shift=0.5 #shift all times by a certain amount. Implemented to plot data at half past the hour to represent time collected at
    delay_time=0
    background_rate=0 # number of events expected in an any given hour in the event catalogue
    
#    event_directory='/Volumes/arc_01/taylorsa/triggers/IQ'
#    event_directory='/Volumes/arc_01/taylorsa/catalogues/final_catalogue/crevassing'
#    v2_event_directory='/Volumes/arc_01/taylorsa/triggers/IQ'
#    v3_event_directory='/Volumes/arc_01/taylorsa/catalogues/final_catalogue/local'

#    event_directory='/Volumes/arc_02/taylorsa/xcorr_detections/t2_v2'
#    v2_event_directory='/Volumes/arc_02/taylorsa/xcorr_detections/t39_v2'

    event_directory='/home/sam/EVENTS_IT3/TYPE_A/4/'
#    v2_event_directory='/Volumes/arc_01/taylorsa/catalogues/crevassing_cat_v2/sustained'

#    event_directory='/Volumes/arc_02/taylorsa/internal_events/'
#    v2_event_directory='/Volumes/arc_02/taylorsa/marginal_events'
    
    rainfall_times='/Volumes/arc_01/'
    rainfall_data='/Volumes/arc_01/taylorsa/rainfall_data/rainfall_hours_data' 
#    temperature_times='/Volumes/arc_01/taylorsa/rainfall_data/temp_hours'
#    temperature_data='/Volumes/arc_01/taylorsa/rainfall_data/temp_data'
#    lake_level_data='/Volumes/arc_01/FIELD_DATA/TASMAN_LAKE/tasman_lake_level.csv'
    
#    diurnal_event_statistics(start_date,end_date,event_directory,delay_time)
    
    # lake level data
    
#    lake_times, lake_roc, lake_levels = lake_level(lake_level_data, year, start_date, end_date)
    
    # GPS processing and plotting addition...
                                                   
    GPS1_path='/Volumes/arc_01/taylorsa/velocity_data/run1/tas_position_arc5.txt'
    GPS2_path='/Volumes/arc_01/taylorsa/velocity_data/run1/tas_position_arc4.txt'
#    GPS1_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c1_data_3hr_30samp'
#    GPS2_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c3_data_3hr_30samp'
#    GPS1_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c1_data_6hr_30samp'
#    GPS2_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c3_data_6hr_30samp'
    GPS1_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c1_data_24hr_120samp'
    GPS2_vpath='/Volumes/arc_01/taylorsa/velocity_data/run2/c3_data_24hr_120samp'
    GPS_in=[GPS1_path, GPS2_path]
    GPS_vdat=[GPS1_vpath, GPS2_vpath]
    GPS_out, GPS_time, GPS_vel, GPS_vel_time, GPS_vel_unc = GPS_parse(GPS_in, GPS_vdat, start_date, end_date)
    strains=strain_rate(GPS_out, GPS_vel)
    
#    print GPS_vel_unc
    
#    vel_95p_mean = [0, 0]
#    vel_95p_std = [0, 0]
#    vel_95p_mean[0] = np.nanmean(GPS_vel_unc[0])
#    vel_95p_mean[1] = np.nanmean(GPS_vel_unc[1])
#    vel_95p_std[0] = np.nanstd(GPS_vel_unc[0])
#    vel_95p_std[1] = np.nanstd(GPS_vel_unc[1])
#    
#    print vel_95p_mean[0]+2*vel_95p_std[0]
#    print vel_95p_mean[1]+2*vel_95p_std[1]
    
#    sys.exit()
#    
#     # look at direction of horizontal velocity for the two stations
# 
#    vel_dir='/Volumes/arc_01/taylorsa/velocity_data/run2/'
#    
##    veldat=[vel_dir+'c1_dir_24hr_120samp', vel_dir+'c3_dir_24hr_120samp']
#    veldat=[vel_dir+'c1_dir_3hr_30samp', vel_dir+'c3_dir_3hr_30samp']
#        
#    # bastardise a different code to save time (see: trigger_statistics)
#    
#    GPS_vel_out=[[] for GPS in range(len(veldat))]
#    GPS_vel_time_out=[[] for GPS in range(len(veldat))]
#    g=-1
#    for GPS in veldat:
#        g+=1
#        # from HJH
#        v_data=np.genfromtxt(GPS,dtype=None,names=['time','v'])
#        vel_mdatenum=v_data['time']/86400.0
#        vel_pdatenum_obj=[matlab2datetime(ii) for ii in vel_mdatenum]
#        vel_pdatenum=[datetime.datetime.strftime(ii,'%Y-%m-%d %H:%M:%S') for ii in vel_pdatenum_obj]
#        vdatenum=(pl.datestr2num(vel_pdatenum))
#        v = v_data['v']*86400.0
#        # now shape it into my format...
#        for day in range(len(vdatenum)):
#            [odtime, stime]=[int(vdatenum[day]), (vdatenum[day]-(int(vdatenum[day])))]#*86400.0]
#            jdate=((datetime.datetime.fromordinal(odtime)).timetuple()).tm_yday
#            if jdate<start_date: continue
#            if jdate>end_date: continue
#            GPS_vel_time_out[g].append(str(jdate)+'.'+str(stime)[2:])
##            GPS_vel_out[g].append(180-(v[day]*180)) # coordinate system is with "North" in the South (flow) direction
#            GPS_vel_out[g].append(180-(v[day]*180/(2*math.pi)))
#
#    GPS_vel_dir=GPS_vel_out
#    GPS_vel_dir_time_out=GPS_vel_time_out

    # calculate the rainfall in each hour and each day

    dec_julian_days,rainfall_data, julian_days, day_rainfall=chosen_rainfall(start_date,end_date,year,rainfall_times,rainfall_data)
    
    # calculate the number of events occuring in each hour
    
    hourly_event_sum,hourly_event_hours,files,all_times=event_counter(start_date,end_date,event_directory,delay_time)
#    print hourly_event_hours
    v2_hourly_event_sum,v2_hourly_event_hours,v2_files,v2_all_times=event_counter(start_date,end_date,v2_event_directory,delay_time)
#    v3_hourly_event_sum,v3_hourly_event_hours,v3_files,v3_all_times=event_counter(start_date,end_date,v3_event_directory,delay_time)
    
    # remove event catalogue background rate from each hour
    
    for x in range(len(hourly_event_sum)): 
        hourly_event_sum[x]-=background_rate
        if hourly_event_sum[x]<0: hourly_event_sum[x]=0
    
    #rate_of_event=backwards_difference(cumulative_events,1)
   
    #reorgansie times for events to work with rainfall timing
    
    hourly_event_sum_all_hours=time_reconciler(hourly_event_sum,hourly_event_hours,dec_julian_days,start_date,end_date)
    v2_hourly_event_sum_all_hours=time_reconciler(v2_hourly_event_sum,v2_hourly_event_hours,dec_julian_days,start_date,end_date)
#    v3_hourly_event_sum_all_hours=time_reconciler(v3_hourly_event_sum,v3_hourly_event_hours,dec_julian_days,start_date,end_date)
    
    # smooth rainfall, GPS, seismic data...
    
#    smooth_rainfall_data=smooth_data(rainfall_data,12)  
#    smooth_seismic_data=smooth_data(hourly_event_sum_all_hours,12)
#    smooth_v2seismic_data=smooth_data(v2_hourly_event_sum_all_hours,12)
#    smooth_v3seismic_data=smooth_data(v3_hourly_event_sum_all_hours,12)
#    smooth_hvel_lower=smooth_data(GPS_hvel[0],12)
#    smooth_hvel_upper=smooth_data(GPS_hvel[1],12)
#    smooth_hstrain=smooth_data(horstrains, 12)
    
    # hour_shift correction to events and rainfall data (plots at the midpoint of the hour)
    c=0
    dec_julian_days_shifted=[0]*len(dec_julian_days)
    if hour_shift!=0: #if there is a shift to be applied
        for day in dec_julian_days: #this will shift the timing for events and rain (both bin to the start of the hour)
            dec_julian_days_shifted[c]=float("{0:.8f}".format((day-hour_shift/24)))##[:8]
            c+=1
    #observation is made at 0900 for 0800-0900, so plot at 0830
    
    # get position data from GPS lists
    [x1pos,y1pos,z1pos]=[[],[],[]]
    [x2pos,y2pos,z2pos]=[[],[],[]]

    for hour in range(len(GPS_out[0])):
        try:
            [x1, y1, z1]=[GPS_out[0][hour][0], GPS_out[0][hour][1], GPS_out[0][hour][2]]
            x1pos.append(x1)
            y1pos.append(y1)
            z1pos.append(z1)
        except:
            x1pos.append(np.nan)
            y1pos.append(np.nan)
            z1pos.append(np.nan)
            
    for hour in range(len(GPS_out[1])):
        try:
            [x2, y2, z2]=[GPS_out[1][hour][0], GPS_out[1][hour][1], GPS_out[1][hour][2]]
            x2pos.append(x2)
            y2pos.append(y2)
            z2pos.append(z2)
        except:
            x2pos.append(np.nan)
            y2pos.append(np.nan)
            z2pos.append(np.nan)
            
    # de-mean the vertical positions
    
    z1pos=np.asarray(z1pos)
    z2pos=np.asarray(z2pos)
    z1posnum=z1pos[~np.isnan(z1pos)] # remove nan values
    z2posnum=z2pos[~np.isnan(z2pos)]
    
    [z1mean, z2mean]=[np.mean(z1posnum), np.mean(z2posnum)]
    [z1dmean, z2dmean]=[[],[]]
    for pos in z1pos:
        z1dmean.append(pos-z1mean)
    for pos in z2pos:
        z2dmean.append(pos-z2mean)        
            
#    print z1mean
#    print z2mean
        
    # now do a running mean for better positions
    
#    z1pos=runmean(z1dmean,360) # 3 sample running mean
#    z2pos=runmean(z2dmean,360) 

    z1pos=z1dmean
    z2pos=z2dmean
    
#    smooth_uplift_lower=smooth_data(z1pos,12)
#    smooth_uplift_upper=smooth_data(z2pos,12)

#    ### do cross-correlation of a few datasets
#    
    # interpolate velocity data to the 1 hour timeframe
    
    centre_days=[125, 130, 132, 133, 135, 136]
    
    vu_ra=[[] for x in range(len(centre_days))]
    vl_ra=[[] for x in range(len(centre_days))]
    vu_vl=[[] for x in range(len(centre_days))]
    cr_ra=[[] for x in range(len(centre_days))]
    vu_cr=[[] for x in range(len(centre_days))]
    vl_cr=[[] for x in range(len(centre_days))]
    ra_st=[[] for x in range(len(centre_days))]
    st_cr=[[] for x in range(len(centre_days))]
    
    d=-1
    for day in centre_days:
        d+=1
        interp_start=day-1
        interp_end=day+1
        interpx=np.linspace(interp_start,interp_end,(interp_end-interp_start)*24+1)
        
        interpx,int_GPS_vel_lower=interp1d(GPS_vel_time[0], GPS_vel[0], interpx)
        interpx,int_GPS_vel_upper=interp1d(GPS_vel_time[1], GPS_vel[1], interpx)
        interpx,int_strain=interp1d(GPS_vel_time[0], strains, interpx)
        
        # smooth other data to the 24 hour running mean
        
        smooth_rainfall_data=smooth_data(rainfall_data[dec_julian_days.index(interp_start):dec_julian_days.index(interp_end)+1],3)
        smooth_hourly_event_sum_all_hours=smooth_data(hourly_event_sum_all_hours[dec_julian_days.index(interp_start):dec_julian_days.index(interp_end)+1],3)    
        
        # normalise all inputs
        
#        int_GPS_vel_lower = [float(i)/max(int_GPS_vel_lower) for i in int_GPS_vel_lower]
#        int_GPS_vel_upper = [float(i)/max(int_GPS_vel_upper) for i in int_GPS_vel_upper]
#        int_strain = [float(i)/max(int_strain) for i in int_strain]
#        smooth_rainfall_data = [float(i)/max(smooth_rainfall_data) for i in smooth_rainfall_data]
#        smooth_hourly_event_sum_all_hours = [float(i)/max(smooth_hourly_event_sum_all_hours) for i in smooth_hourly_event_sum_all_hours]

#        int_GPS_vel_lower = [(float(i)-np.mean(int_GPS_vel_lower))/np.std(int_GPS_vel_lower) for i in int_GPS_vel_lower]
#        int_GPS_vel_upper = [(float(i)-np.mean(int_GPS_vel_upper)) /np.std(int_GPS_vel_upper)for i in int_GPS_vel_upper]
#        int_strain = [(float(i)-np.mean(int_strain))/np.std(int_strain) for i in int_strain]
#        smooth_rainfall_data = [(float(i)-np.mean(smooth_rainfall_data))/np.std(smooth_rainfall_data) for i in smooth_rainfall_data]
#        smooth_hourly_event_sum_all_hours = [(float(i)-np.mean(smooth_hourly_event_sum_all_hours))/np.std(smooth_hourly_event_sum_all_hours) for i in smooth_hourly_event_sum_all_hours]
        
#        print int_GPS_vel_lower
#        print '\n'
#        print int_GPS_vel_upper
#        print '\n'
#        print int_strain
#        print '\n'
#        print smooth_rainfall_data
#        print '\n'
#        print smooth_hourly_event_sum_all_hours
        
        [int_GPS_vel_upper_sqsum, int_GPS_vel_lower_sqsum, int_strain_sqsum, smooth_rainfall_data_sqsum, smooth_hourly_event_sum_all_hours_sqsum]=[0,0,0,0,0]
        for i in int_GPS_vel_upper: int_GPS_vel_upper_sqsum+=i**2
        for i in int_GPS_vel_lower: int_GPS_vel_lower_sqsum+=i**2
        for i in int_strain: int_strain_sqsum+=i**2
        for i in smooth_rainfall_data: smooth_rainfall_data_sqsum+=i**2
        for i in smooth_hourly_event_sum_all_hours: smooth_hourly_event_sum_all_hours_sqsum+=i**2   
        
        # cross correlate data
        
        vu_ra[d]=[xcorrlag(int_GPS_vel_upper,smooth_rainfall_data,1/24.0)[0],xcorrlag(int_GPS_vel_upper,smooth_rainfall_data,1/24.0)[1]/math.sqrt(int_GPS_vel_upper_sqsum*smooth_rainfall_data_sqsum)]
        vl_ra[d]=[xcorrlag(int_GPS_vel_lower,smooth_rainfall_data,1/24.0)[0],xcorrlag(int_GPS_vel_lower,smooth_rainfall_data,1/24.0)[1]/math.sqrt(int_GPS_vel_lower_sqsum*smooth_rainfall_data_sqsum)]
        vu_vl[d]=[xcorrlag(int_GPS_vel_lower,int_GPS_vel_upper,1/24.0)[0],xcorrlag(int_GPS_vel_lower,int_GPS_vel_upper,1/24.0)[1]/math.sqrt(int_GPS_vel_upper_sqsum*int_GPS_vel_lower_sqsum)]
        cr_ra[d]=[xcorrlag(smooth_hourly_event_sum_all_hours,smooth_rainfall_data,1/24.0)[0],xcorrlag(smooth_hourly_event_sum_all_hours,smooth_rainfall_data,1/24.0)[1]/math.sqrt(smooth_hourly_event_sum_all_hours_sqsum*smooth_rainfall_data_sqsum)]
        vu_cr[d]=[xcorrlag(int_GPS_vel_upper,smooth_hourly_event_sum_all_hours,1/24.0)[0],xcorrlag(int_GPS_vel_upper,smooth_hourly_event_sum_all_hours,1/24.0)[1]/math.sqrt(smooth_hourly_event_sum_all_hours_sqsum*int_GPS_vel_upper_sqsum)]
        vl_cr[d]=[xcorrlag(int_GPS_vel_lower,smooth_hourly_event_sum_all_hours,1/24.0)[0],xcorrlag(int_GPS_vel_lower,smooth_hourly_event_sum_all_hours,1/24.0)[1]/math.sqrt(smooth_hourly_event_sum_all_hours_sqsum*int_GPS_vel_lower_sqsum)]
        ra_st[d]=[xcorrlag(smooth_rainfall_data,int_strain,1/24.0)[0],xcorrlag(smooth_rainfall_data,int_strain,1/24.0)[1]/math.sqrt(smooth_rainfall_data_sqsum*int_strain_sqsum)]
        st_cr[d]=[xcorrlag(smooth_hourly_event_sum_all_hours,int_strain,1/24.0)[0],xcorrlag(smooth_hourly_event_sum_all_hours,int_strain,1/24.0)[1]/math.sqrt(smooth_hourly_event_sum_all_hours_sqsum*int_strain_sqsum)]
    
#    print vu_vl
#    print vu_ra
#    print vl_ra
#    print cr_ra
#    print vu_cr
#    print vl_cr
#    print ra_st
#    print st_cr
#    
#    import sys
#    sys.exit()
    
#    smoothed_rainfall_data=smooth_data(rainfall_data,3)
#    smoothed_hourly_event_sum_all_hours=smooth_data(hourly_event_sum_all_hours,3)    
    
    ### plot hourly events against hourly rainfall ( individual event category ) ###
    
#    import matplotlib
#    matplotlib.use('Qt4Agg')  
    
    #axis & grid definition
    
    # unsmoothed:
        
    rain_major_ticks=np.arange(0,25,5)
#    rain_minor_ticks=np.arange(0,roundtocell(int(max(rainfall_data)),5),1)
    velocity_major_ticks=np.arange(0,3.5,1)
#    velocity_minor_ticks=np.arange(-1,4,0.2)  
#    uplift_major_ticks=np.arange(-2,2,1)
    events_major_ticks=np.arange(0,60,15) # crevssing
#    events_minor_ticks=np.arange(0,roundtocell(max(hourly_event_sum_all_hours),5)+5,1)   
    v2events_major_ticks=np.arange(0,20,5) # margin crevassing
#    v2events_minor_ticks=np.arange(0,roundtocell(max(v2_hourly_event_sum_all_hours),2)+2,1)
#    v3events_major_ticks=np.arange(0,120,30) #locals
#    v3events_major_ticks=np.arange(0,roundtocell(max(v3_hourly_event_sum_all_hours),2)+2,1)
#    lake_level_major_ticks=np.arange(0, 3, 0.5)
#    strain_major_ticks=np.arange(-0.0002, 0.0010, 0.0002)
    strain_major_ticks=np.arange(-0.0005, 0.0015, 0.0005)
    
    # smoothed:
    
#    rain_major_ticks=np.arange(0,roundtocell(int(max(smooth_rainfall_data)),2),2)
#    rain_minor_ticks=np.arange(0,roundtocell(int(max(smooth_rainfall_data)),5),1)
#    velocity_major_ticks=np.arange(-0.5,2,0.5)
#    velocity_minor_ticks=np.arange(-0.5,2,0.1)  
#    uplift_major_ticks=np.arange(-0.5,1.5,0.5)
#    uplift_minor_ticks=np.arange(-0.5,1.5,0.1),
#    events_major_ticks=np.arange(0,roundtocell(max(smooth_seismic_data),2)+2,2)
#    events_minor_ticks=np.arange(0,roundtocell(max(smooth_seismic_data),5)+2,1)   
#    v2events_major_ticks=np.arange(0,roundtocell(max(smooth_v2seismic_data),2)+2,1)
#    v2events_minor_ticks=np.arange(0,roundtocell(max(smooth_v2seismic_data),2)+2,1)
#    v3events_major_ticks=np.arange(0,roundtocell(max(smooth_v3seismic_data),2)+2,5)
#    v3events_major_ticks=np.arange(0,roundtocell(max(smooth_v3seismic_data),2)+2,1)
    
    start_date=122
    decimal_days_major_ticks=np.arange(start_date,end_date+2,3)
    decimal_days_minor_ticks=np.arange(start_date,end_date+2,1)

#    decimal_days_major_ticks=np.arange(start_date,end_date,1)
#    decimal_days_minor_ticks=np.arange(start_date,end_date,1/24.0)
    
    #plotting...
#    bar_width=0.2 # use this for all-time plot
    bar_width=0.039 # use this for zooming in
#    bar_width=0.05
#    line_width=1.5
    
#    plt.figure(figsize=(7.48,7.48)) # 190 mm width (full page)
    plt.figure(figsize=(3.74,7.48))

#    plt.figure(figsize=(20,15.98))
    plt.subplot(611)
    ax1=plt.gca()    
    
    ax1.bar(dec_julian_days_shifted,rainfall_data,lw=0, width=bar_width, color = colors[0])
#    ax1.plot(dec_julian_days_shifted,smooth_rainfall_data,color='b',lw=2*line_width)
    ax1.set_ylabel('Rain rate\n(mm hr$^{-1}$)',fontsize=10, labelpad=20)
    ax1.set_xticks(decimal_days_major_ticks)
    ax1.set_xticks(decimal_days_minor_ticks,minor='True')
    ax1.set_xlim(start_date,end_date)
    ax1.set_yticks(rain_major_ticks)
#    ax1.set_yticks(rain_minor_ticks,minor='True')
    ax1.set_ylim(0,roundtocell(max(rainfall_data),5))
#    ax1.set_ylim(0,roundtocell(max(smooth_rainfall_data),5))
    ax1.tick_params(axis='both',which='major',labelsize=10)
    ax1.tick_params(axis='x',which='minor',labelsize=None)
    ax1.get_yaxis().set_label_coords(-0.1,0.5)
    ax1.spines['top'].set_visible(True)
    ax1.spines['right'].set_visible(True)
    ax1.yaxis.set_ticks_position('both')
    ax1.xaxis.set_ticks_position('both')
    ax1.tick_params(which = 'both', length = 3)
    ax1.tick_params(which = 'both', width = 1.1)
#    ax1.grid(which='major',alpha=0.5)
#    ax1.grid(which='minor',alpha=0.25)
    
#    plt.title('Crevassing',fontsize=16,fontweight='bold',y=1.03)
    
    plt.subplot(613, sharex=ax1)
    ax2=plt.gca()

    ax2.plot(GPS_vel_time[1],GPS_vel[1], lw=std_linewidth/2, color = colors[2], linestyle='--')    
    ax2.plot(GPS_vel_time[0],GPS_vel[0], lw=std_linewidth/2, color = colors[1])
#    ax2.plot(dec_julian_days_shifted,smooth_hvel_lower,color='green', lw=line_width)
#    ax2.plot(dec_julian_days_shifted,smooth_hvel_upper,color='purple', lw=line_width)
    ax2.set_xticks(decimal_days_major_ticks)
    ax2.set_xticks(decimal_days_minor_ticks,minor='True')
    ax2.set_ylim(0,3.5)
#    ax2.set_ylim(-0.5,2)
    ax2.set_yticks(velocity_major_ticks)
#    ax2.set_yticks(velocity_minor_ticks, minor='True')
    ax2.set_ylabel('Horizontal\nvelocity\n(m d$^{-1}$)', fontsize=10, labelpad=20)
    ax2.set_xlim(start_date,end_date)
    ax2.tick_params(axis='both',which='major',labelsize=10)
    ax2.tick_params(axis='x',which='minor',labelsize=None)
    ax2.get_yaxis().set_label_coords(-0.1,0.5)
    ax2.spines['top'].set_visible(True)
    ax2.spines['right'].set_visible(True)
    ax2.yaxis.set_ticks_position('both')
    ax2.xaxis.set_ticks_position('both')
    ax2.tick_params(which = 'both', length = 3)
    ax2.tick_params(which = 'both', width = 1.1) 
#    ax2.grid(which='major',alpha=0.5)
    
#    plt.subplot(614, sharex=ax1)
#    ax0=plt.gca()
#    
#    
#    ax0.plot(GPS_vel_dir_time_out[0],GPS_vel_dir[0],color='green', lw=line_width)
#    ax0.plot(GPS_vel_dir_time_out[1],GPS_vel_dir[1],color='purple', lw=line_width)
#    ax0.set_xticks(decimal_days_major_ticks)
#    ax0.set_xticks(decimal_days_minor_ticks,minor='True')
#    ax0.set_ylabel('Velocity\nBearing \n(Degrees)', fontsize=10, labelpad=15)
#    ax0.set_xlim(start_date,end_date)
#    ax0.set_yticks(np.arange(150, 250, 5))
##    ax0.set_ylim()
#    ax0.tick_params(axis='both',which='major',labelsize=10)
#    ax0.tick_params(axis='x',which='minor',labelsize=None)
#    ax0.get_yaxis().set_label_coords(-0.1,0.5)  
##    ax0.grid(which='major',alpha=0.5)
    
    plt.subplot(612, sharex=ax1)
    ax3=plt.gca()
    
    ax3.plot(GPS_time[1],z2pos, lw=std_linewidth/2, color=colors[2], linestyle='--')    
    ax3.plot(GPS_time[0],z1pos, lw=std_linewidth/2, color=colors[1])
#    ax3.plot(dec_julian_days_shifted,smooth_uplift_lower,color='green', lw=line_width)
#    ax3.plot(dec_julian_days_shifted,smooth_uplift_upper,color='purple', lw=line_width)
#    ax3.set_xticks(decimal_days_major_ticks)
#    ax3.set_xticks(decimal_days_minor_ticks,minor='True')
#    ax3.set_ylim(-2,2)
    ax3.set_yticks(np.arange(-0.5,3,0.25))
#    ax3.set_yticks(uplift_minor_ticks, minor='True')
    ax3.set_ylabel('Detrended\nelevation\n(m)', fontsize=10, labelpad=20)
#    ax3.set_xlim(start_date,end_date)
    ax3.tick_params(axis='both',which='major',labelsize=10)
    ax3.tick_params(axis='x',which='minor',labelsize=None)
    ax3.get_yaxis().set_label_coords(-0.1,0.5)
    ax3.spines['top'].set_visible(True)
    ax3.spines['right'].set_visible(True)
    ax3.yaxis.set_ticks_position('both')
    ax3.xaxis.set_ticks_position('both')
    ax3.tick_params(which = 'both', length = 3)
    ax3.tick_params(which = 'both', width = 1.1)
#    ax3.grid(which='major',alpha=0.5)
    
#    ax1.set_xticklabels([]) 
    
    plt.subplot(614, sharex=ax1)
    ax4=plt.gca()
    
    ax4.plot(GPS_vel_time[0],strains, lw=std_linewidth, color=colors[3])
#    ax4.plot(dec_julian_days_shifted,smooth_hstrain,color='orange', lw=line_width)
    ax4.set_xticks(decimal_days_major_ticks)
    ax4.set_xticks(decimal_days_minor_ticks,minor='True')
    ax4.set_ylabel('Strain rate\n(s$^{-1}$)', fontsize=10, labelpad=20)
    ax4.set_xlim(start_date,end_date)
    ax4.set_yticks(strain_major_ticks)
    ax4.set_ylim(-0.0005, 0.0013)
    ax4.tick_params(axis='both',which='major',labelsize=10)
    ax4.tick_params(axis='x',which='minor',labelsize=None)
    ax4.get_yaxis().set_label_coords(-0.1,0.5)  
    ax4.spines['top'].set_visible(True)
    ax4.spines['right'].set_visible(True)
    ax4.yaxis.set_ticks_position('both')
    ax4.xaxis.set_ticks_position('both')
    ax4.tick_params(which = 'both', length = 3)
    ax4.tick_params(which = 'both', width = 1.1)
#    ax4.grid(which='major',alpha=0.5)

#    ax1.set_xticklabels([])      
    
    plt.subplot(615, sharex=ax1)
    ax6=plt.gca()
    
    ax6.bar(dec_julian_days_shifted,hourly_event_sum_all_hours,lw=0, width=bar_width, color=colors[4])
#    ax6.plot(dec_julian_days_shifted,smooth_seismic_data,color='r',lw=2*line_width)
#    ax6.set_ylabel('"Spaced"\nCrevassing\n(count hr$^{-1}$)',fontsize=10, labelpad=25)
    ax6.set_ylabel('Type A\nseismicity\n(count hr$^{-1}$)',fontsize=10, labelpad=20)
    ax6.set_xticks(decimal_days_major_ticks)
    ax6.set_xticks(decimal_days_minor_ticks,minor='True')
    ax6.set_xlim(start_date,end_date)
    ax6.set_yticks(events_major_ticks)
#    ax6.set_yticks(events_minor_ticks,minor='True')
    ax6.set_ylim(0,45)
#    ax6.set_ylim(0,roundtocell(max(smooth_seismic_data),2)+2)
    ax6.tick_params(axis='both',which='major',labelsize=10)
    ax6.tick_params(axis='x',which='minor',labelsize=None)
    ax6.get_yaxis().set_label_coords(-0.1,0.5)
    ax6.spines['top'].set_visible(True)
    ax6.spines['right'].set_visible(True)
    ax6.yaxis.set_ticks_position('both')
    ax6.xaxis.set_ticks_position('both')
    ax6.tick_params(which = 'both', length = 3)
    ax6.tick_params(which = 'both', width = 1.1)
#    ax6.grid(which='major',alpha=0.5)
#    ax6.grid(which='minor',alpha=0.25)
  
    ax1.set_xticklabels([])
 
    plt.subplot(616)#, sharex=ax1)
    ax7=plt.gca()
   
    ax7.bar(dec_julian_days_shifted,v2_hourly_event_sum_all_hours,lw=0, width=bar_width, color=colors[4])    
#    ax7.plot(dec_julian_days_shifted,smooth_v2seismic_data,color='r',lw=2*line_width)
#    ax7.set_ylabel('"Sustained"\nCrevassing\n(count hr$^{-1}$)',fontsize=10, labelpad=25)
    ax7.set_ylabel('Type B\nseismicity\n(count hr$^{-1}$)',fontsize=10, labelpad=15)
    ax7.set_xticks(decimal_days_major_ticks)
    ax7.set_xticks(decimal_days_minor_ticks,minor='True')
    ax7.set_xlim(start_date,end_date)
    ax7.set_yticks(v2events_major_ticks)
    ax7.tick_params(axis='both',which='major',labelsize=10)
#    ax7.set_yticks(v2events_major_ticks)
#    ax7.set_yticks(v2events_minor_ticks,minor='True')
    ax7.set_ylim(0,15)
#    ax7.set_ylim(0,roundtocell(max(smooth_v2seismic_data),2)+2)

    ax7.set_xlabel('Day of year (2016)', fontsize=10, labelpad=5)
    ax7.set_xticks(decimal_days_major_ticks)
    ax7.set_xticks(decimal_days_minor_ticks,minor='True')
    ax7.tick_params(axis='x',which='major',labelsize=10)
    ax7.spines['top'].set_visible(True)
    ax7.spines['right'].set_visible(True)
    ax7.yaxis.set_ticks_position('both')
    ax7.xaxis.set_ticks_position('both')
    ax7.tick_params(which = 'both', length = 3)
    ax7.tick_params(which = 'both', width = 1.1)

#    ax7.tick_params(axis='both',which='y',labelsize=10)
##    ax7.tick_params(axis='both',which='x',labelsize=10)
#    ax7.tick_params(axis='x',which='minor',labelsize=None)
##    ax7.get_xaxis().set_label_coords(0.5,-0.11)
    ax7.get_yaxis().set_label_coords(-0.1,0.5)
#    ax7.grid(which='major',alpha=0.5)
#    ax7.grid(which='minor',alpha=0.25)
    
    ax6.set_xlim(122,139)
    ax7.set_xlim(122,139)
    
#    ax1.set_xticklabels([])
#    
#    plt.subplot(616)#, sharex=ax1)
##    
#    ax8=plt.gca()
#    ax8.plot(dec_julian_days, lake_levels, color='orange', lw=line_width)
#    ax8.set_ylabel('Lake Level\n(m)', fontsize=10, labelpad=25)
#    ax8.set_ylim(728.5,732)
#    ax8.set_yticks(np.arange(728.5, 732, 1))
#    ax8.set_xlabel('Day of Year (2016)', fontsize=16, labelpad=15)
#    ax8.set_xticks(decimal_days_major_ticks)
#    ax8.set_xticks(decimal_days_minor_ticks,minor='True')
##    ax8.set_xlim(start_date,end_date)
#    ax8.tick_params(axis='x',which='major',labelsize=10)
#    ax8.tick_params(axis='x',which='major',labelsize=10)
#    ax8.tick_params(axis='x',which='minor',labelsize=None)
##    ax8.grid(which='major',alpha=0.5)
#    ax8.get_yaxis().set_label_coords(-0.07,0.5)   
#    ax8.tick_params(axis='y',which='major',labelsize=10)
#
#    ax9=ax8.twinx()
#    ax9.plot(dec_julian_days,lake_roc,color='cyan',lw=line_width)
#    ax9.set_ylabel('Lake Level\nRate of Change\n(m s$^{-1}$)',fontsize=10, labelpad=50, rotation=270)
#    ax9.set_yticks(lake_level_major_ticks) 
#    ax9.set_ylim(-1, 2.5)
#    ax9.tick_params(axis='y',which='major',labelsize=10)
#    ax9.set_xlim(start_date,end_date)
#    
    plt.show()
##    plt.savefig('/Users/home/taylorsa/FIGURE_U_WANT.pdf', format='pdf', dpi=300)
#    
#    # plot the catalogue events in ascending order of number per hour
#    
#    event_statistics=sorted(hourly_event_sum_all_hours)
#    
#    print hourly_event_sum_all_hours
#    
#    # remove 0 values
#    
#    event_statistics=filter(lambda a: a != 0.0, event_statistics)
##    
##    # find the distribution of event counts
##    
#    event_distribution=[0 for x in range(int(max(event_statistics))+1)]
#    for count in event_statistics:
#        event_distribution[int(count)]+=1
#    
#    # find the continuous count of events:
#    
#    # contiuous count:
#    c=1
#    count=1
#    cumulative_events=[]
#    cumulative_events.append(event_statistics[0]) #pre-set the first entry so subsequent entries can add it to their value
#    for entry in event_statistics: #look at every entry in the list of event counts per hour
#        cumulative_events.append(event_statistics[c]+cumulative_events[c-1]) #add the event count for a given hour to all of those preceeding
#        c+=1 #increment the position in the list
#        if c==len(event_statistics): #when the index reaches the end of the list terminate the loop
#            break
#        count+=event_statistics[c] #increment the number of events by the number in that hour
#    
#    # find the first and second derivatives of this (dx = 1 hour)
#    cumulative_events_1D,cumulative_events_2D=central_difference(cumulative_events,1)
#        
#    fig=plt.figure()
#    ax = fig.add_subplot(111)    
#    ax.plot(range(len(cumulative_events)),cumulative_events,lw=2)
#    ax.plot(range(len(cumulative_events)),cumulative_events_1D,lw=1)
#    ax.plot(range(len(cumulative_events)),cumulative_events_2D,lw=0.5)
#
#    print np.mean(event_statistics)
#    print np.std(event_statistics)
#    print scipy.stats.mode(event_statistics)
#    print np.median(event_statistics)
#
#    plt.close('all')
#    fig=plt.figure(figsize=(8.27,8.27))
#    ax=fig.add_subplot(111)
##    ax.plot(range(len(event_statistics)),event_statistics)
#    ax.plot(range(len(event_distribution)),event_distribution, lw=3, color='red')
#    ax.set_xlabel('Event Count (count hr$^{-1}$)', fontsize=14, labelpad=15)
#    ax.set_ylabel('Event Count Occurrence', fontsize=14, labelpad=15)
#    ax.plot([1 for x in range(351)],range(351), lw=1.5, color='black')
#    ax.plot([5 for x in range(351)],range(351), lw=1.5, color='black')
#    ax.set_xticks(np.arange(0,36,5))
#    ax.set_xticks(np.arange(0,36,1),minor='True')
#    ax.set_xlim(0,36)
#    ax.set_yticks(np.arange(0,350,50))
#    ax.set_yticks(np.arange(0,350,10),minor='True')
#    ax.set_ylim(0,350)
#    ax.tick_params(axis='both',which='major',labelsize=12)
#    ax.tick_params(axis='x',which='minor',labelsize=None)
##    ax.grid(which='major',alpha=0.5)
##    ax.grid(which='minor',alpha=0.25)
#    plt.show()  
#    plt.savefig('/Users/home/taylorsa/FIGURE_U_WANT.pdf', format='pdf', dpi=300)
#    
#
#run_event_statistics()    
    
    
    
def run_event_statistics_windows(outdir, event_directory,window_start_date,window_end_date):
    # creates a plot of the 3 rain temporal views (derivatives) vs the seismicity of each window
    
    import csv 
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    import scipy
    
    year=2016
    hour_shift=0.5 #shift all times by a certain amount. Implemented to plot data at half past the hour to represent time collected at
    delay_time=7
    background_rate=0 # number of events expected in an any given hour in the event catalogue
    
    rainfall_times='/Volumes/arc_01/taylorsa/rainfall_data/rainfall_hours'
    rainfall_data='/Volumes/arc_01/taylorsa/rainfall_data/rainfall_hours_data' 
    temperature_times='/Volumes/arc_01/taylorsa/rainfall_data/temp_hours'
    temperature_data='/Volumes/arc_01/taylorsa/rainfall_data/temp_data'
    
    # calculate temperature for each hour

#    T_dec_julian_days, temp_data, T_julian_days, day_temp=chosen_rainfall(window_start_date,window_end_date,year,temperature_times,temperature_data)
    
    # GPS processing
                                                       
#    GPS1_path='/Volumes/arc_01/taylorsa/velocity_data/tas_position_arc5_1.txt'
#    GPS2_path='/Volumes/arc_01/taylorsa/velocity_data/tas_position_arc4_1.txt'
#    GPS_in=[GPS1_path, GPS2_path]
#    GPS_out=GPS_parse(GPS_in)
#    GPS_pos, GPS_vel=GPS_VELOCITY(GPS_out, window_start_date, window_end_date)
#    strains=strain_rate(GPS_pos, GPS_vel)

    # calculate the rainfall in each hour and each day

    dec_julian_days,rainfall_data, julian_days, day_rainfall=chosen_rainfall(window_start_date,window_end_date,year,rainfall_times,rainfall_data)
    
    # calculate the number of events occuring in each hour
    
    hourly_event_sum,hourly_event_hours,files,all_times=event_counter(window_start_date,window_end_date,event_directory,delay_time)
    
    # remove event catalogue background rate from each hour
    
    for x in range(len(hourly_event_sum)): 
        hourly_event_sum[x]-=background_rate
        if hourly_event_sum[x]<0: hourly_event_sum[x]=0
   
    #reorgansie times for events to work with rainfall timing
    
    hourly_event_sum_all_hours=time_reconciler(hourly_event_sum,hourly_event_hours,dec_julian_days,window_start_date,window_end_date)
    
    # smooth rainfall data to 24 hour means
    
    smooth_rainfall_data=smooth_data(rainfall_data,12)  
    
    # sum hourly rainfall into the cumulative amount
    
    cumulative_rainfall=accumulate_data(rainfall_data)
    
    # find the "acceleration" of rainfall
    
    acceleration_rainfall=central_difference(rainfall_data,1/24.0)[0]
    
    # hour_shift correction to events and rainfall data (plots at the midpoint of the hour)
    c=0
    dec_julian_days_shifted=[0]*len(dec_julian_days)
    if hour_shift!=0: #if there is a shift to be applied
        for day in dec_julian_days: #this will shift the timing for events and rain (both bin to the start of the hour)
            dec_julian_days_shifted[c]=float("{0:.8f}".format((day-hour_shift/24)))##[:8]
            c+=1
    #observation is made at 0900 for 0800-0900, so plot at 0830
    
    ### plot hourly events against hourly rainfall ( individual event category ) ###
    
    import matplotlib
    matplotlib.use('Qt4Agg')  
    
    #axis & grid definition

#    temp_major_ticks=np.arange(0,int(max(temp_data)+max(temp_data)/10.0),int(math.ceil(max(temp_data)/5)))
#    temp_minor_ticks=np.arange(0,int(max(temp_data)+max(temp_data)/10.0),int(math.ceil(max(temp_data)/10)))
    
    rain_major_ticks=np.arange(0,roundtocell(int(max(rainfall_data)),5)+5,5)
    rain_minor_ticks=np.arange(0,roundtocell(int(max(rainfall_data)),5)+5,1)
    cum_rain_major_ticks=np.arange(0,roundtocell(int(max(cumulative_rainfall)),20)+20,20)
    cum_rain_minor_ticks=np.arange(0,roundtocell(int(max(cumulative_rainfall)),20)+20,10)
#    acc_rain_major_ticks=np.arange(int(min(acceleration_rainfall)),int(max(acceleration_rainfall)+max(acceleration_rainfall)/10.0),int(math.ceil(max(acceleration_rainfall)+abs(min(acceleration_rainfall)))/5))
#    acc_rain_minor_ticks=np.arange(int(min(acceleration_rainfall)),int(max(acceleration_rainfall)+max(acceleration_rainfall)/10.0),int(math.ceil(max(acceleration_rainfall)+abs(min(acceleration_rainfall)))/10))
    acc_rain_major_ticks=np.arange(-roundtocell(int(max([max(acceleration_rainfall),abs(min(acceleration_rainfall))])),60)-60,roundtocell(int(max([max(acceleration_rainfall),abs(min(acceleration_rainfall))])),60)+60,60)
    acc_rain_minor_ticks=np.arange(-roundtocell(int(max([max(acceleration_rainfall),abs(min(acceleration_rainfall))])),60)-60,roundtocell(int(max([max(acceleration_rainfall),abs(min(acceleration_rainfall))])),60)+60,30)

    events_major_ticks=np.arange(0,roundtocell(max(hourly_event_sum_all_hours),5)+5,5)
    events_minor_ticks=np.arange(0,roundtocell(max(hourly_event_sum_all_hours),5)+5,1)
    
    decimal_days_major_ticks=np.arange(window_start_date,window_end_date,1)
    decimal_days_minor_ticks=np.arange(window_start_date,window_end_date,1/24.0)
    
    #plotting...
#    bar_width=0.5 # use this for all-time plot
    bar_width=0.02 # use this for zooming in
    line_width=1

    plt.figure(figsize=(12,9))
    plt.subplot(511) # 5 rows ( 1 temp, 3 rain, 1 seis )
    ax1=plt.gca()
    
    plt.title('Weather Relationships with Glacier Seismicity Peaks',fontsize=16,fontweight='bold',y=1.03)
    
#    ax1.plot(dec_julian_days_shifted,temp_data,color='purple',lw=2)
    ax1.set_ylabel('T (deg. C)',fontsize=14)
    ax1.set_xticks(decimal_days_major_ticks)
    ax1.set_xticks(decimal_days_minor_ticks,minor='True')
    ax1.set_xlim(window_start_date,window_end_date)
#    ax1.set_yticks(temp_major_ticks)
#    ax1.set_yticks(temp_minor_ticks,minor='True')
#    ax1.set_ylim(0,max(temp_data))
#    ax1.tick_params(axis='y',which='major',labelsize=14)
#    ax1.tick_params(axis='x',which='major',labelsize=None)
#    ax1.tick_params(axis='x',which='minor',labelsize=None)
#    ax1.get_xaxis().set_label_coords(0.5,-0.11)
#    ax1.get_yaxis().set_label_coords(-0.04,0.5)
#    ax1.grid(which='major',alpha=0.5)
#    ax1.grid(which='minor',alpha=0.25)
#    ax1.set_xticklabels([])
    
    plt.subplot(513, sharex=ax1)
    ax3=plt.gca()    

    ax3.bar(dec_julian_days_shifted,smooth_rainfall_data,color='b',width=bar_width,lw=line_width)
    ax3.set_ylabel(r'$(\frac{dR}{dt})$',fontsize=14)
    ax3.set_xticks(decimal_days_major_ticks)
    ax3.set_xticks(decimal_days_minor_ticks,minor='True')
    ax3.set_xlim(window_start_date,window_end_date)
    ax3.set_yticks(rain_major_ticks)
    ax3.set_yticks(rain_minor_ticks,minor='True')
    ax3.set_ylim(0,roundtocell(max(rainfall_data),5)+5)
    ax3.tick_params(axis='y',which='major',labelsize=10)
    ax3.tick_params(axis='x',which='major',labelsize=None)
    ax3.tick_params(axis='x',which='minor',labelsize=None)
    ax3.get_yaxis().set_label_coords(-0.04,0.5)
    ax3.grid(which='major',alpha=0.5)
    ax3.grid(which='minor',alpha=0.25)
    ax3.set_xticklabels([])
    
    plt.subplot(514, sharex=ax1)
    ax4=plt.gca()    

    ax4.bar(dec_julian_days_shifted,acceleration_rainfall,color='b',width=bar_width,lw=line_width)
    ax4.set_ylabel(r'$(\frac{d^2R}{dt^2})$',fontsize=14)
    ax4.set_xticks(decimal_days_major_ticks)
    ax4.set_xticks(decimal_days_minor_ticks,minor='True')
    ax4.set_xlim(window_start_date,window_end_date)
    ax4.set_yticks(acc_rain_major_ticks)
    ax4.set_yticks(acc_rain_minor_ticks,minor='True')
    ax4.set_ylim(-roundtocell(max([abs(min(acceleration_rainfall)),max(acceleration_rainfall)]),60)-60,roundtocell(max([abs(min(acceleration_rainfall)),max(acceleration_rainfall)]),60)+60)
    ax4.tick_params(axis='y',which='major',labelsize=10)
    ax4.tick_params(axis='x',which='major',labelsize=None)
    ax4.tick_params(axis='x',which='minor',labelsize=None)
    ax4.get_yaxis().set_label_coords(-0.04,0.5)
    ax4.grid(which='major',alpha=0.5)
    ax4.grid(which='minor',alpha=0.25)
    ax4.set_xticklabels([])
    
    plt.subplot(512, sharex=ax1)
    ax2=plt.gca()    

    ax2.bar(dec_julian_days_shifted,cumulative_rainfall,color='b',width=bar_width,lw=line_width)
    ax2.set_ylabel(' R (mm) ',fontsize=14)
    ax2.set_xticks(decimal_days_major_ticks)
    ax2.set_xticks(decimal_days_minor_ticks,minor='True')
    ax2.set_xlim(window_start_date,window_end_date)
    ax2.set_yticks(cum_rain_major_ticks)
    ax2.set_yticks(cum_rain_minor_ticks,minor='True')
    ax2.set_ylim(0,roundtocell(max(cumulative_rainfall),20)+20)
    ax2.tick_params(axis='y',which='major',labelsize=10)
    ax2.tick_params(axis='x',which='major',labelsize=None)
    ax2.tick_params(axis='x',which='minor',labelsize=None)
    ax2.get_yaxis().set_label_coords(-0.04,0.5)
    ax2.grid(which='major',alpha=0.5)
    ax2.grid(which='minor',alpha=0.25)
    ax2.set_xticklabels([])
    
    plt.subplot(515, sharex=ax1)
    ax5=plt.gca()
    
    ax5.bar(dec_julian_days_shifted,hourly_event_sum_all_hours,color='r',width=bar_width,lw=line_width)
    ax5.set_xlabel('Julian Day (2016)', fontsize=14)
#    ax5.set_ylabel(event_directory.split("/")[-2]+' events',fontsize=14)
    ax5.set_ylabel('Event Count',fontsize=14)
    ax5.set_xticks(decimal_days_major_ticks)
    ax5.set_xticks(decimal_days_minor_ticks,minor='True')
    ax5.set_xlim(window_start_date,window_end_date)
    ax5.set_yticks(events_major_ticks)
    ax5.set_yticks(events_minor_ticks,minor='True')
    ax5.set_ylim(0,roundtocell(max(hourly_event_sum_all_hours),5)+5)
    ax5.tick_params(axis='both',which='major',labelsize=10)
    ax5.tick_params(axis='x',which='minor',labelsize=None)
    ax5.get_xaxis().set_label_coords(0.5,-0.11)
    ax5.get_yaxis().set_label_coords(-0.04,0.5)
    ax5.grid(which='major',alpha=0.5)
    ax5.grid(which='minor',alpha=0.25)
    ax5.set_xticklabels(range(window_start_date,window_end_date))
    
    ax5.tick_params(axis='x',which='major',labelsize=10)
    ax5.tick_params(axis='x',which='minor',labelsize=None)
    
    plt.savefig(outdir+'/rain_relationship_smoothed.png',format='png',dpi=150)
    plt.close('all')
    
    # run the drainage model for this window
    drainage_model(outdir, window_start_date, window_end_date,dec_julian_days_shifted, rainfall_data, dec_julian_days_shifted, hourly_event_sum_all_hours)
    plt.close('all')

## run the above function for many windows
#
#import csv
#import math
#import os
#import matplotlib.pyplot as plt
#    
#window_directory='/Volumes/arc_02/taylorsa/windows/'
#windows=os.listdir(window_directory)
#for window in windows:
#    print window
#    if window in ['1','6','8']: 
#        continue
#    event_directory=window_directory+window+'/'
#    with open(event_directory+'time_windows.csv','r') as csvfile:
#        row_reader=csv.reader(csvfile, delimiter=',')
#        for row in row_reader: # takes the LAST row (i.e. most recent windowing for the given peak ( assuming threshold doesn't change between window limits ))
#            window_start_date=int(math.floor(float(row[0])))
#            window_end_date=int(math.ceil(float(row[-1])))
#    outdir=window_directory+window
#    run_event_statistics_windows(outdir,event_directory,window_start_date,window_end_date)
#    plt.close('all')

def run_event_statistics_catalogue(event_directory,subplot):
    
    #this runs the stand-alone event statistics script. Otherwise the script
    #acts as a number of functions to be called for more advanced processes
    #(this version is for subplotting many event catalogues)
               
    import csv 
    import matplotlib.pyplot as plt
    import numpy as np
    import math
    
    start_date=112
    end_date=201
    year=2016
    hour_shift=0.5 #shift all times by a certain amount. Implemented to plot data at half past the hour to represent time collected at
    delay_time=0
    rainfall_times='/Volumes/arc_01/taylorsa/rainfall_data/rainfall_hours'
    rainfall_data='/Volumes/arc_01/taylorsa/rainfall_data/rainfall_hours_data' 
    
#    diurnal_event_statistics(start_date,end_date,event_directory,delay_time)
    
    dec_julian_days,rainfall_data, julian_days, day_rainfall=chosen_rainfall(start_date,end_date,year,rainfall_times,rainfall_data)
    
    hourly_event_sum,hourly_event_hours,files,all_times=event_counter(start_date,end_date,event_directory,delay_time)
    
    #rate_of_event=backwards_difference(cumulative_events,1)
   
    #reorgansie times for events to work with rainfall timing
    hourly_event_sum_all_hours=time_reconciler(hourly_event_sum,hourly_event_hours,dec_julian_days,start_date,end_date)
    
    # hour_shift correction to events and rainfall data (plots rainfall at the midpoint of the hour)
    c=0
    dec_julian_days_shifted=[0]*len(dec_julian_days)
    if hour_shift!=0: #if there is a shift to be applied
        for day in dec_julian_days: #this will shift the timing for events and rain (both bin to the start of the hour)
            dec_julian_days_shifted[c]=float("{0:.8f}".format((day-hour_shift/24)))##[:8]
            c+=1
    #observation is made at 0900 for 0800-0900, so plot at 0830

    ### plot hourly events against hourly rainfall ( many event categories ) ###
    
    import matplotlib
    matplotlib.use('Qt4Agg')  
    
    #axis & grid definition
    
    rain_major_ticks=np.arange(0,20,5)
#    rain_minor_ticks=np.arange(0,20,2.5)
#
#    events_major_ticks=np.arange(0,250,50)
#    events_minor_ticks=np.arange(0,250,25)
#    
#    decimal_days_major_ticks=np.arange(start_date,end_date,10)
#    decimal_days_minor_ticks=np.arange(start_date,end_date,1)

#    events_major_ticks=np.arange(0,roundtocell(max(hourly_event_sum_all_hours),20)+20,(roundtocell(max(hourly_event_sum_all_hours),20)+20)/5.0)

    if event_directory.split("/")[-2]=='noise':
        event_major_ticks=np.arange(0,120,30)
    elif event_directory.split("/")[-2]=='earthquake':
        event_major_ticks=np.arange(0,4,1)
    elif event_directory.split("/")[-2]=='crevassing':
        event_major_ticks=np.arange(0,40,10)
    elif event_directory.split("/")[-2]=='calving':        
        event_major_ticks=np.arange(0,8,2)
    elif event_directory.split("/")[-2]=='icequake':
        event_major_ticks=np.arange(0,8,2)
    elif event_directory.split("/")[-2]=='helicopter':
        event_major_ticks=np.arange(0,4,1)
    elif event_directory.split("/")[-2]=='lightning':
        event_major_ticks=np.arange(0,20,5)
    else:
        event_major_ticks=np.arange(0,200,50)
    
#    decimal_days_major_ticks=np.arange(start_date,end_date,1)
#    decimal_days_minor_ticks=np.arange(start_date,end_date,1/24.0)

    decimal_days_major_ticks=np.arange(start_date,end_date,10)
    decimal_days_minor_ticks=np.arange(start_date,end_date,1)
    
    bar_width=0.2
#    bar_width=0.02 # use this for zooming in
    line_width=1

    plt.subplot(811)
    ax1=plt.gca()

    # rainfall (daily)

    ax1.bar(dec_julian_days_shifted,rainfall_data,color='b',width=bar_width,lw=line_width)
    ax1.set_ylabel('Rain Rate\n(mm hr$^{-1}$)',fontsize=8)
    ax1.set_xticks(decimal_days_major_ticks)
    ax1.set_xticks(decimal_days_minor_ticks,minor='True')
    ax1.set_xlim(start_date,end_date)
    ax1.set_yticks(rain_major_ticks)
#    ax1.set_yticks(rain_minor_ticks,minor='True')
    ax1.set_ylim(0,roundtocell(int(max(rainfall_data)),5)+5,5)
    ax1.tick_params(axis='y',which='major',labelsize=8)
    ax1.tick_params(axis='x',which='major',labelsize=None)
    ax1.tick_params(axis='x',which='minor',labelsize=None)
    ax1.get_yaxis().set_label_coords(-0.07,0.5)
#    ax1.grid(which='major',alpha=0.5)
#    ax1.grid(which='minor',alpha=0.25)
    ax1.set_xticklabels([])

    # plot seismic data at subplot position
    
    plt.subplot(subplot)
    ax2=plt.gca()
    
    ax2.bar(dec_julian_days_shifted,hourly_event_sum_all_hours,color='r',width=bar_width,lw=line_width)
    ytitle=event_directory.split("/")[-2]
    Ytitle=ytitle[:1].upper() + ytitle[1:]
    ax2.set_ylabel(Ytitle+'\n(count hr$^{-1}$)',fontsize=8)
#    ax2.set_ylabel('Event Count',fontsize=14)
    ax2.set_xticks(decimal_days_major_ticks)
    ax2.set_xticks(decimal_days_minor_ticks,minor='True')
    ax2.set_xlim(start_date,end_date)
    ax2.set_yticks(event_major_ticks)
#    ax2.set_yticks(events_minor_ticks,minor='True')
    ax2.tick_params(axis='y',which='major',labelsize=8)
    ax2.tick_params(axis='x',which='major',labelsize=None)   
    ax2.tick_params(axis='x',which='minor',labelsize=None)
    ax2.get_yaxis().set_label_coords(-0.07,0.5)
#    ax2.grid(which='major',alpha=0.5)
#    ax2.grid(which='minor',alpha=0.25)
    
    if str(subplot)[-1]=='8':
        ax2.tick_params(axis='both',which='major',labelsize=8) 
        ax2.tick_params(axis='x',which='minor',labelsize=None)
    else:
        ax2.set_xticklabels([])


## run this section for many event catalogues
##
## NOTE: need to manually change the number of rows in the subplot
#
#import matplotlib.pyplot as plt
#
##catalogue_directory='/Volumes/arc_01/taylorsa/catalogues/master_catalogue/'
##catalogue_subdirectories=['local','crevassing']
#catalogue_directory='/Volumes/arc_01/taylorsa/triggers/sta_locals/'
#catalogue_subdirectories=['TSNM1','TSNM2','TSNM3']
#
##catalogue_directory='/Volumes/arc_01/taylorsa/catalogues/plotting_master_cat/'
##catalogue_subdirectories=['noise','earthquake','crevassing','calving','icequake','helicopter','lightning']
#
#si=1   
#
#for subdir in catalogue_subdirectories:
#    print subdir
#    si+=1
#    event_directory=catalogue_directory+subdir+'/'
#    subplot=int(str(len(catalogue_subdirectories)+1)+str(1)+str(si))
#    run_event_statistics_catalogue(event_directory, subplot)
#    
#  
#plt.gcf()
#plt.xlabel('Day of Year (2016)',fontsize=16,labelpad=15)
#    
#plt.show()

# run this section for a single event catalogue
