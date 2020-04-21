import obspy
import os
import matplotlib.pyplot as plt

# create a list of mseed files to load

mseedloc = '/home/samto/EVENTS_IT3/TYPE_D/3/'
mseed=os.listdir(mseedloc)
mseed.sort()
# plot all such files        
        
templates=[]
for entry in mseed:
    print(entry)
    st=obspy.read(mseedloc + entry)
    stog=obspy.Stream()
    for tr in st:
        if tr.stats.station not in ['TSNM1','TSNM2','TSNM3','TSNC4']:
            stog.append(tr)
    st=stog
    st.normalize
    fig=plt.plot(show=False)
    plt.close()
    st.plot(fig=fig, show=False)
    st.spectrogram()
#    matplotlib.use("Qt4Agg")     
#    mngr = plt.get_current_fig_manager()
#    mngr.window.setGeometry(3840,0,1080,1900)
    plt.show()
    
#    print('Template?')
#    ans=input('[Y/n]: ')
#    if ans=='Y':
#        print(entry+' will be saved for template generation\n')
#        templates.append(entry)
#    plt.close('all')

print(templates)


