import json
import math
import matplotlib.pyplot as plt
import datetime
import numpy as np
import requests
data_colors = ['black','blue','forestgreen','red','forestgreen']
fit_colors  =['grey','aqua','lightgreen','pink','lightgreen']
var_to_plot = 'positive'
url = 'https://covidtracking.com/api/us/daily'
r = requests.get(url)
d = r.json()

date_array=[]
num_pos_array=[]
for r in d:
    date_int = r["date"]
    year = math.floor(date_int/10000)
    month_day = date_int % 10000
    month = math.floor(month_day/100)
    day = month_day%100
    dt = datetime.date(year=year,month=month,day=day)
    date_array.append(dt)
    num_pos_array.append(r[var_to_plot])

ts_all= {'US' : {'date':date_array,
                var_to_plot:num_pos_array}}

state_list = ['CA','WA','NY']
states = {'CA':0,'WA':1,'NY':2}


url ='https://covidtracking.com/api/states/daily'
r = requests.get(url)
d = r.json()
#file ='C:/Users/Mears/Dropbox/docs/covid/daily_states.json'
#with open(file) as f:
 #   d = json.load(f)

for state2 in state_list:
    date_array=[]
    num_pos_array=[]
    for r in d:
        state = r["state"]
        #print(state2,state)
        if state2 == state:
            date_int = r["date"]
            year = math.floor(date_int/10000)
            month_day = date_int % 10000
            month = math.floor(month_day/100)
            day = month_day%100
            #print(year,month,day)
            dt = datetime.date(year=year,month=month,day=day)
            date_array.append(dt)
            num_pos_array.append(r[var_to_plot])
    ts_all[state2] = {"date" : date_array,
                      var_to_plot : num_pos_array}

                
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, title='COVID-19 Stats from covidtracking.com', 
                     xlabel='Date', 
                     ylabel=f'Number of {var_to_plot}s')

axs=[]
for i,key in enumerate(ts_all.keys()):
    ts = ts_all[key]
    a = ax.plot(ts['date'],ts[var_to_plot],label=key,color=data_colors[i])
    axs.append(a)

plt.legend()

fig2 = plt.figure(figsize=(12, 9))
xlim = [datetime.date(year=2020,month=3,day=1),datetime.date(year=2020,month=4,day=1)]
ax2 = fig2.add_subplot(111, title='COVID-19 Stats from covidtracking.com', xlabel='Date', 
       ylabel=f"Number of {var_to_plot}s",xlim=xlim)

axs=[]
for i,key in enumerate(ts_all.keys()):
    ts = ts_all[key]
    a = ax2.semilogy(ts['date'],ts[var_to_plot],label=key,color=data_colors[i])
    dt_arr=[]
    for date in ts['date']:
        dt = date-datetime.date(year=2020,month=3,day=1)
        print(i,dt)
        dt_arr.append(dt.days)
    
    dt_arr = np.array(dt_arr)
    ts_var_to_plot = np.array(ts[var_to_plot],dtype=np.float)
    ok = np.isfinite(ts_var_to_plot)
    fit=np.polyfit(dt_arr[ok], np.log(ts_var_to_plot[ok]), 1, w=np.sqrt(ts_var_to_plot[ok]))
    start_day = np.min(np.array(dt_arr))
    print(i,start_day)
    a =[]
    b =[]
    for j in np.arange(0,30):
        a.append(start_day+j)
        b.append(datetime.date(year=2020,month=3,day=1)+datetime.timedelta(days=float(start_day+j)))
    yfit=np.exp(fit[1])*np.exp(fit[0]*np.array(a))
    axs.append(a)
    afit=ax2.semilogy(b,yfit,color=fit_colors[i])

plt.legend()
plt.grid()
plt.show()

print()