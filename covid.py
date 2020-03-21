import json
import math
import matplotlib.pyplot as plt

state_list = ['CA','WA','NY']
states = {'CA':0,'WA':1,'NY':2}

date_array = [[],[],[]]
num_pos_array = [[],[],[]]

state_dicts = [[],[],[]]

file ='C:/Users/mears/Desktop/daily_all.json'
with open(file) as f:
    d = json.load(f)

q = enumerate(state_list)
for i,state2 in enumerate(state_list):
    dt = state_dicts[i]
    for record in d:
        state = record["state"]
        if state2 == state:
            state_dicts[i].append(record)
                
for i,d in enumerate(state_dicts):
    print(i)
    for r in d:
        date_int = r["date"]
        year = math.floor(date_int/10000)
        month_day = date_int % 10000
        month = math.floor(month_day/100)
        day = month_day%100
        print(year,month,day)

        date_array[i].append(r["date"])
        num_pos_array[i].append(r["positive"])

    


print()