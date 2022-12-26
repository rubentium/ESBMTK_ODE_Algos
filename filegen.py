import os
import numpy as np
import pandas as pd
from scipy import interpolate
from IPA import IPA
import polymc as pmc

# Variables (ie user defined input)
time = 'time'
raw_data = 'data'
min_step = int
data = pd.read_csv(r'file-path', usecols=[time, raw_data])
xval_raw = np.linspace(float, float, int) # interval bounds and number of timesteps
coef = 1 # default 1

# Interpolation and interval patitioning
xval = sorted(list(xval_raw) + list(data[time]))
interp = interpolate.interp1d(data[time], data[raw_data])
ipa = IPA()
ipa.parter(xval, interp(xval), min_step)
lister = ipa.lister

def func_val(timestamp, interp = interp):
      '''interpolation function'''
      return float(interp(timestamp))

def files(data, lister=lister, time=time, ydata=raw_data):
      '''Turn the output list of IPA into a folder of
      csv files for PolyMC to be able to work with and 
      then calls PolyMC to output the stepsize for each file'''

      part_size = len(lister)
      list_data = data.values.tolist()
      list_time = data[time].tolist()
      for tup in lister:
            first = tup[0]
            second = tup[1]
            if first not in list_time:
                  list_data.append( [first, func_val(first)] )
            if second not in list_time:
                  list_data.append( [second, func_val(second)] )

      list_data.sort(key=lambda x: x[0])
      list_of_list_data = []
      for tup in lister:
            first = tup[0]
            second = tup[1]
            subinterval_start = list_data.index( [first, func_val(first)] )
            subinterval_end = list_data.index( [second, func_val(second)] )
            list_of_list_data.append(list_data[subinterval_start: subinterval_end+1]) if subinterval_end != len(list_data)-1 else list_of_list_data.append(list_data[subinterval_start: -1])

      # create the folder and the files
      directory = './IPA/'
      try:
        if not os.path.exists(directory):
            os.makedirs(directory)
      except OSError:
            print ('Error: Creating directory. ' +  directory)

      for  i in range(part_size):
            sublist = list_of_list_data[i]
            subdata = pd.DataFrame(sublist)
            subdata.to_csv(directory + f'IPA-f{i}.csv', index=False, header=[time, ydata])

      # calculating stepsize using PolyMC
      polymc_list = []
      for filename in os.listdir(directory):
            file = os.path.join(directory, filename)
            polymc_data = pd.read_csv(file, usecols=[0, 1])
            return_meta = pmc.meta(polymc_data)
            return_hull = pmc.ConvexHull(return_meta[0]).hull
            return_funcgen = pmc.funcgen(return_hull)
            return__mc_area = pmc.mc_area(
            polymc_data[time][0],
            polymc_data[time][len(polymc_data[time]) - 1],
            return_meta[1],
            return_meta[2],
            return_funcgen,
            )
            polymc_list.append((filename, round(pmc.step_size(return__mc_area, coef), 2)))
            
      return polymc_list


print(files(data))