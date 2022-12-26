import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
# import pandas as pd

# time = 'Time [year]'
# raw_data = 'C [mol/yr]'
# data = pd.read_csv(r'C:\Users\ruben\OneDrive\Документы\Python Scripts\test\ISA92a-scenario.csv', usecols=[time, raw_data])
# xval_raw = np.linspace(0, 2788, 2788) # change the number of timesteps
# xval = sorted(list(xval_raw) + list(data[time]))
# interp = np.interp(xval, data[time], data[raw_data])

class IPA:
    def __init__(self) -> None:
        self.lister = []

    def parter(self, interval: np.ndarray, interp: np.array, min_size: int) -> None:
        if len(interval) < min_size:
            raise ValueError('something went wrong, last section cant be an interval')
        elif len(interval) == min_size:
            # if the inputted interval is of the same length as 
            # a minimum interval then we take the whole thing as
            # an interval
            self.lister.append( (interval[0], interval[-1]) )
        else:
            # starts from 0
            pivot_index = 0
            # interpolation value at 0
            f_pivot = interp[pivot_index]
            # index where the minumum interval ends
            min_interval_end = pivot_index + min_size-1
            # interval that we take into mean, goes from the second
            # to last index in the minimum interval
            mean_interval = interp[pivot_index + 1: min_interval_end]
            # the mean itself
            mean = stats.tmean(mean_interval)
            # slope of the line from the first to mean
            slope = (mean-f_pivot)/(min_interval_end - pivot_index)
            # the function of the form y = m(x - x0) + f(x0)
            L = lambda x: slope*(x - pivot_index) + f_pivot
            # the standard devation of the minimum inteval
            stdev = stats.tstd(mean_interval)
            # walk one-by-one starting from the pivot point
            walker_stepper = 1
            while min_interval_end + walker_stepper < len(interval):
                # starting from the first index after the last index of
                # the original interval, it calculates the function value for
                # L and adds and subtracts the standard devation of the interval
                smaller = L(min_interval_end + walker_stepper)-stdev
                larger = L(min_interval_end + walker_stepper)+stdev
                # interpolation value after at the specific index
                func_val = interp[min_interval_end + walker_stepper]
                # if the next point after the interval is in the range,
                # nothing happens as it's now in the interval and the
                # walker_stepper continues increasing
                temporary_pivot = min_interval_end+walker_stepper
                if smaller <= func_val <= larger:
                    if len(interval[temporary_pivot:]) < min_size:
                        self.lister.append( (interval[pivot_index], interval[-1]) )
                        break
                else:
                    # once we hit a point not in the range, we set it to be
                    # a temporary (potentially new) pivot and consider the
                    # points after it
                    # if we cannot fit a new interval after the temporary pivot
                    # we just ignore it and make the whole section into a one
                    # large interval
                    if len(interval[temporary_pivot:]) < min_size:
                        self.lister.append( (interval[pivot_index], interval[-1]) )
                        break
                    else:
                        # if we can fit an interval after pivot, we check if the first-removed
                        # average is in the original range of function +- standard devation
                        mean = stats.tmean(interp[temporary_pivot+1: temporary_pivot+min_size])
                        new_smaller = L(temporary_pivot + min_size)-stdev
                        new_larger = L(temporary_pivot + min_size)+stdev
                        if new_smaller <= mean <= new_larger:
                            # we make the walker_stepper 0 as it as the new
                            # min_interval_end has been updated
                            walker_stepper = 0
                            min_interval_end = temporary_pivot + min_size
                        else:
                            # if it's not in the range then we create the interval and plug it in
                            # the function again just with the original interval removed
                            self.lister.append( (interval[pivot_index], interval[temporary_pivot-1]) )
                            self.parter(interval[temporary_pivot:], interp[temporary_pivot:], min_size)
                            # breaks as this interval is complete
                            break
                walker_stepper += 1


# P = IPA()
# P.parter(xval, interp, 200) # define minimum interval size
# plt.plot(xval, interp)
# print(P.lister)
# for e in list(P.lister):
#     plt.axvline(x=e[0], color='red')

# plt.axvline(x=list(list(P.lister))[-1][1], color='red')
# plt.xlabel(time)
# plt.ylabel(raw_data)
# plt.title('title-of-the-data')
# plt.show()
