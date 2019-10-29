import pandas as pd
import time as time
import pytz
import numpy as np
import CSUPVLib as pv

TIMEZONE = 'America/Denver'
TZ = pytz.timezone(TIMEZONE)
UTC = pytz.utc

#source this to wherever you store the file.
CSUinput = pd.read_csv('../FCweatherData/CSUweatherData.csv', header=0)

#adjusting index for correct timestamp that is "timezone-aware" in Mountain Standard Time (MST)
CSUinput['timestamp'] = CSUinput['date'] + ' ' + CSUinput['time']
CSUinput['timestamp'] = pd.to_datetime(CSUinput.timestamp, format="%m/%d/%Y %H:%M")
CSUinput['utcTS'] = CSUinput['timestamp'] + pd.Timedelta(hours=7) #convert TZ-naive timestamp to UTC time,  but still naive of any timezone consideration
CSUinput['utcAware'] = CSUinput.utcTS.dt.tz_localize('UTC') #gives timezone to the UTC timestamp
localTZ = pd.DatetimeIndex(CSUinput.utcAware)
CSUinput['local_time'] = localTZ.tz_convert(tz=TZ)
CSUinput.index = CSUinput['local_time'] #setting local Mountain Standard Time datetime index as the index
CSUinput = CSUinput.drop(columns=['date', 'time', 'utcTS', 'utcAware', 'timestamp', 'local_time']) #get rid of extra columns






downsampled = CSUinput.resample('5T').fillna('bfill')




resampled=downsampled.resample('15T', label='right', closed='right').mean()


resampled = resampled.drop(columns=['direction'])

resampled = resampled[resampled.index.year >2014]

resampled['DOY'] = resampled.index.dayofyear
resampled['year'] = resampled.index.year
resampled['hour'] = resampled.index.hour
resampled['weekDay'] = ((resampled.index.dayofweek) // 5 == 1).astype(int)
resampled['dayOfWeek'] = resampled.index.dayofweek
resampled['timestamp'] = resampled.index
resampled.index = resampled.index.tz_convert(tz=None)
resampled.to_csv("../FCweatherData/15minWeather.csv",  index_label='idx')




#here I just add another blank year's worth of 15-timestamps for future use.
# It's useful because the PVL model runs off of timestamp and calculates solar position from it.
start = finalResample.index[-1] + pd.Timedelta(minutes=15)
end = finalResample.index[-1] + pd.Timedelta(weeks=52)
addedDTS = pd.date_range(start=start, end=end, freq='15T', tz=TZ)
fullDateRange = pd.DataFrame(index=addedDTS)
fullDateRange.temp = np.nan
fullDateRange.windSpeed = np.nan
fullDateRange.irradiance = np.nan
finalResample = finalResample.append(fullDateRange)
finalResample.index = finalResample.index.tz_convert(tz=TZ)
finalResample['DOY'] = finalResample.index.dayofyear
finalResample['year'] = finalResample.index.year
finalResample['hour'] = finalResample.index.hour
finalResample['dayOfWeek'] = (((finalResample.index.dayofweek) // 5 == 1).astype(int))
finalResample['timestamp'] = finalResample.index
#to wherever you want to save the file.
finalResample.to_csv('../Datasets/weather/fullWeather.csv', index_label='idx')




#if you want to see how to gather only "daytime" or "nighttime" timestamps, this is it below
#basically I run pv.azimuth (first function of the PVL model) can capture day time stamps as when the zenith of the sun
#is less than 89 degrees, meaning it is above horizon. timestamps right around sunrise/sunset will be problematic because the
#cosine of the zenith angle is used in some calculations and it creates abnormal values as zenith approaches 90 deg.
#because cos(90deg) approaches zero.



#create DF for solar positioning from CSU PVL
df = pv.azimuth(finalResample, LAT_DEGREES, LONG_DEGREES)
df.index = pd.DatetimeIndex(df.index).tz_convert(tz=TZ)

dayTimestamps = df[df.zenith <= 89]  #sun above horizon
dayTimestamps.to_csv("../Datasets/weather/dayTimestamps.csv", index_label='idx')

nightTimestamps = df[df.zenith > 92] #sun well below horizon
nightTimestamps.to_csv("../Datasets/weather/nightTimestamps.csv", index_label='idx')

dayWeather15min = finalResample[finalResample.index.isin(dayTimestamps.index)]
dayWeather15min['DOY'] = dayWeather15min.index.dayofyear
dayWeather15min['year'] = dayWeather15min.index.year
dayWeather15min['hour'] = dayWeather15min.index.hour
dayWeather15min['dayOfWeek'] = (((dayWeather15min.index.dayofweek) // 5 == 1).astype(int))
dayWeather15min.to_csv('../Datasets/weather/15minDayWeather.csv', index_label='idx')

nightWeather15min = finalResample[finalResample.index.isin(nightTimestamps.index)]
nightWeather15min['DOY'] = nightWeather15min.index.dayofyear
nightWeather15min['year'] = nightWeather15min.index.year
nightWeather15min['hour'] = nightWeather15min.index.hour
nightWeather15min['dayOfWeek'] = (((nightWeather15min.index.dayofweek) // 5 == 1).astype(int))
nightWeather15min.to_csv('../Datasets/weather/15minNightWeather.csv', index_label='idx')

