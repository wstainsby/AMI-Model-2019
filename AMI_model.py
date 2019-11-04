import pandas as pd
import datetime as dt
import numpy as np
import pandas.io.sql as sqlio
import psycopg2 as p2
import time as time
import pytz
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.dates as mdates
import CSU_PVSTEM as pv
import csv as csv
from matplotlib.ticker import PercentFormatter


TIMEZONE = 'America/Denver'
TZ = pytz.timezone(TIMEZONE)
UTC = pytz.utc

LAT_DEGREES = 40.59541667  # Set for Fort Collins (all systems)
LONG_DEGREES = -105.0844   # Same
ALT = 1525                 # elevation
UNAME = #need to fill in with access login
PASSWORD = #need to fill in with access login
dayNow = time.strftime("%m-%d-%Y")

# set up DB connection
connStr = #connection to the database
print('Connected to DB')
conn = p2.connect(connStr)

batteryList = ['14871','14873','22653','25248','30545','30607','34323','44250','46852','46852','50009','74382',
              '78468','84603','88923','89343','92874'] #exclude prem's with batteries installed


resRateCodes = ['A100', 'A115', 'A125', 'A130', 'A140', 'B140', 'E100', 'E110',
                'E115', 'E120', 'E125', 'E130', 'E135', 'E140']

DGassets = pd.read_csv('inputs/...) #a list of DG assets
DGassets['dateInstalled'] = pd.to_datetime(DGassets.installDate, infer_datetime_format=True)
DGassets = DGassets.sort_values('premiseID', ascending=True)
earlyAssets = DGassets[DGassets['dateInstalled'].dt.year < 2016]
DGassets = DGassets[DGassets['dateInstalled'].dt.year > 2015]
DGassets = DGassets[DGassets['dateInstalled'] < '04-01-2019']
DGassets = DGassets.sort_values('dateInstalled', ascending=True)
DGassets = DGassets.sort_values('premiseID', ascending=True)
DGassets = DGassets[(DGassets.premiseRateCode.isin(resRateCodes)) | (DGassets.rateClass.isin(resRateCodes))]
DGassets = DGassets[~DGassets.premiseID.isin(batteryList)]
DGassets = DGassets[~DGassets.tilt.isnull()]
DGassets = DGassets[~DGassets.azimuth.isnull()]


'''range of timestamps that occur when sun is at least 1 degree above horizon according to the PVL model'''
dayTimestamps = pd.read_csv('inputs/dayTimestamps.csv', index_col='idx')
dayTimestamps.index = pd.to_datetime(dayTimestamps.index, utc=True, infer_datetime_format=True)
dayTimestamps.index = pd.DatetimeIndex(dayTimestamps.index).tz_convert(tz=TZ)

Daylight15minWeather = pd.read_csv('inputs/15minDayWeather.csv', index_col='idx')
Daylight15minWeather.index = pd.to_datetime(Daylight15minWeather.index, utc=True, infer_datetime_format=True)
Daylight15minWeather.index = pd.DatetimeIndex(Daylight15minWeather.index).tz_convert(tz=TZ)
Daylight15minWeather = Daylight15minWeather[Daylight15minWeather.index.year > 2014]

fullWeather = pd.read_csv('inputs/fullWeather.csv', index_col='idx')
fullWeather.index = pd.to_datetime(fullWeather.index, utc=True, infer_datetime_format=True)
fullWeather.index = pd.DatetimeIndex(fullWeather.index).tz_convert(tz=TZ)
fullWeather = fullWeather[fullWeather.index.year > 2014]
fullWeather = fullWeather[fullWeather.index < '04-01-2019']

'''Snow loss factors with tilt and percetage generation loss associated as "loss factor".'''
tilt_degree = [5, 10, 15, 20, 25, 30, 46]
loss_factor = [0.925,	0.83,	0.735,	0.64,	0.545,	0.45,	0.26]
a = zip(loss_factor, tilt_degree)
lossTable = pd.DataFrame(a, columns=['lossFactor', 'tilt'])

snow = pd.read_csv('inputs/FCsnow.csv')
idx = pd.DatetimeIndex(snow.date).tz_localize(tz=TZ)
snow = snow.set_index(idx)
snow = snow.fillna(value=0)
snowHourly = snow.resample('H').ffill()

A = time.time()


'''uncomment this section if you want to recreate the JSON comparable timestamp file from a new set of parameters.
holiday list for each year [ NY day, memorial day, july 4, labor day, thanksgiving (thurs+fri), xmas eve, xmas, 26th]'''

# holidayDict = {}
# holidayDict[2015] = [1, 145, 185, 250, 330, 331, 358, 359, 360]
# holidayDict[2016] = [1, 151, 186, 248, 329, 330, 359, 360, 361]
# holidayDict[2017] = [1, 149, 185, 247, 327, 328, 358, 359, 360]
# holidayDict[2018] = [1, 148, 185, 246, 326, 327, 358, 359, 360]
# holidayDict[2019] = [1]
# hourDelta = 4
#
# compTimestamps = pd.DataFrame(index=Daylight15minWeather.index, columns=['listTS'])
# compTimestamps = compTimestamps[compTimestamps.index.year > 2014]
#
# for timestamp, timeData in Daylight15minWeather.iterrows():
#
#     #make day of year TS's as a separate set from all three components
#     shoulderDaysBegin = (365-15 + timeData['DOY']) <= (Daylight15minWeather['DOY'])
#     shoulderDaysEnd = (365-timeData['DOY']) <= (15-Daylight15minWeather['DOY'])
#     midYearDays = ((timeData['DOY'] >= (Daylight15minWeather['DOY'] - 15)) & (timeData['DOY'] <= (Daylight15minWeather['DOY'] + 15)))
#     dayOfYearDays =((timeData['DOY'] >= (Daylight15minWeather['DOY'] - 15)) & (timeData['DOY'] <= (Daylight15minWeather['DOY'] + 15)))
#     allShoulderDays = np.logical_or(shoulderDaysBegin, shoulderDaysEnd)
#     combinedTS = np.logical_or(allShoulderDays, midYearDays)
#
#     dayOfWeekTS = (timeData['dayOfWeek'] == Daylight15minWeather['dayOfWeek'])
#
#     matchingTS = pd.DataFrame(Daylight15minWeather.loc[dayOfWeekTS & combinedTS])
#
#     tsHours = (timeData['hour'] >= matchingTS['hour'] - hourDelta) & \
#               (timeData['hour'] <= (hourDelta + matchingTS['hour']))
#     tsHours_before = (24 + timeData['hour'] - hourDelta) <= matchingTS['hour']
#     tsHours_after = (24 - timeData['hour'] <= (hourDelta - matchingTS['hour']))
#     tsHours_mid = (timeData['hour'] >= matchingTS['hour'] - hourDelta) & (
#             timeData['hour'] <= matchingTS['hour'] + hourDelta)
#     matchingTS = matchingTS.loc[tsHours_before | tsHours_after | tsHours_mid]
#
#     tempHours = (timeData['temp'] >= (matchingTS['temp'] - (0.3*matchingTS['temp'].std()))) &\
#                 (timeData['temp'] <= (matchingTS['temp'] + (0.3*matchingTS['temp'].std())))
#
#     irrHours = (timeData['irradiance'] >= (matchingTS['irradiance'] - (0.4*matchingTS['irradiance'].std()))) & \
#                (timeData['irradiance'] <= (matchingTS['irradiance'] + (0.4*matchingTS['irradiance'].std())))
#
#     matchingTS = matchingTS.loc[tempHours & irrHours]
#
#     matchingTS = matchingTS.loc[matchingTS.index != timestamp]
#
#     # remove major american holidays
#     for year, holiday in holidayDict.items():
#         try:
#             holidaysInYear = matchingTS['year'] == np.logical_and(matchingTS['year'], matchingTS['DOY'].isin(holiday))
#             matchingTS = matchingTS[~holidaysInYear]
#             matchingTS = matchingTS.index
#         except:
#             continue
#
#     listHours = list(matchingTS)
#     #print(len(listHours))
#     compTimestamps.loc[timestamp, 'listTS'] = listHours
#     #make day of year days as a separate set from all three components
#
# compTimestamps.to_json("../Datasets/Comparable/compTS_final.json")

compTS_3 = pd.read_json('inputs/compTS.json')
compTS_3 = pd.DataFrame(compTS_3)
compTS_3.index = pd.to_datetime(compTS_3.index, utc=True, infer_datetime_format=True)
compTS_3.index = compTS_3.index.tz_localize(tz='UTC')
compTS_3.index = compTS_3.index.tz_convert(tz=TZ)
compTS_3 = compTS_3.loc[compTS_3.index.dropna()]

B = time.time()
avgLength= []
print('{:.2f} seconds to load dictionaries'.format(B-A))
for key, data in compTS_3.iterrows():
    a=len(data.listTS)
    avgLength.append(a)

print(np.mean(avgLength))


'''AMI Query, splitting data into pre- and post-PV dataframes'''

SQLFMT1 = '''
SELECT *, "timestamp" AT TIME ZONE 'UTC' as "timeUTC"
FROM "AMIReadings"
WHERE
    "premiseID" = '{}'

ORDER BY "timestamp" ASC
'''


#dfAZ = pv.azimuth(Daylight15minWeather, LAT_DEGREES, LONG_DEGREES)
dfAZ = pv.azimuth(fullWeather, LAT_DEGREES, LONG_DEGREES)
df1 = pv.haurwitz(dfAZ)  # clearsky GHI
df1['GHI'] = fullWeather['irradiance']
#df1['GHI'] = Daylight15minWeather['irradiance']
df2 = pv.calcErbs(df1)  # Erbs derives DNI and DHI from inputed GHI
energyInput = pd.merge(pd.merge(dfAZ, df2[['GHI', 'DNI', 'DHI']], left_index=True, right_index=True),
                       fullWeather[['windSpeed', 'temp', 'irradiance']], left_index=True, right_index=True)

outList = []
outputSums = []

day_offset = 20 #buffer around PV install date

print("Done with PVL")

with open('../output/outputSums15minDay.csv', "w") as oFile:
    csvW = csv.writer(oFile)
    csvW.writerow(['premID', 'AMI', 'PVL', 'systemCapacity', 'tilt', 'azimuth', 'installDate', 'rateCode',
                   'estErrorMedian', 'estErrorMean', 'WHerrorMedian', 'WHerrorMean','tsPercentage', 'numOutputOverages', 'percentIrregular'])

    '''iterate through each premise in DG assets'''
    for prem, premData in DGassets.iterrows():

        plt.close()

        irregularList = []
        systemCapacity = premData.capacitykW

        premID = str(premData.premiseID)

        tilt = premData.tilt
        azimuth = premData.azimuth
        rateClass = premData.rateClass
        rateCode = premData.premiseRateCode
        ACrating = premData.maxAC



        #calculating my two buffer dates for the queries below.
        installDate = premData.dateInstalled
        installDate = pd.to_datetime(installDate, utc=None)
        installDate = installDate.tz_localize(tz=TZ)

        query1 = SQLFMT1.format(premID)
        allAMI = sqlio.read_sql_query(query1, conn)

        if len(allAMI) < 200: #discards premises with less 1 week operation time
            continue

        idx = allAMI.timestamp.apply(lambda d: d.replace(tzinfo=None))
        allAMI.index = idx
        allAMI = allAMI.drop(columns=['premiseID', 'isActual', 'timeUTC'])
        allAMI.index= pd.DatetimeIndex(allAMI.index).tz_localize(tz=TZ, ambiguous='NaT')
        #allAMI = allAMI[allAMI.index.isin(dayTimestamps.index)]
        prePV = allAMI[allAMI.index < installDate - dt.timedelta(days=day_offset)]
        postPV = allAMI[allAMI.index >= installDate + dt.timedelta(days=day_offset)]
        originalAMIindex = len(postPV.index)

        del allAMI

        postPV['numCompHours'] = 0
        postPV['medDelBeforePV'] = 0
        postPV['meanDelBeforePV'] = 0
        postPV['RecBeforePV'] = 0
        postPV['compDelAfterPV'] = np.nan

        print("PV query completed")

        if len(prePV) == 0 or len(postPV) == 0:
            continue

        timestampsForPrem = compTS_3[compTS_3.index.isin(postPV.index)]

        for ts, epochDF in timestampsForPrem.iterrows():

            df = pd.to_datetime(epochDF.listTS, utc=True, unit='ms')
            df = df.tz_convert(tz=TZ)

            tsBeforePV = prePV[prePV.index.isin(df)]

            # filling out the AMI data table with the AMI records when more than 5 comparable TS
            postPV.loc[ts, 'numCompHours'] = len(tsBeforePV)
            postPV.loc[ts, 'medDelBeforePV'] = tsBeforePV.deliveredKWH.median()
            postPV.loc[ts, 'meanDelBeforePV'] = tsBeforePV.deliveredKWH.mean()
            postPV.loc[ts, 'RecBeforePV'] = tsBeforePV.receivedKWH.max()

            compHoursAfterPV = postPV[postPV.index.isin(df)] #using this for when delivered postPV power > deliv prePV power
            postPV.loc[ts, 'compDelAfterPV'] = compHoursAfterPV.deliveredKWH.median()

            location = postPV.index.get_loc(ts)
            if len(tsBeforePV) < 3:
                postPV.loc[ts, 'numCompHours'] = postPV.iloc[location-1].numCompHours
                postPV.loc[ts, 'medDelBeforePV'] = postPV.iloc[location-1].medDelBeforePV
                postPV.loc[ts, 'meanDelBeforePV'] = postPV.iloc[location-1].meanDelBeforePV
                postPV.loc[ts, 'RecBeforePV'] = postPV.iloc[location-1].RecBeforePV

        print('{} TS under threshold'.format(len(postPV[postPV['numCompHours'] < 5])))


        '''estimating generation cycling through the equations, starting at gen=0 '''

        postPV['estOutput'] = 0    # default condition

        postPV.loc[postPV['receivedKWH'] > 0, 'estOutput'] = postPV['receivedKWH']  #fourth best estimation

        postPV.loc[postPV['medDelBeforePV'] > postPV['compDelAfterPV'], 'estOutput'] = \
               postPV['medDelBeforePV'] - postPV['compDelAfterPV'] + postPV['receivedKWH']    #third best estimation condition

        postPV.loc[postPV['meanDelBeforePV'] > postPV['deliveredKWH'], 'estOutput'] = \
            postPV['meanDelBeforePV'] - postPV['deliveredKWH'] + postPV['receivedKWH']        #second best estimation condition

        postPV.loc[postPV['medDelBeforePV'] > postPV['deliveredKWH'], 'estOutput'] = \
            postPV['medDelBeforePV'] - postPV['deliveredKWH']+ postPV['receivedKWH']    #best estimation condition

        postPV = postPV[postPV.receivedKWH < systemCapacity]
        postPV = postPV[postPV.estOutput < systemCapacity]

        print("finished comparing timestamps")

        #just want to run PVL on the interval that I'm actually using
        PVLinput = energyInput[energyInput.index.isin(postPV.index)]
        PVLinput = PVLinput[PVLinput.index > (installDate + dt.timedelta(days=day_offset))]
        df3 = pv.calcEnergyByCapacity(PVLinput, tilt, azimuth, systemCapacity)

        if pd.notnull(ACrating):
            df3.loc[(df3['Eest'] > ACrating), 'Eest'] = ACrating

        estimated = pd.DataFrame()
        estimated['PVL'] = df3.Eest / 4
        estimated['AMI'] = postPV.estOutput
        estimated.loc[estimated['PVL'] < 0, 'PVL'] = 0

        estimatedHourly = estimated.resample('H', closed='right', label='right').sum()
        estimatedHourly = estimatedHourly.join(snowHourly['snow'], how='left')
        estimatedHourly = estimatedHourly.join(snowHourly['snwd'], how='left')
        estimatedHourly['tilt'] = tilt
        estimatedHourly['timestamp'] = estimatedHourly.index

        '''snow fall correction based on tilt'''

        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[6].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[6].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[5].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[5].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[4].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[4].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[3].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[3].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[2].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[2].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[1].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[1].lossFactor) * estimatedHourly['PVL']
        estimatedHourly.loc[((estimatedHourly['tilt'] < lossTable.loc[0].tilt) & (
                    (estimatedHourly['snow'] >= 1) | (estimatedHourly['snwd'] >= 3))),
                            'PVL'] = (1 - lossTable.loc[0].lossFactor) * estimatedHourly['PVL']

        estimatedHourly['hourlyError'] = (2 * (estimatedHourly.AMI - estimatedHourly.PVL) /
                                          (estimatedHourly.PVL + estimatedHourly.AMI))

        maxGen = estimatedHourly[['AMI', 'PVL']].values.max()

        estimatedHourly['weightedError'] = ((estimatedHourly[['AMI', 'PVL']].max(axis=1)) / maxGen) *\
                                           estimatedHourly.hourlyError

        f1, ax = plt.subplots(figsize=(12.8, 5), constrained_layout=True)
        ax.yaxis.grid("True", linestyle='--', zorder=0, linewidth=0.4, color='darkgrey')
        a2 = ax.plot(estimatedHourly['PVL'], label='PV-STEM Model', lw=1.2,zorder=5)
        a3 = ax.plot(estimatedHourly['AMI'], label='AMI Model',lw=1.2, zorder=5)
        plt.ylabel('Hourly Generation (kWh)', fontsize=20)
        plt.xticks(rotation=45, fontsize=15)
        myFmt = mdates.DateFormatter('%m-%d-%y')
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_major_locator(mdates.DayLocator(interval=1, tz=TZ))
        ax.tick_params(axis='both', which='major', labelsize=15)
        plt.legend(loc='best', fontsize=20)
        plt.title('{} kW system installed {}'.format(systemCapacity, installDate.strftime('%m/%Y')),
                  fontweight='bold', fontsize=20)
        # plotfile = '../Figures/prem{}_Hourly.png'.format(premID)
        # plt.savefig(plotfile, bbox_inches="tight")



        estimatedHourlyShift = estimatedHourly.shift(-5)
        estimatedDay = estimatedHourlyShift.resample('D').sum()
        estimatedDay = estimatedDay[['PVL', 'AMI']]
        estimatedDay['premiseID'] = premID
        estimatedDay['dailyError'] = (2*(estimatedDay.AMI - estimatedDay.PVL) /
                                        (estimatedDay.PVL + estimatedDay.AMI))

        PVLsum = estimatedDay.PVL.sum()
        AMIsum = estimatedDay.AMI.sum()

        irregular = postPV[(postPV.deliveredKWH > postPV.medDelBeforePV) & (postPV.medDelBeforePV >0)]
        numIrregular = len(irregular) / len(postPV)
        toWrite = (premID, AMIsum, PVLsum, systemCapacity, tilt,
                azimuth, installDate, rateClass, estimatedDay.dailyError.median(), estimatedDay.dailyError.mean(),
                   estimatedHourly.weightedError.median(), estimatedHourly.weightedError.mean(), (len(postPV)/originalAMIindex),
                   len(postPV[postPV.estOutput > systemCapacity]), numIrregular)

        csvW.writerow(toWrite)
        oFile.flush()


        f2, ax = plt.subplots(figsize=(12.8, 5), constrained_layout=True)
        #plt.title('{} kW system installed {}'.format(systemCapacity, installDate.strftime('%m/%Y')),
                 # fontweight='bold', fontsize=20)
        ax.yaxis.grid("True", linestyle='--', zorder=0, linewidth=0.4, color='darkgrey')
        ax.plot(estimatedDay['PVL'], label='PVL Model', lw=1.2, zorder=5)
        ax.plot(estimatedDay['AMI'], label='AMI Model', lw=1.2, zorder=5)
        ax.tick_params(axis='both', which='major', labelsize=15)
        # if ACrating < systemCapacity:
        #     ax.plot(estimatedDay['PVL_maxAC'], label='maxAC', lw=1)

        plt.ylabel('Daily Generation (kWh)', fontsize=20)
        plt.xticks(rotation=45, fontsize=14)
        myFmt = mdates.DateFormatter('%b-%y')
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_major_locator(mdates.MonthLocator(interval=2))
        plt.legend(loc='best', fontsize=15)
        plotfile = '../Figures/prem{}_Day.png'.format(premID)
        plt.savefig(plotfile, bbox_inches="tight")

        estimatedDay.to_csv('../output/{}.csv'.format(premID), index_label='tmp')
        estimatedHourly.to_csv('../output/{}.csv'.format(premID), index_label='idx')


        del estimatedDay
        del estimated
        del estimatedHourly

        plt.close()
        plt.cla()
        plt.clf()

        print('premise {}'.format(premID))
