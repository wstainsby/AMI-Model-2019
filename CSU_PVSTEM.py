import numpy as np
import pandas as pd

MODULE_EFFICIENCY=0.162
INVERTER_EFFICIENCY=0.96
SYSTEM_EFFICIENCY=(MODULE_EFFICIENCY * INVERTER_EFFICIENCY)


ESC = 1367              # Solar Constant (W per m2)
ALBEDO = 0.0            # Ground reflectance (typically 0.2 - see documentation)
DELT = 3.0              # DELT, a1 & b1 depend on panel type - see documentation
a1 = -3.56              # All are metrics for module temperature estimates
b1 = -0.075             # Glass/cell/polymer sheet with open rack assumed
GHI_LIM = 1.0           # Multiplier for GHI limit (GHI > EaH * GHI_LIM)
DNI_LIM = 1.0           # Multiplier for DNI limit (DNI > Ea * DNI_LIM)
TEMP_COEFF = -0.0043
MODULE_AREA = 1.63*0.95    # 1.63 meters squared * 95% cell area:panel area
MODULE_WATTAGE = 250       #250 watts per module


def sind(angleD):
    return np.sin(np.radians(angleD))

def cosd(angleD):
    return np.cos(np.radians(angleD))

def tand(angleD):
    return np.tan(np.radians(angleD))

def acosd(cos):
    return np.degrees(np.arccos(cos))



def calcTimevals(df):
    df['dayOfYear'] = df.index.dayofyear
    df['hour'] = df.index.hour
    df['minute'] = df.index.minute
    df['second'] = df.index.second
    df['tzOffsetHours'] = df.index.map(lambda x: x.utcoffset().total_seconds() / 3600.0)
    df['numDays'] = 365
    df.loc[df.index.is_leap_year, 'numDays'] = 366

#don't use this function

# def calcAzimuthNOAA(df, lat, lon):
#     calcTimevals(df)
#     df['lat'] = lat
#     df['lon'] = lon
#
#     df['gamma'] = (2 * np.pi) / df['numDays'] * (df['dayOfYear'] - 1 + (df['hour'] - 12) / 24)
#     df['eqtime'] = 229.18 * ((0.000075 + 0.001868 * np.cos(df['gamma']) - 0.032077 * np.sin(df['gamma'])
#                               - (0.014615 * np.cos(2 * df['gamma'])) - (0.040849 * np.sin(2 * df['gamma']))))
#     df['declR'] = 0.006918 - 0.399912 * np.cos(df['gamma']) + 0.070257 * np.sin(df['gamma']) \
#                   - 0.006758 * np.cos(2 * df['gamma']) + 0.000907 * np.sin(2 * df['gamma']) \
#                   - 0.002697 * np.cos(3 * df['gamma']) + 0.00148 * np.sin(3 * df['gamma'])
#     df['timeOffsetMinutes'] = df['eqtime'] + (4 * lon) - (60 * df['tzOffsetHours'])
#     df['tstMinutes'] = (df['hour'] * 60) + df['minute'] + (df['second'] / 60) + df['timeOffsetMinutes']
#     df['solarHourAngleD'] =  (df['tstMinutes'] / 4) - 180
#     df['cos_zen'] = sind(lat) * (np.sin(df['declR'])) + (cosd(lat) * np.cos(df['declR']) * cosd(df['solarHourAngleD']))
#     df['phiD'] = acosd(df['cos_zen'])
#
#
#     #WJS 5/5/2019
#     df['180_cosTheta'] = (np.sin(df['declR']) - (sind(lat) * cosd(df['phiD'])))    / (cosd(lat) * sind(df['phiD']))
#
#     df['thetaRadians'] = np.arccos(df['180_cosTheta'])
#     df['thetaDegrees'] = np.degrees(df['thetaRadians'])
#
#
#     df['answer'] = 180 - np.degrees(np.arccos(df['180_cosTheta']))
#
#     #original code
#     #df['cosTheta'] = (np.sin(df['declR']) - (sind(lat) * cosd(df['phiD']))) / (cosd(lat) * sind(df['phiD']))
#     #df['thetaD'] = 360 - acosd(df['cosTheta'])
#     #df['thetaD'] = acosd(df['cosTheta'])
#
#     df.loc[df['solarHourAngleD'] < 0, 'thetaD'] = acosd(df['cosTheta'])
#
#     df['declD'] = np.degrees(df['declR'])
#     df['timeOffsetHours'] = df['timeOffsetMinutes'] / 60.0
#     df['Tsolar'] = df['hour'] + df['timeOffsetHours']
#
#     df['zenith'] = df['phiD']
#     df['azimuth'] = df['thetaD']
#     df['sol_azD'] = df['azimuth']
#     return df

def calcEqt(df):
    df['Eqt'] = 0
    df.loc[(df['dayOfYear'] <= 106), 'Eqt'] = -14.2 * np.sin(np.pi * (df['dayOfYear'] + 7) / 111)
    df.loc[(df['dayOfYear'] > 106) & (df['dayOfYear'] <= 166), 'Eqt'] = 4.0 * np.sin(np.pi * (df['dayOfYear'] - 106) / 59)
    df.loc[(df['dayOfYear'] > 166) & (df['dayOfYear'] <= 246), 'Eqt'] = -6.5 * np.sin(np.pi * (df['dayOfYear'] - 166) / 80)
    df.loc[(df['dayOfYear'] > 246) & (df['dayOfYear'] <= 366), 'Eqt'] = 16.4 * np.sin(np.pi * (df['dayOfYear'] - 247) / 113)  # why is this 247 and not 246?

def azimuth(df, lat, lon):
    retDF = pd.DataFrame(index=df.index)
    calcTimevals(retDF)
    calcEqt(retDF)

    retDF['gamma'] = (2 * np.pi) / retDF['numDays'] * (retDF['dayOfYear'] - 1 + (retDF['hour'] - 12) / 24)
    retDF['eqtNOAA'] = 229.18 * ((0.000075 + 0.001868 * np.cos(retDF['gamma']) - 0.032077 * np.sin(retDF['gamma'])
                              - (0.014615 * np.cos(2 * retDF['gamma'])) - (0.040849 * np.sin(2 * retDF['gamma']))))


    retDF['timeOffsetMinutes'] = retDF['eqtNOAA'] + (4 * lon) - (60 * retDF['tzOffsetHours'])


    retDF['Tsolar'] = retDF['hour'] + retDF['timeOffsetMinutes'] / 60 + \
                      retDF['minute'] / 60 + retDF['second'] / 3600
    retDF['hr_angleD'] = -1 *(np.degrees((np.pi / 12.0) * (12-retDF['Tsolar'])))

    retDF['dec_angleD'] = np.degrees(23.45 * np.pi / 180) * np.sin(2 * np.pi * (284 + retDF['dayOfYear']) / 365)

    retDF['cos_zen'] = (sind(lat) * sind(retDF['dec_angleD']) +
                        cosd(lat) * cosd(retDF['dec_angleD']) * cosd(retDF['hr_angleD']))

    retDF['zen_angleD'] = acosd(retDF['cos_zen'])

    # Apply the numpy arctan2 function:

    retDF['y'] = sind(retDF['hr_angleD'])
    retDF['x'] = (cosd(retDF['hr_angleD']) * sind(lat)) - (tand(retDF['dec_angleD']) * cosd(lat))
    retDF['RES'] = np.arctan2(retDF['y'], retDF['x'])
    retDF['RES_D'] = np.degrees(retDF['RES'])

    retDF['sol_azD'] = (retDF['RES_D'] + 180)

    retDF['zenith'] = retDF['zen_angleD']
    retDF['azimuth'] = retDF['sol_azD']
    return retDF


    # Estimate DNI and DHI given GHI and cos_zen and time values provided in df Or Dict


def calcErbs(df, altitude=0, GHICol='GHI', cosZenCol='cos_zen'):
    retDF = pd.DataFrame(index=df.index)
    calcTimevals(retDF)

    retDF['GHI'] = df[GHICol]             # the purpose of these is to determine DHI/DNI split from DAS GHI
    retDF['cos_zen'] = df[cosZenCol]

    retDF['b'] = (2 * np.pi) / retDF['numDays'] * retDF['dayOfYear']
    retDF['Ea'] = (ESC * (1.00011 + 0.034221 * np.cos(retDF['b']) + 0.00128 * np.sin(retDF['b']) +
                          0.000719 * np.cos(2.0 * retDF['b']) + 0.0000077 * np.sin(2.0 * retDF['b'])))

    retDF['EaH'] = 0.0
    retDF.loc[retDF['cos_zen'] > 0.0, 'EaH'] = retDF['Ea'] * retDF['cos_zen']   # "sun up" only

    retDF['GHIMax'] = retDF['EaH'] * GHI_LIM                                    # setting maximum GHI value
    retDF.loc[retDF['GHI'] > retDF['GHIMax'], 'GHI'] = retDF['GHIMax']          # limiting GHI to EaH

    retDF['Kt'] = 0.0
    retDF.loc[retDF['cos_zen'] > 0.0, 'Kt'] = retDF['GHI'] / retDF['EaH']       # "sun up" only

    retDF['Kd'] = (0.9511 - 0.1604 * retDF['Kt'] + 4.388 * retDF['Kt'] ** 2
                   - 16.638 * retDF['Kt'] ** 3.0 + 12.336 * retDF['Kt'] ** 4)

    retDF.loc[retDF['Kt'] < 0.22, 'Kd'] = 1 - 0.09 * retDF['Kt']
    retDF.loc[retDF['Kt'] > 0.80, 'Kd'] = 0.165

    retDF['DHI'] = retDF['GHI'] * retDF['Kd']

    retDF['DNI'] = (retDF['GHI'] - retDF['DHI']) / retDF['cos_zen']

    retDF['DNIMax'] = retDF['Ea'] * DNI_LIM                                # setting maximum DNI value
    retDF.loc[retDF['DNI'] > retDF['DNIMax'], 'DNI'] = retDF['DNIMax']     # limiting DNI to Ea
    return retDF






def calcEnergy(df, tiltD, azimuthD, cellArea, systemCapacity, testTemp=25, losses=0.087, systemEfficiency=SYSTEM_EFFICIENCY,
               zenithCol='zenith', solAzCol='sol_azD',
               ghiCol='GHI', dniCol='DNI', dhiCol='DHI',
               windSpeedCol='windSpeed', tempCol='temp', TEMP_COEFF = -0.0043
               ):
    retDF = pd.DataFrame(index=df.index)
    retDF['zenith'] = df[zenithCol]
    retDF['sol_azD'] = df[solAzCol]
    retDF['GHI'] = df[ghiCol]
    retDF['DNI'] = df[dniCol]
    retDF['DHI'] = df[dhiCol]
    retDF['windSpeed'] = df[windSpeedCol]
    retDF['temp'] = df[tempCol]


    retDF['AOID'] = acosd(cosd(retDF['zenith']) * cosd(tiltD) +
                          (sind(retDF['zenith']) * sind(tiltD) * cosd(retDF['sol_azD'] - azimuthD)))



    retDF['Eb'] = retDF['DNI'] * cosd(retDF['AOID'])
    retDF['Ed'] = (retDF['DHI'] * (1 + cosd(tiltD)) / 2 + retDF['GHI'] * (0.012 * retDF['zenith'] - 0.04) * (1 - cosd(tiltD)) / 2)
    retDF['Epoa'] = retDF['Eb'] + retDF['Ed']

    retDF['Tmod'] = retDF['Epoa'] * np.exp(a1 + b1 * retDF['windSpeed']) + retDF['temp']
    retDF['Tcell'] = retDF['Tmod'] + (retDF['Epoa'] / 1000) * DELT
    retDF['Tcorr'] = 1.0 + TEMP_COEFF * (retDF['Tcell'] - testTemp)

    retDF['Eest'] = (retDF['Epoa'] * retDF['Tcorr'] * (1.0 - losses) * systemEfficiency * cellArea / 1000)

    retDF.loc[retDF['Eest'] > systemCapacity, 'Eest'] = systemCapacity

    return retDF

def calcEnergyByCapacity(df, tiltD, azimuthD, capacity, testTemp=25, losses=0.087, systemEfficiency=SYSTEM_EFFICIENCY,
               zenithCol='zenith', solAzCol='sol_azD',
               ghiCol='GHI', dniCol='DNI', dhiCol='DHI',
               windSpeedCol='windSpeed', tempCol='temp', TEMP_COEFF = -0.0043):

    cellArea = capacity * 1000 * (MODULE_AREA/MODULE_WATTAGE)

    return calcEnergy(df, tiltD, azimuthD, cellArea, capacity, testTemp, losses, systemEfficiency,
                      zenithCol, solAzCol, ghiCol, dniCol, dhiCol, windSpeedCol, tempCol, TEMP_COEFF
                      )


'''
The following functions are derived from the Sandia Laboratory's PV_Lib Toolbox
 https://pvpmc.sandia.gov/applications/pv_lib-toolbox/
'''

#
# Estimate GHI from the cos_zen value passed in df
#

def haurwitz(df, cosZenCol='cos_zen'):
    '''
    Determine clear sky GHI from Haurwitz model.

    Implements the Haurwitz clear sky model for global horizontal
    irradiance (GHI) as presented in [1, 2]. A report on clear
    sky models found the Haurwitz model to have the best performance
    in terms of average monthly error among models which require only
    zenith angle [3].

    Parameters
    ----------
    apparent_zenith : Series
        The apparent (refraction corrected) sun zenith angle
        in degrees.

    Returns
    -------
    pd.DataFrame
    The modeled global horizonal irradiance in W/m^2 provided
    by the Haurwitz clear-sky model.

    Initial implementation of this algorithm by Matthew Reno.

    References
    ----------

    [1] B. Haurwitz, "Insolation in Relation to Cloudiness and Cloud
     Density," Journal of Meteorology, vol. 2, pp. 154-166, 1945.

    [2] B. Haurwitz, "Insolation in Relation to Cloud Type," Journal of
     Meteorology, vol. 3, pp. 123-124, 1946.

    [3] M. Reno, C. Hansen, and J. Stein, "Global Horizontal Irradiance Clear
     Sky Models: Implementation and Analysis", Sandia National
     Laboratories, SAND2012-2389, 2012
    '''
    retDF = pd.DataFrame(index=df.index)
    retDF['cos_zen'] = df[cosZenCol]
    retDF['GHI'] = 0.0
    retDF.loc[retDF['cos_zen'] > 0, 'GHI'] = 1098.0 * retDF['cos_zen'] * np.exp(-0.059 / retDF['cos_zen'])
    return retDF
