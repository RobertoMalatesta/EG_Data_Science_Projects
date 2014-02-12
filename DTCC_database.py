# -*- coding: utf-8 -*-
"""
This program downloads all of a year's interest rate product data from the Depository Trust
and Clearing Corporation. The data contain prices, quantities, and characteristics for individual and block
trades. The program compiles a list of trading days, downloads the .zip files from the DTCC,
unzips these files into .csv files, and appends all the .csv files to make a flat database.

@author: Hamed Faquiryan
"""

# The following lines import all the modules we'll need
import urllib2
import numpy
import datetime
import calendar
import zipfile
import os


year = 2013 #Year of current data
m_lim = 11  #Current month
m_low = 3   #First month of data for the first year of data


# This all the DTCC website's info as well as the filepath for the downloads
baseurl = 'https://kgc0418-tdw-data-0.s3.amazonaws.com/slices/CUMULATIVE_EQUITIES_' + str(year) + '_'
ext = '.zip'
filepath = '/Users/Hamed/Desktop/dtcc_equities//'
rates = 'EQUITIES.csv'

months = []
x = numpy.arange(m_low,m_lim+1,1)
for it in range(len(x)):
    months.append(str(x[it]))
    if len(months[it]) < 2:
        months[it] = '0' + months[it]
del x

## Now that all our .zip files are downloaded, we're goin to extract the .csv files in each
## of them and rename them by year_month_day

for month in range(len(months)):
    days_lim = calendar.monthrange(year, int(months[month]))[1]
    days = []
    for day in numpy.arange(1, days_lim + 1, 1):
        if len(str(day)) < 2 and datetime.date(year, int(months[month]), day).weekday() <= 4:
            days.append('0' + str(day))
        elif datetime.date(year, int(months[month]), day).weekday() <= 4:
            days.append(str(day))
    for close in range(len(days)):
        try:
            doi = urllib2.urlopen(baseurl + months[month] + '_' + days[close] + ext)
            with open(filepath + months[month] + '_' + days[close] + ext, 'wb') as code:
                code.write(doi.read())
            x = open(filepath + months[month] + '_' + days[close] + ext, 'rb')
            zfile = zipfile.ZipFile(x)
            zfile.extract(rates, filepath)
            os.rename(filepath + rates, filepath + months[month] + '_' + days[close] + '.csv')
            #print months[month] + '_' + days[close] + ext + ' has been saved as a .csv file.'
        except urllib2.URLError:
            print 'Either this date was a holiday, the URL was invalid, or something has gone terribly wrong ...'
    print 'All trading days in ' + calendar.month_name[month + 1] + ' have been saved as .csv files.'

# Now we consolidate all the .csv files into one .csv file
dtcc_rates = open(filepath + 'dtcc_equities.csv', 'a')

for line in open(filepath + '02_28.csv'):
    dtcc_rates.write(line)

for month in range(len(months)):
    days_lim = calendar.monthrange(year, int(months[month]))[1]
    days = []
    for day in numpy.arange(1, days_lim + 1, 1):
        if len(str(day)) < 2 and datetime.date(year, int(months[month]), day).weekday() <= 4:
            days.append('0' + str(day))
        elif datetime.date(year, int(months[month]), day).weekday() <= 4:
            days.append(str(day))
    for close in range(len(days)):
        try:
            f = open(filepath + months[month] + '_' + days[close] + '.csv')
            f.next()
            for line in f:
                dtcc_rates.write(line)
            f.close()
        except IOError:
            break
        print "We're all done, friend."
dtcc_rates.close()
            




