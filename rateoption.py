import math as ma
import numpy as np
from statistics import NormalDist
from scipy.stats import norm
import pandas as pd
import datetime as dt
from matplotlib import pyplot as plt
import QuantLib as Ql
import calendar as cal
import dateutil.relativedelta
from model import hw1f_model as hw1f

curve_zero = [0.019999452,
              0.019999452,
              0.022624443,
              0.025524245,
              0.028602631,
              0.029838992,
              0.031748738,
              0.033907276,
              0.039328792,
              0.041775496,
              0.043108492,
              0.043986181,
              0.044609803,
              0.045252978,
              0.045678818,
              0.046117978,
              0.046572489,
              0.046810916,
              0.046917368,
              0.045489601]

curve_date = [0,
              0.002739726,
              0.079452055,
              0.164383562,
              0.249315068,
              0.498630137,
              0.750684932,
              1.016438356,
              2.002739726,
              3.010958904,
              4.008219178,
              5.005479452,
              6.005479452,
              7.008219178,
              8.008219178,
              9.021917808,
              10.0109589,
              12.0109589,
              15.01643836,
              20.02191781]

today_date = Ql.Date(1, 2, 2010)
Mat = 10
Notional = 100
sigma = 0.16550020
K = 0.0463552554

# ###
effective_date = Ql.Date(today_date.serialNumber() + 1)
calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)
while not calendar.isBusinessDay(effective_date):
    effective_date = Ql.Date(effective_date.serialNumber() + 1)

optionMatD = dt.datetime(effective_date.year(), effective_date.month(), effective_date.dayOfMonth()) + dateutil.relativedelta.relativedelta(months=12*Mat)
optionMatD = Ql.Date(optionMatD.day, optionMatD.month, optionMatD.year)

freq = Ql.Period('3M')
convention = Ql.ModifiedFollowing
term_date_convention = Ql.ModifiedFollowing
my_rule = Ql.DateGeneration.Backward
end_month = False
schedule = Ql.Schedule(effective_date, optionMatD, freq, calendar, convention, term_date_convention, my_rule, end_month)
# Unadj_couponDate == schdule.__getitem__()
# Unadj_couponDate == fixRateStart
# fixRateStart == schdule.__getitem__()
FixingD = np.zeros(4*Mat)
fixRateStart = np.zeros(4*Mat)
fixRateEnd = np.zeros(4*Mat)
for i in range(4*Mat):
    FixingD[i] = (hw1f.my_working_day(schedule.__getitem__(i), 'KOREA', -1) - today_date)/365
    fixRateStart[i] = (schedule.__getitem__(i) - today_date)/365
    fixRateEnd[i] = (schedule.__getitem__(i+1) - today_date)/365

fixPer = fixRateEnd-fixRateStart

methods = Ql.LinearInterpolation(curve_date, curve_zero)
df1 = np.zeros(4*Mat)
df2 = np.zeros(4*Mat)
fixRate = np.zeros(4*Mat)
PaymentD = np.zeros(4*Mat)
DayContFraction = np.zeros(4*Mat)
for i in range(4*Mat):
    df1[i] = ma.exp(-methods(fixRateStart[i], allowExtrapolation=True)*fixRateStart[i])
    df2[i] = ma.exp(-methods(fixRateEnd[i], allowExtrapolation=True)*fixRateEnd[i])
    fixRate[i] = ((df1[i] - df2[i])/df2[i])/fixPer[i]
    PaymentD[i] = (schedule.__getitem__(i+1) - today_date)/365
    DayContFraction[i] = (schedule.__getitem__(i+1) - schedule.__getitem__(i))/365


normal_d = np.zeros(4*Mat)
d1 = np.zeros(4*Mat)
d2 = np.zeros(4*Mat)
normal_cap = np.zeros(4*Mat)
DF = np.zeros(4*Mat)
Caplet = np.zeros(4*Mat)
for i in range(4*Mat):
    normal_d[i] = (fixRate[i] - K)/(sigma*ma.sqrt(FixingD[i]))
    normal_cap[i] = (fixRate[i] - K)*norm.cdf(normal_d[i])+sigma*ma.sqrt(FixingD[i])*norm.pdf(normal_d[i])
    DF[i] = ma.exp(-methods(PaymentD[i], allowExtrapolation=True)*PaymentD[i])
    Caplet[i] = Notional*normal_cap[i]*DayContFraction[i]*DF[i]

ON_DF = ma.exp(-methods(curve_date[1], allowExtrapolation=True)*curve_date[1])
price = sum(Caplet)/ON_DF
print(price)

# Normal Swaption
today_date = Ql.Date(1, 2, 2010)
Mat = 2
swapMat = 3
Notional = 100
sigma = 0.176
X = 0.05

effective_date = Ql.Date(today_date.serialNumber() + 1)
calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)
while not calendar.isBusinessDay(effective_date):
    effective_date = Ql.Date(effective_date.serialNumber() + 1)

optionMatD = dt.datetime(effective_date.year(), effective_date.month(), effective_date.dayOfMonth()) + dateutil.relativedelta.relativedelta(months=12*Mat)
optionMatD = Ql.Date(optionMatD.day, optionMatD.month, optionMatD.year)

term_date = Ql.Date(optionMatD.dayOfMonth(), optionMatD.month(), optionMatD.year() + swapMat)
freq = Ql.Period('3M')
convention = Ql.ModifiedFollowing
term_date_convention = Ql.ModifiedFollowing
my_rule = Ql.DateGeneration.Backward
end_month = False
schedule = Ql.Schedule(optionMatD, term_date, freq, calendar, convention, term_date_convention, my_rule, end_month)
# Unadj_couponDate == schdule.__getitem__()

PaymentD = np.zeros(4*swapMat)
DayContFraction = np.zeros(4*swapMat)
for i in range(4*swapMat):
    DayContFraction[i] = (schedule.__getitem__(i+1) - schedule.__getitem__(i))/365
    PaymentD[i] = (schedule.__getitem__(i+1) - today_date)/365

alpha = (optionMatD-today_date)/365
methods = Ql.LinearInterpolation(curve_date, curve_zero)
df_a = ma.exp(-methods(alpha, allowExtrapolation=True)*alpha)
df_ = np.zeros(4*swapMat)
for i in range(4*swapMat):
    df_[i] = ma.exp(-methods(PaymentD[i])*PaymentD[i])
df_b = df_[df_.size-1]

temp = sum(DayContFraction*df_)
swapR = (df_a-df_b)/temp

d = (swapR-X)/(sigma*ma.sqrt(alpha))
price = ((swapR-X)*norm.cdf(d)+sigma*ma.sqrt(alpha)*norm.pdf(d))*temp*Notional
print(price)

test_date0 = Ql.Date(10, 12, 2020)
print('test_date is ', test_date0)

test_date = hw1f.my_working_day(test_date0, 'KOREA', 1)
print('test_date is ', test_date)

test_date = hw1f.my_working_day(test_date0, 'KOREA', 2)
print('test_date is ', test_date)

test_date0 = Ql.Date(10, 12, 2020)
print('test_date is ', test_date0)

test_date = hw1f.my_working_day(test_date0, 'KOREA', -1)
print('test_date is ', test_date)

test_date = hw1f.my_working_day(test_date0, 'KOREA', -2)
print('test_date is ', test_date)

test_date = hw1f.my_working_day(test_date0, 'KOREA', -5)
print('test_date is ', test_date)

print(schedule.__getitem__(1))
print(schedule.previousDate(schedule.__getitem__(1)))