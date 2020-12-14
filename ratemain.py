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

# curve_zero = [0.0199994520748,
#               0.0199994520748,
#               0.0226244425935,
#               0.0255242446688,
#               0.0286026309564,
#               0.0298389921414,
#               0.031748737524,
#               0.0339072761895,
#               0.0393287919992,
#               0.0417754957289,
#               0.0431084917179,
#               0.0439861811485,
#               0.0446098027925,
#               0.0452529775785,
#               0.0456788183665,
#               0.0461179776996,
#               0.0465724887031,
#               0.0468109160698,
#               0.0469173680078,
#               0.045489601316]

curve_zero = [0.00459997101394344,
              0.00459997101394344,
              0.00531051622861427,
              0.00595292309511829,
              0.00657271601951631,
              0.00688178568597534,
              0.00703533341029391,
              0.00718717686537302,
              0.00776453736367852,
              0.00834215241195633,
              0.00917268978120266,
              0.0100321798532742,
              0.0105895716739209,
              0.0109701894430125,
              0.0112756010709692,
              0.011635223255438,
              0.0118920955393558,
              0.0120978357905755,
              0.0124872120683835,
              0.0120927281586671,
              0.0115636198036743,
]

# curve_date = [0.00,
#               0.0027397260274,
#               0.0794520547945,
#               0.1643835616438,
#               0.2493150684932,
#               0.4986301369863,
#               0.7506849315068,
#               1.0164383561644,
#               2.0027397260274,
#               3.0109589041096,
#               4.0082191780822,
#               5.0054794520548,
#               6.0054794520548,
#               7.0082191780822,
#               8.0082191780822,
#               9.0219178082192,
#               10.0109589041096,
#               12.0109589041096,
#               15.0164383561644,
#               20.0219178082192]

curve_date = [0,
              0.00273972602739726,
              0.0931506849315069,
              0.172602739726027,
              0.249315068493151,
              0.501369863013699,
              0.753424657534247,
              1.0027397260274,
              1.5013698630137,
              2.0027397260274,
              3.00821917808219,
              4.0054794520548,
              5.0054794520548,
              6.0054794520548,
              7.0054794520548,
              8.01369863013699,
              9.01095890410959,
              10.0082191780822,
              12.0109589041096,
              15.013698630137,
              20.0191780821918]

swaptionMat = [0.25, 0.25, 0.25, 0.25, 0.25, 0.25, 0.25,
               0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50,
               1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
               2.00, 2.00, 2.00, 2.00, 2.00, 2.00, 2.00,
               3.00, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00,
               4.00, 4.00, 4.00, 4.00, 4.00, 4.00, 4.00,
               5.00, 5.00, 5.00, 5.00, 5.00, 5.00, 5.00,
               7.00, 7.00, 7.00, 7.00, 7.00, 7.00, 7.00,
               10.00, 10.00, 10.00, 10.00, 10.00, 10.00, 10.00]

swapMat = [1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10,
           1, 2, 3, 4, 5, 7, 10]

swapR = [0.00761531417955446, 0.00873129252736225, 0.0095959667488482, 0.0103776694381925, 0.0108572267535491, 0.0114874057720390, 0.0122232334109845,
         0.00812627942943411, 0.00918076141703711, 0.0100339742373390, 0.0107269752201200, 0.0111426112379323, 0.0117375074257768, 0.0124045626904274,
         0.00943478192056441, 0.01017641967244130, 0.0109962444031612, 0.0114361844742637, 0.0117091207613959, 0.0122390554739941, 0.0127585594529357,
         0.01092580717105860, 0.01178879102050540, 0.0121179509528559, 0.0122939769763345, 0.0124520867818963, 0.0128915823491173, 0.0132866795939864,
         0.01265981334898270, 0.01272424227069890, 0.0127610327997699, 0.0128454598129341, 0.0130993554270435, 0.0133434013161059, 0.0133035790633999,
         0.01282423519471640, 0.01284794631207350, 0.0129344659425773, 0.0132504032256154, 0.0133766853370088, 0.0135705692904096, 0.0131210290898562,
         0.01287121973581180, 0.01298894193709630, 0.0133703980280466, 0.0135275000392293, 0.0136035047229035, 0.0138264167504430, 0.0128814968024174,
         0.01414709044193680, 0.01404167533568780, 0.0140392158819286, 0.0140791827096224, 0.0141790282560767, 0.0132012473225838, 0.0123576976783184,
         0.01424258406941620, 0.01442642327383080, 0.0132238452242891, 0.0125704800676329, 0.0121098314925947, 0.0115976075451230, 0.0110667964101808]

today_date = Ql.Date(7, 12, 2020)
effective_date = Ql.Date(today_date.serialNumber() + 1)

calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)
while not calendar.isBusinessDay(effective_date):
    effective_date = Ql.Date(effective_date.serialNumber() + 1)

notional = 100
mat = 10
K = 0.0463552554
# a = 0.0948
a = -0.0021
# sigma = 0.01413
sigma = 0.0043
X = 0.046355

term_date = Ql.Date(today_date.dayOfMonth(), today_date.month(), today_date.year() + mat)
freq = Ql.Period('3M')
conv = Ql.ModifiedFollowing
term_date_conv = Ql.ModifiedFollowing
my_rule = Ql.DateGeneration.Backward
end_month = False
schedule = Ql.Schedule(effective_date, term_date, freq, calendar, conv, term_date_conv, my_rule, end_month)

PaymentD = np.ones((1, 4 * mat))
for i in range(40):
    PaymentD[0, i] = (schedule.__getitem__(i + 1) - today_date) / 365

effective_day = (effective_date - today_date) / 365

methods = Ql.LinearInterpolation(curve_date, curve_zero)

n = PaymentD.size

sum_ = 0
for i in range(n):
    if i == 0:
        tau_ = PaymentD[0, i] - effective_day
        BondP = ma.exp(-methods(effective_day)*effective_day)
        BondN = ma.exp(-methods(PaymentD[0, i]) * PaymentD[0, i])
        B = (1-ma.exp(-a*tau_))/a
        sigmaP = sigma*ma.sqrt((1-ma.exp(-2*a*effective_day))/2/a)*B
        h = ma.log(BondN*(1+X*tau_)/BondP)/sigmaP+sigmaP/2
        sum_ = sum_ + BondP*norm.cdf(-h+sigmaP) - (1+X*tau_)*BondN*norm.cdf(-h)
    else:
        tau_ = PaymentD[0, i] - PaymentD[0, i-1]
        BondP = ma.exp(-methods(PaymentD[0, i-1]) * PaymentD[0, i-1])
        BondN = ma.exp(-methods(PaymentD[0, i]) * PaymentD[0, i])
        B = (1-ma.exp(-a*tau_))/a
        sigmaP = sigma * ma.sqrt((1 - ma.exp(-2 * a * PaymentD[0, i-1])) / 2 / a) * B
        h = ma.log(BondN * (1 + X * tau_) / BondP) / sigmaP + sigmaP / 2
        sum_ = sum_ + BondP*norm.cdf(-h+sigmaP) - (1+X*tau_)*BondN*norm.cdf(-h)

Price = notional*sum_
print(Price)

effective_date
optionMatD = dt.datetime(effective_date.year(), effective_date.month(), effective_date.dayOfMonth()) + dateutil.relativedelta.relativedelta(months=12*swaptionMat[0])
optionMatD = Ql.Date(optionMatD.day, optionMatD.month, optionMatD.year)

term_date = Ql.Date(optionMatD.dayOfMonth(), optionMatD.month(), optionMatD.year() + swapMat[0])
freq = Ql.Period('3M')
conv = Ql.ModifiedFollowing
term_date_conv = Ql.ModifiedFollowing
my_rule = Ql.DateGeneration.Backward
end_month = False

swaption_schedule = Ql.Schedule(optionMatD, term_date, freq, calendar, conv, term_date_conv, my_rule, end_month)

n = 4*swapMat[0]
swaption_PaymentD = np.zeros(n)
daycontfrac = np.zeros(n)

for i in range(n):
    swaption_PaymentD[i] = (swaption_schedule.__getitem__(i+1) - today_date)/365
    daycontfrac[i] = (swaption_schedule.__getitem__(i+1) - swaption_schedule.__getitem__(i)) / 365

alpha = (optionMatD - today_date)/365
C = swapR[0]*daycontfrac
C[n-1] = 1+C[n-1]

PM_0_t = ma.exp(-methods(alpha)*alpha)
dtt = 0.00001
DF0 = ma.exp(-methods(alpha)*alpha)
DF1 = ma.exp(-methods(alpha+dtt)*(alpha+dtt))
InstFwd = -ma.log(DF1/DF0)/dtt

PM_0_T = np.zeros(n)
B_t_T = np.zeros(n)

for i in range(n):
    PM_0_T[i] = ma.exp(-methods(swaption_PaymentD[i])*swaption_PaymentD[i])
    B_t_T[i] = (1-ma.exp(-a*(swaption_PaymentD[i]-alpha)))/a

# Newton Raphson Method

xT_ = 0
dx = 0
F = 100
err = 10e-15

while F > err:
    xT_ = xT_ + dx
    F = 0
    FPrime = 0
    for k in range(n):
        TempV = C[k]*PM_0_T[k]/PM_0_t*ma.exp(B_t_T[k]*InstFwd-pow(sigma, 2)/(4*a)*(1-ma.exp(-2*a*alpha))*pow(B_t_T[k], 2)-B_t_T[k]*xT_)
        F = F + TempV
        FPrime = FPrime + TempV*(-B_t_T[k])

    F = F-1
    dx = -F/FPrime

# #############

X_Bond_P = np.zeros(n)
for i in range(n):
    X_Bond_P[i] = PM_0_T[i]/PM_0_t*ma.exp(-pow(sigma, 2)/(4*a)*(1-ma.exp(-2*a*alpha))*pow(B_t_T[i], 2)-B_t_T[i]*xT_+B_t_T[i]*InstFwd)

swapion_sum_ = 0
x_t = 0
t = 0
# HW_1F_BondP
MktBond_T = ma.exp(-methods(alpha, allowExtrapolation=True)*alpha)
MktBond_t = ma.exp(-methods(t, allowExtrapolation=True)*t)
B = (1-ma.exp(-a*(alpha - t)))/a
BondPrice_t_T_ = MktBond_T/MktBond_t*ma.exp(-pow(sigma, 2)/(4*a)*(1-ma.exp(-2*a*t))*pow(B, 2)-pow(sigma, 2)/(2*pow(a, 2))*pow((1-ma.exp(-a*t)), 2)*B - B*x_t)
BondPrice_t_T = hw1f.hw1f_bond_p(curve_date, curve_zero, a, sigma, 0, alpha, x_t)
# ####

for i in range(n):
    BondPrice_t_S = hw1f.hw1f_bond_p(curve_date, curve_zero, a, sigma, 0, swaption_PaymentD[i], x_t)
    sigmaP = sigma*ma.sqrt((1-ma.exp(-2*a*alpha))/(2*a))*B_t_T[i]
    h = ma.log(BondPrice_t_S/(BondPrice_t_T*X_Bond_P[i]))/sigmaP + sigmaP/2
    ZBP = X_Bond_P[i]*BondPrice_t_T*norm.cdf(-h+sigmaP) - BondPrice_t_S*norm.cdf(-h)
    swapion_sum_ = swapion_sum_+C[i]*ZBP

PayerSwaptionPrice = notional*swapion_sum_

print(swaption_schedule.__getitem__(0))