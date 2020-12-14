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


def my_working_day(ref_date, my_calendar, n):
    if my_calendar == 'KOREA':
        calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)

    if n > 0:
        for i in range(n):
            ref_date = Ql.Date(ref_date.serialNumber() + 1)
            while not calendar.isBusinessDay(ref_date):
                ref_date = Ql.Date(ref_date.serialNumber() + 1)
    elif n < 0:
        for i in range(-n):
            ref_date = Ql.Date(ref_date.serialNumber() - 1)
            while not calendar.isBusinessDay(ref_date):
                ref_date = Ql.Date(ref_date.serialNumber() - 1)
    else:
        ref_date = ref_date

    return ref_date


def hw1f_bond_p(curve_date, curve_zero, a, sigma, t, mat, x):
    methods = Ql.LinearInterpolation(curve_date, curve_zero)
    mkt_bond_mat = ma.exp(-methods(mat, allowExtrapolation=True) * mat)
    mkt_bond_t = ma.exp(-methods(t, allowExtrapolation=True) * t)
    b = (1 - ma.exp(-a * (mat - t))) / a
    bond_price_t_mat = mkt_bond_mat/mkt_bond_t*ma.exp(-pow(sigma, 2)/(4*a)*(1-ma.exp(-2*a*t))*pow(b, 2)-pow(sigma, 2)/(2*pow(a, 2))*pow((1-ma.exp(-a*t)), 2)*b - b*x)
    return bond_price_t_mat
