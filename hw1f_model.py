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


def normal_swaption(curve_zero, curve_date, today_date, mat, swapmat, notional, sigma, x):
    effective_date = Ql.Date(today_date.serialNumber() + 1)
    calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)
    while not calendar.isBusinessDay(effective_date):
        effective_date = Ql.Date(effective_date.serialNumber() + 1)

    option_mat_d = dt.datetime(effective_date.year(), effective_date.month(), effective_date.dayOfMonth()) + dateutil.relativedelta.relativedelta(months=12 * mat)
    option_mat_d = Ql.Date(option_mat_d.day, option_mat_d.month, option_mat_d.year)

    term_date = Ql.Date(option_mat_d.dayOfMonth(), option_mat_d.month(), option_mat_d.year() + swapmat)
    freq = Ql.Period('3M')
    convention = Ql.ModifiedFollowing
    term_date_convention = Ql.ModifiedFollowing
    my_rule = Ql.DateGeneration.Backward
    end_month = False
    schedule = Ql.Schedule(option_mat_d, term_date, freq, calendar, convention, term_date_convention, my_rule, end_month)

    payment_d = np.zeros(4 * swapmat)
    day_cont_fraction = np.zeros(4 * swapmat)
    for i in range(4 * swapmat):
        day_cont_fraction[i] = (schedule.__getitem__(i + 1) - schedule.__getitem__(i)) / 365
        payment_d[i] = (schedule.__getitem__(i + 1) - today_date) / 365

    alpha = (option_mat_d - today_date) / 365
    methods = Ql.LinearInterpolation(curve_date, curve_zero)
    df_a = ma.exp(-methods(alpha, allowExtrapolation=True) * alpha)
    df_ = np.zeros(4 * swapmat)
    for i in range(4 * swapmat):
        df_[i] = ma.exp(-methods(payment_d[i]) * payment_d[i])
    df_b = df_[df_.size - 1]

    temp = sum(day_cont_fraction * df_)
    swap_r = (df_a - df_b) / temp

    d = (swap_r - x) / (sigma * ma.sqrt(alpha))
    price = ((swap_r - x) * norm.cdf(d) + sigma * ma.sqrt(alpha) * norm.pdf(d)) * temp * notional
    return price


def normal_cap(curve_zero, curve_date, today_date, mat, notional, sigma, k):
    effective_date = Ql.Date(today_date.serialNumber() + 1)
    calendar = Ql.SouthKorea(Ql.SouthKorea.KRX)
    while not calendar.isBusinessDay(effective_date):
        effective_date = Ql.Date(effective_date.serialNumber() + 1)

    option_mat_d = dt.datetime(effective_date.year(), effective_date.month(), effective_date.dayOfMonth()) + dateutil.relativedelta.relativedelta(months=12 * mat)
    option_mat_d = Ql.Date(option_mat_d.day, option_mat_d.month, option_mat_d.year)

    freq = Ql.Period('3M')
    convention = Ql.ModifiedFollowing
    term_date_convention = Ql.ModifiedFollowing
    my_rule = Ql.DateGeneration.Backward
    end_month = False
    schedule = Ql.Schedule(effective_date, option_mat_d, freq, calendar, convention, term_date_convention, my_rule, end_month)
    # Unadj_couponDate == schdule.__getitem__()
    # Unadj_couponDate == fixRateStart
    # fixRateStart == schdule.__getitem__()
    fixing_d = np.zeros(4 * mat)
    fix_rates_start = np.zeros(4 * mat)
    fix_rate_end = np.zeros(4 * mat)
    for i in range(4 * mat):
        fixing_d[i] = (my_working_day(schedule.__getitem__(i), 'KOREA', -1) - today_date) / 365
        fix_rates_start[i] = (schedule.__getitem__(i) - today_date) / 365
        fix_rate_end[i] = (schedule.__getitem__(i + 1) - today_date) / 365

    fix_per = fix_rate_end - fix_rates_start

    methods = Ql.LinearInterpolation(curve_date, curve_zero)
    df1 = np.zeros(4 * mat)
    df2 = np.zeros(4 * mat)
    fix_rate = np.zeros(4 * mat)
    payment_d = np.zeros(4 * mat)
    day_cont_fraction = np.zeros(4 * mat)
    for i in range(4 * mat):
        df1[i] = ma.exp(-methods(fix_rates_start[i], allowExtrapolation=True) * fix_rates_start[i])
        df2[i] = ma.exp(-methods(fix_rate_end[i], allowExtrapolation=True) * fix_rate_end[i])
        fix_rate[i] = ((df1[i] - df2[i]) / df2[i]) / fix_per[i]
        payment_d[i] = (schedule.__getitem__(i + 1) - today_date) / 365
        day_cont_fraction[i] = (schedule.__getitem__(i + 1) - schedule.__getitem__(i)) / 365

    normal_d = np.zeros(4 * mat)
    d1 = np.zeros(4 * mat)
    d2 = np.zeros(4 * mat)
    normal_caplet = np.zeros(4 * mat)
    df = np.zeros(4 * mat)
    caplet = np.zeros(4 * mat)
    for i in range(4 * mat):
        normal_d[i] = (fix_rate[i] - k) / (sigma * ma.sqrt(fixing_d[i]))
        normal_caplet[i] = (fix_rate[i] - k) * norm.cdf(normal_d[i]) + sigma * ma.sqrt(fixing_d[i]) * norm.pdf(normal_d[i])
        df[i] = ma.exp(-methods(payment_d[i], allowExtrapolation=True) * payment_d[i])
        caplet[i] = notional * normal_caplet[i] * day_cont_fraction[i] * df[i]

    on_df = ma.exp(-methods(curve_date[1], allowExtrapolation=True) * curve_date[1])
    price = sum(caplet) / on_df
    return price
