from QuantLib import *
import matplotlib.pyplot as plt
import numpy as np
# from scipy.integrate import simps, cumtrapz, romb
# matplotlib inline
import math

# option parameters
strike_price = 110.0
payoff = PlainVanillaPayoff(Option.Call, strike_price)
# option data
maturity_date = Date(15, 1, 2016)
spot_price = 127.62
strike_price = 130
volatility = 0.20
# the historical vols for a year
dividend_rate = 0.0163
option_type = Option.Call
risk_free_rate = 0.001
day_count = Actual365Fixed()
calendar = UnitedStates()
calculation_date = Date(8, 5, 2015)
Settings.instance().evaluationDate = calculation_date

# construct the European Option
payoff = PlainVanillaPayoff(option_type, strike_price)
exercise = EuropeanExercise(maturity_date)
european_option = VanillaOption(payoff, exercise)

v0 = volatility*volatility
# spot variance
kappa = 0.1
theta = v0
sigma = 0.1
rho = -0.75
spot_handle = QuoteHandle(SimpleQuote(spot_price))

flat_ts = YieldTermStructureHandle(FlatForward(calculation_date, risk_free_rate, day_count))
dividend_yield = YieldTermStructureHandle(FlatForward(calculation_date, dividend_rate, day_count))
heston_process = HestonProcess(flat_ts, dividend_yield,spot_handle, v0, kappa,theta, sigma, rho)

engine = AnalyticHestonEngine(HestonModel(heston_process),0.01, 1000)
european_option.setPricingEngine(engine)
h_price = european_option.NPV()
print("The Heston model price is",h_price)

flat_vol_ts = BlackVolTermStructureHandle(BlackConstantVol(calculation_date, calendar,volatility, day_count))
bsm_process = BlackScholesMertonProcess(spot_handle, dividend_yield,flat_ts, flat_vol_ts)
european_option.setPricingEngine(AnalyticEuropeanEngine(bsm_process))
bs_price = european_option.NPV()
print("The Black-Scholes-Merton model price is ", bs_price)