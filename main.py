import math as ma
import QuantLib as Ql
import numpy as np
import datetime as dt
import rv_ran4_generator as ran4
import openpyxl
from statistics import NormalDist
from matplotlib import pyplot as plt

wb = openpyxl.load_workbook("/Users/tablo/PycharmProjects/Equity/My2Calc200922_SofrCurve.xlsm", data_only=True)
ws = wb['SpreadOption']


effective_date = Ql.Date(25,9,2020)
termination_date = Ql.Date(25,9,2023)
tenor = Ql.Period(6, Ql.Months)
# calendar = Ql.UnitedStates()
calendar = Ql.SouthKorea()
business_convention = Ql.Following
termination_biz_convention = Ql.Following
date_gen = Ql.DateGeneration.Forward
end_of_month = False

schedule = Ql.Schedule(effective_date,termination_date,tenor,calendar,business_convention,termination_biz_convention,date_gen,end_of_month)
total_Days = schedule[6] - effective_date
number_simulation = 10
number_exercise = 6

S0 = ws.cell(1,17).value
BP = S0
strike = np.zeros(number_exercise)
rf = np.zeros(number_exercise)
vol = np.zeros(number_exercise)
div = np.zeros(number_exercise)

 = 0.03
 = 0.30
 = 0.01



dc = 365
dt = 1/dc

St = np.zeros((total_Days+1,number_simulation))

for i in range(number_simulation):
    St[0,i] = S0 / BP * 100

seed = 1
seed = ran4.rv_ran4_seed(seed)
sum_payoff = 0

for i in range(number_simulation):
    previous_days = 0
    for k in range(number_exercise):
        if k == 0:

            exercise_days = schedule[k+1] - effective_date
        else:
            exercise_days = schedule[k+1] - schedule[k]

        for j in range(exercise_days):
            rand = ran4.rv_ran4_generator(seed)
            seed = ran4.rv_ran4_seed(seed)

            if k == 0:
                St[j+1, i] = St[j, i]*ma.exp((rf - div-0.5*pow(vol, 2))*dt+vol*ma.sqrt(dt)*NormalDist().inv_cdf(rand))
            else:
                St[j + previous_days + 1, i] = St[j + previous_days, i] * ma.exp(
                    (rf - div - 0.5 * pow(vol, 2)) * dt + vol * ma.sqrt(dt) * NormalDist().inv_cdf(rand))

        if k != number_exercise-1:
            if k == 0:
                if St[j + 1, i] > strike:
                    payoff = max(St[j + 1, i] - strike, 0)
                    sum_payoff = sum_payoff + payoff
                    break
                else:
                    payoff = 0
            else:
                if St[j + previous_days + 1, i] > strike:
                    payoff = max(St[j + previous_days + 1, i]-strike,0)
                    sum_payoff = sum_payoff + payoff
                    break
                else:
                    payoff = 0
        else:
            if St[j + previous_days + 1, i] > strike:
                payoff = max(St[j + previous_days + 1, i] - strike, 0)
                sum_payoff = sum_payoff + payoff
            else:
                payoff = 0

        previous_days = previous_days + exercise_days


plt.plot(St)
plt.show()
