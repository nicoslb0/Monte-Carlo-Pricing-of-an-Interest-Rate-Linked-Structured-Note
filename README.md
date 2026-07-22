# Monte-Carlo-Pricing-of-an-Interest-Rate-Linked-Structured-Note

You are working as a junior quantitative analyst at a hedge fund specialising in fixed￾income and structured products. The fund is currently evaluating the potential 
investment in a structured note whose coupon payments depend on the evolution of 
benchmark interest rates. In particular, the note’s coupon payments are linked to two 
benchmark interest rates, denoted BIR1 and BIR2, while discounting will be based on a 
separate discounting interest rate (DIR). As part of the analysis team, you have been asked 
to value this instrument using simulation methods and to assess its attractiveness as an 
investment opportunity. Your task is to simulate the relevant interest-rate dynamics at 
monthly frequency, estimate the fair value of the instrument using Monte Carlo 
techniques, and provide a clear investment recommendation. The analysis should be 
written as a coherent report, supported by tables, figures, and numerical results where 
appropriate, explaining both the methodology used and the economic intuition behind 
your conclusions. 
 
Structured Note Characteristics 
The hedge fund is analysing a structured note whose coupon payments depend on the 
evolution of two benchmark interest rates. The main characteristics of the instrument are 
as follows: 
• Issue date: assumed to be today (valuation at issuance) 
• Maturity: to be chosen by the student 
• Coupon payments: annual 
• Underlying benchmark rates: 
o BIR1: Benchmark Interest Rate 1 
o BIR2: Benchmark Interest Rate 2 
• Discounting rate: 
o DIR: Discounting Interest Rate 
• Frequency of interest-rate simulation: monthly 
At each coupon date 𝑡, the coupon rate of the note is determined according to the 
following rule: 
𝑐𝑡 = min⁡{𝐶,  max⁡[𝐹,  𝛼(𝐵𝐼𝑅1𝑡 − 𝐵𝐼𝑅2𝑡) + 𝛽𝐵𝐼𝑅1𝑡
]} 
 
where:  • 𝐹denotes the floor rate, • 𝐶denotes the cap rate,  • 𝛼 and 𝛽are coefficients chosen by the student. 

The coupon payment at each date is determined by applying the above rate to the chosen 
face value of the note. At maturity, the investor receives the final coupon payment together 
with the repayment of the face value. 

Students must simulate the paths of BIR1, BIR2, and DIR at monthly frequency, 
construct the coupon payments implied by the rule above, discount all cash flows using 
DIR, and estimate the fair value of the structured note using Monte Carlo simulation. 
Students are required to clearly justify their modelling choices, including the selected 
maturity, face value, floor (F), cap (C), and coefficients (𝛼, 𝛽). These choices should 
be economically reasonable and explained in the context of interest-rate dynamics and 
structured product design. 
 
Modelling Framework 
In order to simulate the evolution of the benchmark interest rates (BIR1, BIR2, and DIR), 
students should employ a continuous-time stochastic interest rate model. In particular, 
students may use either the Vasicek model or the Cox–Ingersoll–Ross (CIR) model [or 
both], depending on which specification they consider most appropriate for the 
behaviour of the simulated rates.

 
Tasks 
Task 1 – Baseline Valuation and Investment Recommendation 
Assume that the structured note is purchased at issuance. Using Monte Carlo simulation, 
simulate the paths of BIR1, BIR2, and DIR at monthly frequency over the chosen maturity. 
Based on the simulated paths, construct the annual coupon payments implied by the 
coupon rule and compute the present value of all future cash flows using the simulated 
DIR. 
Estimate the fair value of the structured note and present the distribution of simulated 
prices, including appropriate summary statistics and confidence intervals. Based on your 
results, provide a clear investment recommendation, explaining whether the note 
appears attractive relative to its face value. 
 
Task 2 – Alternative Discounting Rate 
Repeat the valuation exercise using an alternative discounting interest rate (DIR). 
Compare the resulting fair value with the one obtained in Task 1 and discuss how and 
why the valuation changes. In your discussion, explain the role that the discount rate plays 
in the pricing of fixed-income instruments. 
 
Task 3 – Shorter Maturity Scenario 
Re-evaluate the structured note assuming a maturity equal to half (or approximately 
half) of the originally chosen maturity. Keeping the same modelling framework and 
parameters, simulate the interest-rate paths again and compute the new fair value of the 
note. 
Compare the results with those obtained in Task 1 and explain how the change in maturity 
affects the valuation and the investment decision. 
 
