### Title Slide
- Hello and thanks for this opportunity to present our work, titled:
  **Inferring incidence rate ratios from cross-sectional odds ratios
  using individual-based modelling - applied to depression & harmful drinking**
### Background & Motivation
- First, a bit of background
- **Major depression** is characterized by low mood and reduced enjoyment of daily activities
  - Depending on the specific definition used and population studied overall depression prevalence typically ranges from 1 to 10%
- **Harmful drinking** reflects patterns of alcohol consumption that meaningfully increase risks of related outcomes, including dependence
  - Definitions of harmful drinking are perhaps even more varied (than depression) but overall prevalence typically ranges from 2 to 20%
- Many studies have identified **associations** between depression & drinking including odds ratios ranging from 1.0 to 4.2
- If depression **causes** drinking, we could imagine two possible mechanisms:
  - the onset rate of drinking could increase while depressed and/or
  - the recovery rate from drinking could decrease while depressed (as compared to when not depressed)
  - and these effectively represent **incidence rate ratios**
- The **figure** here illustrates a simple conceptual framework, where individuals may transition in/out of depression (horizontally) and drinking (vertically), but while depressed, on the right hand side, drinking transitions are modified by these two rate ratios
- Now, most epidemiological studies are cross sectional, and estimate ORs using *prevalence data*
- But prevalence reflects a **balance of incidence & duration **and here individuals can recover from both the exposure (depression) and outcome (drinking)
- So this begs the question ...
### Research Question & Objectives
- **What can ORs (actually) tell us about these onset & recovery IRRs (in this context)?**
- (Objectives)
- Our **first objective** focuses on a **reference case**, and seeks to** **characterize the relationship between ORs and**
  - IRR of drinking **onset** while depressed (red arrow)
  - IRR of drinking **recovery** while depressed (green arrow)
- Our **second objective** is then a sensitivity analysis to determine whether this relationship (focusing on the onset effect) varies by context, specifically with:
  - base rates of **depression** onset and recovery, and
  - base rates of **drinking** onset and recovery
### Methods: Individual-Based Simulation Model
- To address these objectives, we built an individual-based simulation model
- The model represents an open population, where all individuals enter the modelled population at age 10, with no history of depression or drinking, and exit at age 60
- From age 10, individuals can experience onset of **depression and/or drinking** at **base rates** derived from the literature, shown here per 100 person-years
- Finally, we can apply "effects" of depression on drinking as IRRs, namely:
  - faster drinking onset while depressed (shown in red)
  - slower drinking recovery while depressed (shown in green)
### Methods: Simulations & Odds Ratios
- Regarding simulations & calculation of ORs
- For the **reference case** (Objective 1), we analyzed the model under the default base rates described above, varying
  - the onset IRR from 1 to 8, and
  - the recovery IRR from 1 to 1/8
- For reference, in the null case (both IRR = 1) the default rates yield approximately
  - 3.8% overall prevalence of depression,
  - 5.4% overall prevalence of drinking, and
  - 0.2% overall prevalence of depressed *and* drinking
- For the **sensitivity analysis** (Objective 2), we further varied:
  - base rates of onset & recovery by factors of 0.5, 1.0, and 1.5
  - first for depression, and then for drinking
- For **calculating ORs**:
  - we used the entire population of 10,000 modelled individuals at equilibrium
  - and we report mean and 95% intervals across 41 stochastic model runs
### Result {1.a} OR can underestimate onset IRR by factor of 4+
- Now, on to results
- Here we plot the OR (y-axis) versus IRR onset (x-axis), in log-base-2 scale, with a dotted line of equality
- We can see that:
  - First, when IRR = 1, we also have OR approximately 1 so there is minimal "offset" bias around the null
  - However, as IRR increases, the OR tends to underestimate IRR by a factor of approximately 4
- For example, if we have OR of 2, we could **look-up** to find the corresponding IRR of about 5 (though the confidence interval is quite wide)
- The main **implication** of this result it that the underlying effect of drinking onset while depressed may actually be approximately 4 times what the OR suggests
  - This is perhaps surprising, since we are often warned that OR "overestimate" associations
  - But here, the underestimation derives from that balance of onset & recovery of both exposure & outcome which acts to continually **erode **the degree of association in the prevalence data
### Result {1.b} OR hardly influenced by recovery IRR
- Now considering Objective {1.b} where depression reduces the rate of drinking *recovery* while depressed
- In this case, we expect a greater OR as the IRR *decreases* below 1, so the dotted line is flipped to slope down and the x-axis ends at 1.
- We can see that OR is **hardly influenced** by the recovery IRR, so that achieving OR above approximately 1.5 is almost impossible, even when depressed individuals recover from drinking at 1/8 the rate of non-depressed individuals
- The **main implication** is that it would be very hard to actually identify the recovery effect from the OR here
  - However if we observe OR > 1.5, we likely *need* to have some onset effect
### Result {2.a} ...
- Now moving to Objective {2} -- our sensitivity analysis
- As a reminder, here we focus again on the *onset effect*
- (first slide)
- First, we vary the **depression onset rate**, shown here as lighter blue with higher onset rate
- We can see that increasing depression onset rate actually has minimal influence on the OR-IRR relationship
- (second slide)
- Next, we vary the **depression recovery rate**, again with lighter colours indicating higher recovery rate
- Here we see that the OR-IRR mismatch gets worse as recovery rates increase (shorter episodes)
- In other words, OR better approximates IRR when episodes of depression (the exposure) are *longer*
### Result {2.b} ...
- Turning next to rates of drinking onset and recovery (our outcome)
- (first slide)
- Varying the **drinking onset rate**, we find a similar result as with depression onset rate: minimal influence on the OR-IRR mismatch
- (next slide)
- However, varying the **drinking recovery rate**, we now see that higher recovery rates actually *reduce* the bias (opposite to depression)
- In other words, OR better approximates IRR when episodes of drinking (the outcome) are *shorter*
### Why do recovery rates influence OR more than onset rates?
- You may be wondering: Why do *recovery* rates seem to influence OR more than *onset* rates?
- The **short answer** is that recovery rates are generally much larger than onset rates (in this context) -- for example:
  - Rates of depression **onset** are *at most 10% per year,*
  - whereas **recovery** rates are likely *at least 50% per year*
  - So, we can pretty much always say that rates of depression recovery will be much greater than rates of onset
- **Mathematically**, we derived the analytic expression for OR, ignoring "age effects" due to all individuals entering the model *not* depressed & *not* drinking at age 10
- If all the onset rates (the blue and purple terms) are "very small" we can simplify this to the second expression
- We can then see how:
  - If depression recovery rate (b, orange) is much greater than drinking recovery rate (d, yellow): the expression tends towards 1
  - If the opposite is true (b is much smaller than d): the expression tends towards the ratio of IRRs (O / R) which we might call "unbiased"
### Summary & Conclusions
- So, to conclude:
- Our first insight is that OR > 1 in this "2x2" system can **mechanistically** derive from some combination of onset and/or recovery effects (which we conceptualize as IRRs)
- However, ORs **substantially underestimate** these IRRs
  - Specifically, we derived the expression here, from which we can see that bias towards the null comes from:
      - faster recovery from depression (our exposure)
      - slower recovery from drinking (our outcome)
  - Now, our results here have used **homogeneous rates** for all individuals, but we confirmed our findings are qualitatively unchanged when individuals have different (heterogeneous) base rates of depression and drinking onset and recovery
- The major implication of these findings is that assuming OR approximates these IRRs would underestimate the projected impact of depression interventions on drinking (and similarly in other applications)
- Finally, two limitations are that we have not considered effects of age on base rates, and we have not examined weighted sampling when calculating ORs
### References
- With that, here are our references, and ...
### Thanks
- Thank you for your attention, and thanks to our supporting institutions and funders
- We have code on github if folks are interested at the link here
- Over to the moderators for questions ...
