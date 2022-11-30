Psyc 575 Final Analyses
================
Chelsey Pan + Mengzhao Yan
2022-12-01

``` r
# Load packages
library(tidyverse)
library(haven)
library(here)
library(performance)
library(lme4)
library(lmerTest)
library(modelsummary)
library(brms)
library(sjPlot)
theme_set(theme_classic() +
    theme(panel.grid.major.y = element_line(color = "grey92")))

# Load data
df <- read_dta('psyc575finalupdated.dta')
```

Reformat data wide to long

``` r
# Reformat from wide to long
df_long <- df %>%
  select(-c(`_merge1`, `_merge2`, wave1:wave3)) %>%
  pivot_longer(c(-ID),
               names_to = c('.value', 'wave'),
               names_pattern = "(age|female|race|marital|edulevel|deptot|adltotb|socsuptot|weight_adj|weight_sel)([1-3])",
               names_transform = list(wave = as.integer)) %>%
  mutate(wave = wave - 1) # Convert wave from 1-3 to 0-2
```

## Variable Summary

-   `ID`: participant id

**Demographics**

-   `age`: mean age1 = 69, mean age2 = 73, mean age3 = 68 Since NSHAP
    recruits new research participants in each round, the mean age does
    not increase across time points.
-   `race`: participant race. 1= non-Hispanic White, 2 = non-Hispanic,
    Black, 3 = Hispanic, 4 = Other
-   `female`: participant gender, binary-coded. 1 = Female, 0 = Male
-   `marital`: marital status. 1= Married, 0 = Unmarried
-   `educlevel`: education level. 1 = \>12 years, 0 = \<or=12 years

**Scales**

-   `deptot`: aggregate depression score on an 11-item scale. Each
    question is scored from 1 (rarely or none of the time) - 4 (most of
    the time).
-   `adltotb`: total number of reported functional/physical difficulties
    out of 7 questions about activities of daily living.
-   `socsuptot`: aggregate social support score on a 6-item scale,
    scored from 1 (hardly ever or never) to 3 (often)

Data exploration

``` r
# Spread of predictors
df_long %>%
    select(wave, age, female, race, deptot, adltotb, socsuptot) %>%
    psych::pairs.panels(jiggle = TRUE, factor = 0.5, ellipses = FALSE,
                        cex.cor = 1, cex = 0.5)
```

![](final_files/figure-gfm/exploration-1.png)<!-- -->

``` r
# Analysis of attrition
# Add complete/incomplete variable
df_comp <- df %>%
    # Compute summaries by rows
    rowwise() %>%
    # First compute the number of missing occasions
    mutate(nmis_deptot = sum(is.na(c_across(deptot1:deptot3))),
           # Complete only when nmis_read = 0
           complete = if_else(nmis_deptot == 0, "complete", "incomplete")) %>%
    ungroup()
# Compare the differences
datasummary((deptot1 + socsuptot1 + adltotb1 + age1 + female1 + race1 + marital1 + edulevel1) ~
              complete * (Mean + SD), data = df_comp)
```

|            | complete / Mean | complete / SD | incomplete / Mean | incomplete / SD |
|:-----------|----------------:|--------------:|------------------:|----------------:|
| deptot1    |           15.84 |          4.95 |             17.29 |            5.44 |
| socsuptot1 |           12.91 |          3.70 |             11.29 |            4.35 |
| adltotb1   |            0.61 |          1.31 |              1.27 |            1.87 |
| age1       |           66.69 |          6.81 |             72.10 |            7.93 |
| female1    |            0.54 |          0.50 |              0.49 |            0.50 |
| race1      |            1.44 |          0.77 |              1.45 |            0.77 |
| marital1   |            0.67 |          0.47 |              0.52 |            0.50 |
| edulevel1  |            0.58 |          0.49 |              0.43 |            0.49 |

It looks like there is a noticeable difference in the mean age of
participants at wave 1 who didn’t complete all 3 waves. It’s possible
that some of the participants didn’t complete all 3 waves due to passing
away, however this is hard to ascertain given that not all participants
in the study began at wave 1. The other prominent difference was a
higher average amount of physical difficulties among the participants
who didn’t complete all 3 waves.

ICC

``` r
# ICC for outcome - depression symptoms score
m0_dep <- lmer(deptot ~ (1 | ID), data = df_long)
performance::icc(m0_dep)
```

    ## # Intraclass Correlation Coefficient
    ## 
    ##     Adjusted ICC: 0.570
    ##   Unadjusted ICC: 0.570

``` r
# ICC for predictor 1: physical health difficulties 
m0_adl <- lmer(adltotb ~ (1 | ID), data = df_long)
performance::icc(m0_adl)
```

    ## # Intraclass Correlation Coefficient
    ## 
    ##     Adjusted ICC: 0.546
    ##   Unadjusted ICC: 0.546

``` r
# ICC for predictor 2: social support score
m0_socsup <- lmer(socsuptot ~ (1 | ID), data = df_long)
performance::icc(m0_socsup)
```

    ## # Intraclass Correlation Coefficient
    ## 
    ##     Adjusted ICC: 0.522
    ##   Unadjusted ICC: 0.522

Variability at the individual level accounts for about 57% of the total
variability in total depression symptom scores. The ICC for the physical
health score predictor was 0.546, while the ICC for the social support
predictor was 0.522.

Bonus: ICC of outcome using brms

``` r
# Run unconditional model predicting depression scores
m0 <- brm(deptot ~ (1 | ID), data = df_long,
          seed = 123,
          file = 'dep_icc')

# Get summary
summary(m0)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (1 | ID) 
    ##    Data: df_long (Number of observations: 10578) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6069) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     4.05      0.06     3.92     4.17 1.00      808     1656
    ## 
    ## Population-Level Effects: 
    ##           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept    16.70      0.07    16.57    16.83 1.00     1801     2578
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.52      0.04     3.45     3.59 1.00     1088     2075
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# Obtain ICC
performance::icc(m0)
```

    ## # Intraclass Correlation Coefficient
    ## 
    ##     Adjusted ICC: 0.569
    ##   Unadjusted ICC: 0.569

``` r
# Get conditional ICC
m_cs <- brm(deptot ~ 0 + factor(wave) + (1 | ID),
    data = df_long,
    seed = 123,
    file = 'dep_icc_cs')

# Get summary
summary(m_cs)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ 0 + factor(wave) + (1 | ID) 
    ##    Data: df_long (Number of observations: 10578) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6069) 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sd(Intercept)     4.04      0.06     3.92     4.17 1.00      828     1499
    ## 
    ## Population-Level Effects: 
    ##             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## factorwave0    16.58      0.09    16.40    16.76 1.00     2453     2572
    ## factorwave1    16.30      0.09    16.13    16.47 1.00     2237     2418
    ## factorwave2    17.01      0.08    16.86    17.16 1.00     2555     3125
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.50      0.04     3.43     3.58 1.00     1075     2653
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

``` r
# Obtain ICC
performance::icc(m_cs)
```

    ## # Intraclass Correlation Coefficient
    ## 
    ##     Adjusted ICC: 0.571
    ##   Unadjusted ICC: 0.569

Using brms, you get about the same ICC.

Separate time-varying predictors into within-person and between-person
levels

``` r
df_long <- df_long %>%
  group_by(ID) %>%
  mutate(across(c(adltotb, socsuptot, age),
         list("pm" = ~ mean(., na.rm = TRUE),
                       "pmc" = ~ . - mean(., na.rm = TRUE)))) %>%
  mutate(race = as.factor(race)) %>%
  mutate(age_w0 = first(age))
```

**Model Equations**

**Lvl 1:**

$$
\text{deptot}_{ti} = \beta_{0i} + \beta_{1i} adltotbpmc_{ti} + e_{ti}
$$

**Lvl 2:**

$$
\beta_{0i} = \gamma_{00} + \gamma_{01} adltotbpm_i + \gamma_{02} race_i + \gamma_{03} adltotbpm_i \times
 race_i + u_{0i}
$$

$$\beta_{1i} = \gamma_{10} + \gamma_{11} race_i + u_{1i}$$

## Preliminary analysis

Base model: Is there an association between reported physical
difficulties and depressive symptoms across the 3 waves, and does this
interact with participant race (treated as a lvl 2 predictor in this
case)

``` r
# Model with just functional health difficulties
m1_prelim <- brm(deptot ~ (adltotb_pm + adltotb_pmc) * race + (adltotb_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_race')

# Get model summary
summary(m1_prelim)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (adltotb_pm + adltotb_pmc) * race + (adltotb_pmc | ID) 
    ##    Data: df_long (Number of observations: 10543) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6052) 
    ##                            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## sd(Intercept)                  3.67      0.06     3.56     3.79 1.00      922
    ## sd(adltotb_pmc)                0.93      0.08     0.77     1.08 1.00     1084
    ## cor(Intercept,adltotb_pmc)     0.25      0.07     0.11     0.39 1.00     2097
    ##                            Tail_ESS
    ## sd(Intercept)                  1912
    ## sd(adltotb_pmc)                1960
    ## cor(Intercept,adltotb_pmc)     2119
    ## 
    ## Population-Level Effects: 
    ##                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept            15.38      0.08    15.22    15.54 1.00     1655     2357
    ## adltotb_pm            1.35      0.05     1.25     1.45 1.00     1625     2155
    ## adltotb_pmc           0.60      0.07     0.45     0.74 1.00     2527     2834
    ## race3                 0.70      0.22     0.27     1.11 1.00     1618     1894
    ## race2                 0.61      0.19     0.23     0.99 1.00     1452     2261
    ## race4                 0.42      0.37    -0.32     1.13 1.00     1940     2664
    ## adltotb_pm:race3     -0.16      0.12    -0.39     0.07 1.00     1823     2185
    ## adltotb_pm:race2     -0.28      0.10    -0.47    -0.08 1.00     1658     2344
    ## adltotb_pm:race4     -0.15      0.24    -0.62     0.31 1.00     2301     2433
    ## adltotb_pmc:race3    -0.20      0.19    -0.57     0.16 1.00     2564     2715
    ## adltotb_pmc:race2    -0.16      0.15    -0.46     0.13 1.00     2812     2871
    ## adltotb_pmc:race4     0.21      0.46    -0.68     1.11 1.00     4627     3162
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.34      0.04     3.26     3.42 1.00      821     2074
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Table summarizing model results

``` r
# Summarize the model results
msummary(m1_prelim,
         estimate = c("{estimate} [{conf.low}, {conf.high}]"),
         statistic = NULL,  # suppress the extra rows for SEs
    shape = effect + term ~ model,
    title = 'Table 1: Model coefficients')
```

|        |                                    |          Model 1          |
|:-------|:-----------------------------------|:-------------------------:|
| fixed  | b_Intercept                        | 15.377 \[15.217, 15.539\] |
|        | b_adltotb_pm                       |  1.352 \[1.249, 1.448\]   |
|        | b_adltotb_pmc                      |  0.595 \[0.453, 0.741\]   |
|        | b_race3                            |  0.696 \[0.268, 1.107\]   |
|        | b_race2                            |  0.607 \[0.231, 0.991\]   |
|        | b_race4                            |  0.416 \[-0.317, 1.135\]  |
|        | b_adltotb_pm × race3               | -0.156 \[-0.388, 0.073\]  |
|        | b_adltotb_pm × race2               | -0.277 \[-0.470, -0.081\] |
|        | b_adltotb_pm × race4               | -0.141 \[-0.615, 0.315\]  |
|        | b_adltotb_pmc × race3              | -0.207 \[-0.573, 0.162\]  |
|        | b_adltotb_pmc × race2              | -0.159 \[-0.462, 0.127\]  |
|        | b_adltotb_pmc × race4              |  0.213 \[-0.677, 1.110\]  |
|        | sigma                              |  3.339 \[3.265, 3.419\]   |
| random | sd_ID\_\_Intercept                 |  3.673 \[3.558, 3.788\]   |
|        | sd_ID\_\_adltotb_pmc               |  0.928 \[0.771, 1.085\]   |
|        | cor_ID\_\_Intercept\_\_adltotb_pmc |  0.251 \[0.111, 0.388\]   |
|        | Num.Obs.                           |           10543           |
|        | R2                                 |           0.600           |
|        | R2 Adj.                            |           0.350           |
|        | R2 Marg.                           |           0.128           |
|        | ELPD                               |         -30337.4          |
|        | ELPD s.e.                          |           97.8            |
|        | LOOIC                              |          60674.8          |
|        | LOOIC s.e.                         |           195.6           |
|        | WAIC                               |          59612.3          |
|        | RMSE                               |           2.60            |
|        | r2.adjusted.marginal               |           0.123           |

Table 1: Model coefficients

Figures showing association between main predictor and outcome

``` r
# Figure showing association between adltotb_pm and depression symptoms,
# split by race
figures[4]
```

    ## $`adltotb_pm:race`

![](final_files/figure-gfm/figures-1.png)<!-- -->

``` r
# Figure showing association between adltotb_pmc and depression symptoms
# split by race
figures[5]
```

    ## $`adltotb_pmc:race`

![](final_files/figure-gfm/figures-2.png)<!-- -->

Writeup:

We observed evidence for a positive association between depression
symptoms and physical difficulties. The 68% plausible range showed that
for most white non-Hispanic participants, a 1 unit increase in reported
physical difficulties was associated with \[0.53, 0.67\] increase in
predicted depression symptoms. We also see a significant positive
association across white non-Hispanic participants as well, with one
unit increase in reported physical difficulties associated with \[1.3,
1.4\] increase in predicted depression symptoms. We also observe
evidence for a positive association between race and depression
symptoms, except among participants who reported their race as “other.”
Based on the 95% CIs, we observe a evidence for an interaction only
between reported physical difficulties and participant race among
non-Hispanic Black participants, where reported physical difficulties
negatively interacts with participant race (68% plausible range \[-0.38,
-0.18\]).

## Followup analyses

Version 2 of prelim: rerun with age added as a lvl 2 covariate

``` r
# Model with functional health difficulties and age
m1a_prelim <- brm(deptot ~ (adltotb_pm + adltotb_pmc) * race + age_w0 + (adltotb_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_race_agel2')

# Get model summary
summary(m1a_prelim)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (adltotb_pm + adltotb_pmc) * race + age + (adltotb_pmc | ID) 
    ##    Data: df_long (Number of observations: 10543) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6052) 
    ##                            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## sd(Intercept)                  3.63      0.06     3.52     3.75 1.01      730
    ## sd(adltotb_pmc)                0.91      0.08     0.75     1.07 1.00     1187
    ## cor(Intercept,adltotb_pmc)     0.25      0.07     0.11     0.39 1.00     2170
    ##                            Tail_ESS
    ## sd(Intercept)                  1321
    ## sd(adltotb_pmc)                2186
    ## cor(Intercept,adltotb_pmc)     2107
    ## 
    ## Population-Level Effects: 
    ##                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept            17.58      0.38    16.83    18.35 1.00     1927     2697
    ## adltotb_pm            1.38      0.05     1.29     1.48 1.00     1470     2170
    ## adltotb_pmc           0.63      0.07     0.48     0.77 1.00     2619     2663
    ## race3                 0.62      0.21     0.21     1.03 1.00     1620     2380
    ## race2                 0.56      0.19     0.18     0.93 1.00     1275     2310
    ## race4                 0.31      0.37    -0.40     1.06 1.00     1593     2244
    ## age                  -0.03      0.01    -0.04    -0.02 1.00     1818     2525
    ## adltotb_pm:race3     -0.15      0.11    -0.37     0.07 1.00     1570     2021
    ## adltotb_pm:race2     -0.29      0.10    -0.48    -0.10 1.01     1296     1937
    ## adltotb_pm:race4     -0.16      0.24    -0.62     0.31 1.00     1832     2787
    ## adltotb_pmc:race3    -0.23      0.18    -0.57     0.11 1.00     2707     3159
    ## adltotb_pmc:race2    -0.18      0.15    -0.48     0.11 1.00     2793     3161
    ## adltotb_pmc:race4     0.19      0.44    -0.68     1.05 1.00     4166     3274
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.35      0.04     3.28     3.44 1.01      695     1132
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Model 1: looking at the relationship between physical difficulties and
depression without race

``` r
# Physical difficulties and age
m1 <- brm(deptot ~ adltotb_pm + adltotb_pmc + age_pm + age_pmc + (adltotb_pmc + age_pmc| ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_age')

# Get summary
summary(m1)
```

    ## Warning: There were 4 divergent transitions after warmup. Increasing adapt_delta
    ## above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-
    ## transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ adltotb_pm + adltotb_pmc + age_pm + age_pmc + (adltotb_pmc + age_pmc | ID) 
    ##    Data: df_long (Number of observations: 10578) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6069) 
    ##                            Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## sd(Intercept)                  3.65      0.06     3.54     3.77 1.00      857
    ## sd(adltotb_pmc)                0.93      0.08     0.77     1.09 1.01      923
    ## sd(age_pmc)                    0.05      0.03     0.00     0.12 1.01      253
    ## cor(Intercept,adltotb_pmc)     0.22      0.07     0.08     0.35 1.00     1946
    ## cor(Intercept,age_pmc)         0.33      0.33    -0.49     0.90 1.00     1788
    ## cor(adltotb_pmc,age_pmc)      -0.14      0.42    -0.82     0.75 1.01     1466
    ##                            Tail_ESS
    ## sd(Intercept)                  1544
    ## sd(adltotb_pmc)                1681
    ## sd(age_pmc)                     321
    ## cor(Intercept,adltotb_pmc)     2297
    ## cor(Intercept,age_pmc)         1665
    ## cor(adltotb_pmc,age_pmc)       1999
    ## 
    ## Population-Level Effects: 
    ##             Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept      19.31      0.43    18.46    20.12 1.00     1562     1966
    ## adltotb_pm      1.33      0.04     1.25     1.41 1.00     1553     2330
    ## adltotb_pmc     0.49      0.06     0.37     0.61 1.00     3185     3010
    ## age_pm         -0.06      0.01    -0.07    -0.04 1.00     1400     1873
    ## age_pmc         0.04      0.01     0.01     0.06 1.00     4985     2872
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.32      0.04     3.24     3.40 1.00      559     1177
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Model 2: predicting reported depression symptoms from both physical
difficulties and social support

``` r
# Physical difficulties, social support, and age
m2 <- brm(deptot ~ adltotb_pm + adltotb_pmc + age_pm + age_pmc + socsuptot_pm + socsuptot_pmc + (adltotb_pmc + age_pmc + socsuptot_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_socsup')

# Get summary
summary(m2)
```

    ## Warning: Parts of the model have not converged (some Rhats are > 1.05). Be
    ## careful when analysing the results! We recommend running more iterations and/or
    ## setting stronger priors.

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ adltotb_pm + adltotb_pmc + age_pm + age_pmc + socsuptot_pm + socsuptot_pmc + (adltotb_pmc + age_pmc + socsuptot_pmc | ID) 
    ##    Data: df_long (Number of observations: 10578) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6069) 
    ##                                Estimate Est.Error l-95% CI u-95% CI Rhat
    ## sd(Intercept)                      3.55      0.06     3.43     3.66 1.01
    ## sd(adltotb_pmc)                    0.88      0.08     0.72     1.04 1.01
    ## sd(age_pmc)                        0.04      0.03     0.00     0.10 1.02
    ## sd(socsuptot_pmc)                  0.30      0.04     0.23     0.37 1.01
    ## cor(Intercept,adltotb_pmc)         0.23      0.07     0.08     0.37 1.00
    ## cor(Intercept,age_pmc)             0.17      0.34    -0.60     0.80 1.00
    ## cor(adltotb_pmc,age_pmc)          -0.14      0.41    -0.84     0.74 1.00
    ## cor(Intercept,socsuptot_pmc)      -0.22      0.08    -0.38    -0.05 1.00
    ## cor(adltotb_pmc,socsuptot_pmc)    -0.07      0.17    -0.41     0.26 1.01
    ## cor(age_pmc,socsuptot_pmc)        -0.15      0.41    -0.86     0.73 1.08
    ##                                Bulk_ESS Tail_ESS
    ## sd(Intercept)                       838     1678
    ## sd(adltotb_pmc)                    1080     2238
    ## sd(age_pmc)                         322      669
    ## sd(socsuptot_pmc)                   708     1706
    ## cor(Intercept,adltotb_pmc)         2895     3069
    ## cor(Intercept,age_pmc)             5200     2432
    ## cor(adltotb_pmc,age_pmc)           2528     2463
    ## cor(Intercept,socsuptot_pmc)       3078     2851
    ## cor(adltotb_pmc,socsuptot_pmc)      619     1443
    ## cor(age_pmc,socsuptot_pmc)           39      175
    ## 
    ## Population-Level Effects: 
    ##               Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept        22.66      0.46    21.78    23.55 1.00     2974     2986
    ## adltotb_pm        1.21      0.04     1.14     1.29 1.00     2803     2934
    ## adltotb_pmc       0.49      0.06     0.37     0.60 1.00     5256     3389
    ## age_pm           -0.05      0.01    -0.07    -0.04 1.00     2898     2890
    ## age_pmc           0.02      0.01    -0.00     0.04 1.00     7886     3369
    ## socsuptot_pm     -0.28      0.01    -0.31    -0.26 1.00     3088     3018
    ## socsuptot_pmc    -0.09      0.02    -0.13    -0.05 1.00     6204     3773
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.21      0.04     3.13     3.30 1.01      625     1187
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Model 3: predicting reported depression symptoms from physical
difficulties and social support, including an interaction between the
two

``` r
# Physical difficulties x social support, and age
m3 <- brm(deptot ~ (adltotb_pm + adltotb_pmc) * (socsuptot_pm + socsuptot_pmc) + age_pm + age_pmc  + (adltotb_pmc + age_pmc + socsuptot_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_socsup_int')

# Get summary
summary(m3)
```

    ## Warning: Parts of the model have not converged (some Rhats are > 1.05). Be
    ## careful when analysing the results! We recommend running more iterations and/or
    ## setting stronger priors.

    ## Warning: There were 1 divergent transitions after warmup. Increasing adapt_delta
    ## above 0.8 may help. See http://mc-stan.org/misc/warnings.html#divergent-
    ## transitions-after-warmup

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (adltotb_pm + adltotb_pmc) * (socsuptot_pm + socsuptot_pmc) + age_pm + age_pmc + (adltotb_pmc + age_pmc + socsuptot_pmc | ID) 
    ##    Data: df_long (Number of observations: 10578) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 6069) 
    ##                                Estimate Est.Error l-95% CI u-95% CI Rhat
    ## sd(Intercept)                      3.54      0.06     3.43     3.66 1.01
    ## sd(adltotb_pmc)                    0.87      0.08     0.70     1.03 1.00
    ## sd(age_pmc)                        0.04      0.03     0.00     0.11 1.01
    ## sd(socsuptot_pmc)                  0.30      0.04     0.23     0.38 1.00
    ## cor(Intercept,adltotb_pmc)         0.24      0.07     0.09     0.38 1.00
    ## cor(Intercept,age_pmc)             0.17      0.34    -0.58     0.79 1.00
    ## cor(adltotb_pmc,age_pmc)          -0.16      0.41    -0.82     0.70 1.01
    ## cor(Intercept,socsuptot_pmc)      -0.21      0.09    -0.39    -0.05 1.00
    ## cor(adltotb_pmc,socsuptot_pmc)    -0.02      0.17    -0.38     0.30 1.01
    ## cor(age_pmc,socsuptot_pmc)        -0.20      0.38    -0.86     0.58 1.06
    ##                                Bulk_ESS Tail_ESS
    ## sd(Intercept)                       406     1040
    ## sd(adltotb_pmc)                     878     1431
    ## sd(age_pmc)                         358      342
    ## sd(socsuptot_pmc)                   520      688
    ## cor(Intercept,adltotb_pmc)         1919     2171
    ## cor(Intercept,age_pmc)             2811     2194
    ## cor(adltotb_pmc,age_pmc)           1260     2086
    ## cor(Intercept,socsuptot_pmc)       1410     1046
    ## cor(adltotb_pmc,socsuptot_pmc)      298      704
    ## cor(age_pmc,socsuptot_pmc)           50      115
    ## 
    ## Population-Level Effects: 
    ##                           Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS
    ## Intercept                    22.50      0.47    21.57    23.43 1.00     1508
    ## adltotb_pm                    1.36      0.10     1.16     1.56 1.00     1101
    ## adltotb_pmc                   0.09      0.22    -0.32     0.52 1.00     1468
    ## socsuptot_pm                 -0.27      0.02    -0.30    -0.24 1.00     1201
    ## socsuptot_pmc                -0.08      0.02    -0.13    -0.03 1.00     2804
    ## age_pm                       -0.05      0.01    -0.07    -0.04 1.00     1377
    ## age_pmc                       0.02      0.01    -0.00     0.04 1.00     3943
    ## adltotb_pm:socsuptot_pm      -0.01      0.01    -0.03     0.00 1.00     1199
    ## adltotb_pm:socsuptot_pmc     -0.01      0.02    -0.04     0.02 1.00     2731
    ## adltotb_pmc:socsuptot_pm      0.03      0.02    -0.00     0.07 1.00     1535
    ## adltotb_pmc:socsuptot_pmc    -0.05      0.03    -0.11     0.00 1.00     1610
    ##                           Tail_ESS
    ## Intercept                     2099
    ## adltotb_pm                    1909
    ## adltotb_pmc                   2173
    ## socsuptot_pm                  2164
    ## socsuptot_pmc                 2455
    ## age_pm                        2265
    ## age_pmc                       2742
    ## adltotb_pm:socsuptot_pm       2067
    ## adltotb_pm:socsuptot_pmc      2714
    ## adltotb_pmc:socsuptot_pm      2436
    ## adltotb_pmc:socsuptot_pmc     1888
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.22      0.05     3.12     3.31 1.01      339      745
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

**Interacting each variable of interest with race**

We ran into convergence issues with models 1-3, which may be related to
the fact that age as a level 1 variable is confounded with time. In
these following analyses, age is instead treated as a level 2 variable,
with age_w0 equal to initial age at the first wave.

Model 4a: Looking at interaction between physical difficulties and race,
with social support and age (lvl2) included

``` r
m4a <- brm(deptot ~ (adltotb_pm + adltotb_pmc) * race + (socsuptot_pm + socsuptot_pmc) + age_w0 + (adltotb_pmc + socsuptot_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_race_int_all')

# Get summary
summary(m4a)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (adltotb_pm + adltotb_pmc) * race + (socsuptot_pm + socsuptot_pmc) + age_w0 + (adltotb_pmc + socsuptot_pmc | ID) 
    ##    Data: df_long (Number of observations: 6829) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 2993) 
    ##                                Estimate Est.Error l-95% CI u-95% CI Rhat
    ## sd(Intercept)                      3.33      0.07     3.20     3.47 1.00
    ## sd(adltotb_pmc)                    0.87      0.08     0.71     1.04 1.00
    ## sd(socsuptot_pmc)                  0.33      0.04     0.26     0.40 1.00
    ## cor(Intercept,adltotb_pmc)         0.21      0.08     0.07     0.37 1.00
    ## cor(Intercept,socsuptot_pmc)      -0.21      0.08    -0.37    -0.06 1.00
    ## cor(adltotb_pmc,socsuptot_pmc)     0.02      0.17    -0.34     0.35 1.01
    ##                                Bulk_ESS Tail_ESS
    ## sd(Intercept)                      1072     2160
    ## sd(adltotb_pmc)                    1050     1968
    ## sd(socsuptot_pmc)                   649     1655
    ## cor(Intercept,adltotb_pmc)         2262     2588
    ## cor(Intercept,socsuptot_pmc)       2694     2543
    ## cor(adltotb_pmc,socsuptot_pmc)      436      795
    ## 
    ## Population-Level Effects: 
    ##                   Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept            20.11      0.80    18.58    21.70 1.00     1525     2548
    ## adltotb_pm            1.18      0.06     1.06     1.31 1.00     1688     2409
    ## adltotb_pmc           0.51      0.08     0.36     0.67 1.00     2764     2755
    ## race3                 0.43      0.29    -0.13     1.00 1.00     1417     2218
    ## race2                 0.23      0.25    -0.28     0.71 1.00     1531     2811
    ## race4                 0.69      0.57    -0.44     1.85 1.00     1655     2350
    ## socsuptot_pm         -0.29      0.02    -0.34    -0.25 1.00     1806     2435
    ## socsuptot_pmc        -0.10      0.02    -0.14    -0.05 1.00     4278     3253
    ## age_w0               -0.02      0.01    -0.04    -0.00 1.00     1633     2094
    ## adltotb_pm:race3     -0.03      0.15    -0.32     0.26 1.00     1620     2116
    ## adltotb_pm:race2     -0.21      0.12    -0.44     0.03 1.00     1667     2286
    ## adltotb_pm:race4     -0.48      0.33    -1.11     0.16 1.00     1943     2658
    ## adltotb_pmc:race3    -0.16      0.20    -0.55     0.22 1.00     3108     2986
    ## adltotb_pmc:race2    -0.16      0.16    -0.48     0.15 1.00     2683     3318
    ## adltotb_pmc:race4     0.28      0.47    -0.64     1.20 1.00     5134     2956
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.19      0.05     3.10     3.29 1.00      796     1664
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Model 4b: Looking at interaction between social support and race, with
physical difficulties and age (lvl2) included

``` r
m4b <- brm(deptot ~ (adltotb_pm + adltotb_pmc) + (socsuptot_pm + socsuptot_pmc) * race + age_w0 + (adltotb_pmc + socsuptot_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_socsup_race_int_all')

# Get summary
summary(m4b)
```

    ##  Family: gaussian 
    ##   Links: mu = identity; sigma = identity 
    ## Formula: deptot ~ (adltotb_pm + adltotb_pmc) + (socsuptot_pm + socsuptot_pmc) * race + age_w0 + (adltotb_pmc + socsuptot_pmc | ID) 
    ##    Data: df_long (Number of observations: 6829) 
    ##   Draws: 4 chains, each with iter = 2000; warmup = 1000; thin = 1;
    ##          total post-warmup draws = 4000
    ## 
    ## Group-Level Effects: 
    ## ~ID (Number of levels: 2993) 
    ##                                Estimate Est.Error l-95% CI u-95% CI Rhat
    ## sd(Intercept)                      3.33      0.07     3.19     3.46 1.00
    ## sd(adltotb_pmc)                    0.87      0.08     0.70     1.04 1.00
    ## sd(socsuptot_pmc)                  0.33      0.04     0.26     0.41 1.00
    ## cor(Intercept,adltotb_pmc)         0.22      0.08     0.07     0.37 1.00
    ## cor(Intercept,socsuptot_pmc)      -0.21      0.08    -0.36    -0.05 1.00
    ## cor(adltotb_pmc,socsuptot_pmc)     0.01      0.17    -0.35     0.32 1.01
    ##                                Bulk_ESS Tail_ESS
    ## sd(Intercept)                       991     2077
    ## sd(adltotb_pmc)                    1080     2005
    ## sd(socsuptot_pmc)                   869     1758
    ## cor(Intercept,adltotb_pmc)         2901     2880
    ## cor(Intercept,socsuptot_pmc)       2782     3337
    ## cor(adltotb_pmc,socsuptot_pmc)      579     1102
    ## 
    ## Population-Level Effects: 
    ##                     Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## Intercept              20.24      0.85    18.55    21.92 1.00     2829     3042
    ## adltotb_pm              1.13      0.05     1.04     1.23 1.00     2698     3014
    ## adltotb_pmc             0.47      0.06     0.34     0.59 1.00     5336     3275
    ## socsuptot_pm           -0.31      0.03    -0.36    -0.26 1.00     2982     3256
    ## socsuptot_pmc          -0.12      0.03    -0.17    -0.07 1.00     5412     3414
    ## race3                   0.68      0.83    -0.94     2.29 1.00     2466     2852
    ## race2                  -1.54      0.66    -2.86    -0.24 1.00     2529     2733
    ## race4                   3.91      1.79     0.46     7.42 1.00     2658     2780
    ## age_w0                 -0.02      0.01    -0.04     0.00 1.00     2906     2753
    ## socsuptot_pm:race3     -0.02      0.07    -0.17     0.11 1.00     2476     2856
    ## socsuptot_pm:race2      0.14      0.06     0.03     0.25 1.00     2568     2758
    ## socsuptot_pm:race4     -0.32      0.15    -0.61    -0.03 1.00     2630     3018
    ## socsuptot_pmc:race3     0.02      0.07    -0.12     0.16 1.00     6065     3465
    ## socsuptot_pmc:race2     0.10      0.06    -0.01     0.21 1.00     6091     3610
    ## socsuptot_pmc:race4     0.11      0.13    -0.15     0.37 1.00     6077     3061
    ## 
    ## Family Specific Parameters: 
    ##       Estimate Est.Error l-95% CI u-95% CI Rhat Bulk_ESS Tail_ESS
    ## sigma     3.19      0.05     3.10     3.29 1.01      793     1713
    ## 
    ## Draws were sampled using sampling(NUTS). For each parameter, Bulk_ESS
    ## and Tail_ESS are effective sample size measures, and Rhat is the potential
    ## scale reduction factor on split chains (at convergence, Rhat = 1).

Table summarizing models 1-3

``` r
# Summarize the model results
msummary(list(FuncDiff = m1, 
              FD_SocSup = m2, 
              FDxSocSup = m3),
         estimate = c("{estimate} [{conf.low}, {conf.high}]"),
         statistic = NULL,  # suppress the extra rows for SEs
    shape = effect + term ~ model,
    title = 'Table 1: Model coefficients')
```

|        |                                        |         FuncDiff          |         FD_SocSup         |         FDxSocSup         |
|:-------|:---------------------------------------|:-------------------------:|:-------------------------:|:-------------------------:|
| fixed  | b_Intercept                            | 19.321 \[18.456, 20.123\] | 22.664 \[21.782, 23.554\] | 22.504 \[21.573, 23.426\] |
|        | b_adltotb_pm                           |  1.328 \[1.246, 1.406\]   |  1.213 \[1.141, 1.288\]   |  1.354 \[1.156, 1.561\]   |
|        | b_adltotb_pmc                          |  0.492 \[0.370, 0.610\]   |  0.486 \[0.369, 0.601\]   |  0.093 \[-0.318, 0.524\]  |
|        | b_age_pm                               | -0.055 \[-0.067, -0.043\] | -0.054 \[-0.066, -0.042\] | -0.054 \[-0.065, -0.042\] |
|        | b_age_pmc                              |  0.036 \[0.015, 0.058\]   |  0.021 \[-0.002, 0.043\]  |  0.021 \[-0.002, 0.043\]  |
|        | sigma                                  |  3.324 \[3.242, 3.400\]   |  3.215 \[3.127, 3.300\]   |  3.217 \[3.124, 3.308\]   |
|        | b_socsuptot_pm                         |                           | -0.283 \[-0.310, -0.255\] | -0.270 \[-0.300, -0.237\] |
|        | b_socsuptot_pmc                        |                           | -0.088 \[-0.130, -0.047\] | -0.079 \[-0.128, -0.031\] |
|        | b_adltotb_pm × socsuptot_pm            |                           |                           | -0.014 \[-0.031, 0.004\]  |
|        | b_adltotb_pm × socsuptot_pmc           |                           |                           | -0.014 \[-0.044, 0.016\]  |
|        | b_adltotb_pmc × socsuptot_pm           |                           |                           |  0.033 \[-0.002, 0.068\]  |
|        | b_adltotb_pmc × socsuptot_pmc          |                           |                           | -0.054 \[-0.113, 0.005\]  |
| random | sd_ID\_\_Intercept                     |  3.655 \[3.539, 3.771\]   |  3.552 \[3.434, 3.664\]   |  3.543 \[3.426, 3.663\]   |
|        | sd_ID\_\_adltotb_pmc                   |  0.933 \[0.766, 1.094\]   |  0.882 \[0.724, 1.043\]   |  0.871 \[0.700, 1.032\]   |
|        | sd_ID\_\_age_pmc                       |  0.040 \[0.003, 0.120\]   |  0.031 \[0.001, 0.104\]   |  0.030 \[0.002, 0.105\]   |
|        | cor_ID\_\_Intercept\_\_adltotb_pmc     |  0.217 \[0.081, 0.355\]   |  0.230 \[0.081, 0.374\]   |  0.237 \[0.093, 0.381\]   |
|        | cor_ID\_\_Intercept\_\_age_pmc         |  0.339 \[-0.489, 0.897\]  |  0.182 \[-0.596, 0.802\]  |  0.183 \[-0.581, 0.789\]  |
|        | cor_ID\_\_adltotb_pmc\_\_age_pmc       | -0.179 \[-0.824, 0.748\]  | -0.187 \[-0.835, 0.740\]  | -0.208 \[-0.823, 0.698\]  |
|        | sd_ID\_\_socsuptot_pmc                 |                           |  0.302 \[0.228, 0.371\]   |  0.305 \[0.229, 0.375\]   |
|        | cor_ID\_\_Intercept\_\_socsuptot_pmc   |                           | -0.216 \[-0.379, -0.054\] | -0.211 \[-0.387, -0.051\] |
|        | cor_ID\_\_adltotb_pmc\_\_socsuptot_pmc |                           | -0.069 \[-0.407, 0.257\]  | -0.013 \[-0.378, 0.298\]  |
|        | cor_ID\_\_age_pmc\_\_socsuptot_pmc     |                           | -0.176 \[-0.858, 0.727\]  | -0.211 \[-0.864, 0.583\]  |
|        | Num.Obs.                               |           10578           |           10578           |           10578           |
|        | R2                                     |           0.603           |           0.629           |           0.628           |
|        | R2 Adj.                                |           0.353           |           0.380           |           0.379           |
|        | R2 Marg.                               |           0.135           |           0.175           |           0.175           |
|        | ELPD                                   |         -30400.8          |         -30205.6          |         -30214.2          |
|        | ELPD s.e.                              |           97.7            |           97.2            |           97.4            |
|        | LOOIC                                  |          60801.6          |          60411.2          |          60428.5          |
|        | LOOIC s.e.                             |           195.4           |           194.4           |           194.8           |
|        | WAIC                                   |          59721.9          |          59215.0          |          59233.1          |
|        | RMSE                                   |           2.59            |           2.46            |           2.46            |
|        | r2.adjusted.marginal                   |           0.133           |           0.175           |           0.174           |

Table 1: Model coefficients

Table summarizing models 4a and 4b

``` r
# Coefficients for models with variables interacting with race
msummary(list(fd_race = m4a,
              ss_race = m4b),
         estimate = c("{estimate} [{conf.low}, {conf.high}]"),
         statistic = NULL,  # suppress the extra rows for SEs
    shape = effect + term ~ model,
    title = 'Table 2: Model coefficients for race interactions')
```

|        |                                        |          fd_race          |          ss_race          |
|:-------|:---------------------------------------|:-------------------------:|:-------------------------:|
| fixed  | b_Intercept                            | 20.104 \[18.578, 21.703\] | 20.239 \[18.551, 21.925\] |
|        | b_adltotb_pm                           |  1.184 \[1.061, 1.307\]   |  1.131 \[1.038, 1.232\]   |
|        | b_adltotb_pmc                          |  0.515 \[0.365, 0.671\]   |  0.465 \[0.343, 0.588\]   |
|        | b_race3                                |  0.426 \[-0.134, 0.996\]  |  0.680 \[-0.944, 2.293\]  |
|        | b_race2                                |  0.235 \[-0.281, 0.706\]  | -1.535 \[-2.858, -0.238\] |
|        | b_race4                                |  0.681 \[-0.440, 1.846\]  |  3.887 \[0.464, 7.416\]   |
|        | b_socsuptot_pm                         | -0.295 \[-0.339, -0.251\] | -0.308 \[-0.360, -0.256\] |
|        | b_socsuptot_pmc                        | -0.096 \[-0.139, -0.052\] | -0.121 \[-0.174, -0.067\] |
|        | b_age_w0                               | -0.020 \[-0.039, -0.001\] | -0.019 \[-0.038, 0.001\]  |
|        | b_adltotb_pm × race3                   | -0.032 \[-0.321, 0.258\]  |                           |
|        | b_adltotb_pm × race2                   | -0.205 \[-0.443, 0.031\]  |                           |
|        | b_adltotb_pm × race4                   | -0.486 \[-1.112, 0.162\]  |                           |
|        | b_adltotb_pmc × race3                  | -0.163 \[-0.548, 0.222\]  |                           |
|        | b_adltotb_pmc × race2                  | -0.161 \[-0.475, 0.148\]  |                           |
|        | b_adltotb_pmc × race4                  |  0.278 \[-0.636, 1.205\]  |                           |
|        | sigma                                  |  3.194 \[3.105, 3.285\]   |  3.190 \[3.099, 3.285\]   |
|        | b_socsuptot_pm × race3                 |                           | -0.024 \[-0.165, 0.109\]  |
|        | b_socsuptot_pm × race2                 |                           |  0.140 \[0.029, 0.251\]   |
|        | b_socsuptot_pm × race4                 |                           | -0.319 \[-0.608, -0.026\] |
|        | b_socsuptot_pmc × race3                |                           |  0.016 \[-0.121, 0.160\]  |
|        | b_socsuptot_pmc × race2                |                           |  0.100 \[-0.010, 0.207\]  |
|        | b_socsuptot_pmc × race4                |                           |  0.113 \[-0.146, 0.371\]  |
| random | sd_ID\_\_Intercept                     |  3.330 \[3.198, 3.466\]   |  3.328 \[3.190, 3.459\]   |
|        | sd_ID\_\_adltotb_pmc                   |  0.870 \[0.705, 1.038\]   |  0.874 \[0.702, 1.037\]   |
|        | sd_ID\_\_socsuptot_pmc                 |  0.335 \[0.261, 0.405\]   |  0.335 \[0.257, 0.410\]   |
|        | cor_ID\_\_Intercept\_\_adltotb_pmc     |  0.213 \[0.067, 0.365\]   |  0.220 \[0.069, 0.369\]   |
|        | cor_ID\_\_Intercept\_\_socsuptot_pmc   | -0.214 \[-0.367, -0.060\] | -0.204 \[-0.359, -0.050\] |
|        | cor_ID\_\_adltotb_pmc\_\_socsuptot_pmc |  0.023 \[-0.336, 0.351\]  |  0.016 \[-0.345, 0.319\]  |
|        | Num.Obs.                               |           6829            |           6829            |
|        | R2                                     |           0.613           |           0.613           |
|        | R2 Adj.                                |           0.393           |           0.391           |
|        | R2 Marg.                               |           0.166           |           0.168           |
|        | ELPD                                   |         -19186.1          |         -19188.2          |
|        | ELPD s.e.                              |           79.4            |           79.5            |
|        | LOOIC                                  |          38372.2          |          38376.4          |
|        | LOOIC s.e.                             |           158.9           |           159.1           |
|        | WAIC                                   |          37824.4          |          37814.3          |
|        | RMSE                                   |           2.55            |           2.55            |
|        | r2.adjusted.marginal                   |           0.159           |           0.161           |

Table 2: Model coefficients for race interactions

Create figures

``` r
# Figures for model 1
plot(
    conditional_effects(m1),
    points = TRUE,
    point_args = list(
        size = 0.5, width = 0.05, height = 0.1)
)
```

![](final_files/figure-gfm/unnamed-chunk-9-1.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-2.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-3.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-4.png)<!-- -->

``` r
# Figures for model 2
plot(
    conditional_effects(m2),
    points = TRUE,
    point_args = list(
        size = 0.5, width = 0.05, height = 0.1)
)
```

![](final_files/figure-gfm/unnamed-chunk-9-5.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-6.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-7.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-8.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-9.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-10.png)<!-- -->

``` r
# Figures for model 3
plot(
    conditional_effects(m3),
    points = TRUE,
    point_args = list(
        size = 0.5, width = 0.05, height = 0.1)
)
```

![](final_files/figure-gfm/unnamed-chunk-9-11.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-12.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-13.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-14.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-15.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-16.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-17.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-18.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-19.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-20.png)<!-- -->

``` r
# Figures for model 4a
plot(
    conditional_effects(m4a),
    points = TRUE,
    point_args = list(
        size = 0.5, width = 0.05, height = 0.1)
)
```

![](final_files/figure-gfm/unnamed-chunk-9-21.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-22.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-23.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-24.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-25.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-26.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-27.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-28.png)<!-- -->

``` r
# Figures for model 4b
plot(
    conditional_effects(m4b),
    points = TRUE,
    point_args = list(
        size = 0.5, width = 0.05, height = 0.1)
)
```

![](final_files/figure-gfm/unnamed-chunk-9-29.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-30.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-31.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-32.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-33.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-34.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-35.png)<!-- -->![](final_files/figure-gfm/unnamed-chunk-9-36.png)<!-- -->
