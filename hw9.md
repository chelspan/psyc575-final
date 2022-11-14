HW 9
================
Chelsey Pan + Mengzhao Yan
2022-11-13

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

![](hw9_files/figure-gfm/exploration-1.png)<!-- -->

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
  mutate(across(c(adltotb, socsuptot),
         list("pm" = ~ mean(., na.rm = TRUE),
                       "pmc" = ~ . - mean(., na.rm = TRUE)))) %>%
  mutate(race = as.factor(race))
```

**Model Equations**

**Lvl 1:**

$$
\text{deptot}_{ti} = \beta_{0i} + \beta_{1i} \text{adltotb_pmc}_{ti} + e_{ti}
$$

**Lvl 2:**

$$
\beta_{0i} = \gamma_{00} + \gamma_{01} \text{adltotb_pm}_i + \gamma_{02} \text{race}_i + \gamma_{03} \text{adltotb_pm}_i \times
 \text{race}_i + u_{0i}
$$

$$\beta_{1i} = \gamma_{10} + \gamma_{11}\text{race}_i + u_{1i}$$

## Preliminary analysis

Base model: Is there an association between reported physical
difficulties and depressive symptoms across the 3 waves, and does this
interact with participant race (treated as a lvl 2 predictor in this
case)

``` r
# Model with just functional health difficulties
m1 <- brm(deptot ~ (adltotb_pm + adltotb_pmc) * race + (adltotb_pmc | ID),
          data = df_long,
          seed = 123,
          file = 'dep_physicalhealth_race')

# Get model summary
summary(m1)
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
msummary(m1,
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

![](hw9_files/figure-gfm/figures-1.png)<!-- -->

``` r
# Figure showing association between adltotb_pmc and depression symptoms
# split by race
figures[5]
```

    ## $`adltotb_pmc:race`

![](hw9_files/figure-gfm/figures-2.png)<!-- -->

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
