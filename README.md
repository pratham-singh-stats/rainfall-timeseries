# 🌧️ Modeling Rainfall Patterns with El Niño Influence Using Time Series

![R](https://img.shields.io/badge/Language-R-276DC3?style=flat&logo=r)
![Status](https://img.shields.io/badge/Status-Completed-brightgreen)
![Domain](https://img.shields.io/badge/Domain-Environmental%20Statistics-2E75B6)

## 📌 Overview
This project investigates long-term Indian rainfall patterns (1950–2021) and quantifies the influence of the **El Niño–Southern Oscillation (ENSO)** using advanced time series modeling. The goal is to produce reliable 12-month-ahead forecasts applicable to **agricultural planning and water resource management**.

---

## 🎯 Objectives
- Decompose 70+ years of monthly rainfall data into trend, seasonal, and residual components
- Test for stationarity and identify appropriate model orders
- Compare four forecasting models: ARIMA, SARIMA, SARIMAX, and ETS
- Quantify the impact of the El Niño (ONI) index as an exogenous variable
- Evaluate forecast accuracy using AIC, BIC, and RMSE

---

## 📊 Dataset
| Attribute | Details |
|-----------|---------|
| Source | IMD (rainfall), NOAA (ONI index) |
| Period | January 1950 – December 2021 |
| Observations | ~864 monthly records |
| Variables | Monthly rainfall (mm), El Niño ONI index |

---

## 🔬 Methodology

### 1. Exploratory Data Analysis
- Visualized raw time series with LOESS trend overlay
- Plotted average monthly rainfall to confirm monsoon seasonality (Jun–Sep peak)
- Scatter plot of ONI index vs rainfall to assess ENSO correlation

### 2. STL Decomposition
- Separated series into **Trend + Seasonal + Remainder** components
- Confirmed strong seasonal structure (strength > 0.85)

### 3. Stationarity Testing
- Applied **Augmented Dickey-Fuller (ADF) test**
- Result: Series stationary after seasonal differencing

### 4. Model Comparison

| Model | AIC | BIC | Test RMSE |
|-------|-----|-----|-----------|
| ARIMA (auto-selected) | — | — | Baseline |
| SARIMA(1,1,1)(0,1,1)[12] | — | — | ↓ ~10% |
| **SARIMAX + ONI Index** | **Best** | **Best** | **↓ 15–25%** |
| ETS | — | — | ↓ ~8% |

> ✅ **SARIMAX** achieved the best forecast accuracy, confirming El Niño as a statistically significant predictor of rainfall variation.

### 5. Residual Diagnostics
- Ljung-Box test: residuals confirmed as white noise (p > 0.05)
- No significant autocorrelation remaining in residuals

---

## 📈 Key Findings
- Digital ENSO signals (ONI index) explain a meaningful portion of inter-annual rainfall variability
- SARIMAX reduced RMSE by **15–25%** compared to the baseline ARIMA model
- Peak rainfall months (July–August) showed highest forecast uncertainty
- Model outputs can directly inform **kharif crop planning** and **reservoir management**

---

## 🛠️ Tech Stack
- **Language:** R
- **Key Packages:** `forecast`, `tseries`, `tidyverse`, `ggplot2`, `lubridate`, `gridExtra`

---

## 🚀 How to Run
```r
# 1. Install required packages
install.packages(c("tidyverse", "lubridate", "tseries",
                   "forecast", "ggplot2", "gridExtra"))

# 2. Open rainfall_analysis.R in RStudio

# 3. Run the full script (Ctrl+Shift+Enter)
#    or source it:
source("rainfall_analysis.R")
```

---

## 📁 Repository Structure
```
01_rainfall_timeseries/
│
├── rainfall_analysis.R     # Full analysis script (EDA → Modeling → Forecasting)
├── README.md               # This file
└── plots/                  # Output plots (generated on run)
```

---

## 🔗 References
- India Meteorological Department: https://www.imd.gov.in
- NOAA Climate Prediction Center (ONI): https://www.cpc.ncep.noaa.gov
- Hyndman & Athanasopoulos (2021), *Forecasting: Principles and Practice* (3rd ed.)
