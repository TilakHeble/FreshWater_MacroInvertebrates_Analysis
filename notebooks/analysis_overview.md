# Freshwater Macroinvertebrates Analysis – Overview

## 📌 Project Summary

This project analyses long-term freshwater macroinvertebrate data from England to understand biodiversity trends over time. The dataset spans multiple decades and includes complex ecological observations collected using different methodologies.

---

## ⚠️ Key Challenge

The data presents several real-world challenges:

* Mixed data types (categorical abundance bands vs numeric counts)
* Missing and inconsistent observations
* Changes in sampling methodology over time
* Duplicate samples and irregular recording

---

## 🧹 Approach

A structured data pipeline was developed to:

* Clean and standardise raw datasets
* Classify observations into **bin-based** and **count-based** data
* Integrate multiple data sources (sites, metrics, taxa)
* Build modelling-ready datasets

---

## 📊 Modelling Strategy

Multiple statistical approaches were applied:

* Presence/Absence (Binomial GAM)
* Ordered categorical models
* Censored Poisson models
* Censored Normal models

These approaches were used to account for uncertainty and mixed data formats.

---

## 📈 Key Insights

* Evidence of long-term changes in macroinvertebrate abundance
* Differences between seasonal patterns (spring vs autumn)
* Importance of handling data uncertainty in ecological modelling

---

## 🧠 Key Learning

This project demonstrates the importance of:

* Robust data cleaning for real-world datasets
* Careful handling of mixed data types
* Using appropriate statistical models for uncertain data

---

## 🔗 Full Analysis

For complete details, see:

* `scripts/` for full pipeline
* `outputs/report.pdf` for full report

