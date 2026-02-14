# Polynomial NARX Model

This project focuses on system identification and modeling using real input-output experimental data. The implementation is developed in MATLAB and uses provided dataset iddata-05.mat for parameter estimation and model validation.

Implementation of a polynomial Nonlinear ARX (NARX) model for black-box dynamic system modeling, including one-step-ahead prediction, free-run simulation, and MSE-based model evaluation.

---

##  Project Overview

This project focuses on identifying a nonlinear dynamic system using a polynomial NARX (Nonlinear AutoRegressive with eXogenous input) structure.

The model is developed using measured input-output data and aims to approximate system behavior without prior knowledge of the internal dynamics (black-box modeling).

---

##  Model Description

The NARX model structure is defined as:

y(k) = f(y(k-1), ..., y(k-na), u(k-1), ..., u(k-nb)) + e(k)

Where:
- y(k) = system output
- u(k) = system input
- na, nb = model orders
- m = polynomial degree
- e(k) = measurement noise

The polynomial degree controls model complexity and is selected using validation data to avoid underfitting and overfitting.

---

## Technologies
- MATLAB

---

##  Features

- Polynomial basis expansion
- Configurable model orders (na, nb)
- Adjustable polynomial degree (m)
- One-Step-Ahead Prediction mode
- Free-Run Simulation mode
- Performance evaluation using Mean Squared Error (MSE)
- Model tuning and overfitting analysis

---

##  Evaluation Metrics

Model performance is evaluated using:

- **Prediction MSE** – short-term accuracy
- **Simulation MSE** – long-term accuracy

The optimal polynomial degree is selected based on validation MSE.

---

## 

