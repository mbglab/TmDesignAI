# TmDesignAI
Intelligent Design of *Escherichia coli* Terminators

This repository contains two separate models related to *E. coli* terminators:

1. **Terminator Strength Prediction Model**: A model designed to predict the strength of terminators based on given DNA sequence data.
2. **Terminator Generation Model**: A model aimed at generating new terminator sequences.

Each model is stored in its respective directory.

---

## Directory Structure

```
📂 Terminator_Strength_Prediction_Model/
   ├── features.txt  
   ├── seq606_RNAfold.txt   
   ├── SVR_model.py     

📂 Terminator_Generation_Model/
   ├── data/             
   ├── utils/
   ├── models.py     
   ├── wgan_gp_gene.py 
```

---

## Terminator Strength Prediction Model

### Description
This model is trained to predict the strength of terminators based on provided sequence data.

### Files
- `features.txt`: a dataset containing terminator sequences, features, and corresponding terminator strength values.
- `seq606_RNAfold.txt`: the file containing minimum free energy (MFE) values computed using RNAfold.
- `SVR_model.py`: the main script for data processing, feature extraction, and model training.

---

## Terminator Generation Model

### Description
This model is designed to generate novel terminator sequences.

### Files
- `data/`: a directory containing necessary datasets for training and validation.
- `utils/`: a directory containing utility scripts for data preprocessing, PyTorch operations, and visualization.
- `models.py`: a Python script defining the generator and discriminator models used in WGAN-GP.
- `wgan_gp_gene.py`: the main script for training and evaluating the WGAN-GP model.
---


## Dependencies

Both models may require specific Python libraries. Install dependencies using:

```bash
pip install -r requirements.txt
