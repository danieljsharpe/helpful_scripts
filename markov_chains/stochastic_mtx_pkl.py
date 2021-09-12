# very simple script to pickle a toy DTMC

import numpy as np

T = np.array([[0.50,0.20,0.15,0.15,0.00],
              [0.15,0.75,0.10,0.00,0.00],
              [0.20,0.10,0.45,0.10,0.15],
              [0.05,0.00,0.25,0.70,0.00],
              [0.00,0.00,0.10,0.00,0.90]],dtype=float)

T.dump("transnmtx.pkl")
