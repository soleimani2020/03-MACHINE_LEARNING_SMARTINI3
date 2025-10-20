import os
import numpy as np
from itertools import islice

def get_enthalpy(log_file="log_Enthelpy.txt", line_number=6, value_column=1):
    if not os.path.isfile(log_file):
        raise FileNotFoundError(f"Log file '{log_file}' not found.")

    with open(log_file) as f:
        line = next(islice(f, line_number - 1, line_number))
        values = line.split()
        return float(values[value_column])

enthalpy = get_enthalpy()
print(enthalpy)
