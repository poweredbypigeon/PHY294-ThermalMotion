import pandas as pd 
import os 
import uncertainties 
import math

def get_x2_series(df):
    x_series = df.iloc[:, 0]
    x_series -= x_series[0]
    x_series = x_series * x_series 

    return x_series

def get_distances(df):
    distances = []
    x_series = df.iloc[:, 0]
    y_series = df.iloc[:, 1]
    for i in range(1, len(x_series)):
        distances.append(math.sqrt((x_series[i] - x_series[i-1])**2 + (y_series[i] - y_series[i-1])**2))
    return distances

SOURCE  = "C:\\Users\\compl\\Documents\\Engineering Science\\Y2S2\\PHY294H1\\Labs\\L3+4\\good_data"

filenames = os.listdir(SOURCE)
# everything in good-data folder

# Load the data, skipping the header row

x2_avg = [] # average x^2 values

distances = [] # distances between the two points

s = None 
s_total = None 

for filename in filenames:
    df = pd.read_csv(os.path.join(SOURCE, filename), delimiter="\t", skiprows=2)
    # print("Number of rows", len(df), "for", filename)
    print(filename)
    s = get_x2_series(df)
    distances.extend(get_distances(df))
    if s_total is None:
        s_total = s
    else:
        s_total += s

s_total = s_total / len(filenames)

print(s_total)

print(distances)
""" 
def pre_process(df):
    x_series = df.iloc[:, 0]
    x_series -= x_series[0]
    x_series = x_series * x_series 

    y_series = df.iloc[:, 1]
    y_series -= y_series[0]
    y_series = y_series * y_series

    return x_series.tolist(), y_series.tolist()

# Extract columns into lists
# xs = x_series.tolist()
# ys = y_series.tolist()




# print(xs[:5], ys[:5])  # Print first few values to check

 """