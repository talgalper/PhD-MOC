import pandas as pd
import os
import shutil
from tqdm import tqdm

structures = os.listdir("/home/ubuntu/Desktop/pocketminer/results/txt_results/")

if os.path.isdir("/home/ubuntu/Desktop/pocketminer/results/small_structures"):
    print("Destination folder found, moving small structures...")
else:
    print("destination folder does not exist, creating now...")
    print("Moving small structures...")
    os.mkdir("/home/ubuntu/Desktop/pocketminer/results/small_structures")

count = 0

for i in tqdm(structures):
    data = pd.read_table(f"/home/ubuntu/Desktop/pocketminer/results/txt_results/{i}", header=None)
    length = len(data)
    if length < 80:
        shutil.move(src = f"/home/ubuntu/Desktop/pocketminer/results/txt_results/{i}", dst = f"/home/ubuntu/Desktop/pocketminer/results/small_structures/{i}")
        count = count + 1

print(f"{count} small structres (<80 aa) moved")