import os

# pip install snakeviz
os.system("py -3.7 -O -m cProfile -o rendeer.prof rendeer.py")
os.system("snakeviz rendeer.prof")
