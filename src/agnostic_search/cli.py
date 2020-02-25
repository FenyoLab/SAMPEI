import argparse

parser = argparse.ArgumentParser()
parser.add_argument("MGF Query File", type=str)
parser.add_argument("MGF Target File", type=str)
parser.add_argument("Max Peaks per Target", type=int)
parser.add_argument("MZ Error", type=int)
parser.add_argument("MZ Error Type", type=str, choices=["ppm"])
parser.add_argument("ID File", type=str)
parser.add_argument("Output Directory", type=str)
parser.add_argument("filter", type=bool)

args = parser.parse_args()
print(args)
