import argparse

from src.agnostic_search.driver import main as driver


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("mgfQueryFile", type=str)
    parser.add_argument("mgfTargetFile", type=str)
    parser.add_argument("IdFile", type=str)
    parser.add_argument("-p", "--max-peaks-per-scan", type=int, default=20)
    parser.add_argument("-e", "--mz-error", type=int, default=20)
    parser.add_argument(
        "-t", "--mz-error-type", type=str, choices=["ppm"], default="ppm"
    )
    parser.add_argument("-O", "--output-directory", type=str, default="output")
    parser.add_argument("-f", "--filter", action="store_true")

    args = parser.parse_args()
    driver(args)


if __name__ == "__main__":
    main()
