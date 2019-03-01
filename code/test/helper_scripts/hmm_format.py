import sys


def main():
    with sys.stdin as reader:
        print(reader.readline().strip())
        for line in reader:
            tokens = line.strip().split('\t')
            for i, t in enumerate(tokens):
                if i < 2:
                    print(t, end="\t")
                else:
                    print(f"{float(t):.6}", end="\t")
            print()


if __name__ == "__main__":
    main()
