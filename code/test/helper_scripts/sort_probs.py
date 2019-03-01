import sys


def main():
    with sys.stdin as reader:
        sorting = []
        for line in reader:
            tokens = line.strip().split('\t')
            if sorting == []:
                keys = [(t.split(':')[0], i) for i, t in enumerate(tokens[2:])]
                inds = [2 + x[1] for x in sorted(keys, key=lambda k:k[0])]
                sorting = [0, 1] + inds

            print("\t".join([tokens[s] for s in sorting]))


if __name__ == "__main__":
    main()
