


def main():
    parser = argparse.ArgumentParser(
        description="Recursively aggregate a single resolution cooler file into a multi-resolution file.")
    parser.add_argument(
        "cooler_file",
        help="Cooler file",
        metavar="COOLER_PATH")
    parser.add_argument(
        "--out", "-o",
        help="Output multires file")
    parser.add_argument(
        "--n_cpus", "-n",
        help="Number of cpus to use in process pool (Default=1, i.e. no pool)",
        default=1,
        type=check_ncpus)
    parser.add_argument(
        "--chunk-size", "-c",
        help="Chunk size",
        default=int(10e6),
        type=int)

    parser.add_argument('--balance', dest='balance', action='store_true')
    parser.add_argument('--no-balance', dest='balance', action='store_false')
    parser.set_defaults(balance=False)

    args = vars(parser.parse_args())



if __name__ == '__main__':
    main()
