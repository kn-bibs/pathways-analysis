#!/usr/bin/python3.6
import sys
from command_line import CLI


def run(argv):
    args = CLI().parse(argv[1:])
    print(args)
    return args.method.run(args.experiment)


if __name__ == '__main__':
    run(sys.argv)
