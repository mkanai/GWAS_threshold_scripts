#!/usr/bin/env python

import sys

def main():
    argvs = sys.argv
    argc = len(argvs)
    start = int(argvs[1])
    end = int(argvs[2])
    n = (end - start + 1) // int(argvs[3])
    if n <= 1:
        print start, end
        return

    for i in range(start, end, n):
        print i, i+n-1

if __name__ == '__main__':
    main()
