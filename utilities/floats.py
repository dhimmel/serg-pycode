from __future__ import division

import math

def sequence(start, stop, step):
    """
    Adapted from
    carlosvega:
    This is my solution to get ranges with float steps.
    Using this function it's not necessary to import numpy, nor install it.
    I'm pretty sure that it could be improved and optimized. Feel free to do it and post it here.

    http://stackoverflow.com/questions/477486/python-decimal-range-step-value/20549652#20549652
    """
    digits = int(round(math.log(10000, 10)))+1 #get number of digits
    magnitude = 10**digits
    stop = int(magnitude * stop) #convert from
    step = int(magnitude * step) #0.1 to 10 (e.g.)

    start = 10 ** digits * start

    data = list()   #create array

    #calc number of iterations
    end_loop = int((stop-start) // step)

    acc = start

    for i in xrange(0, end_loop + 1):
        data.append(acc / magnitude)
        acc += step

    return data


if __name__ == '__main__':
    print sequence(1, 2.1, 0.01)
    print sequence(1, 2.1, 0.1)
    print sequence(0, 1.1, 0.1)
    print sequence(-1, 0.1, 0.1)
    print sequence(0.0, 5.0, 1.0)
