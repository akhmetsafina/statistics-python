import random
import math

def getRandomVector(k):
    return [random.random() for i in xrange(k)]

def sampleMoment(b, k):
    return 1. * sum(list(map(lambda x: x ** k, b))) / len(b)

def sampleCentralMoment(b, k):
    return sampleMoment(b, k) - sampleMoment(b, 1)

def sampleVariance(b):
    return sampleMoment(b, 2) - (sampleMoment(b, 1) ** 2)

def assymetry(b):
    return sampleCentralMoment(b, 3) / (sampleCentralMoment(b, 2) ** 1.5)

def assymetryPrime(b):
    n = len(b)
    return assymetry(b) * math.sqrt(n * (n - 1) * 1.) / (n - 2.)

def varianceOfAssymetryPrime(b):
    n = len(b)
    return 6. * n * (n - 1.) * 1. / (n - 2.) / (n + 1.) / (n + 3.)

def Giri(b):
    n = len(b)
    ans = 0
    avg = sampleMoment(b, 1)
    for i in xrange(n):
        ans += abs(b[i] - avg)
    ans /= n
    ans /= math.sqrt(sampleVariance(b))
    return ans

def meanValueGiri(b):
    n = len(b)
    return math.sqrt(2. / math.pi) * (1 + 2. / (8. * n - 9))

def varianceGiri(b):
    n = len(b)
    return 1. / n * ((1 - 3. / math.pi) - 1. / (4. * math.pi * n))

def GiriStarStatistic(b):
    return (Giri(b) - meanValueGiri(b)) / (math.sqrt(varianceGiri(b)))

def kolmogorovStatistic(c):
    #print "!!!"
    b = c[:]
    b.sort()
    maxVal = -100000000
    for i in xrange(len(b)):
    #    print i, len(b), b[i], abs(((i + 1.) / len(b)) - b[i]), abs(b[i] - (i * 1.) / len(b))
        maxVal = max(maxVal, max(abs((i + 1.) / len(b) - b[i]), abs(b[i] - (i * 1.) / len(b))))
    #print maxVal
    return maxVal * math.sqrt(len(b))

def omegaSquareStatistic(c):
    b = c[:]
    b.sort()
    n = len(b)
    ans = 0.
    ans += 1. / (12 * n)
    for i in xrange(n):
        ans += (b[i] - (2. * (i + 1) - 1.) / (2. * n)) ** 2
    return ans

def assymetryStatistic(b):
    return assymetryPrime(b) / math.sqrt(varianceOfAssymetryPrime(b))

def getColumnWithNumber(col, firstN):
    a = []
    with open('input.txt') as f:
        a = f.readlines()
    c = []
    for i in xrange(1, min(len(a), firstN + 1)):
        c.append(float(((a[i].split())[col - 1]).replace(',', '.')))
    return c

def simulateWithReversedFunction(reversedFunction, b):
    return list(map(reversedFunction, b))

def empericDistributionFunction(b, x):
    cnt = 0
    for i in b:
        if x >= i:
            cnt += 1
    return cnt * 1. / len(b)

def SmirnovDPlus(x, y):
    maxv = -5
    z = sorted(x)
    for i in xrange(len(z)):
        maxv = max(maxv, (i + 1) * 1. / len(z) - empericDistributionFunction(y, z[i]))
    return maxv

def SmirnovDMinus(x, y):
    maxv = -5
    z = sorted(y)
    for i in xrange(len(z)):
        maxv = max(maxv, (i + 1) * 1. / len(z) - empericDistributionFunction(x, z[i]))
    return maxv

def oneSidedSmirnovStatistic(x, y):
    return math.sqrt(len(x) * len(y) * 1. / (len(x) + len(y))) * SmirnovDPlus(x, y)

def SmirnovD(x, y):
    smirnovDPlus = SmirnovDPlus(x, y)
    smirnovDMinus = SmirnovDMinus(x, y)
    return max(smirnovDPlus, smirnovDMinus)

def SmirnovStatistic(x, y):
    smirnovD = SmirnovD(x, y)
    someSqrt = math.sqrt(len(x) * len(y) * 1. / (len(x) + len(y)))
    return smirnovD * someSqrt

def MannWhitney(x, y):
    cnt = 0
    for i in x:
        for j in y:
            if i < j:
                cnt += 1
    return cnt

def MannWhitneyExpectation(x, y):
    return len(x) * len(y) * 1. / 2

def MannWhitneyVariance(x, y):
    return len(y) * len(x) * 1. * (1 + len(x) + len(y)) / 12

def MannWhitneyStatistic(x, y):
    return (MannWhitney(x, y) - MannWhitneyExpectation(x, y)) * 1. / math.sqrt(MannWhitneyVariance(x, y))

col = 15
sampleSize = 100

c = getColumnWithNumber(col, sampleSize)
c = simulateWithReversedFunction(lambda x: math.sqrt(x), c)

#print kolmogorovStatistic(c)
#print omegaSquareStatistic(c)
#print sampleCentralMoment([1, 2, 3], 2)
#print assymetryStatistic(c)
#print GiriStarStatistic(c)


a = getColumnWithNumber(1, sampleSize)
#c = getColumnWithNumber(col - 1, sampleSize)
b = list(map(lambda x, y, z: max(x, y, z), getColumnWithNumber(col + 1, sampleSize), getColumnWithNumber(col + 2, sampleSize), getColumnWithNumber(col + 3, sampleSize)))



print SmirnovStatistic(a, b)
print SmirnovStatistic(a, c)
print MannWhitneyStatistic(a, b)
print oneSidedSmirnovStatistic(a, c)