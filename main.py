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

def median(b):
    z = sorted(b)
    if (len(z) % 2 == 1):
        return z[len(z) / 2]
    else:
        return (z[len(z) / 2] + z[len(z) / 2 - 1]) * 1. / 2

def WalshMedian(b):
    a = [(b[i] + b[j]) * 1. / 2 for i in xrange(len(b)) for j in xrange(i, len(b))]
    return median(a);

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

# for \theta = 0.5
def signsStatistic(b):
    return (sum([1 for i in b if i > 0.5]) - len(b) * 1. / 2) / math.sqrt(len(b) * 1. / 4)

# for \theta = 0.5 (suggested by Willcockson).
def signRanksStatistic(b):
    z = sorted(b, key = lambda x: abs(x - 0.5))
    t = sum([(i + 1) for i in xrange(len(z)) if z[i] > 0.5])
    return (t - (len(b) * (len(b) + 1) * 1. / 4)) / math.sqrt(len(b) * (len(b) + 1) * (2 * len(b) + 1) * 1. / 24)

def OrlovStatistic(b):
    sum = 0
    for i in xrange(len(b)):
        sum += ((1 - empericDistributionFunction(b, b[i]) - empericDistributionFunction(b, -b[i])) ** 2)
    return sum

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

def oneSidedSmirnovSignificancyLevel(x, y):
    return math.exp(-2 * oneSidedSmirnovStatistic(x, y))

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

#len(x) must be greater than 1.
def correctedSampleVariance(x):
    return sampleVariance(x) * 1.  * len(x) / (len(x) - 1)

def totalSampleVariance(x, y):
    return (sampleVariance(x) * len(x) + sampleVariance(y) * len(y)) * 1. / (len(x) + len(y) - 2)

def StudentStatistic(x, y):
    return math.sqrt(len(x) * len(y) * 1. / (len(x) + len(y))) * (sampleMoment(y, 1) - sampleMoment(x, 1)) * 1. / totalSampleVariance(x, y)

def chiSquareStatisticUniformIsHomogenious(x, y, N):
    chiSquare = 0
    for j in xrange(N):
        temp = generalFrequencyOfHitsInGapsUniform(x, y, N, j)
        chiSquare += ((numberOfHitsInGapsUniform(x, N, j) - len(x) * temp) ** 2) * 1. / (len(x) * temp)
        chiSquare += ((numberOfHitsInGapsUniform(y, N, j) - len(y) * temp) ** 2) * 1. / (len(y) * temp)
    return chiSquare

def chiSquareStatisticIsUniform(x, N):
    chiSquare = 0
    for j in xrange(N):
        temp = numberOfHitsInGapsUniform(x, N, j)
        chiSquare += ((temp - (len(x) * 1. / N)) ** 2) * 1. / (len(x) * 1. / N)
    return chiSquare

def numberOfHitsInGapsUniform(x, N, j):
    num = 0
    for i in xrange(len(x)):
        if ((x[i] > j * 1. / N) and (x[i] <= (j + 1) * 1. / N)):
            num += 1
    return num

def generalFrequencyOfHitsInGapsUniform(x, y, N, j):
    return (numberOfHitsInGapsUniform(x, N, j) + numberOfHitsInGapsUniform(y, N, j)) / (len(x) * 1. + len(y) * 1.)

col = 15
sampleSize = 100

x = getColumnWithNumber(col, sampleSize)
y = getColumnWithNumber(35, sampleSize)

print chiSquareStatisticUniformIsHomogenious(x, y, 7)
print chiSquareStatisticIsUniform(x, 7)
