#!/usr/bin/env python2.7

import math
import copy
import numbers

# Implemented:
# Initialization from float (or int, doesn't matter), zero initialization value is treated specially.
# add, + (__add__, __radd__), addition with zero on either side.
# mul, * (__mul__, __rmul__), multiplication with zero on either side. 
# div, / (__div__, __rdiv__), division of zero is taken into account.

def isRegularNum(arg):
    if isinstance(arg, numbers.Integral) or isinstance(arg, float):
        return True
    else:
        return False




class bvvBigNum:
    """Numbers with large exponents"""
    
    def ln(self):
        return math.log(self.Mant) + math.log(10) * self.Exp

    def lg(self):
        return math.log10(self.Mant) + self.Exp

    def isSameType(self, obj):
        if obj.__class__.__name__ == self.__class__.__name__:
            return True
        else:
            return False

    def Convert(self, obj):
        if self.isSameType(obj):
            return obj
        elif isRegularNum(obj):
            return bvvBigNum(obj)
            
    def __init__(self, number):
        if number == 0.0:
            self.Exp = 0
            self.Mant = 0.0
        else:
            self.Exp = int(math.floor(math.log10(number))) # + 1 -- if you want the format 0.xxx * 10^yyy instead of x.xx * 10^yyy
            # print "Exp: ", self.Exp
            if self.Exp < -323:
                # 10**(-324) is floating point zero:
                self.Exp = 0
                self.Mant = 0
            else:
                self.Mant = float(number) / 10**int(self.Exp)
       
    def string(self):
        return str(self.Mant) + " " + str(self.Exp)
    
    def add(self, other0):
        other = self.Convert(other0)

        if self.Mant == 0.0:
            self.Mant = other.Mant
            self.Exp = other.Exp
            return

        if other.Mant == 0.0: return self

        if self.Exp < other.Exp:
           NumSmallestExponent = self
           NumLargestExponent = other
        else:
           NumSmallestExponent = other
           NumLargestExponent = self
        ExpDiff = NumLargestExponent.Exp - NumSmallestExponent.Exp

        #if ExpDiff > 10:
        #    print "ExpDiff: ", self.string(), other.string()

        if ExpDiff > 250:
            NewMantissa = NumLargestExponent.Mant
        else:
            NewMantissa = NumLargestExponent.Mant +  NumSmallestExponent.Mant / math.pow(10, ExpDiff)
        # NewMantissa may become too small (< 1) or too large (>= 10).
        NumNewMantissa = bvvBigNum(NewMantissa)
        self.Mant = NumNewMantissa.Mant
        self.Exp = NumLargestExponent.Exp + NumNewMantissa.Exp
    
    def __add__(self, other):
        result = copy.deepcopy(self)
        result.add(other)
        return result


    def __radd__(self, other):
        result = self.Convert(other)
        result.add(self)
        return result

    def mul(self, other0):
        other = self.Convert(other0)
        if self.Mant == 0: return
        if other.Mant == 0:
            self.Mant = 0
            return
            
        NewMantissa = self.Mant * other.Mant
        # NewMantissa may become too small (< 1) or too large (>= 10).
        NumNewMantissa = bvvBigNum(NewMantissa)
        self.Mant = NumNewMantissa.Mant
        self.Exp = self.Exp + other.Exp + NumNewMantissa.Exp
        

    def __mul__(self, other):
        result = copy.deepcopy(self)
        result.mul(other)
        return result

    def __rmul__(self, other):
        result = self.Convert(other)
        result.mul(self)
        return result

    
    def div(self, other0):
        other = self.Convert(other0)
        if self.Mant == 0:
            return
        NewMantissa = self.Mant / other.Mant
        # NewMantissa may become too small (< 1) or too large (>= 10).
        NumNewMantissa = bvvBigNum(NewMantissa)
        self.Mant = NumNewMantissa.Mant
        self.Exp = self.Exp - other.Exp + NumNewMantissa.Exp

    def __div__(self, other):
        result = copy.deepcopy(self)
        result.div(other)
        return result

    def __rdiv__(self, other):
        result = self.Convert(other)
        result.div(self)
        return result
        
    def float(self):
        return self.Mant * 10**self.Exp
    

################################
# A module test until the end: #
################################
#Num1 = 101.0
#Num2 = 5.0
#
#x1 = bvvBigNum(Num1)
#x2 = bvvBigNum(Num2)
#print "The contents before is ", x1.string(), x2.string()
#
#x1.div(Num2)
#print "Div in place: ", x1.string()
#
#print "Div should be: ", Num1 / Num2

#print "The contents of addition on the fly ", (Num2 + x1).string()
#print "Addition should be ", Num1 + Num2

#x2.add(x1)
#print "Adding in place: ", x2.string()

#print x1.ln(), math.log(Num1)
#print x1.lg(), math.log10(Num1)
