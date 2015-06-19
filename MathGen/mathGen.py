#complexity may be multiplicative comp1*comp2*...+complexity
#or additive comp1+comp2+...+complexity

import numpy
import random

class math_obj:
    def __init__(self):
        pass
    
    def complexity(self):
        raise NotImplementedError("complexity not implemented")

    def evaluate(self, **kwargs):
        for key in kwargs:
            command = key + " = kwargs[\'" + key + "\']"
            exec(command)
        return eval(str(self))

class value(math_obj):
    def __init__(self, value, complexity=1):
        self._value = value
        self._complexity = complexity
    
    def __str__(self):
        return str(self._value)
        
    def complexity(self):
        return self._complexity
        
class function(math_obj):
    def __init__(self, func, arglist, complexity=1.0, iscompmulty = False):
        self._function = func
        self._arglist = arglist
        self._complexity = complexity
        self._iscompmulty = iscompmulty
        
    def __str__(self):
        output = self._function+"("
        first = True
        for arg in self._arglist:
            if first:
                first=False
            else:
                output = output + ","
            output = output + str(arg)
        output = output + ")"
        return output
        
    def complexity(self):
        if self._iscompmulty:
            comp = 1.0
            for arg in self._arglist:
                comp *= arg.complexity()
            return comp + self._complexity
        else:
            comp = self._complexity
            for arg in self._arglist:
                comp += arg.complexity()
            return comp

class operator(math_obj):
    def __init__(self, oper, pre, post, complexity=1.0, iscompmulty = False):
        self._operator = oper
        self._pre = pre
        self._post = post
        self._complexity = complexity
        self._iscompmulty = iscompmulty
        
    def __str__(self):
        return "(" + str(self._pre) + self._operator + str(self._post) + ")"
        
    def complexity(self):
        if self._iscompmulty:
            return self._pre.complexity()*self._post.complexity() + self._complexity
        else:
            return self._pre.complexity() + self._post.complexity() + self._complexity

class equation:
    def __init__(self, operators, functions, variables, value_complexity=1.0, seed=None):
        #operators, functions, varables should be dicts with the value (complexity, ismulti, numargs), (complexity, ismulti) or complexity.
        self.equation = None
        self.target = None
        self.operators = operators
        self.functions = functions
        self.variables = variables
        self.value_complexity = value_complexity
        
        self._max_depth = 20
        
        random.seed(seed)
        
    def random_equation(self, operator_prob, function_prob, variable_prob, random_mean, random_deviation):
        if operator_prob+function_prob+variable_prob > 1.0:
            raise ValueError("The summ input probibilitys must be lessthan or equal to one.")
        if operator_prob == 0.0:
            oper_val = -1.0
        else:
            oper_val = operator_prob
        if function_prob == 0.0:
            func_val = -1.0
        else:
            func_val = operator_prob + function_prob
        if variable_prob == 0.0:
            var_val = -1.0
        else:
            var_val = operator_prob + function_prob + variable_prob
        return self._set_random_element(oper_val, func_val, var_val, random_mean, random_deviation, 0)
        
    def _set_random_element(self, oper_val, func_val, var_val, random_mean, random_deviation, depth):
        element = None
        number = random.uniform(0.0, 1.0)
        
        if depth > self._max_depth:
            if number < 0.5:
                variable = random.choice(list(self.variables.keys()))
                element = value(variable, self.variables[variable])
            else:
                element = value(random.gauss(random_mean, random_deviation), self.value_complexity)
        else:
            if number < oper_val:
                oper = random.choice(list(self.operators.keys()))
                element = operator(oper, self._set_random_element(oper_val, func_val, var_val, random_mean, random_deviation, depth+1), self._set_random_element(oper_val, func_val, var_val, random_mean, random_deviation, depth+1), self.operators[oper][0], self.operators[oper][1])
            elif number < func_val:
                func = random.choice(list(self.functions.keys()))
                args = []
                for i in range(self.functions[func][2]):
                    args.append(self._set_random_element(oper_val, func_val, var_val, random_mean, random_deviation, depth+1))
                element = function(func, args, self.functions[func][0], self.functions[func][1])
            elif number < var_val:
                variable = random.choice(list(self.variables.keys()))
                element = value(variable, self.variables[variable])
            else:
                element = value(random.gauss(random_mean, random_deviation), self.value_complexity)
        return element
    
    def find_error(self):
        error = numpy.sum(numpy.power(self.equation.evaluate() - self.target, 2))
        return error

test = equation({"+":(1, False), "-":(1, False), "*":(2, False)}, {"numpy.cos":(3, False, 1), "numpy.sin":(3, False, 1), "numpy.power":(3, True, 2)}, {"2":1, "numpy.e":1, "numpy.pi":1}, 1)
for i in range(3):
    equ = test.random_equation(0.2, 0.2, 0.1, 10.0, 5)
    print(str(equ), equ.complexity())
    print(eval(str(equ)))
