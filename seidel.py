from random import randrange
import numpy
import copy

#https://www.pythonpool.com/gaussian-elimination-python/
def gauss(constraints):
    matrix_copy = constraints.copy()
    n = len(matrix_copy)
    x = numpy.zeros(n)
    for i in range(n):
        for j in range(i+1, n):
            ratio = matrix_copy[j][i]/matrix_copy[i][i]
            for k in range(n+1):
                matrix_copy[j][k] = matrix_copy[j][k] - ratio * matrix_copy[i][k]
    x[n-1] = matrix_copy[n-1][n]/matrix_copy[n-1][n-1]
    for i in range(n-2,-1,-1):
        x[i] = matrix_copy[i][n]
    for j in range(i+1,n):
        x[i] = x[i] - matrix_copy[i][j]*x[j]
    x[i] = x[i]/matrix_copy[i][i]
    return x.tolist()

def seidel(objective_function: list, constraints: list, signs:list):
    #ONE VARIABLE
    if len(constraints[0]) == 2:
        #Treat all inequalities in Ax ≤ b as equalities and solve for x in each.
        #Return the point out of these that maximizes c.
        results = []
        for const in constraints:
            try:
                results.append(((const[-1]/const[0]) * objective_function[1]))
            except ZeroDivisionError:
                pass
        maximum = max(results)
        final_result_x1 = (constraint[-1]-(constraint[1]*maximum))/constraint[0]
        return final_result_x1
    
    #ONE CONSTRAINT
    if len(constraints) == 1:
        x2 = constraints[0][-1]/constraints[0][1]
        return [0.0,x2]

    #AS MANY CONSTRAINTS AS VARIABLES
    if len(constraints[0])-1 == len(constraints):
        #Solve the linear system Ax = b using Gaussian elimination and return the resulting point.
        matrix = copy.deepcopy(constraints)
        try:
            point = gauss(matrix)
        except ZeroDivisionError:
            return None
        opt = 0
        for i in range(len(objective_function)):
            opt += point[i]*objective_function[i]
        if all(i >= 0 for i in point):
            return [point]
        else:
            return None

    #DIFFERENT NUMBER OF CONSTRAINTS AND VARIABLES
    else:
        x = None
        while x == None:
            h = randrange(0,len(constraints))
            constraint = constraints[h]
            constraints_modified = copy.deepcopy(constraints)
            constraints_modified.pop(h)
            x = seidel(objective_function, constraints_modified, signs)
            if x is not None:
                break
        equation = 0
        for i in range(len(objective_function)):
            equation += constraint[i]*x[0][i]
        if (signs[h] == '<=') and (equation <= constraint[-1]):
            return x
        elif (signs[h] == '>=') and (equation >= constraint[-1]):
            return x
        elif (signs[h] == '=') and (equation == constraint[-1]):
            return x
        else:
            signs[h] = '='
            seidel(objective_function, constraints, signs)

# #TEST 1
# constraints = [[-1.0,1.0,0.0],[1.0,1.0,8.0],[1.0,1.0,6.0]]
# #either <=, >= or =
# signs = ['<=', '<=', '<=']
# objective_function = [0,1]

# #TEST 2
# constraints = [[-3.0,1.0,4.0],[2.0,1.0,4.0],[6.0,2.0,24.0],[-2.0,1.0,0.0]]
# #either '<=', '>=' or '='
# signs = ['<=', '<=', '<=', '<=']
# objective_function = [-1,1]

constraints = [[1.0,1.0,0.0],[-1.0,-1.0,2.0]]
#either '<=', '>=' or '='
signs = ['<=', '<=']
objective_function = [-1,1]


#Fix for recursion error
result = None
while result is None:
    try:
        result = seidel(objective_function, constraints, signs)
    except RecursionError:
        pass
print('Result: ' + str(result))
