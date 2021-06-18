import sys

from Public.UploadAlgorithm import uploadAlgorithm
from Public.SaveApproximationSet import saveApproximationSet

if __name__ == '__main__':
    if (str(sys.argv[1]) == '--help'):
        f = open('../README.txt',"r")
        contents = f.read()
        f.close()
        print(contents)
    else:
        if (len(sys.argv) != 8):
            sys.exit("Incorrect number of arguments. Use: main.py --help")
        algorithm = str(sys.argv[1])
        H1 = int(sys.argv[2])
        H2 = int(sys.argv[3])
        problem = str(sys.argv[4])
        m = int(sys.argv[5])
        max_generations = int(sys.argv[6])
        runs = int(sys.argv[7])
        
        valid_algorithms = ['NSGA-III', 'NSGA-III-RSE', 'NSGA-III-GAE', 'NSGA-III-COU', 
                            'NSGA-III-PT', 'NSGA-III-MPT', 'NSGA-III-GPT', 'NSGA-III-KRA']
        if (algorithm not in valid_algorithms):
            sys.exit("Invalid MOEA name. Use: main.py --help")
        if (H1 <= 0):
            sys.exit("Invalid value for the number of divisions per objective function for the boundary layer. Use: main.py --help")
        if (H2 < 0):
            sys.exit("Invalid value for the number of divisions per objective function for the inner layer. Use: main.py --help")
        valid_problems = ['DTLZ1', 'DTLZ2', 'DTLZ3', 'DTLZ4', 'DTLZ5', 'DTLZ6', 'DTLZ7', 
                          'DTLZ1_MINUS', 'DTLZ2_MINUS', 'DTLZ3_MINUS', 'DTLZ4_MINUS', 
                          'DTLZ5_MINUS', 'DTLZ6_MINUS', 'DTLZ7_MINUS']
        if (problem not in valid_problems):
            sys.exit("Invalid MOP name. Use: main.py --help")
        if (m < 2):
            sys.exit("Invalid value for the number of objective functions. Use: main.py --help")
        if (max_generations < 0):
            sys.exit("Invalid value for the maximum number of generations. Use: main.py --help")
        if (runs <= 0):
            sys.exit("Invalid value for the number of independent runs. Use: main.py --help")
        
        if (algorithm == 'NSGA-III-RSE'):
            flag = 1
        elif (algorithm == 'NSGA-III-GAE'):
            flag = 2
        elif (algorithm == 'NSGA-III-COU'):
            flag = 3
        elif (algorithm == 'NSGA-III-PT'):
            flag = 4
        elif (algorithm == 'NSGA-III-MPT'):
            flag = 5
        elif (algorithm == 'NSGA-III-GPT'):
            flag = 6
        elif (algorithm == 'NSGA-III-KRA'):
            flag = 7
        for i in range(1, runs+1):
            print('Algorithm:', algorithm, '| H1:', H1, '| H2:', H2, '| Problem:', problem, 
                  '| Objectives:', m, '| Generations:', max_generations, '| Run:', i)
            if (algorithm == 'NSGA-III'):
                main = uploadAlgorithm(algorithm)
                P = main(H1, H2, problem, m, max_generations)
            else:
                main = uploadAlgorithm('NSGA-III-K')
                P = main(H1, H2, flag, problem, m, max_generations)
            saveApproximationSet(P.obj, algorithm, problem, i, "save_all")
